#  This script processes ITS2 amplicons sequenced using an Illumina Miseq sequencer
#  to characterize the community of parasitic strongyle nematodes infecting horses, as described in:
#
#  Poissant J, Gavriliuc S, Bellaw J, Redman EM, Avramenko RW, Robinson D, 
#  Workentine ML, Shury TK, Jenkins EJ, McLoughlin PD, Nielsen MK, Gilleard JS. 
#  A repeatable and quantitative DNA metabarcoding assay to characterize mixed 
#  strongyle infections in horses. International Journal for Parasitology.
# 
#  *Please acknowledge the above-mentioned article when reusing this code or associated data. 
#  The following article should also be cited when re-using the curated equine parasitic nematodes ITS2 
#  database provided with this script: 
#
#  Workentine ML, Chen R, Zhu S, Gavriliuc S, Shaw N, de Rijke J, Redman EM, Avramenko RW, Wit J, 
#  Poissant J, Gilleard JS (2020) A database for ITS2 sequences from nematodes. BMC Genetics. 21:74.
#
#   This code implements a modifed verion of the DADA2 pipeline available at: 
#   https://benjjneb.github.io/dada2/ITS_workflow.html
#   https://benjjneb.github.io/dada2/tutorial.html
#
#   The following R packages are needed: DADA2, ShortRead, Biostrings
#
#   The following softwares are also required: cutadapt, ITSx
#   (https://cutadapt.readthedocs.io/en/stable/installation.html)
#   (https://microbiology.se/software/itsx/)
#
# OVERVIEW:
#   1. Remove primers using Cutadapt
#   2. Quality filtering: remove poor quality sequences
#   3. Build error models: learn error rates from sequence quality scores
#   4. Denoise: remove artefacts from PCR and sequencing
#   5. Merging: combining the forward and reverse read. By default DADA2 merges 
#      paired reads using a value for the maximum number of mismatches. However, 
#      the ITS2 gene has a variable number of bases and therefore filtering out 
#      reads that have more than a set amount of mismatches will bias against
#      longer ITS2 variants (which are more likely to have more mismatches). Therefore,
#      we set the max num. of mismatches to a proportion of the length of
#      the overlapping region (1.5%).
#   6. Chimera removal: using the pooled method
#   7. Removal flanking regions: using ITSx
#   8. Taxonomy assignment: using the RDP classifier 
#

# Initialization
rm(list=ls())
pkgs <- c("dada2", "ShortRead", "Biostrings", "dplyr")
lapply(pkgs, library, character.only=T, quiet=T, warn.conflicts=F)

# Provide path to folder where fastq sequence files are located 
path <- "replace this with path to folder containing sequences"

# Provide path to folder where cutadapt is installed
cutadapt <- "replace this with path to folder containing cutadapt"

#list sequences files present in your specified folder
list.files(path)

# Create character vector of locations + names of forward seqs
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
# Same for reverse
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE)) 

# NC1 and NC2 primers
FWD <- "ACGTCTGGTTCAGGGTTGTT"
REV <- "TTAGTTTCTTTTCCTCCGCT"

# Function that returns all orientations (forward, reverse, their complements
#   and reverse complements) of primers
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna),
               Reverse = Biostrings::reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}

# Get all of our orientation
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Specify location + name of reads filtered for ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# Filter out all sequences containing Ns
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN=0, compress = FALSE,
              multithread = FALSE)

#Function to count the number of primer hits
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Count hits
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

system(paste(cutadapt, "--version")) #  ensure that cutadapt is working

# This will be the directory containing reads with primers remove
path.cut <- file.path(path, paste0("cutadapt"))
# Make the directory
if(!dir.exists(path.cut)) dir.create(path.cut)
# Specify directory + names of reads with primers removed
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Specify cutadapt parameters:
#   - g stands for front of the read,
#   - a stands for the adaptor (primer) to remove
#   - Capitalized parameters correspond to the R2 read
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run cutadapt from r: n is the number of rounds to search for and 
#   remove a primer, set to 5 in case primers appear more than once.
#   o refers to output name
#   p indicates that outputs consist of paired-end reads
for (i in seq_along(fnFs)) {
  system(paste(cutadapt, R1.flags, R2.flags, "-n", 5, "-o", fnFs.cut[i], 
               "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}

# Now count primer hits in our sequences, should all be 0
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Ensure the pipeline has the same results every time
set.seed(1)

# Input the path of sequences ("path" variable), 
# DADA2 formatted taxonomy reference FASTA file ("db_path"),
# allowed overlap mismatch proportion ("overlap_mismatch_proportion"),
# and confidence threshold 
path <- "seqs/cutadapt"
db_path <- "Equine_Parasitic_Nematodes_Reference_Taxonomy_Database_v1.fasta"
overlap_mismatch_proportion <- 0.015
taxonomy_confidence_threshold <- 0

list.files(path)

# Extracting directory + filenames of primerless sequences
cutFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names by removing absolute path from file name
#   and cutting off sample name using Illumina identifier
#   (eg: "_L001_R1_001.fastq" becomes "")
get.sample.name <- function(fname) strsplit(basename(fname), "_L")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

# Plot forward and reverse quality profiles (get an idea of how our reads look)
jpeg('cutFs_quality_profile.jpg') ; plotQualityProfile(cutFs[1:2]) ; dev.off()
# Quality of reverse reads crashes quickly
jpeg('cutRs_quality_profile.jpg') ; plotQualityProfile(cutRs[1:2]) ; dev.off() 

# Create a vector of file paths + names to where the filtered sequences
#   will be stored
filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))
saveRDS(filtFs, "filtFs.rds") ; saveRDS(filtRs, "filtRs.rds")
filtFs <- readRDS("filtFs.rds") ; filtRs <- readRDS("filtRs.rds")

# 1. Conduct quality filtering
#   - maxN: Set maximum number of ambiguous bases to 0 
#     (should already be 0 from primer removal step)
#   - maxEE: vector containing the maximum expected errors for the
#     forward and reverse read of the sample
#   - rm.phix: check for PhiX sequences and remove if found
#   - compress: store filtered sequences in compressed format to save space
#   - truncQ: truncate the remainder of a sequence at the first occurrence
#     of quality score 2 or lower
#   - minLen: reads should be at least 200 base pairs
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN=0, maxEE=c(2,5),
		     rm.phix=T, compress=T, multithread=F, truncQ=2, minLen=200)
saveRDS(out, "out.rds")


# 2. Construct error models (mostly defaults)
# MAX_CONSIST: the number of rounds in self-consistency loop
#   to reach error convergence, increased from 10 to 15
errF <- learnErrors(filtFs, multithread = FALSE, MAX_CONSIST=15)
errR <- learnErrors(filtRs, multithread = FALSE, MAX_CONSIST=15)
saveRDS(errF, "errF.rds") ; saveRDS(errR, "errR.rds")

# Plot error models
# Estimated error rates (black line) should fit with observed error
#   rates (points), and error rates should decrease with higher quality score
jpeg("error-f.jpg") ; plotErrors(errF, nominalQ = TRUE) ; dev.off()
jpeg("error-r.jpg") ; plotErrors(errR, nominalQ = TRUE) ; dev.off()

# 3. Denoising
#   Enable full alignment
errF <- readRDS("errF.rds")
errR <- readRDS("errR.rds")
setDadaOpt(BAND_SIZE=-1)

dadaFsPooled <- dada(filtFs, err = errF, multithread = FALSE, pool = T)
dadaRsPooled <- dada(filtRs, err = errR, multithread = FALSE, pool = T)
saveRDS(dadaFsPooled, "dadaFsPooled.rds") ; saveRDS(dadaRsPooled, "dadaRsPooled.rds")

# 4. Merging
#   Initial merge to compute overlap lengths (set maxMismatch to infinity)
mergers_w_rejects <- mergePairs(dadaFsPooled, filtFs, dadaRsPooled, filtRs, 
                                verbose=TRUE, maxMismatch = Inf)

# Create empty vector to contain merged reads by sample
mergers <- vector("list", length(mergers_w_rejects)) 
names(mergers) <- sample.names #  Name samples

for (i in c(1:length(mergers_w_rejects))) {
  mergers_w_rejects[[i]]$length <- mergers_w_rejects[[i]]$nmatch + mergers_w_rejects[[i]]$nmismatch
  
  mergers_w_rejects[[i]]$percent_mismatch <-
    (mergers_w_rejects[[i]]$nmismatch + mergers_w_rejects[[i]]$nindel)/mergers_w_rejects[[i]]$length*100
  
  # First remove all ASVs below the percent cutoff
  mergers_w_rejects[[i]] <- mergers_w_rejects[[i]][mergers_w_rejects[[i]]$percent_mismatch < percent, ]
  
  # And also only keep sequences with accept = T
  mergers_w_rejects[[i]] <- mergers_w_rejects[[i]][mergers_w_rejects[[i]]$accept == TRUE, ]
  
  
  mergers[[i]] <- mergers_w_rejects[[i]]
  
}
saveRDS(mergers, "mergers.rds")
rm(mergers_w_rejects)
seqtab <- makeSequenceTable(mergers) ; dim(seqtab)

# 5. Remove chimeras using pooled method
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", 
                                    multithread=FALSE, verbose=TRUE)
# Create a table of ASV length distributions
table(nchar(getSequences(seqtab.nochim)))

# Create a table showing how many sequences were retained after 
# each step of the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFsPooled, getN), sapply(dadaRsPooled, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", 
                     "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

#save this file
write.table(track, "trackSequences.csv", quote=F, sep=",")

# Creating fasta file, + merging count and taxonomy tables
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# Loop to create headers for ASV fasta file (>ASV1, >ASV2, ...)
for (i in 1:dim(seqtab.nochim)[2]) 
{
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#	Making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")
pre_extract <- readDNAStringSet("ASVs.fa")

# 6. Remove flanking regions using ITSx, running on ASVs.fa
system("ITSx -i ASVs.fa --cpu 2 --partial 200 -t M --complement F --preserve T")

# Read extracted ITSx sequences into seqtab.nochim
extracted_seqs <- readDNAStringSet("ITSx_out.ITS2.full_and_partial.fasta")

indx <- which(pre_extract@ranges@NAMES %in% extracted_seqs@ranges@NAMES)
seqtab.nochim <- seqtab.nochim[ ,indx]

colnames(seqtab.nochim) <- as.character(extracted_seqs)
asv_seqs_ex <- colnames(seqtab.nochim)

#	Making and writing out a fasta of our final ASV seqs:
asv_headers <- asv_headers[indx]
asv_fasta_ex <- c(rbind(asv_headers, asv_seqs_ex))
write(asv_fasta_ex, "ASVs_extracted.fa")

# 7. Taxonomy using assignTaxonomy (RDP algorithm) with pre-specified database
#   and confidence threshold from above
assignTax <- assignTaxonomy(seqtab.nochim, db_path, multithread=TRUE,
                            minBoot = taxonomy_confidence_threshold, 
                            tryRC=TRUE, outputBootstraps=TRUE)
saveRDS(assignTax, "assignTax.rds")

#	Count table: number of times ASV is observed/sample
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)

# Taxonomy tables: ASV # and corresponding taxonomy
asv_tax <- assignTax$tax
asv_tax_boot <- assignTax$boot
colnames(asv_tax) <- c("domain", "kingdom", "phylum", "class", 
                       "order", "family", "genus", "species")
colnames(asv_tax_boot) <- c("domain.conf", "kingdom.conf", "phylum.conf", 
                            "class.conf", "order.conf", "family.conf", 
                            "genus.conf", "species.conf")
rownames(asv_tax) <- sub(">", "", asv_headers)
rownames(asv_tax_boot) <- sub(">", "", asv_headers)

# Merging taxonomy assignments and bootstrap values
tax_merged <- cbind(asv_tax, asv_tax_boot)

# Rearrange for viewing convenience
tax_merged_rearranged <-data.frame(tax_merged[,c(1,9,2,10,3,11,4,12,
                                                 5,13,6,14,7,15,8,16)])

# Merging ASV counts and new taxonomy table
final_table <- data.frame(cbind(asv_tab, tax_merged_rearranged), check.names=F)

head(final_table)

write.csv(final_table, "ASV_and_taxonomic_assignments_table.csv", 
            quote=F, row.names=T)

# The resulting file contains all ASVs identified by dada2, as well as corresponding taxonomic
# assignments and their confidence. Results from this file can be further processed and analyzed 
# according to specific research goals. 
