maxMismatchDiagnostic <- function(mergers) {
  # To make these plots I will need some extra packages not needed in the rest of the pipeline
  library(dplyr) ; library(ggplot2) ; library(ggpubr)
  
  # Generate table of all samples, first append sample name
  for (i in 1:length(mergers)) {mergers[[i]]$sample <- names(mergers)[[i]]}
  # use dplyr bind rows to create a singular data frame
  mergers.table <- bind_rows(mergers)
  
  # Add some summary statistics to this table, first total length of overlap region
  # equal to the number of mismatches + number of matches. Then also add
  # Length category which refers to overlap length (either 'short' <100bp or 'long' <200bp)
  mergers.table$length <- mergers.table$nmatch + mergers.table$nmismatch
  mergers.table$length_category <- ifelse(mergers.table$length < 100, "<100", "<200")
  # Compute total mismatch as the sum of mismatch + indels
  mergers.table$total_mismatch <- mergers.table$nmismatch + mergers.table$nindel
  # Compute percent mismatch as a percentage of length of overlap region
  mergers.table$percent_mismatch <- mergers.table$total_mismatch/mergers.table$length*100
  # For plotting purposes, assign a "percent_category"
  mergers.table <- mergers.table %>% mutate(percent_category = case_when(
    percent_mismatch == 0 ~ "0%",
    (percent_mismatch > 0 & percent_mismatch <= 0.5) ~ "< 0.5%",
    (percent_mismatch > 0.5 & percent_mismatch <= 1) ~ "< 1%",
    (percent_mismatch > 1 & percent_mismatch <= 1.5) ~ "< 1.5%",
    (percent_mismatch > 1.5 & percent_mismatch <= 2) ~ "< 2%",
    (percent_mismatch > 2 & percent_mismatch <= 2.5) ~ "< 2.5%",
    percent_mismatch > 2.5 ~ "> 2.5%"
  ))
  
  # First plot (panel a): total number of mismatches versus ASV abundance
  plot_a_mismatch_v_abundance <- ggplot(mergers.table, 
                                        aes(x = total_mismatch, y = abundance)) + 
    geom_point()
  
  mismatch.constant.agg <- aggregate(abundance ~ total_mismatch + length_category,
                                     data = mergers.table, FUN = sum)
  
  # Show count information up to 6 mismatches
  mismatch.constant.sum <- data.frame(
    total_mismatch = rep(c(0,1,2,3,4,5,6),2), 
    length_category = c(rep(c("<100"), 7), rep(c("<200"), 7)),
    count = c(cumsum(mismatch.constant.agg %>%
                       filter(length_category == "<100") %>%
                       top_n(7) %>% pull(abundance)),
              cumsum(mismatch.constant.agg %>%
                       filter(length_category == "<200") %>%
                       top_n(7) %>% pull(abundance)))
  )
  
  # Get total abundances of each kind of overlap category (short and long)
  sum.short.overlap <- mergers.table %>% filter(length_category == "<100") %>%
    pull(abundance) %>% sum()
  sum.long.overlap <- mergers.table %>% filter(length_category == "<200") %>%
    pull(abundance) %>% sum()
  
  # Label the length category 
  mismatch.constant.sum$total <- ifelse(mismatch.constant.sum$length_category == "<100", 
                                        sum.short.overlap, 
                                        sum.long.overlap)
  
  # Compute the proportin of reads excluded for each mismatch value + length 
  # category
  mismatch.constant.sum$proportion_excluded <- 
    1 - mismatch.constant.sum$count/mismatch.constant.sum$total
  
  # Second plot (panel b): proportion of reads excluded as a function of total mismatches
  plot_b_mismatch_v_prop_excluded <- ggplot(mismatch.constant.sum,
                                            aes(x = total_mismatch, y = proportion_excluded,
                                                color = length_category)) +
    geom_point()
  
  # Third plot (panel c): percent mismatch versus ASV abundance
  plot_c_mismatch_v_percent_mismatch <- ggplot(mergers.table,
                                               aes(x = percent_mismatch, y = abundance)) +
    geom_point()
  
  mismatch.percent.agg <- aggregate(abundance ~ percent_category + length_category,
                                    data = mergers.table, FUN = sum)
  
  # Order the above table by percent category
  # Make factor to specify order
  mismatch.percent.agg$percent_category <- factor(mismatch.percent.agg$percent_category,
                                                  levels = c("0%","<0.5%","<1%","<1.5%","<2%", "<2.5%", ">2.5%"))
  mismatch.percent.agg <- mismatch.percent.agg[order(mismatch.percent.agg$percent_category), ]
  
  mismatch.percent.sum <- data.frame(
    percent_mismatch = rep(c("0%","<0.5%","<1%","<1.5%","<2%", "<2.5%", ">2.5%"),2), 
    length_category = c(rep(c("<100"), 7), rep(c("<200"), 7)),
    count = c(cumsum(mismatch.percent.agg %>%
                       filter(length_category == "<100") %>%
                       slice_min(n = 7, order_by = percent_category, 
                                 with_ties= F) %>% pull(abundance)),
              cumsum(mismatch.percent.agg %>%
                       filter(length_category == "<200") %>%
                       slice_min(n = 7, order_by = percent_category, 
                                 with_ties = F) %>% pull(abundance)))
  )
  
  mismatch.percent.sum$total <- ifelse(mismatch.percent.sum$length_category == "<100",
                                       sum.short.overlap, sum.long.overlap)
  
  mismatch.percent.sum$proportion_excluded <- 
    1 - mismatch.percent.sum$count/mismatch.percent.sum$total
  
  # Make factor to specify order for plot
  mismatch.percent.sum$percent_mismatch <- factor(mismatch.percent.sum$percent_mismatch,
                                                  levels = c("0%","<0.5%","<1%","<1.5%","<2%", "<2.5%", ">2.5%"))
  
  # Fourth plot (panel d): proportion of reads excluded as a function of percent mismatch
  plot_d_percent_mismatch_v_prop_excluded <- 
    ggplot(mismatch.percent.sum, aes(x = percent_mismatch,
                                     y = proportion_excluded,
                                     color = length_category)) + 
    geom_point()
  
  ggpubr::ggarrange(plot_a_mismatch_v_abundance,
                    plot_b_mismatch_v_prop_excluded,
                    plot_c_mismatch_v_percent_mismatch,
                    plot_d_percent_mismatch_v_prop_excluded,
                    common.legend = TRUE, legend = "bottom")
  
}

maxMismatchDiagnostic(mergers_w_rejects)