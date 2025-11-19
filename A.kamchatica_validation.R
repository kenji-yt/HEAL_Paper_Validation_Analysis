#################
# Load Packages #
#################

library(devtools)
library(data.table)
library(gridExtra)
library(healr)

####################
# Define Functions #
####################

# A function to count entries that are not null in a heal list with a "diff" item (sample wise differences between two lists).
count_non_zero <- function(diff_list){
  
  wa <- lapply(diff_list, function(prog){
    just_diff <- prog$diff[,5:ncol(prog$diff)]
    non_zero <- just_diff[, lapply(.SD, function(x) sum(x != 0, na.rm = T))]
    perc <- non_zero/nrow(prog$diff)
    return(perc)
  })
  
  we <- rbindlist(wa)
  
  sums <- we[, lapply(.SD, sum)]
  return(as.vector(unlist(sums)))
}

# Function to get only the differences over both subgenomes.
get_diff_dt <- function(diff_list){
  
  wa <- lapply(diff_list, function(prog){
    just_diff <- prog$diff[,5:ncol(prog$diff)]
    return(just_diff)
  })
  
  we <- rbindlist(wa)
  
  return(we)
}


#############################################
# Load data, filter and check GC correction #
#############################################

# Load in the heal lists (precounted by HEAL)
wgbs_list <- read_heal_list("/srv/kenlab/kenji/Rs7_wgbs_healr/")
dna_list <- read_heal_list("/export/kenlabsmb/kenji/DNAseq_Stefan_PAPER_reanalysis/back_up_out/DNA_seq_reanalysis/HEAL/healr/healr_list")
subset_50M_dna <- read_heal_list("/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/50M_run/results/healr/healr_list")
subset_20M_dna <- read_heal_list("/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/HEAL/results/healr/healr_list")

# filter
dna_filt <- filter_bins(dna_list)
wgbs_filt <- filter_bins(wgbs_list)
sub_50_filt <- filter_bins(subset_50M_dna)
sub_20_filt <- filter_bins(subset_20M_dna)

# GC correction and plot generation (commented out here )
dn_gc_cor <- correct_gc(dna_filt, n_windows = 10, loess_span = 1, pch = 16, alpha = 0.1, cex = 0.4, ymax = 5 )#, device = "png", output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomic_paper/")
sub_50_gc_cor <- correct_gc(sub_50_filt, n_windows = 10, loess_span = 1, pch = 16, alpha = 0.1, cex = 0.4, ymax = 5 )#, device = "png", output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomic_paper/")
sub_20_gc_cor <- correct_gc(sub_20_filt, n_windows = 10, loess_span = 1, pch = 16, alpha = 0.1, cex = 0.4, ymax = 5 )#, device = "png", output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomic_paper/")

## Correction seems not necessarily desirable.

# Get copy number on the filtered by not GC corrected data
dna_cn <- get_copy_number(dna_filt)
wgbs_cn <- get_copy_number(wgbs_filt)
sub_50_cn <- get_copy_number(sub_50_filt)
sub_20_cn <- get_copy_number(sub_20_filt)

#####################################################
#  Count what % of bins have a CN different from 2  #
#####################################################
 
dna_cn_summary<- summarize_cn(dna_cn)
wgbs_cn_summary <- summarize_cn(wgbs_cn)

norm_perc_vec_dna <- unlist(lapply(dna_cn_summary, function(smp){
  which_norm <- colnames(smp$total_count_table) == "2"
  smp$total_count_table[2, which_norm]
}))
summary(100-norm_perc_vec_dna)

norm_perc_vec_wgbs <- unlist(lapply(wgbs_cn_summary, function(smp){
  which_norm <- colnames(smp$total_count_table) == "2"
  smp$total_count_table[2, which_norm]
}))
summary(100-norm_perc_vec_wgbs)


###########################################
# Error in sample file naming: Correcting #
###########################################
### We found that the sample names for the WGBS and WGS data were not matching.
### This is visible when making plots of CN for samples from WGS vs WGBS:

#### HM_G1 ####
dna_plot <- plot_all_bins(dna_cn, view_samples = c("HM_RS7_G1_1","HM_RS7_G1_2","HM_RS7_G1_3"), plot_cn = T, return_plot = T)
wgbs_plot <- plot_all_bins(wgbs_cn, view_samples = c("HM_RS7_G1_1","HM_RS7_G1_2","HM_RS7_G1_3"), plot_cn = T, return_plot = T)

grid.arrange(dna_plot, wgbs_plot, nrow = 2)
### We can see that the sample IDs have been mixed up. We use this plot to rename the WGS datasets.


# HM_RS7_G1_1 => HM_RS7_G1_3
names(dna_cn$A.halleri$bins)[6] <- "HM_RS7_G1_3"
names(dna_cn$A.halleri$CN)[5] <- "HM_RS7_G1_3"
names(dna_cn$A.lyrata$bins)[6] <- "HM_RS7_G1_3"
names(dna_cn$A.lyrata$CN)[5] <- "HM_RS7_G1_3"
# HM_RS7_G1_2 => HM_RS7_G1_1
names(dna_cn$A.halleri$bins)[7] <- "HM_RS7_G1_1"
names(dna_cn$A.halleri$CN)[6] <- "HM_RS7_G1_1"
names(dna_cn$A.lyrata$bins)[7] <- "HM_RS7_G1_1"
names(dna_cn$A.lyrata$CN)[6] <- "HM_RS7_G1_1"
# HM_RS7_G1_3 => HM_RS7_G1_2
names(dna_cn$A.halleri$bins)[8] <- "HM_RS7_G1_2"
names(dna_cn$A.halleri$CN)[7] <- "HM_RS7_G1_2"
names(dna_cn$A.lyrata$bins)[8] <- "HM_RS7_G1_2"
names(dna_cn$A.lyrata$CN)[7] <- "HM_RS7_G1_2"

# 50M
# HM_RS7_G1_1 => HM_RS7_G1_3
names(sub_50_cn$A.halleri$bins)[6] <- "HM_RS7_G1_3"
names(sub_50_cn$A.halleri$CN)[5] <- "HM_RS7_G1_3"
names(sub_50_cn$A.lyrata$bins)[6] <- "HM_RS7_G1_3"
names(sub_50_cn$A.lyrata$CN)[5] <- "HM_RS7_G1_3"
# HM_RS7_G1_2 => HM_RS7_G1_1
names(sub_50_cn$A.halleri$bins)[7] <- "HM_RS7_G1_1"
names(sub_50_cn$A.halleri$CN)[6] <- "HM_RS7_G1_1"
names(sub_50_cn$A.lyrata$bins)[7] <- "HM_RS7_G1_1"
names(sub_50_cn$A.lyrata$CN)[6] <- "HM_RS7_G1_1"
# HM_RS7_G1_3 => HM_RS7_G1_2
names(sub_50_cn$A.halleri$bins)[8] <- "HM_RS7_G1_2"
names(sub_50_cn$A.halleri$CN)[7] <- "HM_RS7_G1_2"
names(sub_50_cn$A.lyrata$bins)[8] <- "HM_RS7_G1_2"
names(sub_50_cn$A.lyrata$CN)[7] <- "HM_RS7_G1_2"

# 20M
# HM_RS7_G1_1 => HM_RS7_G1_3
names(sub_20_cn$A.halleri$bins)[6] <- "HM_RS7_G1_3"
names(sub_20_cn$A.halleri$CN)[5] <- "HM_RS7_G1_3"
names(sub_20_cn$A.lyrata$bins)[6] <- "HM_RS7_G1_3"
names(sub_20_cn$A.lyrata$CN)[5] <- "HM_RS7_G1_3"
# HM_RS7_G1_2 => HM_RS7_G1_1
names(sub_20_cn$A.halleri$bins)[7] <- "HM_RS7_G1_1"
names(sub_20_cn$A.halleri$CN)[6] <- "HM_RS7_G1_1"
names(sub_20_cn$A.lyrata$bins)[7] <- "HM_RS7_G1_1"
names(sub_20_cn$A.lyrata$CN)[6] <- "HM_RS7_G1_1"
# HM_RS7_G1_3 => HM_RS7_G1_2
names(sub_20_cn$A.halleri$bins)[8] <- "HM_RS7_G1_2"
names(sub_20_cn$A.halleri$CN)[7] <- "HM_RS7_G1_2"
names(sub_20_cn$A.lyrata$bins)[8] <- "HM_RS7_G1_2"
names(sub_20_cn$A.lyrata$CN)[7] <- "HM_RS7_G1_2"

### Check again!
dna_plot <- plot_all_bins(dna_cn, view_samples = c("HM_RS7_G1_1","HM_RS7_G1_2","HM_RS7_G1_3"), plot_cn = T, return_plot = T)
wgbs_plot <- plot_all_bins(wgbs_cn, view_samples = c("HM_RS7_G1_1","HM_RS7_G1_2","HM_RS7_G1_3"), plot_cn = T, return_plot = T)

grid.arrange(dna_plot, wgbs_plot, nrow = 2)
### Now it's ok!

#### HM_G4 ####

dna_plot <- plot_all_bins(dna_cn, view_samples = c("HM_RS7_G4_1","HM_RS7_G4_2","HM_RS7_G4_3"), plot_cn = T, return_plot = T)
wgbs_plot <- plot_all_bins(wgbs_cn, view_samples = c("HM_RS7_G4_1","HM_RS7_G4_2","HM_RS7_G4_3"), plot_cn = T, return_plot = T)

 grid.arrange(dna_plot, wgbs_plot, nrow = 2)
### Looks alright!

#### LL_G1 ####
dna_plot <- plot_all_bins(dna_cn, view_samples = c("LL_RS7_G1_1","LL_RS7_G1_2","LL_RS7_G1_3"), plot_cn = T, return_plot = T)
wgbs_plot <- plot_all_bins(wgbs_cn, view_samples = c("LL_RS7_G1_1","LL_RS7_G1_2","LL_RS7_G1_3"), plot_cn = T, return_plot = T)

grid.arrange(dna_plot, wgbs_plot, nrow = 2)
### Looks alright! 

#### LL_G4 ####
dna_plot <- plot_all_bins(dna_cn, view_samples = c("LL_RS7_G4_1","LL_RS7_G4_2","LL_RS7_G4_3"), plot_cn = T, return_plot = T)
wgbs_plot <- plot_all_bins(wgbs_cn, view_samples = c("LL_RS7_G4_1","LL_RS7_G4_2","LL_RS7_G4_3"), plot_cn = T, return_plot = T)

grid.arrange(dna_plot, wgbs_plot, nrow = 2)
### Looks alright! 


######################
# Remove short spans #
######################
wgbs_ns_cn <- remove_short_spans(wgbs_cn, max_length = 5)
dna_ns_cn <- remove_short_spans(dna_cn, max_length = 5)
sub_50_ns_cn <- remove_short_spans(sub_50_cn, max_length = 5)
sub_20_ns_cn <- remove_short_spans(sub_20_cn, max_length = 5)

###############################
# Read classification metrics #
###############################
### The input files used here were create with the 'qc_eagle_rc.sh' script. 
### Column V2 is the number of input reads and V3 is the percentage of reads kept. 

### wgs 
wgs_rc <- read.table("/export/kenlabsmb/kenji/DNAseq_Stefan_PAPER_reanalysis/back_up_out/DNA_seq_reanalysis/HEAL/results/eagle_rc/summary_rc_metrics.tsv", header = F, skip = 1)
mean(wgs_rc$V2/1000000)
min(wgs_rc$V2/1000000)
max(wgs_rc$V2/1000000)

mean(as.numeric(gsub("%", "", wgs_rc$V3)))
min(as.numeric(gsub("%", "", wgs_rc$V3)))
max(as.numeric(gsub("%", "", wgs_rc$V3)))

### wgbs
wgbs_rc <- read.table("/srv/kenlab/kenji/wgbs_rc_stats.tsv", header = F, skip = 1)
mean(wgbs_rc$V2/1000000)
min(wgbs_rc$V2/1000000)
max(wgbs_rc$V2/1000000)

mean(as.numeric(gsub("%", "", wgbs_rc$V3)))
min(as.numeric(gsub("%", "", wgbs_rc$V3)))
max(as.numeric(gsub("%", "", wgbs_rc$V3)))


######################################################
# check relationship between coverage and difference #
######################################################

#### Histogram of the medians
median_counts_dna <- unlist(lapply(get_sample_stats(dna_cn), function(list){
  list[[1]]
}))
median_counts_50Mdna <- unlist(lapply(get_sample_stats(sub_50_cn), function(list){
  list[[1]]
}))
median_counts_20Mdna <- unlist(lapply(get_sample_stats(sub_20_cn), function(list){
  list[[1]]
}))
median_counts_wgbs <- unlist(lapply(get_sample_stats(wgbs_cn), function(list){
  list[[1]]
}))

par(mfrow=c(4,1))
hist(median_counts_dna/10000*150, xlim = c(000, 50.000), main = "WGS", col = "red4", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)
hist(median_counts_50Mdna/10000*150, xlim = c(000, 50.000), main = "50M", col = "orangered", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)
hist(median_counts_20Mdna/10000*150, xlim = c(000, 50.000), main = "20M", col = "goldenrod1", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)
hist(median_counts_wgbs/10000*150, xlim = c(000, 50.000), main = "WGBS", col = "purple4", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)

#### Percentage of difference to WGS per sample
diff_50M_list <- compare_heal_list(dna_cn, sub_50_cn, n_threads = 12, abs=TRUE)
diff_perc_50M <- count_non_zero(diff_50M_list)

diff_20M_list <- compare_heal_list(dna_cn, sub_20_cn, n_threads = 12, abs=TRUE)
diff_perc_20M <- count_non_zero(diff_20M_list)

diff_wgbs_list <- compare_heal_list(dna_cn, wgbs_cn, n_threads = 12, abs=TRUE)
diff_perc_wgbs <- count_non_zero(diff_wgbs_list)

par(mfrow=c(3,1))
hist(diff_perc_50M*100, xlim=c(0, 0.2*100), col="orangered", main = "50M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_20M*100, xlim=c(0, 0.2*100), col="goldenrod1", main = "20M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_wgbs*100, col="purple4", main = "WGBS", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)

par(mfrow=c(2,1))
hist(diff_perc_50M, xlim=c(0, 0.03))
hist(diff_perc_20M, xlim=c(0, 0.03))


diff_gc_20M_list <- compare_heal_list(gc_dna_cn, gc_sub_20_cn, n_threads = 12, abs=TRUE)
diff_gc_perc_20M <- count_non_zero(diff_gc_20M_list)
diff_gc_50M_list <- compare_heal_list(gc_dna_cn, gc_sub_50_cn, n_threads = 12, abs=TRUE)
diff_gc_perc_50M <- count_non_zero(diff_gc_50M_list)

par(mfrow=c(2,1))
hist(diff_gc_perc_20M, xlim=c(0, 0.03))
hist(diff_perc_50M, xlim=c(0, 0.03))

# Do it with no short spans 
# sub_50_ns_cn <- remove_short_spans(sub_50_cn, max_length = 5)
# sub_20_ns_cn <- remove_short_spans(sub_20_cn, max_length = 5)
# 
# diff_ns_50M_list <- compare_heal_list(dna_ns_cn, sub_50_ns_cn, n_threads = 12, abs=TRUE)
# diff_perc_ns_50M <- count_non_zero(diff_ns_50M_list)
# 
# diff_ns_20M_list <- compare_heal_list(dna_ns_cn, sub_20_ns_cn, n_threads = 12, abs=TRUE)
# diff_perc_ns_20M <- count_non_zero(diff_ns_20M_list)
# 
# par(mfrow=c(4,1))
# hist(diff_perc_ns_50M, xlim=c(0, 0.03))
# hist(diff_perc_50M, xlim=c(0, 0.03))
# hist(diff_perc_ns_20M, xlim=c(0, 0.03))
# hist(diff_perc_20M, xlim=c(0, 0.03))


##################
#    WGBS data   #
##################

median_counts_wgbs <- unlist(lapply(get_sample_stats(wgbs_cn), function(list){
  list[[1]]
}))

par(mfrow=c(4,1))
hist(median_counts_dna/10000*150, xlim = c(000, 50.000), main = "WGS", col = "red4", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)
hist(median_counts_50Mdna/10000*150, xlim = c(000, 50.000), main = "50M", col = "orangered", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)
hist(median_counts_20Mdna/10000*150, xlim = c(000, 50.000), main = "20M", col = "goldenrod1", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)
hist(median_counts_wgbs/10000*150, xlim = c(000, 50.000), main = "WGBS", col = "purple4", xlab = "Coverage (per basepair)", cex.main = 2, cex.lab = 1.2)


diff_wgbs_list <- compare_heal_list(dna_cn, wgbs_cn, n_threads = 12, abs=TRUE)
# diff_wgbs_counts <- unlist(lapply(get_sample_stats(diff_wgbs_list, method = "mean"), function(list){
#   list[[1]]
# }))
diff_perc_wgbs <- count_non_zero(diff_wgbs_list)

# 
# wgbs_no_short <- remove_short_spans(wgbs_cn, max_length = 3)
# diff_ns_wgbs_list <- compare_heal_list(dna_cn, wgbs_no_short, n_threads = 12, abs=TRUE)
# diff_ns_wgbs_counts <- unlist(lapply(get_sample_stats(diff_ns_wgbs_list, method = "mean"), function(list){
#   list[[1]]
# }))
# diff_perc_ns_wgbs <- count_non_zero(diff_ns_wgbs_list)


par(mfrow=c(3,1))
hist(diff_perc_50M*100, xlim=c(0, 0.2*100), col="orangered", main = "50M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_20M*100, xlim=c(0, 0.2*100), col="goldenrod1", main = "20M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_wgbs*100, col="purple4", main = "WGBS", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
#hist(diff_perc_ns_wgbs)#, xlim=c(0,0.35))
# 
# diff_ns_wgbs_list <- compare_heal_list(dna_ns_cn, wgbs_ns_cn, n_threads = 12, abs=TRUE)
# diff_perc_ns_wgbs <- count_non_zero(diff_ns_wgbs_list)
# 
# par(mfrow=c(4,1))
# hist(diff_perc_50M, xlim=c(0, 0.2), col="blue4")
# hist(diff_perc_20M, xlim=c(0, 0.2), col="blue4")
# hist(diff_perc_ns_wgbs, col="blue4")
# hist(diff_perc_wgbs, col="blue4")
### No big difference!


diff_wgbs_list <- compare_heal_list(gc_dna_cn, wgbs_cn, n_threads = 12, abs=TRUE)
# diff_wgbs_counts <- unlist(lapply(get_sample_stats(diff_wgbs_list, method = "mean"), function(list){
#   list[[1]]
# }))
diff_gc_perc_wgbs <- count_non_zero(diff_wgbs_list)

par(mfrow=c(4,1))
hist(diff_perc_50M*100, xlim=c(0, 0.2*100), col="orangered", main = "50M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_20M*100, xlim=c(0, 0.2*100), col="goldenrod1", main = "20M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_wgbs*100, col="purple4", main = "WGBS", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_gc_perc_wgbs*100, col="purple", main = "WGBS", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)

# No difference.
######################################
# Check if the difference is biased! #
######################################
diff_wgbs_no_abs_list <- compare_heal_list(dna_cn, wgbs_cn, n_threads = 12)
count_non_zero(diff_wgbs_no_abs_list)
diffs_dt <- get_diff_dt(diff_wgbs_no_abs_list)

par(mfrow = c(4, 3))
for(i in 1:ncol(diffs_dt)){
  smp_nme <- colnames(diffs_dt)[i]
  n_poz <- sum(diffs_dt[[i]]>0, na.rm = T)
  n_neg <- sum(diffs_dt[[i]]<0, na.rm = T)
  diff=n_poz-n_neg
  hist(diffs_dt[[i]], xlim = c(-2, 2), col = "darkgreen", main = paste(n_neg,";", diff, ";", n_poz), xlab = smp_nme)
}

diff_no_abs_wgbs_counts <- unlist(lapply(get_sample_stats(diff_wgbs_no_abs_list, method = "mean"), function(list){
  list[[1]]
}))

diff_wgbs_no_abs_list$A.halleri$bins <- diff_wgbs_no_abs_list$A.halleri$diff
diff_wgbs_no_abs_list$A.lyrata$bins <- diff_wgbs_no_abs_list$A.lyrata$diff
diff_wgbs_no_abs_list$A.halleri$CN <- diff_wgbs_no_abs_list$A.halleri$diff
diff_wgbs_no_abs_list$A.lyrata$CN <- diff_wgbs_no_abs_list$A.lyrata$diff

plot_all_bins(diff_wgbs_no_abs_list, plot_cn = T, add_bins = F, cn_line_width = 0.9, color_map = c("purple3", "orange2"))

# Seems like there are a lot of short lived ones
waza_short <- remove_short_spans(diff_wgbs_no_abs_list, max_length = 5)
waza_short$A.halleri$diff <- waza_short$A.halleri$CN
waza_short$A.lyrata$diff <- waza_short$A.lyrata$CN

plot_all_bins(waza_short, plot_cn = T, add_bins = F, cn_line_width = 0.9, color_map = c("purple3", "orange2"))
##  visually a big effect but not much actually 
diffs_dt <- get_diff_dt(waza_short)

par(mfrow = c(4, 3))
for(i in 1:ncol(diffs_dt)){
  smp_nme <- colnames(diffs_dt)[i]
  n_poz <- sum(diffs_dt[[i]]>0, na.rm = T)
  n_neg <- sum(diffs_dt[[i]]<0, na.rm = T)
  diff=n_poz-n_neg
  hist(diffs_dt[[i]], xlim = c(-2, 2), col = "darkgreen", main = paste(n_neg,";", diff, ";", n_poz), xlab = smp_nme)
}
## NO BIG EFFECT


#######################################
#   Explore what affects accuracy     #
#######################################

# Coverage influences accuracy
par(mfrow=c(1,1))
plot(median_counts_wgbs/10000*150, diff_perc_ns_wgbs*100, pch = 16, ylim = c(0, 19), xlim = c(5.5, 18),
     xlab = "Coverage", ylab = "Percentage of Different Bins", main = "CN Call Differences Between WGBS and WGS vs Coverage in WGBS")
model <- lm(c(diff_perc_ns_wgbs*100) ~ c(median_counts_wgbs/10000*150))
abline(model, col = "red", lwd = 2)

#text(median_counts_wgbs/10000*150, diff_perc_ns_wgbs*100, labels = names(median_counts_wgbs), pos = 3)

# Visualize differences
diff_wgbs_list$A.halleri$bins <- diff_wgbs_list$A.halleri$diff
diff_wgbs_list$A.lyrata$bins <- diff_wgbs_list$A.lyrata$diff
diff_wgbs_list$A.halleri$CN <- diff_wgbs_list$A.halleri$diff
diff_wgbs_list$A.lyrata$CN <- diff_wgbs_list$A.lyrata$diff

plot_all_bins(diff_wgbs_list, plot_cn = T, add_bins = F, cn_line_width = 0.9, color_map = c("purple3", "orange2"))

for(smp in names(aln)){
  
  par(mfrow=c(1,2))
  boxplot(dna_cn$A.halleri$bins$gc_content ~ diff_wgbs_list$A.halleri$diff[[smp]]==0, col = c("red4", "green4"), xlab = "Same as WGS", ylab = "GC Content (per bin)", main = "GC content")
  boxplot(dna_cn$A.halleri$bins$mappability ~ diff_wgbs_list$A.halleri$diff[[smp]]==0, col = c("red4", "green4"), xlab = "Same as WGS", ylab = "Mappability (per bin)", main = "Mappability")
  wilcox.test(dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])
  
}

wilcox.test(dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])
#t.test(dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])

wilcox.test(dna_cn$A.halleri$bins$mappability[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$mappability[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])
#t.test(dna_cn$A.halleri$bins$mappability[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$mappability[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])

# See relationship between certain variables and differences


#############################
#    Investigate spread     #
#############################

# Get aln 
aln <- get_heal_alignment(dna_ns_cn, genespace_dir = "/srv/kenlab/kenji/re_ana_rs7/HEAL/results/genespace/", n_threads = 10)
wgbs_aln <- get_heal_alignment(wgbs_ns_cn, genespace_dir = "/srv/kenlab/kenji/re_ana_rs7/HEAL/results/genespace/", n_threads = 20)

# get summary:
dna_aln_summary <- summarize_aln(aln)
wgbs_aln_summary <- summarize_aln(wgbs_aln)

length_concordant <- unlist(lapply(dna_aln_summary, function(smp){
  prog <- names(smp)
  wa <- c()
  for(p in prog){
    dt <- smp[[p]][[paste0("total_", p)]]
    a <- dt[,3]=="0:4"
    b <- dt[,3]=="1:3"
    c <- dt[,3]=="2:2"
    d <- dt[,3]=="3:1"
    e <- dt[,3]=="4:0"
    which_concord <- rowSums(cbind(a,b,c,d,e))
    bp_length <- sum(dt[as.logical(which_concord),1])
    wa <- c(wa, bp_length)
  }
  return(sum(wa))
}))

length_concordant/(10000*(nrow(dna_cn$A.halleri$bins) + nrow(dna_cn$A.lyrata$bins)))
# HM_RS7_G4_1
plot_bins(dna_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgs", device = "svg") 
plot_bins(wgbs_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgbs", device = "svg") 
plot_bins(sub_20_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/20M", device = "svg")  

#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgbs", device = "svg")  
plot_alignment(alignment = aln, heal_list = dna_cn, add_bins = "all", view_sample = "HM_RS7_G4_1", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs", device = "svg") 
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all", view_sample = "HM_RS7_G4_1", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs", device = "svg") 

#plot_bins(sub_20_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/20M", device = "svg")  
#plot_bins(sub_50_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/50M", device = "svg") 

# LL_RS7_G4_1
#plot_bins(dna_cn, view_sample = "LL_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgs", device = "svg") 
#plot_bins(wgbs_cn, view_sample = "LL_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgbs", device = "svg")  
plot_alignment(alignment = aln, heal_list = dna_ns_cn, add_bins = "all", view_sample = "LL_RS7_G4_1", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs", device = "svg")  
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_ns_cn, add_bins = "all", view_sample = "LL_RS7_G4_1", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs", device = "svg") 


# LL_RS7_G4_2
#plot_bins(dna_cn, view_sample = "LL_RS7_G4_2", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7)#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgs") 
#plot_bins(wgbs_cn, view_sample = "LL_RS7_G4_2", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7)#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgbs") 
plot_alignment(alignment = aln, heal_list = dna_ns_cn, add_bins = "all", view_sample = "LL_RS7_G4_2", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs", device = "svg")  
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_ns_cn, add_bins = "all", view_sample = "LL_RS7_G4_2", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs", device = "svg")  

# LL_RS7_G4_3
plot_bins(dna_cn, view_sample = "LL_RS7_G4_3", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7)#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgs") 
plot_bins(wgbs_cn, view_sample = "LL_RS7_G4_3", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7)#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgbs") 
plot_alignment(alignment = aln, heal_list = dna_ns_cn, add_bins = "all", view_sample = "LL_RS7_G4_3", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs", device = "svg")  
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_ns_cn, add_bins = "all", view_sample = "LL_RS7_G4_3", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs", device = "svg")  



# HM_RS7_G4_3 
# Low coverage, high accuracy
plot_alignment(alignment = aln, heal_list = dna_cn, add_bins = "all", view_sample = "LL_RS7_G4_1", specific_chr = "chr5", color_map = c("purple3", "orange2"), alpha = 0.1)#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs") 
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all", view_sample = "LL_RS7_G4_1", specific_chr = "chr5", color_map = c("purple3", "orange2"), alpha = 0.1)#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs")

plot_bins(full_dna_cn, view_sample = "LL_RS7_G4_1", specific_chr = "chr5", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, add_DNAcopy = T)#, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgbs") 

# plot_bins(sub_50_cn, view_sample = "HM_RS7_G1_1", specific_chr = "chr4", plot_cn = T)
# plot_bins(sub_20_cn, view_sample = "HM_RS7_G1_1", specific_chr = "chr4", plot_cn = T)

# Get the variance in normalized read counts within copy number regions


#########################
#   check chr4 region   #
#########################
mappa <- dna_cn$A.halleri$bins$mappability[dna_cn$A.halleri$bins$chr=="chr4" & dna_cn$A.halleri$bins$start<18000000]
gc <- dna_cn$A.halleri$bins$gc_content[dna_cn$A.halleri$bins$chr=="chr4" & dna_cn$A.halleri$bins$start<18000000]


ex_axis <- diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0
ex_axis[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0] <- "Same as WGS"
ex_axis[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0] <- "Different from WGS"

all_cat <- factor(c(ex_axis, rep("chr4 : 0-18mbp", length(gc))), levels = c("Different from WGS", "Same as WGS", "chr4 : 0-18mbp"))

par(mfrow=c(1,2), mar = c(10, 6, 4, 3))
boxplot(c(dna_cn$A.halleri$bins$gc_content, gc) ~ all_cat, las = 2, col = c("red4", "green4", "purple4"), xlab="", ylab = "GC Content (per bin)", main = "GC content")
boxplot(c(dna_cn$A.halleri$bins$mappability, mappa) ~ all_cat, las = 2, col = c("red4", "green4", "purple4"), xlab = "", ylab = "Mappability (per bin)", main = "Mappability")


wilcox.test(dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0], gc)

wilcox.test(dna_cn$A.halleri$bins$mappability[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0], mappa)


###################
#    HEAT MAPS    #
###################
plot_heal_heat_map(alignment = aln, xrange=c(0:4), yrange=c(0:4), output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/heatmaps/wgs", device = "svg")
plot_heal_heat_map(alignment = wgbs_aln, , xrange=c(0:4), yrange=c(0:4), output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/heatmaps/wgbs", device = "svg")

plot_all_bins(dna_cn, plot_cn = TRUE, add_bins = F, cn_line_width = 0.9, color_map = c("purple3", "orange2"), output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/plot_all/wgs")
plot_all_bins(wgbs_cn, plot_cn = TRUE, add_bins = F, cn_line_width = 0.9, color_map = c("purple3", "orange2"), output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/plot_all/wgbs")

# 
# 
# 
# 
# 
# dna_plot <- plot_all_bins(dna_cn, view_samples = c("HM_RS7_G4_1","HM_RS7_G4_2","HM_RS7_G4_3"), plot_cn = T, return_plot = T)
# wgbs_plot <- plot_all_bins(wgbs_cn, view_samples = c("HM_RS7_G4_1","HM_RS7_G4_2","HM_RS7_G4_3"), plot_cn = T, return_plot = T)
# 
# grid.arrange(dna_plot, wgbs_plot, nrow = 2)
# 
# 
# plot(median_counts_wgbs, diff_counts)
# text(median_counts_wgbs, diff_counts, labels = names(median_counts_wgbs), pos = 3)
# 
# median_counts_dna <- unlist(lapply(get_sample_stats(dna_cn), function(list){
#   list[[1]]
# }))
# par(mfrow=c(2,1))
# hist(median_counts_dna, xlim = c(300, 2700))
# hist(median_counts_wgbs, xlim = c(300, 2700))
# ratio <- median_counts_dna[order(names(median_counts_dna))]/median_counts_wgbs[order(names(median_counts_wgbs))]
# plot(ratio, diff_counts)
# text(ratio, diff_counts, labels = names(median_counts_wgbs), pos = 3)
# 
# diff <- median_counts_dna[order(names(median_counts_dna))] - median_counts_wgbs[order(names(median_counts_wgbs))]
# plot(diff, diff_counts)
# text(ratio, diff_counts, labels = names(median_counts_wgbs), pos = 3)
# 
# 
# ###################
# #    DNA subset   #
# ###################
# subset_dna <- read_heal_list("/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/HEAL/results/healr/healr_list")
# 
# sub_filt <- filter_bins(subset_dna)
# 
# sub_cn <- get_copy_number(sub_filt)
# 
# plot_bins(sub_cn, view_sample = "LL_RS7_G1_1", specific_chr = "chr7", plot_cn = T)
# plot_bins(dna_cn, view_sample = "LL_RS7_G1_1", specific_chr = "chr7", plot_cn = T)
# plot_bins(wgbs_cn, view_sample = "LL_RS7_G1_1", specific_chr = "chr7", plot_cn = T)
# median_counts_sub <- unlist(lapply(get_sample_stats(sub_cn), function(list){
#   list[[1]]
# }))
# 
# par(mfrow=c(3,1))
# hist(median_counts_dna, xlim = c(300, 2700))
# hist(median_counts_sub, xlim = c(300, 2700))
# hist(median_counts_wgbs, xlim = c(300, 2700))
# 
# diff_sub_list <- compare_heal_list(dna_cn, sub_cn, n_threads = 12, abs=TRUE)
# diff_sub_counts <- unlist(lapply(get_sample_stats(diff_sub_list, method = "mean"), function(list){
#   list[[1]]
# }))
# 
# # For some reason there is a line above???
# plot_bins(diff_list, view_sample = "LL_RS7_G1_1", plot_cn = T, add_bins = T)
# 
# plot(median_counts_sub, diff_counts)
# text(median_counts_sub, diff_counts, labels = names(median_counts_sub), pos = 3)
# 
# plot_diff! 
#   Investigate the distribution of differences
# &
#   If GC or mappability affect differences.
# 
# aln <- get_heal_alignment(sub_cn, genespace_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/HEAL/results/genespace/")