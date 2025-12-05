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
#    Get alignment   #
######################

aln <- get_heal_alignment(dna_cn, genespace_dir = "/srv/kenlab/kenji/re_ana_rs7/HEAL/results/genespace/", n_threads = 10)
wgbs_aln <- get_heal_alignment(wgbs_cn, genespace_dir = "/srv/kenlab/kenji/re_ana_rs7/HEAL/results/genespace/", n_threads = 20)

######################
# Remove short spans #
######################
# wgbs_ns_cn <- remove_short_spans(wgbs_cn, max_length = 5)
# dna_ns_cn <- remove_short_spans(dna_cn, max_length = 5)
# sub_50_ns_cn <- remove_short_spans(sub_50_cn, max_length = 5)
# sub_20_ns_cn <- remove_short_spans(sub_20_cn, max_length = 5)

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

########## Figure 10 ###########
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

########## Figure 10 ###########
par(mfrow=c(3,1))
hist(diff_perc_50M*100, xlim=c(0, 0.2*100), col="orangered", main = "50M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_20M*100, xlim=c(0, 0.2*100), col="goldenrod1", main = "20M", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)
hist(diff_perc_wgbs*100, col="purple4", main = "WGBS", xlab = "Percentage of different bins", cex.main = 2.5, cex.lab = 1.3)


##############################
#   Coverage vs Accuracy     #
##############################

########## Figure 11 ###########
# Coverage influences accuracy
par(mfrow=c(1,1))
plot(median_counts_wgbs/10000*150, diff_perc_wgbs*100, pch = 16, ylim = c(0, 19), xlim = c(5.5, 18),
     xlab = "Coverage", ylab = "Percentage of Different Bins", main = "CN Call Differences Between WGBS and WGS vs Coverage in WGBS")
model <- lm(c(diff_perc_wgbs*100) ~ c(median_counts_wgbs/10000*150))
abline(model, col = "red", lwd = 2)


###########################################################################################
#    Plot counts and copy number for a specific sample and chromosome in all data sets    #
###########################################################################################

########## Figure 12 ###########
# HM_RS7_G4_1
plot_bins(dna_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgs", device = "svg") 
plot_bins(wgbs_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/wgbs", device = "svg") 
plot_bins(sub_20_cn, view_sample = "HM_RS7_G4_1", specific_chr = "chr4", plot_cn = T, color_map = c("purple3", "orange2"), alpha = 0.7, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/bins/20M", device = "svg")  


########## Supplementary figure 6 ###########
# GC content and mappability in matching, different and specific region bins for sample HM_RS7_G1_1
mappa <- dna_cn$A.halleri$bins$mappability[dna_cn$A.halleri$bins$chr=="chr4" & dna_cn$A.halleri$bins$start<18000000]
gc <- dna_cn$A.halleri$bins$gc_content[dna_cn$A.halleri$bins$chr=="chr4" & dna_cn$A.halleri$bins$start<18000000]


ex_axis <- diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0
ex_axis[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0] <- "Same as WGS"
ex_axis[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0] <- "Different from WGS"

all_cat <- factor(c(ex_axis, rep("chr4 : 0-18mbp", length(gc))), levels = c("Different from WGS", "Same as WGS", "chr4 : 0-18mbp"))

par(mfrow=c(1,2), mar = c(10, 6, 4, 3))
boxplot(c(dna_cn$A.halleri$bins$gc_content, gc) ~ all_cat, las = 2, col = c("red4", "green4", "purple4"), xlab="", ylab = "GC Content (per bin)", main = "GC content")
boxplot(c(dna_cn$A.halleri$bins$mappability, mappa) ~ all_cat, las = 2, col = c("red4", "green4", "purple4"), xlab = "", ylab = "Mappability (per bin)", main = "Mappability")

# Statistical tests
wilcox.test(dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])
wilcox.test(dna_cn$A.halleri$bins$mappability[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$mappability[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])

# Check the bias of GC and mappabilitty in same vs different bins for all samples.
for(smp in names(aln)){
  
  par(mfrow=c(1,2))
  boxplot(dna_cn$A.halleri$bins$gc_content ~ diff_wgbs_list$A.halleri$diff[[smp]]==0, col = c("red4", "green4"), xlab = "Same as WGS", ylab = "GC Content (per bin)", main = "GC content")
  boxplot(dna_cn$A.halleri$bins$mappability ~ diff_wgbs_list$A.halleri$diff[[smp]]==0, col = c("red4", "green4"), xlab = "Same as WGS", ylab = "Mappability (per bin)", main = "Mappability")
  wilcox.test(dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1!=0], dna_cn$A.halleri$bins$gc_content[diff_wgbs_list$A.halleri$diff$HM_RS7_G1_1==0])
  
}


##############################
#    Plot some alignments    #
##############################

########## Figure 13 ###########

# LL_RS7_G4_1
plot_alignment(alignment = aln, heal_list = dna_cn, add_bins = "all", view_sample = "LL_RS7_G4_1", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs", device = "svg")  
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all", view_sample = "LL_RS7_G4_1", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs", device = "svg") 

# LL_RS7_G4_2
plot_alignment(alignment = aln, heal_list = dna_cn, add_bins = "all", view_sample = "LL_RS7_G4_2", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs", device = "svg")  
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all", view_sample = "LL_RS7_G4_2", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs", device = "svg")  

# LL_RS7_G4_3
plot_alignment(alignment = aln, heal_list = dna_cn, add_bins = "all", view_sample = "LL_RS7_G4_3", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgs", device = "svg")  
plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all", view_sample = "LL_RS7_G4_3", specific_chr = "chr4", color_map = c("purple3", "orange2"), alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/aln/wgbs", device = "svg")  

###################
#    HEAT MAPS    #
###################

########## Figure 14 ###########

plot_heal_heat_map(alignment = aln, xrange=c(0:4), yrange=c(0:4), output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/heatmaps/wgs", device = "svg")
plot_heal_heat_map(alignment = wgbs_aln, , xrange=c(0:4), yrange=c(0:4), output_dir = "/srv/kenlab/kenji/exploring_rs7_BMC_genomics_paper/figures/heatmaps/wgbs", device = "svg")
