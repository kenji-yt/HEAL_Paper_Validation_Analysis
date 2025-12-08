#################
# Load Packages #
#################

library(healr)
# 
# library(devtools)
# load_all()

####################################
# Load data and rename chromosomes #
####################################

final_list <- read_heal_list("/srv/kenlab/kenji/HE_re_analysis_published/Aled_suecica/HEAL/results/healr/healr_list/")

wa <- list(final_list$A.thaliana, final_list$A.arenosa)
names(wa) <- names(final_list)[c(2,1)]

entries <- 6:13
a_replace <- unique(wa$A.arenosa$bins$chr)
names(a_replace) <- entries
for(i in 1:length(entries)){
  wa$A.arenosa$bins$chr[wa$A.arenosa$bins$chr==a_replace[i]] <- rep(as.numeric(entries[i]),sum(wa$A.arenosa$bins$chr==a_replace[i]))
}
wa$A.arenosa$bins$chr <- as.numeric(wa$A.arenosa$bins$chr)
wa$A.arenosa$bins <- data.table::setkey(wa$A.arenosa$bins, chr, start, end, gc_content)

final_list <- wa


#########################################################
# Try a range of filtering schemes and show the results #
#########################################################

# No filter local 
filt_aled <- final_list
cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
count_thresh_1_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"), output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/no_filter", device = "svg")

# Default local
filt_aled <- filter_bins(final_list, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/default_local/filtering.log")
cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output = T)
default_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
# Some spikes (A.thal chr 1 in all.), finds chr1 all changes for both subgenomes (edge up down). 
# A big gain in chr5 of A.arenosa for four samples (centromer?).
# A loss on chromosome 6 A.arenosa for the other sample. 
# Does not label a gain for chr2 A.arenosa end, but chr2 lower count are in 4 samples of A.thal so it might not be linked...
# Let's plot 042 chr1 and 6 for the supplementaries:
plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"))#, output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/default_local", device = "svg")


# # Higher count threshold = 3 
# filt_aled <- filter_bins(final_list, count_threshold = 3)
# cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
# count_thresh_3_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
# Same but with some more high short spans (two more on chr3 of A.th)

# # Mappability 0.8
# filt_aled <- filter_bins(final_list, mappability_threshold = 0.8, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/mappa_0.8/filtering.log")
# cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
# map_08_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
# More short spans on A.are (all for chr3) but less short spans for A.thal (no more for chr3).
# BUT all a.thal chr5 have a short span in the centre (centromere) which is not a high count but just a +1


# # Mappability 0.95
# filt_aled <- filter_bins(final_list, mappability_threshold = 0.95)
# cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
# map_95_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
# Only two short spans on chr3 of A.th. 
# Weird behaviour on chr7 A.are 044 up down, also on chr6 of same sample does down. 
# The chr2 A.are is showing the gain at the end! 
# THe chr1 of A.thal is now not showing the loss! 
# Also part of chr7 of 051 is lower. 
# Some gain in two samples for A.are chr5 (seen already above).

##### Mappability 0.95, count 3
# filt_aled <- filter_bins(final_list, mappability_threshold = 0.95, count_threshold = 3)
# cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
# map_95_count_3_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
# Almost no short spans execpt 1 on chr 3 of A.thal.
# Same chr1 A.th and chr2 A.are as above.
# Same weird chr 7 A.are and loss on chr6 and loss on chr 7 as above.
# Same chr5 gain for 004 and 051 as above. 

##### Mappability 0.99
# filt_aled <- filter_bins(final_list, mappability_threshold = 0.99, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/99_mapa/filtering.log")
# cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
# map_99_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
# plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"), output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/99_mapa", device = "svg")


# Higher count threshold = 1
filt_aled <- filter_bins(final_list, count_threshold = 1, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/count_1/filtering.log")
cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
count_thresh_1_local <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local", return_plot = T)
# The short span gain of chr1 of A.thal is here. 
# Two other short span gains on chr3 A.thal. 
# The large scale gain on 4 chr5 A.are. 
# Loss on chr6 A.are
plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"))#, output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/count_1", device = "svg")


# Default non local
filt_aled <- filter_bins(final_list)
cn_aled <- get_copy_number(filt_aled, n_threads = 10, full_output =T)
default_global <- plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, return_plot = T)
# This is visibly wrong. Let's visualize it better. 
plot_bins(cn_aled, view_sample = "NAT_X_SYN_004", specific_chr = c("LR999451.1"), plot_cn = T, ylim_max = 5)
plot_bins(cn_aled, view_sample = "NAT_X_SYN_004", specific_chr = c("1"), plot_cn = F, ylim_max = 0.1)


######################################
# Get plot of local vs global median #
######################################
local_medians <- get_sample_stats(cn_aled, method = "local_median")
medians <- get_sample_stats(cn_aled)

count_are <- cn_aled$A.arenosa$bins$NAT_X_SYN_004[cn_aled$A.arenosa$bins$chr=="6"]
count_thal <- cn_aled$A.thaliana$bins$NAT_X_SYN_004[cn_aled$A.thaliana$bins$chr=="1"]

pos_are <- cn_aled$A.arenosa$bins$start[cn_aled$A.arenosa$bins$chr=="6"]
pos_thal <- cn_aled$A.thaliana$bins$start[cn_aled$A.thaliana$bins$chr=="1"]


par(mfrow=c(1,2))
plot(pos_thal, count_thal, ylim=c(0,6000), pch=16, main="A.thaliana Subgenome Chromosome 1", ylab="Counts", xlab = "Position (bp)")
abline(h=medians$NAT_X_SYN_004$A.arenosa, col='red', lwd=1.5)
abline(h=local_medians$NAT_X_SYN_004$A.thaliana, col='green', lwd=3)
plot(pos_are, count_are, ylim=c(0,6000), pch=16, main="A.arenosa Subgenome Chromosome 1", ylab="", xlab = "Position (bp)")
abline(h=medians$NAT_X_SYN_004$A.arenosa, col='red', lwd=1.5)
abline(h=local_medians$NAT_X_SYN_004$A.arenosa, col='green', lwd=3)

# Good one: filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)#, gc_quantile = 70)#, mappability_threshold = 0.95)
# Also good one:filt_aled <- filter_bins(final_list, count_threshold = 1, mappability_threshold = 0.98)#, gc_quantile = 70)#, mappability_threshold = 0.95)
filt_aled <- filter_bins(final_list, mappability_threshold = 0.98, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/mapa_98/filtering.log") # PERFECT!!! :D 

cn_aled <- get_copy_number(gc_corected , n_threads = 10, method = "local", full_output =T)
#cn_aled$A.thaliana$DNAcopy$NAT_X_SYN_042
plot_all_bins(cn_aled, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"), plot_cn = T, method = "local")
plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"), output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/mapa_98")


#### Getting the alignment #### # NON MO PA FER SA! 
final_list <- read_heal_list("/srv/kenlab/kenji/HE_re_analysis_published/Aled_suecica/HEAL/results/healr/healr_list/")

filt_aled <- filter_bins(final_list, mappability_threshold = 0.98, gc_threshold = 1.5)#, gc_quantile = 10)#, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/mapa_98/filtering.log") # PERFECT!!! :D 
gc_corrected <- correct_gc(filt_aled, local_normalize = T, n_windows = 10, loess_span = 1, pch = 16, alpha = 0.3, cex = 0.6, output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/gc_filtering_test", device = "svg")# output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/gc_filtering_test", )
cn_aled <- get_copy_number(gc_corrected, n_threads = 10, method = "local", full_output =T)
plot_all_bins(heal_list = cn_aled, plot_cn = T)
wa <- remove_short_spans(cn_aled, max_length = 3)
plot_all_bins(heal_list = wa, plot_cn = T, method = "local", sample_name_size = 4,  
              view_samples =c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"),
                color_map = c("red2", "blue2"),
              add_bins = T, bin_point_size = 1.3, bin_point_alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/", width = 10, height = 8, device = "svg")

aln <- get_heal_alignment(cn_aled, genespace_dir = "/srv/kenlab/kenji/HE_re_analysis_published/Aled_suecica/HEAL/results/genespace/")

dens <- get_concordant_density(aln, cn_aled, normalize = "local")
correct <- correct_cn_with_density(dens, aln, cn_aled, n_threads = 5)
plot_densities(dens, view_sample = "NAT_X_SYN_042", show_discordant = T, heal_list = cn_aled, alignment = aln, corrected_alignment = correct, normalize = "local")
plot_riparian(alignment = aln, heal_list = cn_aled, genespace_dir = "/srv/kenlab/kenji/HE_re_analysis_published/Aled_suecica/HEAL/results/genespace/")
 
# Rename!
wa <- list(cn_aled$A.thaliana, cn_aled$A.arenosa)
names(wa) <- names(cn_aled)[c(2,1)]

entries <- 6:13
a_replace <- unique(wa$A.arenosa$CN$chr)
names(a_replace) <- entries
for(i in 1:length(entries)){
  wa$A.arenosa$CN$chr[wa$A.arenosa$CN$chr==a_replace[i]] <- rep(entries[i],sum(wa$A.arenosa$CN$chr==a_replace[i]))
  wa$A.arenosa$bins$chr[wa$A.arenosa$bins$chr==a_replace[i]] <- rep(entries[i],sum(wa$A.arenosa$bins$chr==a_replace[i]))
  
  for(smp in names(aln)){
    aln[[smp]]$chr_A.arenosa[aln[[smp]]$chr_A.arenosa==a_replace[i]] <- rep(entries[i],sum(aln[[smp]]$chr_A.arenosa==a_replace[i]))
  }
}

plot_alignment(heal_list = wa, alignment = aln, view_sample = "NAT_X_SYN_042", specific_chr = 1, add_bins = "all", method = "local")
plot_alignment(heal_list = wa, alignment = aln, view_sample = "NAT_X_SYN_051", specific_chr = 12, add_bins = "all", method = "local")

plot_riparian(alignment = aln, heal_list = wa, genespace_dir = )

plot_all_bins(wa, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"),
              plot_cn = T, method = "local", sample_name_size = 4, color_map = c("red2", "blue2"),
              add_bins = T, bin_point_size = 1.3, bin_point_alpha = 0.1, output_dir = "/srv/kenlab/kenji/HE_re_analysis_published/Aled_suecica/", width = 10, height = 8)



plot_all_bins(wa, view_samples = c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"),
              plot_cn = T, method = "local", sample_name_size = 4, color_map = c("red2", "green3"),
              add_bins = T, bin_point_size = 1.3, bin_point_alpha = 0.1, output_dir = "/srv/kenlab/kenji/HE_re_analysis_published/Aled_suecica/poster_plot", width = 10, height = 8, cn_line_width = 1.5)

# 
# 
# plot_bins(cn_aled, view_sample = "NAT_X_SYN_054", specific_chr = "2", plot_cn = T, add_bins = T, method = "local", alpha = 0.8)
# 
# cn_sum <- summarize_cn(cn_aled, n_threads = 10)
# write_cn_summary(cn_sum, output_dir = "/srv/kenlab/kenji/summary")
# aln_aled <- get_heal_alignment(cn_aled, "/srv/kenlab/kenji/HE_re_analysis_published/Aled_suecica/HEAL/results/genespace/", n_threads = 10)
# plot_bins(cn_aled, view_sample = "NAT_X_SYN_004", plot_cn = T, add_bins = T, method = "local")
# plot_bins(cn_aled, view_sample = "NAT_X_SYN_054", plot_cn = T, add_bins = T, method = "local")
# plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, add_bins = T, method = "local", specific_chr = "1")
# plot_bins(cn_aled, output_dir = "/srv/kenlab/kenji/plotobitombo", plot_cn = T, add_bins = T, method = "local")
# plot_alignment(cn_aled, aln_aled, view_sample = "NAT_X_SYN_042", method = "local", add_bins = "both", color_map = c("orange", "purple"))
# 
# dens_aled <- get_concordant_density(aln_aled, cn_aled, n_threads = 5, normalize = "local")
# plot_densities(dens_aled, view_sample = "NAT_X_SYN_042", heal_list = cn_aled, normalize = "local")
# corec <- correct_cn_with_density(dens_aled, aln_aled, cn_aled)
# 
# plot_alignment(cn_aled, corec, view_sample = "NAT_X_SYN_042", method = "local", add_bins = "both", color_map = c("orange", "purple"))
