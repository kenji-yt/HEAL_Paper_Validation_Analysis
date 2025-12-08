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

########## Supplementary figure 2. A. ###########
# No filter local 
filt_aled <- final_list
cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)
plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"), output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/no_filter", device = "svg")

########## Supplementary figure 2. B. ###########
# mappability = 0.99, gc = 2, counts = 2
filt_aled <- filter_bins(final_list, mappability_threshold = 0.99, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/default_local/filtering.log")
cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output = T)
plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"), output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/default_local", device = "svg")

########## Supplementary figure 2. B. ###########
# mappability = 0.8, gc = 1, counts = 1
filt_aled <- filter_bins(final_list, mappability_threshold = 0.8, gc_threshold = 1, count_threshold = 1, log_file = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/stringent/filtering.log")
cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output = T)
plot_bins(cn_aled, view_sample = "NAT_X_SYN_042", plot_cn = T, specific_chr = c("1", "6"), method = "local", color_map = c("red2", "blue2"), output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/diff_filter/stringent", device = "svg")

######################################
# Get plot of local vs global median #
######################################

# Set to default filter
filt_aled <- filter_bins(final_list)
cn_aled <- get_copy_number(filt_aled, n_threads = 10, method = "local", full_output =T)

# Get the medians over each subgenome for each sample
local_medians <- get_sample_stats(cn_aled, method = "local_median")
medians <- get_sample_stats(cn_aled)

# Get the read counts for chromosomes 1 and 6
count_are <- cn_aled$A.arenosa$bins$NAT_X_SYN_004[cn_aled$A.arenosa$bins$chr=="6"]
count_thal <- cn_aled$A.thaliana$bins$NAT_X_SYN_004[cn_aled$A.thaliana$bins$chr=="1"]

# Get the bin positions for chromosomes 1 and 6
pos_are <- cn_aled$A.arenosa$bins$start[cn_aled$A.arenosa$bins$chr=="6"]
pos_thal <- cn_aled$A.thaliana$bins$start[cn_aled$A.thaliana$bins$chr=="1"]

########## Figure 2 ###########
# Plot it 
par(mfrow=c(1,2))
plot(pos_thal, count_thal, ylim=c(0,6000), pch=16, main="A.thaliana Subgenome Chromosome 1", ylab="Counts", xlab = "Position (bp)")
abline(h=medians$NAT_X_SYN_004$A.arenosa, col='red', lwd=1.5)
abline(h=local_medians$NAT_X_SYN_004$A.thaliana, col='green', lwd=3)
plot(pos_are, count_are, ylim=c(0,6000), pch=16, main="A.arenosa Subgenome Chromosome 1", ylab="", xlab = "Position (bp)")
abline(h=medians$NAT_X_SYN_004$A.arenosa, col='red', lwd=1.5)
abline(h=local_medians$NAT_X_SYN_004$A.arenosa, col='green', lwd=3)


#####################################
# Settings which match Nibau et al. #
#####################################

# Filtering
filt_aled <- filter_bins(final_list, mappability_threshold = 0.98, gc_threshold = 1.5)

########## Figure 3 ###########
# GC correction 
gc_corrected <- correct_gc(filt_aled, local_normalize = T, n_windows = 10, loess_span = 1, pch = 16, alpha = 0.3, cex = 0.6, output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/gc_filtering_test", device = "svg", output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/gc_filtering_test", )

# Copy Number
cn_aled <- get_copy_number(gc_corrected, n_threads = 10, method = "local", full_output =T)

# Remove short spans
wa <- remove_short_spans(cn_aled, max_length = 3)

########## Figure 4. B. ###########
# Plot all bins
plot_all_bins(heal_list = wa, plot_cn = T, method = "local", sample_name_size = 4,  
              view_samples =c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"),
                color_map = c("red2", "blue2"),
              add_bins = T, bin_point_size = 1.3, bin_point_alpha = 0.1, output_dir = "/srv/kenlab/kenji/exploring_aled_BMC_genomics_paper/", width = 10, height = 8, device = "svg")