# Assumes the working directory is 'HEAL_Paper_Validation_Analysis' 

#################
# Load Packages #
#################

library(healr)

###############
# Set threads #
###############

nThreads <- 10 # Edit to fit your resources 

####################################
# Load data and rename chromosomes #
####################################

final_list <- read_heal_list("data/A.suecica/healr_list/")

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


#######################################
#      Filtering & GC correction      #
#######################################

# Filtering
filt_aled <- filter_bins(final_list, 
                         mappability_threshold = 0.98, gc_threshold = 1.5)

########## Supplementary Figure 2 ###########
# GC correction 
dir.create("Supp_figure_2")
gc_corrected <- correct_gc(filt_aled, local_normalize = TRUE, n_windows = 10,
                           loess_span = 1, pch = 16, alpha = 0.3, cex = 0.6,
                           output_dir = "Supp_figure_2", device = "pdf")


######################################
# Get plot of local vs global median #
######################################

# Set to default filter
filt_aled <- filter_bins(final_list)
cn_aled <- get_copy_number(filt_aled, n_threads = nThreads,
                           method = "local", full_output =T)

# Get the medians over each subgenome for each sample
local_medians <- get_sample_stats(cn_aled, method = "local_median")
medians <- get_sample_stats(cn_aled)

# Get the read counts for chromosomes 1 and 6
count_are <- cn_aled$A.arenosa$bins$NAT_X_SYN_004[cn_aled$A.arenosa$bins$chr=="6"]
count_thal <- cn_aled$A.thaliana$bins$NAT_X_SYN_004[cn_aled$A.thaliana$bins$chr=="1"]

# Get the bin positions for chromosomes 1 and 6
pos_are <- cn_aled$A.arenosa$bins$start[cn_aled$A.arenosa$bins$chr=="6"]
pos_thal <- cn_aled$A.thaliana$bins$start[cn_aled$A.thaliana$bins$chr=="1"]

########## Supplementary Figure 3 ###########

par(mfrow=c(1,2))
plot(pos_thal, count_thal, ylim=c(0,6000), pch=16, main="A.thaliana Subgenome Chromosome 1", ylab="Counts", xlab = "Position (bp)")
abline(h=medians$NAT_X_SYN_004$A.arenosa, col='red', lwd=1.5)
abline(h=local_medians$NAT_X_SYN_004$A.thaliana, col='green', lwd=3)
plot(pos_are, count_are, ylim=c(0,6000), pch=16, main="A.arenosa Subgenome Chromosome 1", ylab="", xlab = "Position (bp)")
abline(h=medians$NAT_X_SYN_004$A.arenosa, col='red', lwd=1.5)
abline(h=local_medians$NAT_X_SYN_004$A.arenosa, col='green', lwd=3)


#####################
# Infer copy number #
#####################

# Copy Number
cn_aled <- get_copy_number(gc_corrected, n_threads = nThreads,
                           method = "local", full_output =T)

# Remove short spans
wa <- remove_short_spans(cn_aled, max_length = 3)


#########################
#     Plot all bins     #
#########################

########## Figure 2. B. ###########
# Plot all bins
dir.create("figure2")
plot_all_bins(heal_list = wa, plot_cn = T, method = "local", sample_name_size = 4,  
              view_samples =c("NAT_X_SYN_004", "NAT_X_SYN_051","NAT_X_SYN_042","NAT_X_SYN_054", "NAT_X_SYN_044"),
                color_map = c("red2", "blue2"),
              add_bins = T, bin_point_size = 1.3, bin_point_alpha = 0.1,
              output_dir = "figure2", width = 10, height = 8, 
              device = "png")