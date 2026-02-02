# Assumes the working directory is 'HEAL_Paper_Validation_Analysis' 

#################
# Load Packages #
#################

library(healr)

###############
# Set threads #
###############

nThreads <- 10 # Edit to fit your resources 

#####################
# Load data, filter #
#####################

# Load in the heal lists (precounted by HEAL)
wgbs_list <- read_heal_list("data/A.kamchatica/healr_list")

# filter
wgbs_filt <- filter_bins(wgbs_list)

#####################
# Infer copy number #
#####################

# Get copy number on the filtered by not GC corrected data
wgbs_cn <- get_copy_number(wgbs_filt, n_threads = nThreads)

######################
#    Get alignment   #
######################

wgbs_aln <- get_heal_alignment(wgbs_cn, genespace_dir = "data/A.kamchatica/syntenicHits/", n_threads = nThreads)

#########################
#    Plot alignments    #
#########################

########## Figure 3.A ########## 

plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all",
               color_map = c("purple4", "orange3"), alpha = 0.1, specific_chr = "chr7",
               view_sample = "LL_RS2L_G1_3", output_dir = "figure3", device = "png",
               n_threads = nThreads) 

###################
#    HEAT MAPS    #
###################

########## Figure 3.B & Supplementary Figure 7 ###########

plot_heal_heat_map(alignment = wgbs_aln, , xrange=c(0:4), yrange=c(0:4),
                   output_dir = "figure3", device = "png")


################################
#    Genome Wide CN heatmap    #
################################

########## Figure 3.C ###########

plot_cn_heat(wgbs_cn, chr_labels = TRUE, subgenome_labels = TRUE, width = 12,
             height = 3, chr_limits = TRUE,  device = "png", output_dir = "figure3")