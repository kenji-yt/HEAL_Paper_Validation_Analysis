# Assumes the working directory is 'HEAL_Paper_Validation_Analysis' 

#################
# Load Packages #
#################

library(devtools)
library(data.table)
library(gridExtra)
library(healr)

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
wgbs_cn <- get_copy_number(wgbs_filt)

######################
#    Get alignment   #
######################

wgbs_aln <- get_heal_alignment(wgbs_cn, genespace_dir = "data/A.kamchatica/syntenicHits/", n_threads = 10)

#########################
#    Plot alignments    #
#########################

########## Figure 3.A ########## 

plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all", color_map = c("purple4", "orange3"), alpha = 0.1, output_dir = "figure3", device = "svg") 

###################
#    HEAT MAPS    #
###################

########## Figure 3.B ###########

plot_heal_heat_map(alignment = wgbs_aln, , xrange=c(0:4), yrange=c(0:4), output_dir = "figure3", device = "svg")


################################
#    Genome Wide CN heatmap    #
################################

########## Figure 3.C ###########

plot_cn_heat(wgbs_cn, output_dir = "figure3")