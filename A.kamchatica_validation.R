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

########## Figure 10.A & Supplementary fig.5 ###########

plot_alignment(alignment = wgbs_aln, heal_list = wgbs_cn, add_bins = "all", color_map = c("purple4", "orange3"), alpha = 0.1, output_dir = "figure10", device = "svg") 

###################
#    HEAT MAPS    #
###################

########## Figure 10.B ###########

plot_heal_heat_map(alignment = wgbs_aln, , xrange=c(0:4), yrange=c(0:4), output_dir = "figure10", device = "svg")


###############################
#    Interesting statistics   #
###############################

# How many concordant bins are changes in CN
100*c(sum(wgbs_aln$LL_RS2L_G1_1$status=="concordant") / length(wgbs_aln$LL_RS2L_G1_1$status), 
                                      sum(wgbs_aln$LL_RS2L_G1_2$status=="concordant") / length(wgbs_aln$LL_RS2L_G1_2$status),
                                      sum(wgbs_aln$LL_RS2L_G1_2$status=="concordant") / length(wgbs_aln$LL_RS2L_G1_2$status))

# How many concordant bins are changes in CN
100*(c(1-sum(wgbs_aln$LL_RS2L_G1_1$status=="concordant" & wgbs_aln$LL_RS2L_G1_1$cn_A.halleri == 2) / sum(wgbs_aln$LL_RS2L_G1_1$status=="concordant"), 
  1-sum(wgbs_aln$LL_RS2L_G1_2$status=="concordant" & wgbs_aln$LL_RS2L_G1_2$cn_A.halleri == 2) / sum(wgbs_aln$LL_RS2L_G1_2$status=="concordant"),
  1-sum(wgbs_aln$LL_RS2L_G1_3$status=="concordant" & wgbs_aln$LL_RS2L_G1_3$cn_A.halleri == 2) / sum(wgbs_aln$LL_RS2L_G1_3$status=="concordant")))


# The percentage of anchor pairs set using concordance optimization
100*c(table(wgbs_aln$LL_RS2L_G1_1$method)[["multiple_concordant"]]/length(wgbs_aln$LL_RS2L_G1_1$method),
                            table(wgbs_aln$LL_RS2L_G1_2$method)[["multiple_concordant"]]/length(wgbs_aln$LL_RS2L_G1_2$method),
                            table(wgbs_aln$LL_RS2L_G1_3$method)[["multiple_concordant"]]/length(wgbs_aln$LL_RS2L_G1_3$method))
