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
wheat_raw <- read_heal_list("data/Wheat/healr_list")

# filter
wheat_filt <- filter_bins(wheat_raw, mappability_threshold = 0.85)

# Check GC correction effects. 
correct_gc(wheat_filt, n_windows = 10, loess_span = 1,
           n_threads = nThreads) # Relationship is flat. Ignore.

#####################
# Infer copy number #
#####################

# Get copy number using standard approach 
## Specifying to output the DNAcopy results (full_output = TRUE)

# Get CN with manual scale 
cn_manual <- get_copy_number(wheat_filt, n_threads = nThreads, 
                             manual_scale = c(0.5, 1.25, 2, 2.75, 3.5))

# Remove short spans
short_span_corrected <- remove_short_spans(cn_manual, max_length = 5)

######################
#    Get alignment   #
######################

# Get alignment
aln_wheat <- get_heal_alignment(short_span_corrected,
                                genespace_dir = "data/Wheat/syntenicHits/",
                                n_threads = nThreads)

######################
# Rename chromosomes #
######################

entries <- paste0("Group ", 1:7)
a_replace <- unique(short_span_corrected$A_subgenome$CN$chr)
d_replace <- unique(short_span_corrected$D_subgenome$CN$chr)
names(a_replace) <- names(d_replace) <- entries
for(i in 1:length(entries)){
  short_span_corrected$A_subgenome$CN$chr[short_span_corrected$A_subgenome$CN$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected$A_subgenome$CN$chr==a_replace[i]))
  short_span_corrected$D_subgenome$CN$chr[short_span_corrected$D_subgenome$CN$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected$D_subgenome$CN$chr==d_replace[i]))
  short_span_corrected$A_subgenome$bins$chr[short_span_corrected$A_subgenome$bins$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected$A_subgenome$bins$chr==a_replace[i]))
  short_span_corrected$D_subgenome$bins$chr[short_span_corrected$D_subgenome$bins$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected$D_subgenome$bins$chr==d_replace[i]))
  
  for(smp in names(aln_wheat)){
    aln_wheat[[smp]]$chr_A_subgenome[aln_wheat[[smp]]$chr_A_subgenome==a_replace[i]] <- rep(entries[i],sum(aln_wheat[[smp]]$chr_A_subgenome==a_replace[i]))
    aln_wheat[[smp]]$chr_D_subgenome[aln_wheat[[smp]]$chr_D_subgenome==d_replace[i]] <- rep(entries[i],sum(aln_wheat[[smp]]$chr_D_subgenome==d_replace[i]))
  }
}

##########################################
# Infer copy number without manual scale #
##########################################

# Get copy number with normal scale for comparison
cn_default <- get_copy_number(short_span_corrected, 
                         n_threads = nThreads, full_output = T)


#########################
#    Plot Copy Number   #
#########################

########## Supplementary Figure 4 ###########

plot_bins(cn_default, plot_cn = TRUE, ylim_max = 4.5,
          output_dir = "Supp_figure_4", add_DNAcopy = T,
          view_sample = "sample_H", specific_chr = "Group 2",
          color_map = c("darkred", "darkgreen"), n_threads = nThreads)

plot_bins(cn_default, plot_cn = TRUE, ylim_max = 4.5,
          output_dir = "Supp_figure_4", add_DNAcopy = T,
          view_sample = "sample_J", specific_chr = "Group 3",
          color_map = c("darkred", "darkgreen"), n_threads = nThreads)


########## Supplementary Figure 5 ###########

plot_bins(short_span_corrected, plot_cn = TRUE, ylim_max = 4.5,
          output_dir = "Supp_figure_5",
          view_sample = "sample_H", specific_chr = "Group 2",
          color_map = c("darkred", "darkgreen"), n_threads = nThreads)

plot_bins(short_span_corrected, plot_cn = TRUE, ylim_max = 4.5,
          output_dir = "Supp_figure_5",
          view_sample = "sample_J", specific_chr = "Group 3",
          color_map = c("darkred", "darkgreen"), n_threads = nThreads)


#######################
#    Plot Riparian    #
#######################

########## Figure 4.A & Supplementary Figures 8.A & 9.A ###########
plot_riparian(alignment = aln_wheat, heal_list = short_span_corrected,
              genespace_dir = "data/Wheat/syntenicRegions/",
              theme="dark", legend_text_size = 13, 
              title_size = 20, width = 20, height = 10,
              output_dir = "wheat_riparian_figures", device_vector = "png")