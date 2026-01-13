#################
# Load Packages #
#################

library(healr)

#####################
# Load data, filter #
#####################

# Load in the heal lists (precounted by HEAL)
wheat_raw <- read_heal_list("data/Wheat/healr_list")

# filter
wheat_filt <- filter_bins(wheat_raw, mappability_threshold = 0.85)

# Check GC correction effects. 
correct_gc(wheat_filt, n_windows = 10, loess_span = 1) # Relationship is flat. Ignore.

#####################
# Infer copy number #
#####################

# Get copy number using standard approach 
## Specifying to output the DNAcopy results (full_output = TRUE)
cn_default <- get_copy_number(wheat_filt, n_threads = 1, full_output = T)

########## Supplementary Figure 4 ###########
# Note that chromosome names are different here. We rename later.
plot_bins(cn_default, plot_cn = TRUE, output_dir = "Supp_figure_4",
          add_DNAcopy = T, color_map = c("darkred", "darkgreen"))
          
# Get CN with manual scale 
cn_manual <- get_copy_number(wheat_filt, n_threads = 1, 
                             manual_scale = c(0.5, 1.25, 2, 2.75, 3.5))

# Note that chromosome names are different here. We rename later.
plot_bins(cn_manual, plot_cn = TRUE, output_dir = "Supp_figure_5",
          color_map = c("darkred", "darkgreen"))

# Remove short spans
short_span_corrected <- remove_short_spans(cn_manual, max_length = 5)

# Get alignment
aln_wheat <- get_heal_alignment(short_span_corrected, genespace_dir = "data/Wheat/syntenicHits/", n_threads = 10)

# Rename chromosomes 
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

# Get copy number with normal scale for comparison
cn_IJ <- get_copy_number(short_span_corrected, n_threads = 1, full_output = T)

########## Figure 5 ###########
# Plot without manual scale
plot_bins(cn_IJ, view_sample = "sample_I", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines")
plot_bins(cn_IJ, view_sample = "sample_H", specific_chr = "Group 2", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines")
plot_bins(cn_IJ, view_sample = "sample_J", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines")

########## Figure 6 ###########
# Plot with manual scale
plot_bins(short_span_corrected, view_sample = "sample_I", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = F, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines/corrected")
plot_bins(short_span_corrected, view_sample = "sample_H", specific_chr = "Group 2", plot_cn = T, add_DNAcopy = F, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines/corrected")
plot_bins(short_span_corrected, view_sample = "sample_J", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = F, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines/corrected")

########## Figure 7, 8 and 9 ###########
### Riparian
plot_riparian(alignment = aln_wheat, heal_list = short_span_corrected, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", theme="dark", n_threads = 10, legend_text_size = 13, title_size = 20, width = 20, height = 10, output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/riparian")