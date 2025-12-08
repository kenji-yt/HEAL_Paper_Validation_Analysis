#################
# Load Packages #
#################

library(healr)

############################
# Load Data and Merge runs #
############################

# We ran two HEAL analyses and merge the results here.

# Load data (these are heal lists created by HEAL)
H_300k <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/300k/healr_list")
JI_300k <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/300k/healr_list/")

# Merge 
JI_300k$A_subgenome$bins$sample_H <- H_300k$A_subgenome$bins$sample_H
JI_300k$D_subgenome$bins$sample_H <- H_300k$D_subgenome$bins$sample_H


####################
# Analyze and plot #
####################

# Filter
filt_IJ <- filter_bins(JI_300k, mappability_threshold = 0.85)

########## Supplementary figure 7 ###########
# Check GC effect
correct_gc(filt_IJ, n_windows = 10, loess_span = 1, pch = 16, alpha = 0.3, cex = 0.6, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/")

# Get CN with manual scale 
cn_manual <- get_copy_number(filt_IJ, n_threads = 1, full_output = T,
                             manual_scale = c(0.5, 1.25, 2, 2.75, 3.5))

# Remove short spans
short_span_corrected <- remove_short_spans(cn_manual, max_length = 5)

# Get alignment
aln_H <- get_heal_alignment(short_span_corrected, genespace_dir = "/srv/kenlab/kenji/syntenicHits/", n_threads = 10)

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
  
  for(smp in names(aln_H)){
    aln_H[[smp]]$chr_A_subgenome[aln_H[[smp]]$chr_A_subgenome==a_replace[i]] <- rep(entries[i],sum(aln_H[[smp]]$chr_A_subgenome==a_replace[i]))
    aln_H[[smp]]$chr_D_subgenome[aln_H[[smp]]$chr_D_subgenome==d_replace[i]] <- rep(entries[i],sum(aln_H[[smp]]$chr_D_subgenome==d_replace[i]))
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
plot_riparian(alignment = aln_H, heal_list = short_span_corrected, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", theme="dark", n_threads = 10, legend_text_size = 13, title_size = 20, width = 20, height = 10, output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/riparian")