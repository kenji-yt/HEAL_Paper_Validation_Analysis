library(devtools)
load_all()

# # Function to add red lines to plot based on DNA copy:
# add_DNAcopy_lines <- function(heal_list, sample, progenitor, chromosome, median_val){
#   
#   if(is.null(heal_list[[progenitor]]$DNAcopy)){
#     stop("The heal_list has DNAcopy object. Exiting..")
#   }
#   
#   which_index <- heal_list[[progenitor]]$DNAcopy[[sample]]$output$chrom==chromosome
#   
#   lengths <- heal_list[[progenitor]]$DNAcopy[[sample]]$output$num.mark[which_index]
#   values <- heal_list[[progenitor]]$DNAcopy[[sample]]$output$seg.mean[which_index]/median_val*2
#   line_vec <- inverse.rle(list(lengths = lengths, values = values))
#   
#   lines(1:length(line_vec), line_vec, col="red")
# }

H_300k <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/300k/healr_list")
JI_300k <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/300k/healr_list/")

### Merge 
JI_300k$A_subgenome$bins$sample_H <- H_300k$A_subgenome$bins$sample_H
JI_300k$D_subgenome$bins$sample_H <- H_300k$D_subgenome$bins$sample_H

# Proceed
filt_IJ <- filter_bins(JI_300k, mappability_threshold = 0.85)

# No effect
correct_gc(filt_IJ, n_windows = 10, loess_span = 1, pch = 16, alpha = 0.3, cex = 0.6, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/")

cn_IJ <- get_copy_number(filt_IJ, n_threads = 1, full_output = T)

plot_bins(cn_IJ, view_sample = "sample_I", specific_chr = "NC_053024.1", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3)
plot_bins(cn_IJ, view_sample = "sample_I", specific_chr = "NC_053037.3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3)

plot_bins(cn_IJ, view_sample = "sample_H", specific_chr = "NC_053023.1", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 6)
plot_bins(cn_IJ, view_sample = "sample_H", specific_chr = "NC_053036.3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 6)

plot_bins(cn_IJ, view_sample = "sample_J", specific_chr = "NC_053024.1", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 6)
plot_bins(cn_IJ, view_sample = "sample_J", specific_chr = "NC_053037.3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 6)

# Get manual scale
cn_manual <- get_copy_number(filt_IJ, n_threads = 1, full_output = T,
                             manual_scale = c(0.5, 1.25, 2, 2.75, 3.5))

cn_ns <- remove_short_spans(cn_manual, max_length = 5)

plot_bins(cn_ns, view_sample = "sample_I", specific_chr = "NC_053024.1", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"))
plot_bins(cn_ns, view_sample = "sample_I", specific_chr = "NC_053037.3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"))

plot_bins(cn_ns, view_sample = "sample_H", specific_chr = "NC_053023.1", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"))
plot_bins(cn_ns, view_sample = "sample_H", specific_chr = "NC_053036.3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"))

plot_bins(cn_ns, view_sample = "sample_J", specific_chr = "NC_053024.1", plot_cn = F, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"))
plot_bins(cn_ns, view_sample = "sample_J", specific_chr = "NC_053037.3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"))

# wa <- cn_IJ$A_subgenome$bins[cn_IJ$A_subgenome$bins$chr=="NC_053024.1",]
# we <- cn_IJ$D_subgenome$bins[cn_IJ$D_subgenome$bins$chr=="NC_053037.3",]
# medians <- get_sample_stats(cn_IJ)
# par(mfrow=c(2,2))
# plot(na.omit(wa$sample_I)/medians$sample_I$A_subgenome*2, pch=16, ylim = c(-0.2, 4.2), ylab ="Normalized Read Count", main = "SAMPLE I")
# abline(h=2)
# abline(h = 0, lty = 2)
# abline(h = 1, lty = 2)
# abline(h = 3, lty = 2)
# abline(h = 4, lty = 2)
# add_DNAcopy_lines(cn_IJ, "sample_I", "A_subgenome", chromosome = "NC_053024.1", median_val = medians$sample_I$A_subgenome)
# plot(wa$sample_J/medians$sample_J$A_subgenome*2, pch=16, ylim = c(-0.2, 4.2), ylab ="Normalized Read Count", main = "SAMPLE J")
# abline(h=2)
# abline(h = 0, lty = 2)
# abline(h = 1, lty = 2)
# abline(h = 3, lty = 2)
# abline(h = 4, lty = 2)
# add_DNAcopy_lines(cn_IJ, "sample_J", "A_subgenome", chromosome = "NC_053024.1", median_val = medians$sample_I$A_subgenome)
# plot(na.omit(we$sample_I)/medians$sample_I$A_subgenome*2, pch=16, ylim = c(-0.2, 4.2), ylab ="Normalized Read Count")
# abline(h=2)
# abline(h = 0, lty = 2)
# abline(h = 1, lty = 2)
# abline(h = 3, lty = 2)
# abline(h = 4, lty = 2)
# add_DNAcopy_lines(cn_IJ, "sample_I", "D_subgenome", chromosome = "NC_053037.3", median_val = medians$sample_I$A_subgenome)
# plot(na.omit(we$sample_J)/medians$sample_J$A_subgenome*2, pch=16, ylim = c(-0.2, 4.2), ylab ="Normalized Read Count")
# abline(h=2)
# abline(h = 0, lty = 2)
# abline(h = 1, lty = 2)
# abline(h = 3, lty = 2)
# abline(h = 4, lty = 2)
# add_DNAcopy_lines(cn_IJ, "sample_J", "D_subgenome", chromosome = "NC_053037.3", median_val = medians$sample_I$A_subgenome)

# Set manual scale: 
# cn_IJ <- get_copy_number(filt_IJ, n_threads = 1, full_output = T, manual_scale = c(0.5, 1.3, 2, 2.7, 3.3))
# plot_all_bins(cn_IJ, plot_cn = T)
# 
# abline(h=2000, col="green")
# 
# cn_IJ$A_subgenome$DNAcopy$sample_I$output[cn_IJ$A_subgenome$DNAcopy$sample_I$output$chrom=="NC_053023.1",]
# #plot_all_bins(cn_H, plot_cn = T)
# plot_all_bins(cn_IJ, plot_cn = T, view_samples = "sample_J")
# plot_all_bins(cn_IJ, plot_cn = T, view_samples = "sample_I")

short_span_corrected <- remove_short_spans(cn_manual, max_length = 5)
plot_all_bins(short_span_corrected, plot_cn = T)

aln_H <- get_heal_alignment(short_span_corrected, genespace_dir = "/srv/kenlab/kenji/syntenicHits/", n_threads = 10)

# rename chromosomes 
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

cn_IJ <- get_copy_number(short_span_corrected, n_threads = 1, full_output = T)

# Before manual
plot_bins(cn_IJ, view_sample = "sample_I", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines")
plot_bins(cn_IJ, view_sample = "sample_H", specific_chr = "Group 2", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines")
plot_bins(cn_IJ, view_sample = "sample_J", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = T, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines")

# After manual and short span
plot_bins(short_span_corrected, view_sample = "sample_I", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = F, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines/corrected")
plot_bins(short_span_corrected, view_sample = "sample_H", specific_chr = "Group 2", plot_cn = T, add_DNAcopy = F, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines/corrected")
plot_bins(short_span_corrected, view_sample = "sample_J", specific_chr = "Group 3", plot_cn = T, add_DNAcopy = F, add_bins = T, color_map = c("darkred", "darkgreen"), ylim_max = 4.5, linewidth = 3, width = 6.70, height = 5.23, device = "svg", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/DNAcopy_lines/corrected")

### Riparian
plot_riparian(alignment = aln_H, heal_list = short_span_corrected, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", theme="dark", n_threads = 10, legend_text_size = 13, title_size = 20, width = 20, height = 10, output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/figures/riparian")



# 
# make_all_plots <- function(path_to_healr_list){
#   
#   H_list <- read_heal_list(path_to_healr_list)
# 
#   filt_H <- filter_bins(H_list, mappability_threshold = 0.8)
#   
#   start <- c(0.1, 0.2, 0.3, 0.4)
#   end <- 4-start
#   midlow <- c(1.1, 1.2, 1.3, 1.4)
#   midend <- 4-midlow
#   
#   for(d in 1:length(start)){
#     
#     for(j in 1:length(midlow)){
#       
#       cat(paste("d:", d, "\n"))
#       cat(paste("j:", j, "\n"))
#       
#       manual_scale <- c(start[d], midlow[j], 2, midend[j], end[d])
#       
#       cn_H <- get_copy_number(filt_H, n_threads = 10, manual_scale = manual_scale)
#       
#       spans <- c(1,3,5)
#       
#       for(s in spans){
#         
#         short_span_corrected <- remove_short_spans(cn_H, max_length = s)
#         
#         aln_H <- get_heal_alignment(short_span_corrected, genespace_dir = "/srv/kenlab/kenji/syntenicHits/", n_threads = 10)
#         
#         # rename chromosomes 
#         entries <- paste0("Group_", 1:7)
#         a_replace <- unique(short_span_corrected$A_subgenome$CN$chr)
#         d_replace <- unique(short_span_corrected$D_subgenome$CN$chr)
#         names(a_replace) <- names(d_replace) <- entries
#         for(i in 1:length(entries)){
#           short_span_corrected$A_subgenome$CN$chr[short_span_corrected$A_subgenome$CN$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected$A_subgenome$CN$chr==a_replace[i]))
#           short_span_corrected$D_subgenome$CN$chr[short_span_corrected$D_subgenome$CN$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected$D_subgenome$CN$chr==d_replace[i]))
#           short_span_corrected$A_subgenome$bins$chr[short_span_corrected$A_subgenome$bins$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected$A_subgenome$bins$chr==a_replace[i]))
#           short_span_corrected$D_subgenome$bins$chr[short_span_corrected$D_subgenome$bins$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected$D_subgenome$bins$chr==d_replace[i]))
#           
#           for(smp in names(aln_H)){
#             aln_H[[smp]]$chr_A_subgenome[aln_H[[smp]]$chr_A_subgenome==a_replace[i]] <- rep(entries[i],sum(aln_H[[smp]]$chr_A_subgenome==a_replace[i]))
#             aln_H[[smp]]$chr_D_subgenome[aln_H[[smp]]$chr_D_subgenome==d_replace[i]] <- rep(entries[i],sum(aln_H[[smp]]$chr_D_subgenome==d_replace[i]))
#           }
#         }
#         
#         out_path <- paste0(dirname(path_to_healr_list), "/riparian_", manual_scale[1], "_", manual_scale[2], "_", "_span_filt_", s)
#         
#         plot_riparian(heal_alignment = aln_H, heal_list = short_span_corrected, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", output_dir = out_path, theme="dark")
#         
#         
#       }
#       
#     }
#   }
# # }
# 
# # For the poster!
# plot_riparian(heal_alignment = aln_H, heal_list = short_span_corrected, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", theme="dark", n_threads = 10, legend_text_size = 10, title_size = 10, width = 20, height = 10, output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/test_HB_rip_poster")
# # sample_H for 300k and 500k
# # make_all_plots("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/300k/healr_list/")
# # #make_all_plots("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/500k/healr_list/")
# path_to_healr_list <- "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/300k/healr_list/"
# # sample_I and J for 300k and 500k
# make_all_plots("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/300k/healr_list/")
# make_all_plots("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/500k/healr_list/")
# # 
# # 
# 
# 
# start <- c(0, 0.1, 0.2, 0.3, 0.4)
# end <- 4-start
# midlow <- c(1, 1.1, 1.2, 1.3, 1.4)
# midend <- 4-midlow
# 
# for(i in 1:length(start)){
#   for(j in 1:length(midlow)){
#     
#     manual_scale <- c(start[i], midlow[j], 2, midend[j], end[i])
#     print(c(manual_scale[1],manual_scale[2]))
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# # H
# H_list <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/300k/healr_list/")
# 
# #filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)
# filt_H <- filter_bins(H_list, mappability_threshold = 0.8)#, gc_quantile = 70)#, mappability_threshold = 0.95)
# #gc_H <- correct_gc(filt_H) 
# 
# cn_H <- get_copy_number(filt_H, n_threads = 10, manual_scale = manual_scale, full_output = T)
# 
# #Hp<- plot_all_bins(cn_H, plot_cn = T, add_bins = T, return_plot = T)
# 
# 
# #plot_bins(cn_H, view_sample = "sample_H", specific_chr = "NC_053023.1", add_bins = T, plot_cn = T)
# short_span_corrected <- remove_short_spans(cn_H, max_length = 2)
# #Hsp<- plot_all_bins(short_span_corrected, plot_cn = T, add_bins = T, return_plot = T)
# aln_H <- get_heal_alignment(short_span_corrected, genespace_dir = "/srv/kenlab/kenji/syntenicHits/", n_threads = 10)
# 
# # rename chromosomes 
# entries <- paste0("Group_", 1:7)
# a_replace <- unique(short_span_corrected$A_subgenome$CN$chr)
# d_replace <- unique(short_span_corrected$D_subgenome$CN$chr)
# names(a_replace) <- names(d_replace) <- entries
# for(i in 1:length(entries)){
#   short_span_corrected$A_subgenome$CN$chr[short_span_corrected$A_subgenome$CN$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected$A_subgenome$CN$chr==a_replace[i]))
#   short_span_corrected$D_subgenome$CN$chr[short_span_corrected$D_subgenome$CN$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected$D_subgenome$CN$chr==d_replace[i]))
#   short_span_corrected$A_subgenome$bins$chr[short_span_corrected$A_subgenome$bins$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected$A_subgenome$bins$chr==a_replace[i]))
#   short_span_corrected$D_subgenome$bins$chr[short_span_corrected$D_subgenome$bins$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected$D_subgenome$bins$chr==d_replace[i]))
#   
#   for(smp in names(aln_H)){
#     aln_H[[smp]]$chr_A_subgenome[aln_H[[smp]]$chr_A_subgenome==a_replace[i]] <- rep(entries[i],sum(aln_H[[smp]]$chr_A_subgenome==a_replace[i]))
#     aln_H[[smp]]$chr_D_subgenome[aln_H[[smp]]$chr_D_subgenome==d_replace[i]] <- rep(entries[i],sum(aln_H[[smp]]$chr_D_subgenome==d_replace[i]))
#   }
# }
# 
# plot_bins(cn_H, view_sample = "sample_H", specific_chr = "NC_053025.1", plot_cn = T, add_bins = T, color_map = c("darkred", "darkgreen"))
# plot_bins(cn_H, view_sample = "sample_H", specific_chr = "NC_053039.3", plot_cn = T, add_bins = T, color_map = c("darkred", "darkgreen"))
# 
# plot_bins(short_span_corrected, view_sample = "sample_H", specific_chr = "Group_6", plot_cn = T, add_bins = T, color_map = c("darkred", "darkgreen"))
# plot_riparian(heal_alignment = aln_H, heal_list = short_span_corrected, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/300k_H_rip", theme="dark")
# 
# 
# # H 300k
# H3k_list <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/300k/healr_list/")
# 
# #filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)
# filt3_H <- filter_bins(H3k_list, mappability_threshold = 0.80)#, gc_quantile = 70)#, mappability_threshold = 0.95)
# #gc_H <- correct_gc(filt_H) 
# 
# cn3_H <- get_copy_number(filt3_H, n_threads = 10, manual_scale = manual_scale, full_output = T)
# 
# p3 <- plot_all_bins(cn3_H, plot_cn = T, add_bins = T, return_plot = T)
# 
# # H 200k
# H2k_list <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_H_check/200k/healr_list/")
# 
# #filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)
# filt2_H <- filter_bins(H2k_list, mappability_threshold = 0.80)#, gc_quantile = 70)#, mappability_threshold = 0.95)
# #gc_H <- correct_gc(filt_H) 
# 
# cn2_H <- get_copy_number(filt2_H, n_threads = 10, manual_scale = manual_scale, full_output = T)
# 
# p2 <- plot_all_bins(cn2_H, plot_cn = T, add_bins = T, return_plot = T)
# 
# # 200k is better because of group 7 going to 4 on D. 200k does have a bit more small scale "noisish" looking things. Namely on group 4 and 5 edges loss and in the group 2 D deletion. 
# 
# grid.arrange(p2, p3)
# 
# #### J and I 
# 
# # JI 200k
# JI2k_list <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/200k/healr_list/")
# 
# #filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)
# filt2_JI <- filter_bins(JI2k_list, mappability_threshold = 0.80)#, gc_quantile = 70)#, mappability_threshold = 0.95)
# #gc_H <- correct_gc(filt_H) 
# 
# cn2_JI <- get_copy_number(filt2_JI, n_threads = 10, manual_scale = manual_scale, full_output = T)
# 
# p2J <- plot_all_bins(cn2_JI, plot_cn = T, view_samples = "sample_J", add_bins = T, return_plot = T)
# p2I <- plot_all_bins(cn2_JI, plot_cn = T, view_samples = "sample_I", add_bins = T, return_plot = T)
# 
# # JI 300k
# JI3k_list <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/300k/healr_list/")
# 
# #filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)
# filt3_JI <- filter_bins(JI3k_list, mappability_threshold = 0.80)#, gc_quantile = 70)#, mappability_threshold = 0.95)
# #gc_H <- correct_gc(filt_H) 
# 
# cn3_JI <- get_copy_number(filt3_JI, n_threads = 10, manual_scale = manual_scale, full_output = T)
# 
# p3J <- plot_all_bins(cn3_JI, plot_cn = T, view_samples = "sample_J", add_bins = T, return_plot = T)
# p3I <- plot_all_bins(cn3_JI, plot_cn = T, view_samples = "sample_I", add_bins = T, return_plot = T)
# 
# grid.arrange(p2J, p3J) # Quite messy for both... Remove short spans then check again
# grid.arrange(p2I, p3I)
# 
# short_span_corrected_JI_2 <- remove_short_spans(cn2_JI, max_length = 2)
# 
# short_span_corrected_JI_3 <- remove_short_spans(cn3_JI, max_length = 2)
# 
# 
# p2J_no_s <- plot_all_bins(short_span_corrected_JI_2, plot_cn = T, view_samples = "sample_J", add_bins = T, return_plot = T)
# p2I_no_s <- plot_all_bins(short_span_corrected_JI_2, plot_cn = T, view_samples = "sample_I", add_bins = T, return_plot = T)
# 
# p3J_no_s <- plot_all_bins(short_span_corrected_JI_3, plot_cn = T, view_samples = "sample_J", add_bins = T, return_plot = T)
# p3I_no_s <- plot_all_bins(short_span_corrected_JI_3, plot_cn = T, view_samples = "sample_I", add_bins = T, return_plot = T)
# 
# 
# # JI 500k
# JI5k_list <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/500k/healr_list/")
# 
# #filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)
# filt5_JI <- filter_bins(JI5k_list, mappability_threshold = 0.80)#, gc_quantile = 70)#, mappability_threshold = 0.95)
# #gc_H <- correct_gc(filt_H) 
# 
# cn5_JI <- get_copy_number(filt5_JI, n_threads = 10, manual_scale = manual_scale, full_output = T)
# 
# p5J <- plot_all_bins(cn5_JI, plot_cn = T, view_samples = "sample_J", add_bins = T, return_plot = T)
# p5I <- plot_all_bins(cn5_JI, plot_cn = T, view_samples = "sample_I", add_bins = T, return_plot = T)
# 
# short_span_corrected_JI_5 <- remove_short_spans(cn5_JI, max_length = 5)
# 
# p5I_no_s <- plot_all_bins(short_span_corrected_JI_5, plot_cn = T, view_samples = "sample_I", add_bins = T, return_plot = T)
# p5J_no_s <- plot_all_bins(short_span_corrected_JI_5, plot_cn = T, view_samples = "sample_J", add_bins = T, return_plot = T)
# 
# grid.arrange(p2J_no_s, p2J) # no short s better
# grid.arrange(p2I_no_s, p2I) # no short s better
# grid.arrange(p3J_no_s, p3J) # no short s better
# grid.arrange(p3I_no_s, p3I) # no short s better
# grid.arrange(p5J_no_s, p5J) # no short s better
# grid.arrange(p5I_no_s, p5I) # no short s better
# 
# grid.arrange(p3J_no_s, p2J_no_s) # 300k is better than 200k!
# grid.arrange(p5J_no_s, p3J_no_s) # 300k is better than 500k!
# grid.arrange(p3I_no_s, p2I_no_s) # 300k is better than 200k!
# grid.arrange(p5I_no_s, p3I_no_s) # 500k is better than 300k!
# 
# # Overall 300k is better.
# # 
# # 
# # # I
# # I_list <- read_heal_list("/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/healr_list")
# # 
# # #filt_aled <- filter_bins(final_list, count_threshold = 1, gc_quantile = 70, mappability_threshold = 0.99)
# # filt_I <- filter_bins(I_list)#, gc_quantile = 40, mappability_threshold = 0.92)
# # #gc_I <- correct_gc(filt_I) 
# # 
# # cn_I <- get_copy_number(filt_I, n_threads = 10)#, manual_scale = manual_scale, full_output = T)
# # plot_all_bins(cn_I, plot_cn = T, add_bins = T)
# # plot_bins(cn_I, view_sample = "sample_J", specific_chr = "NC_053041.3", add_bins = T, plot_cn = T)
# # plot_bins(cn_I, view_sample = "sample_J", specific_chr = "NC_053041.3", add_bins = T, plot_cn = F, color_map = c("red", "green"))
# # plot_bins(gc_I, view_sample = "sample_J", specific_chr = "NC_053041.3", add_bins = T, plot_cn = F, color_map = c("red", "green"))
# # plot_bins(cn_I, view_sample = "sample_J", specific_chr = "NC_053041.3", add_bins = T, plot_cn = T) 
# 
# short_span_corrected_I <- short_span_corrected_JI_3
# 
# aln_I <- get_heal_alignment(short_span_corrected_I, genespace_dir = "/srv/kenlab/kenji/syntenicHits/", n_threads = 10)
# 
# 
# # rename chromosomes 
# entries <- paste0("Group_", 1:7)
# a_replace <- unique(short_span_corrected_I$A_subgenome$CN$chr)
# d_replace <- unique(short_span_corrected_I$D_subgenome$CN$chr)
# names(a_replace) <- names(d_replace) <- entries
# for(i in 1:length(entries)){
#   short_span_corrected_I$A_subgenome$CN$chr[short_span_corrected_I$A_subgenome$CN$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$A_subgenome$CN$chr==a_replace[i]))
#   short_span_corrected_I$D_subgenome$CN$chr[short_span_corrected_I$D_subgenome$CN$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$D_subgenome$CN$chr==d_replace[i]))
#   short_span_corrected_I$A_subgenome$bins$chr[short_span_corrected_I$A_subgenome$bins$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$A_subgenome$bins$chr==a_replace[i]))
#   short_span_corrected_I$D_subgenome$bins$chr[short_span_corrected_I$D_subgenome$bins$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$D_subgenome$bins$chr==d_replace[i]))
#   
#   for(smp in names(aln_I)){
#     aln_I[[smp]]$chr_A_subgenome[aln_I[[smp]]$chr_A_subgenome==a_replace[i]] <- rep(entries[i],sum(aln_I[[smp]]$chr_A_subgenome==a_replace[i]))
#     aln_I[[smp]]$chr_D_subgenome[aln_I[[smp]]$chr_D_subgenome==d_replace[i]] <- rep(entries[i],sum(aln_I[[smp]]$chr_D_subgenome==d_replace[i]))
#   }
# }
# 
# 
# plot_riparian(heal_alignment = aln_I, heal_list = short_span_corrected_I, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", theme = "dark", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/300k/300k_JI_rip")
# 
# 
# ### 500k
# 
# short_span_corrected_I <- short_span_corrected_JI_5
# 
# aln_I <- get_heal_alignment(short_span_corrected_I, genespace_dir = "/srv/kenlab/kenji/syntenicHits/", n_threads = 10)
# 
# 
# # rename chromosomes 
# entries <- paste0("Group_", 1:7)
# a_replace <- unique(short_span_corrected_I$A_subgenome$CN$chr)
# d_replace <- unique(short_span_corrected_I$D_subgenome$CN$chr)
# names(a_replace) <- names(d_replace) <- entries
# for(i in 1:length(entries)){
#   short_span_corrected_I$A_subgenome$CN$chr[short_span_corrected_I$A_subgenome$CN$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$A_subgenome$CN$chr==a_replace[i]))
#   short_span_corrected_I$D_subgenome$CN$chr[short_span_corrected_I$D_subgenome$CN$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$D_subgenome$CN$chr==d_replace[i]))
#   short_span_corrected_I$A_subgenome$bins$chr[short_span_corrected_I$A_subgenome$bins$chr==a_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$A_subgenome$bins$chr==a_replace[i]))
#   short_span_corrected_I$D_subgenome$bins$chr[short_span_corrected_I$D_subgenome$bins$chr==d_replace[i]] <- rep(entries[i],sum(short_span_corrected_I$D_subgenome$bins$chr==d_replace[i]))
#   
#   for(smp in names(aln_I)){
#     aln_I[[smp]]$chr_A_subgenome[aln_I[[smp]]$chr_A_subgenome==a_replace[i]] <- rep(entries[i],sum(aln_I[[smp]]$chr_A_subgenome==a_replace[i]))
#     aln_I[[smp]]$chr_D_subgenome[aln_I[[smp]]$chr_D_subgenome==d_replace[i]] <- rep(entries[i],sum(aln_I[[smp]]$chr_D_subgenome==d_replace[i]))
#   }
# }
# 
# 
# plot_riparian(heal_alignment = aln_I, heal_list = short_span_corrected_I, genespace_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/genespace_dir/", theme = "dark", output_dir = "/srv/kenlab/kenji/exploring_wheat_BMC_genomic_paper/wheat_results/sample_J_I_check/500k/500k_JI_rip")
