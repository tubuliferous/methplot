# Turn this into package and publish on github

# *Also consider generating script that gets gene body enrichments
#   Aggregated (marginalized) over all genes - would need to compress
#   or standardize gene ranges - could get quantiles for each gene,
#   e.g. break into 1000 bins for each gene transcript length plus
#   a little flanking region




library(data.table)
library(plyr)
library(ggplot2)
library(dplyr)
library(parallel)
library(grid)
library(reshape2)
library(methods) # Note: 'methods' package required for Rscript
library(intervals)
require(colorout)
library(readr)
library(stringr)
library(devtools)
# install_github("tslumley/xkcdcolors")
library(xkcdcolors)
df <- data.frame(x = rnorm(100), y = rnorm(100), group = as.factor(sample(2, 100, replace = T)))
color_palette <- c(name2color("black"), name2color("radioactive"))
ggplot(df, aes(x,y)) + geom_point(aes(color=group)) + scale_color_manual(values = color_palette)

# Small helper functions
autosomes <- function(this_table){
  autosomes <- paste("chr", 1:22, sep="") %>% c(., 1:22)
  return(this_table[which(as.character(data.frame(this_table)[,"chr"]) %in% autosomes), ])
}
normal_chroms <- function(this_table){
  auto_x_y <- paste("chr", 1:22, sep="") %>% c(., 1:22) %>% c(., paste("chr", c("X", "Y"), sep="")) %>% c(., "X", "Y")
  return(this_table[which(as.character(data.frame(this_table)[,"chr"]) %in% auto_x_y), ])
}
get_grch_gene_table <- function(grch37_table_path){
  grch37_table <- read_delim(grch37_table_path, delim="\t", col_names=F)
  names(grch37_table) <- c("id", "source", "type", "start", "end", "dot", "strand", "dot_2", "info")
  chrs <- str_match(grch37_table$id, "NC_0+(.+)\\.")[,2]
  chrs[which(chrs=="12920")] <- "M"
  chrs[which(chrs=="23")] <- "X"
  chrs[which(chrs=="24")] <- "Y"
  chrs <- paste("chr", chrs, sep="")
  grch37_table <- grch37_table %>%
                  mutate(chr = chrs) %>%
                  filter(type=="gene") %>%
                  mutate(gene_name = str_match(.$info, "Name=(.+?);")[,2])

  out <- grch37_table[grep("pseudo=true", grch37_table$info, invert=TRUE), ]
  return(out)
}
get_meth_table <- function(meth_path){
  this_meth_table <- fread(meth_path, sep="\t") %>%
                       rename(chr=V1, start=V2, end=V3, perc_meth=V4, meth=V5, unmeth=V6)
}
get_ucsc_refgene_table <-function(refgene_path){
  fread(refgene_path, sep="\t") %>%
  rename(chr = chrom,
         start = txStart,
         end = txEnd,
         gene_name = name2) %>%
         mutate(group = basename(refgene_path)) %>%
         return
}
longest_refgene_transcripts <- function(ref_table){
  ref_table %>%
    # some genes are duplicated, so multiple groups must be used:
    group_by(gene_name, chr, strand) %>%
    summarise(tss = min(start),
              tes = max(end)) %>%
    # autosomes %>%
    normal_chroms %>%
    # .[which(!duplicated(.$gene_name)), ] %>%
    return
}
# Workforce functions
get_bin_ranges_bed <- function(bed_table_path, bin_width, range){
  path_basename <- basename(bed_table_path)
  bed_table <- fread(bed_table_path, sep="\t")
  bed_table <- bed_table %>% rename(chr = V1, start = V2, end = V3)
  bed_table$midpoints <- with(bed_table, round((start + end)/2))
  bin_start <- seq(from = -range, to = (range - bin_width), by=bin_width)
  expanded_bed_table <- bed_table[rep(seq_len(nrow(bed_table)), each=length(bin_start)),]
  bin_start_rep <- rep(bin_start, nrow(bed_table))
  expanded_bed_table$bin_start <- bin_start_rep
  expanded_bed_table$bin_end  <- expanded_bed_table$bin_start + bin_width
  expanded_bed_table %>%
    transmute(chr   = chr,
              start = midpoints + bin_start,
              end   = midpoints   + bin_end,
              bin_start = bin_start,
              bin_end   = bin_end,
              group = path_basename) %>%
    return
}
get_bin_ranges_refgene <- function(refgene_path, range, bin_width){
  refgene_table <- get_ucsc_refgene_table(refgene_path)
  collapsed_refgene <- longest_refgene_transcripts(refgene_table)
  bin_start <- seq(from=-range, to=(range - bin_width), by=bin_width)
  bin_end  <- bin_start + bin_width- 1

  # expanded_refgene <- collapsed_refgene[rep(seq_len(nrow(collapsed_refgene)), each=length(bin_start)), ]

  rep_indices <- nrow(collapsed_refgene) %>%
                      seq_len %>%
                      rep(., each=length(bin_start))
  expanded_refgene <- collapsed_refgene %>% slice(rep_indices)

  expanded_refgene$bin_start <- rep(bin_start, nrow(collapsed_refgene))
  expanded_refgene$bin_end <- expanded_refgene$bin_start + bin_width

  # Initialize bin boundaries
  expanded_refgene$start <- 0
  expanded_refgene$end   <- 0

  # Add plus bin boundaries
  plus_rows  <- which(expanded_refgene$strand == "+")
  expanded_refgene[plus_rows, ]$start <- expanded_refgene$tss[plus_rows] + expanded_refgene$bin_start[plus_rows]
  expanded_refgene[plus_rows, ]$end   <- expanded_refgene$tss[plus_rows] + expanded_refgene$bin_end[plus_rows]

  # Add minus row bin boundaries -- expecting the column annotated "tes" to be
  #   the biological TSS for genes on the "-" strand (i.e. annotated TSS
  #   position will always be < annotated TES position)
  minus_rows <- which(expanded_refgene$strand == "-")
  expanded_refgene[minus_rows, ]$start <- expanded_refgene$tes[minus_rows] - expanded_refgene$bin_end[minus_rows]
  expanded_refgene[minus_rows, ]$end   <- expanded_refgene$tes[minus_rows] - expanded_refgene$bin_start[minus_rows]

  expanded_refgene$group <- refgene_table$group[1]
  return(expanded_refgene %>% select(gene_name, chr, strand, start, end, bin_start, bin_end, group))
}
get_binned_perc_meth <- function(ranges, meth_table){
  ranges <- data.table(ranges)
  setkey(ranges, chr, start, end)

  capture_overlaps <- foverlaps(meth_table, ranges, type="any", nomatch=0L)
  output <- capture_overlaps %>%
    group_by(bin_start, bin_end) %>%
    summarise(meth = sum(meth), unmeth = sum(unmeth)) %>%
    mutate(perc_meth = meth/(meth + unmeth), depth = meth + unmeth, group = ranges$group[1]) %>%
    arrange(bin_start) %>%
    return
}
# Higher level functions
get_tss_perc_meth <- function(refgene_path, meth_table){
  if(range < bin_width)stop("Stopped: Range is less than bin_width!")
  this_refgene_table <- get_ucsc_refgene_table(refgene_path)
  bin_ranges_refgene <- get_bin_ranges_refgene(this_refgene_table, range = range, bin_width = bin_width)
  get_binned_perc_meth(bin_ranges_refgene, meth_table) %>% return
}
# Plotting functions
plot_percent_meth <- function(binned_perc_meth_table, manual_colors=FALSE){
  binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)
  this_plot <- ggplot(binned_perc_meth_table, aes(bin_mid, perc_meth)) +
    geom_line(aes(color=group)) +
    ylab("Percent Methylation") +
    xlab("Distance from Center (bp)") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0.0001), limits = c(0,1)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          # legend.position="none",
          legend.key = element_blank(),
          axis.line = element_line(),
          panel.background = element_blank(),
          axis.text.x = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=14, vjust=-1.5),
          axis.text.y = element_text(color="black", size=12),
          axis.title.y = element_text(color="black", size=14, vjust=3),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.border=element_blank(),
          axis.ticks=element_line(size=0.6, color="black"))

  # Set up color palette to recycle when out of colors; relies on
  #   xkcdcolors package
  if(manual_colors == TRUE){
    color_palette <- c(name2color("blue"),
                       name2color("light red"))
    color_df <- data.frame(names = levels(df$group), color_palette = color_palette)
    recycle_color_palette <- as.character(color_df$color_palette)
    names(recycle_color_palette) <- color_df$names
    this_plot <- this_plot + scale_color_manual(values=recycle_color_palette)

    print(recycle_color_palette)

  }
  this_plot %>% return
}

#########################
# Gene body quantile plots
get_gene_quantiles  <- function(trancript_line, quantiles = 100){
  tss <- transcript_line$tss
  tes <- transcript_line$tes
  chr <- transcript_line$chr
  strand <- transcript_line$strand
  transcript_line <- longest_transcripts[1,]
  quantile_step <- ((tes - tss)/quantiles) %>%
    round(., digits=0)

  start <- c(tss, tss + cumsum(rep(quantile_step, quantiles-1)))
  end   <- c(start + quantile_step) + 1

  if(strand == "-"){
    start <- rev(start)
    end <- rev(end)
  }

  bin_start <- 0:(length(start)-1)
  bin_end <-   bin_start + 1
  group <- paste(quantiles, "geneBodyQuantiles", sep="")

  quantile_df <- data.frame(chr, start, end, bin_start, bin_end, group)
  quantile_df %>% return
}
longest_transcripts <- this_refgene_path %>% get_ucsc_refgene_table %>% longest_refgene_transcripts
gene_quantile_ranges <- ddply(longest_transcripts, .(gene_name), get_gene_quantiles)
quantile_perc_meth <- get_binned_perc_meth(gene_quantile_ranges, meth_table)
plot_percent_meth(quantile_perc_meth)

##########################################

this_refgene_path = "/Users/tubuliferous/Dropbox/_for_ben/scripts/refGene_hg19"
meth_table <- get_meth_table("/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/surya_data/k562_wgbs.txt")
tss_ranges <- get_bin_ranges_refgene(this_refgene_path, bin_width = 1000, range = 100000)
perc_meth_tss <- get_binned_perc_meth(tss_ranges, meth_table)
plot_percent_meth(perc_meth_tss)

##########################################

perc_meth_tss <- perc_meth_tss %>% mutate(scaled_depth = depth/min(depth))
head(perc_meth_tss)
melted <- melt(perc_meth_tss, id.vars = c("bin_start", "bin_end", "meth", "unmeth", "group"))
ggplot(melted, aes(x = (bin_start + bin_end)/2, y = as.numeric(value))) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", nrow = 3)

dim(perc_meth_tss)
head(meth_table)
##########################################
meth_table_30x <- fread("/Users/tubuliferous/Dropbox/_for_ben/k562_wgbs_30x.cov")
gata <- fread("/Users/tubuliferous/Dropbox/_for_ben/GATA2/SL63653_GATA2_rep2.cov")
# peaks <- fread("/Users/tubuliferous/Desktop/GATA2/MACS_SL63653_SL63655_10_peaks.bed")
peaks_rep2 <- fread("/Users/tubuliferous/Desktop/GATA2/MACS_SL60583_SL60585_10_peaks.bed")
peaks_rep2 <- peaks_rep2 %>% rename(chr = V1, start = V2, end = V3)
gata_rep2_chip <- fread("/Users/tubuliferous/Dropbox/_for_ben/GATA2/SL60583_GATA2_rep1.cov")

get_bin_ranges_bed <- bed_table_path("/Users/tubuliferous/Desktop/GATA2/enhancer_bed_files/E123_13_Active_Enhancer_1.bed", bin_width, range)

setkey(peaks_rep2, chr, start, end)
wgbs_gata2_overlap <- foverlaps(meth_table_30x, peaks_rep2, type = "any", nomatch = 0)
wgbs_gata2_overlap_filt <- wgbs_gata2_overlap %>% select(chr, i.start, i.end, i.V4, i.V5, V6) %>% rename(start = i.start, end = i.end, perc_meth = i.V4, meth =  i.V5, unmeth = V6)
enhancer_ranges <- get_bin_ranges_bed(bed_table_path = "/Users/tubuliferous/Desktop/GATA2/enhancer_bed_files/E123_16_Weak_Enhancer_1.bed", bin_width = 100, range = 10000)
binned_enh_meth_gata2 <- get_binned_perc_meth(ranges = enhancer_ranges, meth_table = wgbs_gata2_overlap_filt)
plot_percent_meth(binned_enh_meth_gata2)

meth_table_30x <- meth_table_30x %>% rename(chr = V1, start = V2, end = V3)
# gata <- rename(chr = V1, start = V2, end = V3)
peaks <- peaks %>% rename(chr = V1, start = V2, end = V3)

head(peaks)

chrs <- paste("chr", 1:22, sep="")
sub_list <- lapply(chrs, function(this_chr){
  this_meth_table <- meth_table_30x %>% filter(chr == this_chr)
  this_peak <- peaks %>% filter(chr == this_chr)
  this_meth_table[which(this_meth_table$start %in% this_peak$start), ] %>% return
})

setkey(peaks, chr, start, end)
capture_overlaps <- foverlaps(meth_table_30x, peaks, type="any", nomatch=0)
capture_overlaps <- capture_overlaps %>% select(chr, i.start, i.end, i.V4, i.V5, V6) %>% rename(start = i.start, end = i.end, perc_meth = i.V4, meth =  i.V5, unmeth = V6)
binned_overlap_wgbs_perc_meth <- get_binned_perc_meth(ranges = enhancer_ranges, meth_table = capture_overlaps)
binned_overlap_wgbs_perc_meth$group = "wgbs"
binned_enh_meth_gata2$group = "gata2"
merged <- rbind(binned_overlap_wgbs_perc_meth, binned_enh_meth_gata2)
plot_percent_meth(merged)
original_gata <- fread("/Users/tubuliferous/Desktop/GATA2/wgbs_GATA2.cov")
original_gata <- arrange(original_gata, V1, V2)
capture_overlaps <- arrange(capture_overlaps, chr, start)
identical(capture_overlaps$i.start, original_gata$V2)
