library(data.table)
library(weaverfunctions)
library(ggplot)
library(dplyr)
library(devtools)
# install_github("tubuliferous/methplot")
library(methplot)
library(colorout)

ls("package:methplot")

get_ucsc_refgene_table <-function(refgene_path){
  data.table::fread(refgene_path, sep="\t") %>%
  dplyr::rename(chr = chrom,
         start = txStart,
         end = txEnd,
         gene_name = name2) %>%
         dplyr::mutate(group = basename(refgene_path)) %>%
         return
}
longest_refgene_transcripts <- function(ref_table){
  ref_table %>%
    # some genes are duplicated, so multiple groups must be used:
    dplyr::group_by(gene_name, chr, strand) %>%
    dplyr::summarise(tss = min(start),
              tes = max(end)) %>%
    # autosomes %>%
    normal_chroms %>%
    # .[which(!duplicated(.$gene_name)), ] %>%
    return
}
normal_chroms <- function(this_table){
  auto_x_y <- paste("chr", 1:22, sep="") %>% c(., 1:22) %>% c(., paste("chr", c("X", "Y"), sep="")) %>% c(., "X", "Y")
  return(this_table[which(as.character(data.frame(this_table)[,"chr"]) %in% auto_x_y), ])
}
get_bin_ranges_flanking_gene_body <- function(refgene_path, range, bin_width){
  refgene_table <- get_ucsc_refgene_table(refgene_path)
  collapsed_refgene <- longest_refgene_transcripts(refgene_table)

	# Upstream bins -----------------------------------------------------------

  bin_start <- seq(from=-range, to=0 - bin_width, by=bin_width)
  bin_end  <- bin_start + bin_width - 1

  rep_indices <- nrow(collapsed_refgene) %>%
                      seq_len %>%
                      rep(., each=length(bin_start))
  expanded_refgene_upstream <- collapsed_refgene %>% dplyr::slice(rep_indices)

  expanded_refgene_upstream$bin_start <- rep(bin_start, nrow(collapsed_refgene))
  expanded_refgene_upstream$bin_end <- expanded_refgene_upstream$bin_start + bin_width

  # Initialize bin boundaries
  expanded_refgene_upstream$start <- 0
  expanded_refgene_upstream$end   <- 0

  # Add plus bin boundaries
  plus_rows  <- which(expanded_refgene_upstream$strand == "+")
  expanded_refgene_upstream[plus_rows, ]$start <- expanded_refgene_upstream$tss[plus_rows] + expanded_refgene_upstream$bin_start[plus_rows]
  expanded_refgene_upstream[plus_rows, ]$end   <- expanded_refgene_upstream$tss[plus_rows] + expanded_refgene_upstream$bin_end[plus_rows]

  # Add minus row bin boundaries -- expecting the column annotated "tes" to be
  #   the biological TSS for genes on the "-" strand (i.e. annotated TSS
  #   position will always be < annotated TES position)
  minus_rows <- which(expanded_refgene_upstream$strand == "-")
  expanded_refgene_upstream[minus_rows, ]$start <- expanded_refgene_upstream$tes[minus_rows] - expanded_refgene_upstream$bin_end[minus_rows]
  expanded_refgene_upstream[minus_rows, ]$end   <- expanded_refgene_upstream$tes[minus_rows] - expanded_refgene_upstream$bin_start[minus_rows]

  expanded_refgene_upstream$group <- paste(refgene_table$group[1], "_upstream", sep = "")


  # Downstream bins -----------------------------------------------------------

  bin_start <- seq(from = 0, to = range - bin_width, by = bin_width)
  bin_end  <- bin_start + bin_width - 1

  rep_indices <- nrow(collapsed_refgene) %>%
                      seq_len %>%
                      rep(., each=length(bin_start))
  expanded_refgene_downstream <- collapsed_refgene %>% dplyr::slice(rep_indices)

  expanded_refgene_downstream$bin_start <- rep(bin_start, nrow(collapsed_refgene))
  expanded_refgene_downstream$bin_end <- expanded_refgene_downstream$bin_start + bin_width

  # Initialize bin boundaries
  expanded_refgene_downstream$start <- 0
  expanded_refgene_downstream$end   <- 0

  # Add plus bin boundaries
  plus_rows  <- which(expanded_refgene_downstream$strand == "+")
  expanded_refgene_downstream[plus_rows, ]$start <- expanded_refgene_downstream$tes[plus_rows] + expanded_refgene_downstream$bin_start[plus_rows]
  expanded_refgene_downstream[plus_rows, ]$end   <- expanded_refgene_downstream$tes[plus_rows] + expanded_refgene_downstream$bin_end[plus_rows]

  # Add minus row bin boundaries -- expecting the column annotated "tes" to be
  #   the biological TSS for genes on the "-" strand (i.e. annotated TSS
  #   position will always be < annotated TES position)
  minus_rows <- which(expanded_refgene_downstream$strand == "-")
  expanded_refgene_downstream[minus_rows, ]$start <- expanded_refgene_downstream$tss[minus_rows] - expanded_refgene_downstream$bin_end[minus_rows]
  expanded_refgene_downstream[minus_rows, ]$end   <- expanded_refgene_downstream$tss[minus_rows] - expanded_refgene_downstream$bin_start[minus_rows]

  expanded_refgene_downstream$group <- paste(refgene_table$group[1], "_downstream", sep = "")

  expanded_refgene_merged <- rbind(expanded_refgene_upstream, expanded_refgene_downstream)

  return(expanded_refgene_merged %>% select(gene_name, chr, strand, start, end, bin_start, bin_end, group))
}

get_gene_quantiles  <- function(transcript_line, quantiles = 100, quantile_pseudo_width = 100){
  tss <- transcript_line$tss
  tes <- transcript_line$tes
  chr <- transcript_line$chr
  gene_name <- transcript_line$gene_name
  strand <- transcript_line$strand
  quantile_step <- ((tes - tss)/quantiles) %>%
    round(., digits=0)

  start <- c(tss, tss + cumsum(rep(quantile_step, quantiles-1)))
  end   <- c(start + quantile_step) + 1

  if(strand == "-"){
    start <- rev(start)
    end <- rev(end)
  }

  bin_start <- 0:(length(start)-1) * quantile_pseudo_width 
  bin_end <-   (bin_start + quantile_pseudo_width)
  group <- paste(quantiles, "geneBodyQuantiles", sep="")

  quantile_df <- data.frame(gene_name, chr, strand, start, end, bin_start, bin_end, group)
  quantile_df %>% return
}


refgene_path = "/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/surya/bedtool_test/refGene_hg19"
flanks <- get_bin_ranges_flanking_gene_body(refgene_path, 5000, 100)
refgene_table <- get_ucsc_refgene_table(refgene_path)
collapsed_refgene <- longest_refgene_transcripts(refgene_table)
meth_table <- get_meth_table("/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/k562_wgbs.cov")

bin_ranges_refgene <- get_bin_ranges_refgene("/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/surya/bedtool_test/refGene_hg19", 10000, 1000)
tss_perc_meth <- get_tss_perc_meth("/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/surya/bedtool_test/refGene_hg19", meth_table, range = 10000, bin_width = 1000)


quantiles <- collapsed_refgene %>% 
	rowwise() %>% 
	do(get_gene_quantiles(.))

flanks_downstream_row_idx <- which(flanks$group == "refGene_hg19_downstream")
flanks_downstream_row_idx %>% head

flanks$bin_start[flanks_downstream_row_idx] <- flanks$bin_start[flanks_downstream_row_idx] + max(quantiles$bin_end)
flanks$bin_end[flanks_downstream_row_idx] <- flanks$bin_end[flanks_downstream_row_idx] + max(quantiles$bin_end)

bound <- rbind(flanks, quantiles)

ls("package:methplot")

binned <- get_binned_perc_meth(bound, meth_table)

head(binned)

quartz();plot_percent_meth_with_depth(binned)



