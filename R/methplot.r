#' @docType package
#' @name methplot
#' @title Exploratory plots of genomic methylation data.
#' @import dplyr
#' @import data.table

.onAttach <- function(libname, pkgname){
  packageStartupMessage("Welcome to methplot!")
}

# Small helper functions
#' Return subset of autosomes.
#'
#' @family utility functions
#' @param this_table A character.
#' @return data.frame
autosomes <- function(this_table){
  autosomes <- paste("chr", 1:22, sep="") %>% c(., 1:22)
  return(this_table[which(as.character(data.frame(this_table)[,"chr"]) %in% autosomes), ])
}
#' Return subset of autosomes plus X and Y.
#'
#' @family utility functions
#' @param this_table A character.
#' @return data.frame
normal_chroms <- function(this_table){
  auto_x_y <- paste("chr", 1:22, sep="") %>% c(., 1:22) %>% c(., paste("chr", c("X", "Y"), sep="")) %>% c(., "X", "Y")
  return(this_table[which(as.character(data.frame(this_table)[,"chr"]) %in% auto_x_y), ])
}
#' Import and format grch37 gene model file.
#'
#' @family utility functions
#' @param grch37_table_path A character.
#' @return data.frame
get_grch_gene_table <- function(grch37_table_path){
  grch37_table <- readr::read_delim(grch37_table_path, delim="\t", col_names=F)
  names(grch37_table) <- c("id", "source", "type", "start", "end", "dot", "strand", "dot_2", "info")
  chrs <- stringr::str_match(grch37_table$id, "NC_0+(.+)\\.")[,2]
  chrs[which(chrs=="12920")] <- "M"
  chrs[which(chrs=="23")] <- "X"
  chrs[which(chrs=="24")] <- "Y"
  chrs <- paste("chr", chrs, sep="")
  grch37_table <- grch37_table %>%
                  mutate(chr = chrs) %>%
                  filter(type=="gene") %>%
                  mutate(gene_name = str_match(.$info, "Name=(.+?);")[,2])

  out <- grch37_table[grep("pseudo=false", grch37_table$info, invert=TRUE), ]
  return(out)
}
#' Import and format UCSC refgene gene model file.
#'
#' @family utility functions
#' @param refgene_path A character.
#' @return data.frame
get_ucsc_refgene_table <-function(refgene_path){
  fread(refgene_path, sep="\t") %>%
  rename(chr = chrom,
         start = txStart,
         end = txEnd,
         gene_name = name2) %>%
         mutate(group = basename(refgene_path)) %>%
         return
}
#' Import and format methyl sequencing file.
#'
#' @family utility functions
#' @param meth_path A character.
#' @return data.table
get_meth_table <- function(meth_path){
  this_meth_table <- fread(meth_path, sep="\t") %>%
                       rename(chr=V1, start=V2, end=V3, perc_meth=V4, meth=V5, unmeth=V6)
}
#' Get longest possible gene model from refgene.
#'
#' @family utility functions
#' @param ref_table A data.frame or data.table. 
#' @return data.table or data.frame
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

# Workhorse functions
#' Get ranges of bins from BED element centers.
#'
#' @family workorse functions
#' @param bed_table_path A character.
#' @param bin_width A numeric.
#' @param range A numeric.
#' @return data.table
#' @export
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
#' Get bin ranges for gene TSSs 
#'
#' @family workorse functions
#' @param refgene_path A character.
#' @param range A numeric.
#' @param bin_width A numeric.
#' @return data.table
#' @export
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
#' Get percent methylation table by bin.
#'
#' @family workorse functions
#' @param ranges A data.frame.
#' @param meth_table A data.frame. 
#' @return data.table
#' @export
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
#' Get percent methylation table around gene model TSSs.
#'
#' @family high level functions
#' @param refgene_path A character.
#' @param meth_table A data.frame. 
#' @return data.table
#' @export
get_tss_perc_meth <- function(refgene_path, meth_table){
  if(range < bin_width) stop("Stopped: Range is less than bin_width!")
  bin_ranges_refgene <- get_bin_ranges_refgene(refgene_path, range = range, bin_width = bin_width)
  get_binned_perc_meth(bin_ranges_refgene, meth_table) %>% return
}

# Plotting functions
#' Plot binned percent meth table.
#'
#' @family plotting functions
#' @param binned_perc_meth_table A data.frame or data table. 
#' @param manual_colors A logical. 
#' @return data.table
#' @export
plot_percent_meth <- function(binned_perc_meth_table, manual_colors=FALSE){
  binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)
  this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, perc_meth)) +
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

# Unfinished functions -------------------------
#' Get gene quantiles for single gene.
#'
#' @family plotting functions
#' @param trancript_line A data.frame or data table . 
#' @param quantiles A numeric.  
#' @return data.table
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
