% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/methplot.r
\name{get_bin_ranges_refgene}
\alias{get_bin_ranges_refgene}
\title{Get bin ranges for gene TSSs}
\usage{
get_bin_ranges_refgene(refgene_path, range, bin_width)
}
\arguments{
\item{refgene_path}{A character.}

\item{range}{A numeric.}

\item{bin_width}{A numeric.}
}
\value{
data.table
}
\description{
Get bin ranges for gene TSSs
}
\seealso{
Other workorse functions: \code{\link{add_chip_enrich_to_binned_perc_meth}};
  \code{\link{genomic_complement}};
  \code{\link{get_aggregate_bed_enrichment_over_ranges}};
  \code{\link{get_bin_ranges_bed}};
  \code{\link{get_binned_perc_meth_list}};
  \code{\link{get_binned_perc_meth}};
  \code{\link{get_genome_range_overlaps}};
  \code{\link{get_intersect_single_bases}};
  \code{\link{parallel_genomic_complement}}
}

