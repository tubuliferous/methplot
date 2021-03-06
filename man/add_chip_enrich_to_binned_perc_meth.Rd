% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/methplot.r
\name{add_chip_enrich_to_binned_perc_meth}
\alias{add_chip_enrich_to_binned_perc_meth}
\title{Concatenate binned_perc_meth table with aggreegate enrichment table (e.g. from ChIP peaks), addiding dummy columns as required.}
\usage{
add_chip_enrich_to_binned_perc_meth(binned_perc_meth_table, chip_enrich_table)
}
\arguments{
\item{binned_perc_meth_table}{A data.table.}

\item{chip_enrich_table}{A data.table.}
}
\value{
data.table
}
\description{
Concatenate binned_perc_meth table with aggreegate enrichment table (e.g. from ChIP peaks), addiding dummy columns as required.
}
\seealso{
Other workorse functions: \code{\link{genomic_complement}};
  \code{\link{get_aggregate_bed_enrichment_over_ranges}};
  \code{\link{get_bin_ranges_bed}};
  \code{\link{get_bin_ranges_refgene}};
  \code{\link{get_binned_perc_meth_list}};
  \code{\link{get_binned_perc_meth}};
  \code{\link{get_genome_range_overlaps}};
  \code{\link{get_intersect_single_bases}};
  \code{\link{parallel_genomic_complement}}
}

