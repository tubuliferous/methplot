library(devtools)
devtools::install_github("tubuliferous/methplot")
library(methplot)
meth_table <- get_meth_table("/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/surya/bedtool_test/Pol2_rep1.txt")
bin_ranges_refgene <- get_bin_ranges_refgene("/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/surya/bedtool_test/refGene_hg19", 10000, 1000)
tss_perc_meth <- get_tss_perc_meth("/Users/tubuliferous/Dropbox/The Lives (Projects) of Others/Surya/surya/bedtool_test/refGene_hg19", meth_table, range = 10000, bin_width = 1000)

plot_percent_meth(tss_perc_meth)

