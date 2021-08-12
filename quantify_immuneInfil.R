# Deconvoluting immune cells from tumor fraction in tumor samples 
# https://github.com/icbi-lab/immunedeconv
# setwd("/Volumes/target_nbl_ngs/KP/papThyCarc_bulkRNASeq")

install.packages("BisqueRNA")
library(Biobase)
library(BisqueRNA)

# using counts data
cts.df <- cts.df %>%
  column_to_rownames(var = 'geneSymbol')
# prepare bulk data
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(cts.df))


# get markers
markers.to.get <- read.delim('cell_type_proportions/marker.txt', header = T)
de.genes <- read.delim('DE.geneslist.txt', header = T)
markers <- de.genes %>%
  dplyr::select(1,3)

markers <- merge(markers, markers.to.get, by.x = 'Description', by.y = 'geneSymbol')
markers <- markers[,c(1,3,2)]
names(markers)[1] <- 'gene'
markers <- markers[1:3,]

res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, weighted=F, min_gene = 1)

marker.based.estimates <- res$bulk.props
knitr::kable(marker.based.estimates2, digits = 2)

scaled.true.props <- t(scale(t(true.props)))[rownames(marker.based.estimates),]


# Bisque did not work
# trying different method
remotes::install_github("icbi-lab/immunedeconv")
# to run cibersort: you need to register on the cibersort website, obtain a license, and download the CIBERSORT source code.


library(immunedeconv)

# preparing data
load('samples_TPM_matrix.RData')

tpm.df <- tpm.df %>%
  dplyr::select(-1)
tpm.df <- tpm.df[tpm.df$hgnc_symbol != '',]
tpm.df <- tpm.df %>%
  group_by(hgnc_symbol, sample) %>%
  summarize(max.TPM = max(TPM))
tpm.df <- tpm.df %>%
  distinct() %>%
  spread(key = 'sample', value = 'max.TPM')


tpm.df <- tpm.df %>%
  column_to_rownames(var = 'hgnc_symbol')

# using Timer
#res_quantiseq <- immunedeconv::deconvolute(tpm.df, "timer", indications = c('thca', 'thca'))

# using quanTIseq
res_quantiseq <- immunedeconv::deconvolute(tpm.df, "quantiseq")
names(res_quantiseq) <- gsub('.genes.results','', names(res_quantiseq))

nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

library(randomcoloR)
n <- 50
palette <- distinctColorPalette(n)

res_epic <- deconvolute_epic(tpm.df, tumor = TRUE, scale_mrna = FALSE)
res_quant <- deconvolute_quantiseq(tpm.df, arrays = FALSE ,tumor = TRUE, scale_mrna = FALSE)
res_epic %>%
  data.frame() %>%
  rownames_to_column(var = 'cell_type') %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  labs(title = 'Cell fractions calculated using epic') +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = palette) +
  #scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(res_epic))) +
  theme(plot.title = element_text(hjust = 0.5))


res_quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  labs(title = 'Cell fractions calculated using quanTIseq') +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = palette) +
  #scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq))) +
  theme(plot.title = element_text(hjust = 0.5))

# res_mcp_counter %>%
#   gather(sample, score, -cell_type) %>%
#   ggplot(aes(x=sample, y=score, color=cell_type)) +
#   geom_point(size=4) +
#   facet_wrap(~cell_type, scales="free_x", ncol=3) +
#   scale_color_brewer(palette="Paired", guide=FALSE) +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, filename = 'deconvoluted_celltypes_quantiseq.pdf', width = 10, height = 5)

deconvolute(tpm.df, "timer", indications = c('thca', 'thca')) %>%
  map_result_to_celltypes(c("uncharacterized cell"), "timer")



#write.table(tpm.df, file = 'tpm.df_formatted_2samples.txt', sep = '\t', quote = F, col.names = T, row.names = T)
# running standalone EPIC tool: http://epic.gfellerlab.org/ using counts data


TIC <- read.delim('EPIC_results_TumorInfiltratingcells.txt', skip = 4)
BCI <- read.delim('EPIC_results_bloodCirculatingImmuneCells.txt', skip = 4)

TIC %>%
  gather(key = 'cell_type', value = 'fraction', -sampleID) %>%
  ggplot(., aes(x = sampleID, y = fraction, fill = cell_type)) +
  geom_bar(stat='identity') +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = mycolors)


BCI %>%
  gather(key = 'cell_type', value = 'fraction', -sampleID) %>%
  ggplot(., aes(x = sampleID, y = fraction, fill = cell_type)) +
  geom_bar(stat='identity') +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = mycolors)









