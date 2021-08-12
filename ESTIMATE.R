# using ESTIMATE for predicting tumor purity, and the presence of 
# infiltrating stromal/immune cells in tumor tissues using gene expression data. 
# setwd("/Volumes/target_nbl_ngs/KP/papThyCarc_bulkRNASeq")

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

library(estimate)
help(package="estimate")

# input count matrix should be vst or rlog transformed (log2 transformed)
tumorExpr <- paste0("/Volumes/target_nbl_ngs/KP/papThyCarc_bulkRNASeq/", "vstTransformed_count_matrix_ESTIMATE.txt")

#lncRNAgct <- tempfile(pattern="estimate", fileext=".gct")
filterCommonGenes(input.f=tumorExpr, output.f="tum_genes.gct", id="GeneSymbol")
estimateScore("tum_genes.gct", "tum_estimate_score.gct", platform="illumina")                    

plotPurity(scores="tum_estimate_score.gct", samples="s516", platform="illumina") # !!!

# read in estimate scores
estimate.score <- read.delim('tum_estimate_score.gct', skip = 2)

# reshape data
estimate.score <- estimate.score %>%
  gather(key = 'samples', value = 'score', -c(NAME, Description))

# fix sample names
estimate.score$samples <- gsub('^X','', estimate.score$samples)
estimate.score$samples <- gsub('.genes.results','', estimate.score$samples)

# plot
q <- ggplot(estimate.score, aes(samples, score, fill = NAME)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw()
ggsave(q, filename = 'ESTIMATE_score_comparison_tumorSamples.pdf', width = 10, height = 7)





