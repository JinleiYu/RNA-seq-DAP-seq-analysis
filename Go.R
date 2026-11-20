# NOTE: This script is a template. Before running, provide:
#   - `go`: the GO annotation table (e.g., Osj_GO_annotation.txt)
#   - `set1`: the gene list to test (one gene ID per row)
#   - `file`: working directory

setwd(file)
library("clusterProfiler")

#
if (!require("ggplot2")) {
install.packages("ggplot2")
library("ggplot2")
}

if (!require("scales")) {
  install.packages("scales")
  library("scales")
}

#go <- read.delim("annotation file", header = TRUE)
term2gene <- go[, c(2, 1)]  # GO_ID和gene_ID
term2name <- go[, c(2, 3)]  # GO_ID和GO_term

#set1 <- read.table("file")
genes <- set1[, 1]


# 进行GO富集分析
go_result <- enricher(
  gene = genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  minGSSize = 1,
  maxGSSize = 1000
)

#print(dotplot(go_result, showCategory=30))

#
m <- go_result@result
n <- go[,3:4]  
colnames(n) <- c("Description","Aspect")	 
w <- unique(n) 
me <- merge(m, w, all.x=TRUE)
s <- me[me$Aspect=="P", 1:10] 
#
go_result@result <- s

#
print(dotplot(go_result, showCategory=20))

#
write.table(go_result@result, "file", sep = "\t", quote = FALSE, row.names = FALSE)
