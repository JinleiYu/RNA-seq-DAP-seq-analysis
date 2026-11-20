#RNA-seq workflow 
#Raw reads were filtered and aligned to the reference genome with HISAT2 (v2.2.1)
#hisat2 -p 10 --dta -x rice_index -1 <read1.fq.gz> -2 <read2.fq.gz> -S <sample>.sam

#Index built with hisat2
#hisat2-build -p 10 <reference.fa> rice_index

#SAM files were converted to BAM and sorted with SAMtools (v1.21)
#samtools view -@ 10 -f 0x2 -bS -o <sample>.bam <sample>.sam; samtools sort -@ 10 -o <sample>.sorted.bam <sample>.bam

#Transcript assembly and quantification used Stringtie (v3.3.0)
#stringtie -p 20 -B -e -G <annotation.gff3> -o <sample>.gtf <sample>.sorted.bam

#The resulting t_data.ctab files were used for differential expression analysis.

#DESeq2 (v1.42.1) was run in R (v4.3.3) via tximport

#install related packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("factoextra")
BiocManager::install("rdist")

#library(tximport)
library(DESeq2)
library(factoextra)
library(rdist)

# 
setwd(file)
# 
sampleNames<-c("sample", "sample", "sample", "sample", "sample", "sample")
sampleCondition<-c("Condition", "Condition", "Condition", "Condition", "Condition", "Condition")
sampleTable<-data.frame(sampleName = sampleNames,condition = sampleCondition)
##################### 
data <- c(paste("file/", sampleNames, "_t_data.ctab", sep = ""))

read.table(data[1],h=T)->tmp
head(tmp)
tmp[,c("t_name","gene_id")]->tx2gene

#txi<-tximport(data, type = "stringtie", tx2gene = tx2gene)
head(txi)
###DESeq2###
ddsTxi <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ condition)
dds <- DESeq(ddsTxi)
dds

resultsNames(dds)
as.data.frame(counts(dds,normalize=TRUE))->dds.count

sampleTable$sampleName->colnames(dds.count)

head(dds.count)

write.table(dds.count,"file",quote=F,sep="\t")

#
dds.flt <- dds[rowSums(counts(dds))>10,]
vsd <- vst(dds.flt, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
res <- hcut(sampleDists, k =5, stand = TRUE)

# Visualize
fviz_dend(res,
          type="circular",
          rect = TRUE,
          rect_border="cluster",
          rect_lty=2,
          lwd=1,
          rect_fill = T,
          cex = .8,
          color_labels_by_k=T,
          horiz=T
)

#comparison
res <- results(dds, contrast = c("condition", "condition", "condition"))
head(res)
#resordered<-res[order(res$padj),]
write.table(res,"file",quote=F,sep="\t",row.names=T)

#Differentially expressed genes were defined by padj < 0.05 and |log2FoldChange| ≥ 1.


