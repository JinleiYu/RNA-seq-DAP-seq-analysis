# RNA-seq-DAP-seq-analysis
RNA-seq workflow. 
Raw reads were filtered and aligned to the reference genome with HISAT2 (v2.2.1)
hisat2 -p 10 --dta -x rice_index -1 <read1.fq.gz> -2 <read2.fq.gz> -S <sample>.sam

Index built with hisat2
hisat2-build -p 10 <reference.fa> rice_index

SAM files were converted to BAM and sorted with SAMtools (v1.21)
samtools view -@ 10 -f 0x2 -bS -o <sample>.bam <sample>.sam; samtools sort -@ 10 -o <sample>.sorted.bam <sample>.bam

Transcript assembly and quantification used Stringtie (v3.3.0)
stringtie -p 20 -B -e -G <annotation.gff3> -o <sample>.gtf <sample>.sorted.bam

The resulting t_data.ctab files were used for differential expression analysis.

DESeq2 (v1.42.1) was run in R (v4.3.3) via tximport

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

#setwd(file) 
#sampleNames<-c("sample", "sample", "sample", "sample", "sample", "sample")
sampleCondition<-c("xx", "xx", "xx", "xx", "xx", "xx")
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

#dds.flt <- dds[rowSums(counts(dds))>10,]
vsd <- vst(dds.flt, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
res <- hcut(sampleDists, k =5, stand = TRUE)

#Visualize
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
res <- results(dds, contrast = c("condition", "xx", "xx"))
head(res)
#resordered<-res[order(res$padj),]
write.table(res,"file",quote=F,sep="\t",row.names=T)

Differentially expressed genes were defined by padj < 0.05 and |log2FoldChange| ≥ 1.

GO

setwd(xx)
library("clusterProfiler")

#if (!require("ggplot2")) {
  install.packages("ggplot2")
  library("ggplot2")
}

if (!require("scales")) {
  install.packages("scales")
  library("scales")
}

#go <- read.delim("file", header = TRUE)
term2gene <- go[, c(2, 1)]  # GO_ID和gene_ID
term2name <- go[, c(2, 3)]  # GO_ID和GO_term

#set1 <- read.table("file")
genes <- set1[, 1]

#go_result <- enricher(
  gene = genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  minGSSize = 1,
  maxGSSize = 1000
)

#print(dotplot(go_result, showCategory=30))

#m <- go_result@result
n <- go[,3:4]  
#colnames(n) <- c("Description","Aspect")	 
w <- unique(n) 
#me <- merge(m, w, all.x=TRUE)

s <- me[me$Aspect=="P", 1:10] 
#go_result@result <- s

#print(dotplot(go_result, showCategory=20))

#write.table(go_result@result, "file", sep = "\t", quote = FALSE, row.names = FALSE)

# DAP-seq workflow.
Raw reads were filtered with fastp (v0.23.4)
-l 39 --dont_eval_duplication --qualified_quality_phred 20 --unqualified_percent_limit 8 --adapter_sequence <adapter1> --adapter_sequence_r2 <adapter2> --fix_mgi_id

Clean reads were aligned with STAR (v2.5.3a).
--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --limitBAMsortRAM 64000000000

Read distribution was assessed with RSeQC (v2.6; default parameters). 

Peaks were called with MACS2 (v2.2.7.1) against the corresponding input DNA libraries 
macs2 callpeak -t {ip_bams} -c {input_bams} -g {GENOME_SIZE} -n {sample} --outdir ... -f BAM --bw 300 --mfold 5 50 --slocal=1000 --llocal=10000 --keep-dup=1 -q 0.05 --fe-cutoff 2.0 --buffer-size 100000

Peak annotation and distribution used bedtools (v2.25.0)
intersectBed -a <peak.bed> -b <genome.bed> -wa –wb

Motif analysis used HOMER (v4.10)
findMotifsGenome.pl <peak_file> <genome.fa> <out_dir>

