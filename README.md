# RNA-seq and DAP-seq analysis scripts

This repository contains the R scripts used for the RNA-seq differential expression analysis, the DAP-seq upstream workflow record, and the GO enrichment analysis reported in the manuscript.

## Files

RNA-seq.R: DESeq2 differential expression analysis based on StringTie quantification. 
DAP-seq.R: Record of the DAP-seq upstream workflow, software versions, and parameters. 
Go.R: GO enrichment analysis using clusterProfiler. 

Note: The scripts are analysis templates. Placeholders such as `file`, `sampleNames`, `go`, `set1`, and `"xx"` need to be replaced with your own paths, sample names, annotation table, gene list, and comparison groups before running. The comments inside each script describe the corresponding steps.

## Software versions

- RNA-seq upstream: HISAT2 v2.2.1, SAMtools v1.21, StringTie v3.3.0.
- Differential expression: R v4.3.3, DESeq2 v1.42.1, tximport.
- DAP-seq upstream: fastp v0.23.4, STAR v2.5.3a, RSeQC v2.6, MACS2 v2.2.7.1, bedtools v2.25.0, HOMER v4.10.
- GO enrichment: clusterProfiler, ggplot2, scales.
  
## Reference genome and annotation

RNA-seq: *Oryza sativa* IRGSP-1.0 genome (Ensembl Plants) and MSU/UGA Rice Genome Annotation Project Release 7 GFF3 annotation.
DAP-seq: MSU/UGA Rice Genome Annotation Project Release 7 genome and annotation.
Go:  `Osj_GO_annotation.txt`, compiled from the Gramene / UniProt-GOA rice annotation (GAF format; downloadable from Gramene or UniProt-GOA).

## Data availability

Raw RNA-seq and DAP-seq reads have been deposited in the Genome Sequence Archive under accessions CRA048947 and CRA048932.

## License

The scripts are released under the MIT License.
