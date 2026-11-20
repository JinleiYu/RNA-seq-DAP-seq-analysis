#DAP-seq workflow.
#Raw reads were filtered with fastp (v0.23.4)
#-l 39 --dont_eval_duplication --qualified_quality_phred 20 --unqualified_percent_limit 8 --adapter_sequence <adapter1> --adapter_sequence_r2 <adapter2> --fix_mgi_id

#Clean reads were aligned with STAR (v2.5.3a).
#--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --limitBAMsortRAM 64000000000

#Read distribution was assessed with RSeQC (v2.6; default parameters). 

#Peaks were called with MACS2 (v2.2.7.1) against the corresponding input DNA libraries 
#macs2 callpeak -t {ip_bams} -c {input_bams} -g {GENOME_SIZE} -n {sample} --outdir ... -f BAM --bw 300 --mfold 5 50 --slocal=1000 --llocal=10000 --keep-dup=1 -q 0.05 --fe-cutoff 2.0 --buffer-size 100000

#Peak annotation and distribution used bedtools (v2.25.0)
#intersectBed -a <peak.bed> -b <genome.bed> -wa –wb

#Motif analysis used HOMER (v4.10)
#findMotifsGenome.pl <peak_file> <genome.fa> <out_dir>
  
  