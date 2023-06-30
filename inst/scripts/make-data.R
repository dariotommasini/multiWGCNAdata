# Autism microarray dataset (Voineagu et al. Nature 2011)
# Expression data.frame of log2 transformed and quantile normalized signal was 
# downloaded from the Gene Expression Omnibus (GEO, accession number GSE28521, 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM706469). 
gset <- getGEO("GSE28521", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6883", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }
autism_data = data.frame(X = rownames(ex), ex[,22:79])
save(autism_data, file = "autism_data.rda")

# The metadata was extracted from the GEO entry names.
autism_metadata = data.frame(Sample = paste0("GSM7064", 12:69), 
                      Disease = c(rep("autism", 16), rep("controls", 16), rep("autism", 13), rep("controls", 13)), 
                      Region = c(rep("FC", 32), rep("TC", 26)))
save(autism_metadata, file = "autism_metadata.rda")

# Astrocyte Ribotag dataset (Itoh et al. PNAS 2018)
# The astrocyte data from Itoh et al. PNAS. 2018, was 
# extracted from GEO (GSE100329) using sra-tools, aligned to mm10 using 
# STAR, and gene counts were retrieved using HTSEQ:

# for accession in SRR_list
#   do
# #Fastq
# fasterq-dump --split-files $accession;
# rm -r $accession
# 
# #STAR
# mkdir ALIGNMENT
# /opt/STAR/bin/Linux_x86_64/STAR --runThreadN 12 --genomeDir /saveddata/shared/Reference_Genomes/STAR_indices_overhang75_mm10 --readFilesIn ${accession}_1.fastq ${accession}_2.fastq  --sjdbGTFfile /saveddata/shared/Reference_Genomes/STAR_indices_overhang75_mm10/genes.gtf --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 10 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFilterMismatchNmax 3 --twopassMode Basic --outFileNamePrefix ALIGNMENT/${accession} --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --quantMode TranscriptomeSAM	
# rm *fastq
# rm ALIGNMENT/*Aligned.toTranscriptome.out.bam
# done
# 
# #HTSEQ
# mkdir HTSEQ/
#   mv ALIGNMENT/*Aligned.sortedByCoord.out.bam ./
#   
#   #parallelize this slow step
#   ls *Aligned.sortedByCoord.out.bam | xargs -I{} -n 1 -P 12 samtools index {} {}.bai
# ls *Aligned.sortedByCoord.out.bam | xargs -I{} -n 1 -P 12 sh -c "~/.local/bin/htseq-count -f bam -m intersection-strict -r pos -i gene_name -s reverse {} /saveddata/shared/Reference_Genomes/STAR_indices_overhang75_mm10/genes.gtf > {}.txt"

# Read in HTSEQ files
Countlist = list.files(pattern = '*.txt')
Countfiles = lapply(Countlist, read.table, header=FALSE, sep="\t",fill=TRUE)
Count <- Reduce(function(x, y) {
  merge(x, y, all=TRUE, by="V1")
}, Countfiles)
rownames(Count)=Count$V1
Count$V1=NULL
colnames(Count)<-Countlist
colnames(Count) <- gsub("_Aligned.sortedByCoord.out.bam.LastLinesRem.txt","",colnames(Count))
# filter count by genes listed in length file so that they are matched
Count_filt <- Count[match(rownames(length),rownames(Count)),]
library(edgeR)
RPKM=rpkm(Count_filt, length$Length)
logRPKM=log2(RPKM+1)
colnames(logRPKM) = c("X", c(paste0("EAE_", 1:20), paste0("healthy_", 1:16)))
datExpr = logRPKM
save(datExpr, file = "astrocyte_data.rda")

# The metadata for the astrocyte Ribotag dataset (Itoh et al. PNAS. 2018) was
# generated using the names of the GEO entries. 
# Astrocyte co-expression networks were generated as described in Tommasini and 
# Fogel, BMC Bioinformatics, 2023. and as outlined in the multiWGCNA vignette 
# astrocyte_map_v2.Rmd. 

metadata = data.frame(Sample = c(paste0("EAE_", 1:20), paste0("healthy_", 1:16)), 
                      Disease = c(rep("EAE", 20), rep("WT", 16)), 
                      Region = c(rep(c("Cbl", "Ctx", "Hippo", "Sc"), each = 5), rep(c("Cbl", "Ctx", "Hippo", "Sc"), each = 4)))
save(metadata, file = "astrocyte_metadata.rda")

# SummarizedExperiment objects for each datset were generated from the 
# expression data.frame and metadata like this: 
mat = datExpr[,-1]
rownames(mat) = datExpr$X
colData = metadata
rownames(colData) = metadata$Sample
rowRanges = datExpr["X"]
rownames(rowRanges) = datExpr$X
colnames(rowRanges) = "gene"
se = SummarizedExperiment(assays=list(counts=as.matrix(mat)),
                          rowData=rowRanges, colData=colData)
save(se, file = "astrocyte_se.rda")

# Tau pathology (rTg4510) Mouse Summarized Experiment
# The tau pathology dataset from Castanho et al. Cell Rep. 2020, extracted from 
# GSE125957 supplementary file "GSE125957_processed_data.csv.gz", and processed as 
# described in methods (https://git.exeter.ac.uk/ic322/ad-mice-rna-seq-cell-reports).
# Thus, values are identical to Castanho et al. 2020. Values correspond to 
# DESeq2 rlog-transformed counts. 

rlog_counts <- t(datExpr)
rlog_counts <- rlog(rlog_counts) # we used rlog() 1) to be consistent with the DESeq2 analysis and 2) because from DESeq2 QC rlog transformation showed more robust results
rlog_counts <- t(rlog_counts)
save(datExpr, rlog_counts, datTraits, file = "exp_data_for_WGCNA.RData")
