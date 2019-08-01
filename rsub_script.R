#rsubread

library(Rsubread)
library(limma)
library(edgeR)

# Index building. 

# Build an index for human chromosome 1. This typically takes ~3 minutes. Index
# files with basename 'chr1' will be generated in your current working directory
buildindex(basename="chr1", reference="hg19_chr1.fa")

# Alignment.

# Perform read alignment for all four libraries and report uniquely mapped reads
# only. This typically takes ~5 minutes. BAM files containing the mapping
# results will begenerated in your current working directory.
targets <- readTargets()
align(index="chr1", readfile1=targets$InputFile, output_file=targets$OutputFile)

# Read  summarization.

# Summarize mapped reads to NCBI RefSeq genes. This will onlytake a few seconds.
# Note that the 'featureCountsfunction' contains built-in RefSeq annotations for
# human and mouse genes. featureCountsreturns an R 'List' object, which includes
# raw readcount for each gene in each library and also annotation information
# such as gene identifiersand gene lengths

fc <- featureCounts(files=targets$OutputFile, annot.inbuilt="hg19")

# # Create a DGEListobject
# x <- DGEList(counts=fc$counts, genes=fc$annotation[, c("GeneID", "Length")])
# 
# # Filtering
# 
# # Only keep in the analysis those genes which had >10 reads per million mapped
# # reads in at least two libraries.
# 
# isexpr <- rowSums(cpm(x) > 10) >= 2
# x <- x[isexpr,]
# 
# #Design matrix.Create a design matrix: