
######################################
####  Workflow for batch analysis
######################################

######################################
####  Load library or install if needed
#######################################
# install.packages("ExomeDepth")
library("ExomeDepth")

#######################################
####  Use the internal dataset 
####  included in ExomeDepth software
#######################################
data(exons.hg19)

#######################################
####  Check the data
#######################################
print(head(exons.hg19))

#######################################
####  Load your data as bam files
####  including test sample & control samples
#######################################
my_bam_files <- c(File1.bam, 
                  File2.bam,
                  File3.bam, 
                  etc.)


#######################################
####  Create count data for autosomal
#######################################
my.counts <- getBamCounts(
  bed.frame=exons.hg19,
  bam.files=my_bam_files,
  include.chr=TRUE,
  referenceFasta="/home/ref/hg19/hg19_fullGenome_sort_numeric.fa")


#######################################
####  Convert counts
####  (GRanges class object) into a data frame
#######################################
ExomeCount.dafr <- as(my.counts[, colnames(my.counts)],
                      'data.frame')

#######################################
####  Defining the test and reference samples 
####  which will alternate
#######################################
ExomeCount.dafr$space <- gsub(as.character(ExomeCount.dafr$space), 
                              pattern = 'chr', 
                              replacement = ' ')

source("ExomeDepthBatch.R")
