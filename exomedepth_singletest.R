######################################
####  Workflow for 1 sample
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
#######################################
ExomeCount.dafr <- as(my.counts[, colnames(my.counts)],
                      'data.frame')

ExomeCount.dafr$space <- gsub(as.character(ExomeCount.dafr$space),
                              pattern = 'chr',
                              replacement = ' ')

print(head(ExomeCount.dafr))

#######################################
####  define test sample
####  This is the bam file that we want to
####  investigate for for CNVs in.
####  Make sure matches column heading from 
####  the previous print command (ExomeCount.dafr).
#######################################
my.test <- my.counts$File1.bam

#######################################
####  define refernce sample
####  the rest of the files to compare against.
####  Again make sure matches column heading from 
####  the print command (ExomeCount.dafr).
#######################################
my.ref.samples <-c(
  "File2.bam",
  "File3.bam", 
  etc.)

my.reference.set<-as.matrix(ExomeCount.dafr[,my.ref.samples])

#######################################
####  optimise choice of aggregate refernce set
#######################################
my.choice<-select.reference.set(
  test.counts=my.test,
  reference.counts=my.reference.set,
  bin.length=(ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
  n.bins.reduced=10000)

#######################################
####  Check which of the reference samples 
####  fit the test sample best.
#######################################
print(my.choice[[1]])

#######################################
####  Construct the refernce set
#######################################
my.matrix <-as.matrix( ExomeCount.dafr[, my.choice$reference.choice, 
                                       drop = FALSE])

my.reference.selected<-apply(X=my.matrix,
                             MAR=1,
                             FUN=sum)

#######################################
####  Fit the beta-binomial model on a data frame
#######################################
all.exons <-new('ExomeDepth',
                test=my.test, 
                reference=my.reference.selected,
                formula='cbind(test, reference)~1')

#######################################
####  Call the CNV by running the 
####  underlying hidden Markov model
#######################################
all.exons<-CallCNVs(x=all.exons,
                    transition.probability=10^-4,
                    chromosome=ExomeCount.dafr$space,
                    start=ExomeCount.dafr$start,
                    end=ExomeCount.dafr$end,
                    name=ExomeCount.dafr$names)

#######################################
####  Check the top 6 CNV calls
#######################################
head(all.exons@CNV.calls)

#######################################
####  Load the set of common CNVs identified 
####  in the internat dataset Conrad et al, Nature 2010
####  similar to dbSNP but for CNVs
#######################################
data(Conrad.hg19) 

head(Conrad.hg19.common.CNVs)

#######################################
####  Annotate CNV calls
####  Check for any CNVs in your list that 
####  may be known and mark as known CNVs
#######################################
all.exons<-AnnotateExtra(x=all.exons, 
                         reference.annotation=Conrad.hg19.common.CNVs, 
                         min.overlap=0.5, 
                         column.name='Conrad.hg19') 

print(head(all.exons@CNV.calls))

#######################################
####  Process the Conrad et al data in the GRanges format
#######################################
exons.hg19.GRanges <- GenomicRanges::GRanges(
                        seqnames=exons.hg19$chromosome,
                        IRanges::IRanges(start=exons.hg19$start,
                                         end=exons.hg19$end),
                        names=exons.hg19$name)

#######################################
####  Annotate CNVs with location, gene names, etc.
#######################################
all.exons <- AnnotateExtra(x=all.exons, 
                           reference.annotation=exons.hg19.GRanges, 
                           min.overlap=0.0001, 
                           column.name='exons.hg19')

all.exons@CNV.calls[3:6,]

#######################################
##### Output CNV calls to .csv
####  using your specified filename
####  to the working directory
#######################################
output.file <- 'exome_calls_File1.csv' 


write.csv(file=output.file,
          x=all.exons@CNV.calls,
          row.names=FALSE)


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
