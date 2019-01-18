#==================run this each new R session==================

#### set working directory ( method works only in RStudio! )
{
  #install.packages("rstudioapi")
  library(rstudioapi)
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
  print( getwd() )
}
####

### Guide:
## outside of R:
# 1) place gff annotation and bam files in folders as desired
#
## in R:
# 0) load in functions
# 1) create bigwig coverage files
# 2) import gff annotation file
# 3) configure and generate coverage graphs

### To Do:
# 1) make an easy way to adjust panel spacing (either identical to ggplot panel or with extra space)
# 2) DONE --> make gene names somehow space or clipping aware
# 3) error check for missing or inconsistent file sets
# 4) improve accuracy inferExperiment 'unknown type' output statistics
# 5) account for CIGAR strings properly

#### 0) load functions in the RCoverage file first ####
# Run entire RCoverage file. 
# Ensure ALL packages are installed.

#### 1) calculate coverage ####

# #basic settings:
# dataDir [NULL] #relative path to the location of the bam files
# outDir  [NULL] #relative path to output folder
# TPM     [T]    #T/F: if TRUE, normalized coverage will also be calculated
# 
# #sequencing experiment settings
# inferExperiments  [T}  #T/F: if TRUE, strand and pair information will be infered automatically
#                        # requires genome annotation file (below) and overrides manual settings
# gffFile    [NULL]      #relative path to gff file
# strandMode [0]         #0/1/2: sets strandedness of reads. 0=unstranded, 1=first read forward, 2=first read reverse
# paired     [F]         #T/F: sets whether reads are considered paired
# 
# #paired-read settings:
# countSinglets           [F]    #T/F: if TRUE, coverage of unpaired reads and filtered reads is included 
# onlyProperPairFragments [T]    #T/F: if TRUE, only reads flagged as proper pairs are considered paired
# filter                  [T]    #T/F: if TRUE, reads are considered paired only if the TLEN field is less than maxTLEN
# maxTLEN                 [600]  #sets maxTLEN cutoff value when filtering 

RCoverage(dataDir="example/bam", outDir="example/coverage", TPM=T, 
          inferExperiments=T, gffFile="example/reference/example.gff", strandMode=0, paired=F,
          countSinglets=F, onlyProperPairFragments=T,
          filter=T, maxTLEN=600)


#### 2) import gff file ####
{ 
  file <- list.files("example/reference","\\.gff$", full.names = T, ignore.case = T)[1]
  message("first identified gff file: ",file)
  
  exonsGR <- import(file)
  exonsGR <- exonsGR[exonsGR$type =="gene"]
  exonsGR <- exonsGR[,c("type","gene_biotype","Name")]
  
  exons <- data.frame(left=start(exonsGR), right=end(exonsGR), strand=strand(exonsGR), locus=exonsGR$Name, type=exonsGR$gene_biotype)
  exons <- exons[exons$strand!="*",]
  exons$strand <- as.character(exons$strand)

  head(exons)
}

#### 3) make plots for publication ####

### plotStacked function:
### Data variables
## target = [NULL] GRanges object defining specific chromosome, strand, and coordinates
#           strand can be either "+", "-", or "*" (for unstranded)
## strains = [NULL] a nested list of group and strain names and files
#           1st level: list of groups that are plotted independently (on different y-scales)
#           2nd level: list of strains that are plotted on the same y-scale, but in different tracks
#           3rd level: character vector of file names for each strain
#               note: files should be bigwig (.bw) files where strand is indicated with names ending with:
#                     "_plus.bw", "_minus.bw" or "_both.bw". Strand will be matched to the file name but 
#                     will fall back to unstranded ("_both.bw") if requested strand is not available.
#           Note: the attributes RCov.strainFontFace and RCov.strandRequest can be set at any level of the list to
#                 override the default settings. Attributes will 'cascade' to apply to lower tiers, but will also be 
#                 overridden by attributes manually set at lower levels.
#
### Appearance variables
## palette = an optional character vector defining the colors to apply to each track (in hex representation).
## height = [3] height of graph in inches; default: 3.
## width = [3] width of graph in inches; default: 3.
#
## labelGroups = [T] either "T" or "F". 
#                If FALSE, strain names will replace the group name on the left of each group.
## binsPerInch = [200] controls number of bins/inch used when averaging coverage. Disable by setting to 0.
#                Note that binning is calculated for the supplied graph width and will not update if the graph is resized
## textScaling = [1] text scaling factor. 
#                Proportionally scales all text, as well as the size of the gene track and group label background

strains <- list(
  "group 1"=list(
    "ex1"=list.files("example/coverage/TPM", pattern="^ex1", ignore.case=T, full.names = T),
    "ex2"=list.files("example/coverage/TPM", pattern="^ex2", ignore.case=T, full.names = T)
  ),
  "group 2"=list(
    "ex2"=list.files("example/coverage/TPM", pattern="^ex2", ignore.case=T, full.names = T)
  )
)

#assign attributes to list elements (optional)
attributes(strains$'group 1'$ex2)$RCov.strainFontFace <- "italic"
attributes(strains$'group 2'$ex2)$RCov.strainFontFace <- "italic"
attributes(strains$'group 2'$ex2)$RCov.strandRequest <- "+" #either of "+", "-", or "*"

#review attribute settings:
rapply(strains,attributes)

#palette=c("grey10","darkred","grey10","darkred")
width <- 4
height <- 2
labelGroups <- T

#### ftsZ ####
rna <- "ftsZ"
exons[exons$locus == rna,]

target <- GRanges(seqname="NC_008533.1", 
                  ranges=IRanges(
                    start = exons[exons$locus == rna,]$left-100, 
                    end = exons[exons$locus == rna,]$right+5000 ), 
                  strand=exons[exons$locus == rna,]$strand )

#quartz for plotting to independent window; mac only:
quartz(width=width,height=height)

plotStacked(target = target, strains = strains, palette = palette, labelGroups=labelGroups, height = height, width = width)

#to save a png:
png(paste0("example/coverage - ",rna,".png"),width=width, height=height, units="in", res=300)
plotStacked(target = target, strains = strains, palette=palette, labelGroups=labelGroups, height = height, width = width)
dev.off()