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
          inferExperiments=T, gffFile="example/reference/mini-NC_000913.3.gff", strandMode=0, paired=F,
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



#### ftsZ ####

#Lets start simple.
# We will define a strain list with one group of two tracks. 
# Tracks grouped together are scaled together, so organization is important.
# Each track is a list of bigwig coverage files for a set of replicates.
# Each group is a list of individual track lists.
# The strain list is just a list of individual group lists.
#
# Note: I designed this with TPM-normalized coverage tracks in mind. Keep in mind that the coverage
#       of all the replicates is averaged... this makes the most sense to me with normalized data.
#       If you want to use counts instead, then simply use the count bigwig files below.

strains <- list(
  "mRNA"=list(
    "WT"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T),
    "pnp"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  )
)

#Set the dimensions of the output graph
# Due to some limitations, various calculations for the final graph are not updated if you re-size it later
# It is best to define the size you want and have everything done correctly to the intended size.
width <- 4
height <- 2

#Next define the region and srand to be graphed
# For convenience, here the ompX locus is used to set the left and right window coordinates and the strand
# this can be done manually of course as well
rna <- "ompX"
exons[exons$locus == rna,]

target <- GRanges(seqname="NC_000913.3", 
                  ranges=IRanges(
                    start = exons[exons$locus == rna,]$left-3300, 
                    end = exons[exons$locus == rna,]$right+500 ), 
                  strand=exons[exons$locus == rna,]$strand )

#for plotting to independent window with a defined size, first make the device
dev.new(width=width,height=height, unit="in", noRStudioGD = T)

#then run plotStacked to draw the graph
plotStacked(target = target, strains = strains, labelGroups=T, height = height, width = width)

#to save a png:
# there can be PDF rendering issues on some platforms (macs, at least) if you choose PDF output
# a decently high-resolution png works well however
png(paste0("example/coverage - ",rna,".png"),width=width, height=height, units="in", res=300)
plotStacked(target = target, strains = strains, palette=palette, labelGroups=labelGroups, height = height, width = width)
dev.off()


#What if we want something more complicated?
# How about three groups of two tracks with manually specified strands
strains <- list(
  "plus"=list(
    "WT"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T),
    "pnp"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  ),
  "minus"=list(
    "WT"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T),
    "pnp"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  ),
  "all"=list(
    "WT"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T),
    "pnp"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  )
)

#Assign attributes to list elements (optional)
# You can set two attributes at any level of the strain list. These will be propogated downwards to
# all sub elements. However, an attribute at the group level will be overwriten by any attribute
# manually set at the track level.
#
# Currently two attributes can be set:
#    RCov.strainFontFace : generally R supports "italic", "bold", "plain", "bold-italic"
#    RCov.strandRequest : set to either of "+", "-", or "*"

attributes(strains$'minus')$RCov.strandRequest <- "-"
attributes(strains$'all')$RCov.strandRequest <- "*"

#might as well set our pnp track names italic as well:
attributes(strains$'plus'$pnp)$RCov.strainFontFace <- "italic"
attributes(strains$'minus'$pnp)$RCov.strainFontFace <- "italic"
attributes(strains$'all'$pnp)$RCov.strainFontFace <- "italic"

#review attribute settings:
rapply(strains, attributes, how="replace")

#set new dimensions
width <- 4
height <- 3

#then run plotStacked to draw the graph
dev.new(width=width,height=height, unit="in", noRStudioGD = T)
plotStacked(target = target, strains = strains, labelGroups=T, height = height, width = width)



#### SecA ####

#secA will be used in the following examples
rna <- "secA"
exons[exons$locus == rna,]

target <- GRanges(seqname="NC_000913.3", 
                  ranges=IRanges(
                    start = exons[exons$locus == rna,]$left-2000, 
                    end = exons[exons$locus == rna,]$right+1700 ), 
                  strand=exons[exons$locus == rna,]$strand )

#define a custom color palette as well
palette=c("grey10","darkred","grey10","darkred")

###Group names don't have to be shown. In this case, strain names will be on the left
# first organize our tracks:
strains <- list(
  "group 1"=list(
    "WT +"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T),
    "WT -"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T)
  ),
  "group 2"=list(
    "pnp +"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T),
    "pnp -"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  )
)

#assign attributes once again
attributes(strains$'group 2')$RCov.strainFontFace <- "italic"
attributes(strains$'group 1'$'WT -')$RCov.strandRequest <- "-"
attributes(strains$'group 2'$'pnp -')$RCov.strandRequest <- "-"

#review attribute settings:
rapply(strains, attributes, how="replace")

#run plotStacked to draw the graph
quartz(width=width,height=height)
plotStacked(target = target, strains = strains, palette = palette, labelGroups=F, height = height, width = width)



###Single tracks per group also works just fine
strains <- list(
  "group 1"=list(
    "WT +"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T)
  ),
  "group 2"=list(
    "WT -"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T)
  ),
  "group 3"=list(
    "pnp +"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  ),
  "group 4"=list(
    "pnp -"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  )
)


#assign attributes once again
attributes(strains$'group 3')$RCov.strainFontFace <- "italic"
attributes(strains$'group 4')$RCov.strainFontFace <- "italic"
attributes(strains$'group 2'$'WT -')$RCov.strandRequest <- "-"
attributes(strains$'group 4'$'pnp -')$RCov.strandRequest <- "-"

#review attribute settings:
rapply(strains, attributes, how="replace")

#run plotStacked to draw the graph
quartz(width=width,height=height)
plotStacked(target = target, strains = strains, palette = palette, labelGroups=F, height = height, width = width)



###One group with many tracks is also fine
strains <- list(
  "group 1"=list(
    "WT +"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T),
    "WT -"=list.files("example/coverage/TPM", pattern="KR10000", ignore.case=T, full.names = T),
    "pnp +"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T),
    "pnp -"=list.files("example/coverage/TPM", pattern="NRD999", ignore.case=T, full.names = T)
  )
)


#assign attributes once again
attributes(strains$'group 1'$'pnp +')$RCov.strainFontFace <- "italic"
attributes(strains$'group 1'$'pnp -')$RCov.strainFontFace <- "italic"
attributes(strains$'group 1'$'WT -')$RCov.strandRequest <- "-"
attributes(strains$'group 1'$'pnp -')$RCov.strandRequest <- "-"

#review attribute settings:
rapply(strains, attributes, how="replace")

#run plotStacked to draw the graph
quartz(width=width,height=height)
plotStacked(target = target, strains = strains, palette = palette, labelGroups=F, height = height, width = width)
