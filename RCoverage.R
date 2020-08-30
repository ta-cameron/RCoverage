#=========================do this once=========================
#run the following line to install packages from bioconductor
#install.packages("BiocManager") ; library("BiocManager")
#install(c("GenomicRanges","rtracklayer","Rsamtools","GenomicFiles","GenomicAlignments"))

#==================run this each new R session==================
packages <- c("rstudioapi", "GenomicRanges","rtracklayer","Rsamtools", "GenomicFiles", "GenomicAlignments",
              "ggplot2", "grid","gridExtra","gtable","egg", "viridis")

#force install new packages
#install.packages(packages,dep=TRUE)

ipak <- function(pkg){ #ipak from stevenworthington / GitHub
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages); rm(packages,ipak)

#inferExperiment
# strandedness / paired determination function
# this can be used independently of the other functions to example a bam file of your choice
# simply give it the file path and the gff file path
inferExperiment <- function(file, gffFile=NULL, sampleSize=200000, substituteChromosomeNames=NULL){
  
  ### import gff file ###
  tryCatch( exonsGR <- import(gffFile),
            error=function(e) stop("Error importing genome reference file: ",gffFile,"\n\t",e[1])
  )
  exonsGR <- exonsGR[exonsGR$type =="gene"]
  exonsGR <- exonsGR[,c("gene_biotype","Name")]
  exonsGRStrands <- split(exonsGR, strand(exonsGR))
  rm(exonsGR)
  
  cat("========= Start of inferExperiment =========\n")
  
  #check dependencies
  tryCatch( library("Rsamtools", "GenomicFiles", "GenomicAlignments"),
            error = function(e) stop("Required packages not installed: ",e[1]) )
  
  #warn on small sample size
  if(sampleSize < 1000){
    cat("\x1b[33mSample size is less than 1000... results may be unreliable.\n\x1b[0m")
    Sys.sleep(2)
  }
  
  #set up bam file sampling
  bam <- BamFile(file, yieldSize=1000000)
  yield <- function(x) {
    readGAlignments(x, param=ScanBamParam(what="flag", flag = scanBamFlag(isUnmappedQuery=F)) )
  }
  map <- identity
  #sample reads
  subBAM <- reduceByYield(bam, yield, map, REDUCEsampler(sampleSize, TRUE))
  
  #replace unwanted chromosome names
  if( exists("substituteChromosomeNames") && !is.null(substituteChromosomeNames) ){
    if( !is.data.frame(substituteChromosomeNames) || ncol(substituteChromosomeNames) != 2 ){
      stop("Please supply a data.frame, where the first column is the query and the second column the replacement")
    }
    
    subBAMseq <- seqlevels(subBAM)
    #substitute genome names
    names(substituteChromosomeNames) <- c("original","replacement")
    for(subID in 1:nrow(substituteChromosomeNames) ){
      subBAMseq[grepl(substituteChromosomeNames$original[subID],subBAMseq)] <- as.character(substituteChromosomeNames$replacement[subID])
    }
    subBAM <- renameSeqlevels(subBAM, subBAMseq)  
  }
  #check for non-unique chromosome names now
  if( length(seqlevels(exonsGRStrands)) != length(unique(seqlevels(exonsGRStrands))) ){
    stop("non-unique chromosome names!",
         "\n\t names: ",toString(idx$chr))
  }
  #check for matching chromosome names now
  if(!all(seqlevels(exonsGRStrands) %in% seqlevels(subBAM)) ){
    stop("chromosome names in annotation file are missing from BAM file:",
         "\n\t annotation: ",toString(seqlevels(exonsGRStrands)),
         "\n\t   BAM file: ",toString(seqlevels(subBAM)) )
  }
  
  #split by read order (first / second)
  #also determine here if reads contain pairs
  subBAMOrder <- split(subBAM, bamFlagTest(subBAM@elementMetadata$flag, "isFirstMateRead"))
  switch( length(names(subBAMOrder)), 
          { names(subBAMOrder) <- c("first") ; paired <- F },
          { names(subBAMOrder) <- c("second","first") ; paired <- T }
  )
  
  #prepare table of results
  strandTable <- data.frame(
    hits=rep(0,4),
    row.names=c("forward","reverse","unknown","total")
  )
  
  for(readOrder in names(subBAMOrder)){
    #readOrder <- names(subBAMOrder)[1]
    #split by strand of read
    subBAMStrands <- split(subBAMOrder[[readOrder]], strand(subBAMOrder[[readOrder]]))
    
    switch(readOrder,
           first={
             forwardCount <- "forward"
             reverseCount <- "reverse"
           },
           second={
             forwardCount <- "reverse"
             reverseCount <- "forward"
           })
    
    strandTable["total","hits"] <- strandTable["total","hits"] + length(subBAMStrands$'+') + length(subBAMStrands$'-')
    strandTable[forwardCount,"hits"] <- strandTable[forwardCount,"hits"] + 
      length( na.omit(findOverlaps(subBAMStrands$'+',exonsGRStrands$'+', ignore.strand=T, select="arbitrary") )) +
      length( na.omit(findOverlaps(subBAMStrands$'-',exonsGRStrands$'-', ignore.strand=T, select="arbitrary") ))
    strandTable[reverseCount,"hits"] <- strandTable[reverseCount,"hits"] +
      length( na.omit(findOverlaps(subBAMStrands$'+',exonsGRStrands$'-', ignore.strand=T, select="arbitrary") )) +
      length( na.omit(findOverlaps(subBAMStrands$'-',exonsGRStrands$'+', ignore.strand=T, select="arbitrary") ))
  }
  #number of reads over no genes
  #this is taken as the total number of reads minus the number of forward and reverse hits
  #technically this can be a negative number if there are many ambiguous reads
  #this is not a measure of ambiguous reads, but reflects the number of reads that could not be considered
  strandTable["unknown","hits"] <- strandTable["total","hits"]-sum(strandTable[c("forward","reverse"),"hits"])
  if(strandTable["unknown","hits"] < 0 ){
    strandTable["unknown","hits"] <- 0
  }
  
  #percent of reads mapping to gene
  strandTable$percent <- round(strandTable$hits/strandTable["total","hits"]*100, 1)
  
  #report some suspicious scenarios
  if(sum(strandTable[c("forward","reverse"),"hits"]) < strandTable["unknown","hits"] ){
    cat("\x1b[33mWarning: More than half of all reads were not assigned to genes. Is the annotation correct?\n\x1b[0m")
    Sys.sleep(5)
  }
  if(sum(strandTable[c("forward","reverse"),"hits"]) < 100){
    cat("\x1b[33mWarning: Less than 100 reads total were assigned to genes... this may not produce reliable results.\n\x1b[0m")
    Sys.sleep(5)
  }
  
  cat("Summary Table:\n")
  print(strandTable)
  
  #report support for 
  forward <- round(strandTable["forward","hits"] / sum(strandTable[c("forward","reverse"),"hits"]) *100,1)
  reverse <- round(strandTable["reverse","hits"] / sum(strandTable[c("forward","reverse"),"hits"]) *100,1)
  cat(sep="", "\nForward: ",forward,"% of reads mapped to genes support forward strand. (",strandTable["forward","percent"],"% map to a gene if forward-stranded)\n")
  cat(sep="", "Reverse: ",reverse,"% of reads mapped to genes support reverse strand. (",strandTable["reverse","percent"],"% map to a gene if reverse-stranded)\n\n")
  
  if( forward == reverse || abs(log(forward/reverse,2)) <= 2 ){
    #unstranded (the cutoff for unstranded accepts 20 - 80% on each strand) 80/20 = ratio of 4, hence log(4,2)==2 as cutoff
    strandMode <- 0
    cat("Final determination:\n\t-----> unstranded (0) <-----\n")
  } else if( log(forward/reverse,2) > 2 ){
    #forward stranded (first read on forward strand)
    strandMode <- 1
    cat("Final determination:\n\t-----> forward stranded (1) <-----\n")
  } else {
    #reverse stranded (first read on reverse strand)
    strandMode <- 2
    cat("Final determination:\n\t-----> reverse stranded (2) <-----\n")
  }
  cat("\n========= End of inferExperiment =========\n\n")
  return(list("paired"=paired, "strandMode"=strandMode))
}  


#RCoverage function
RCoverage <- function(dataDir=NULL, outDir=NULL, TPM=T, 
                      gffFile=NULL, inferExperiments=T, strandMode=0, paired=F,
                      countSinglets=F, onlyProperPairFragments=T,
                      filter=T, maxTLEN=600, substituteChromosomeNames=NULL){
  
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
  # 
  # #optional chromosome name replacement
  # # save a data.frame to the variable "substituteChromosomeNames", where the first column is the query and the second column is the replacement
  # # partial matches will be substituted. Comment out next line to disable replacement.
  # #substituteChromosomeNames <- data.frame(c("query"),c("replacement"))
  
  #error check maxTLEN setting
  if(!is.numeric(maxTLEN) || length(maxTLEN)!=1){
    message("maxTLEN should be a single numeric value... attempting coersion")
    tryCatch({maxTLEN <- as.numeric(maxTLEN[1])},
             error= function(e) {stop("maxTLEN must be numeric.",e)}
    )
    message("maxTLEN: ",maxTLEN)
  }
  
  #add trailing / to data & output directory, if absent
  dataDir <- sub("/$","",dataDir)
  dataDir <- sub("$","/",dataDir)
  outDir <- sub("/$","",outDir)
  outDir <- sub("$","/",outDir)
  
  #make a working space
  dir.create(paste0(outDir,"counts"), recursive=T)
  if(TPM){
    dir.create(paste0(outDir,"TPM"))
  }
  
  #here begins the actual evaluation:
  #capture.output(file=paste0(outDir,"RbedToolsCoverage.log"), split=T, {
    
    ### coverage loop ###
    files <- list.files(path=dataDir, pattern="\\.bam$", ignore.case=T)
    files <- sub("^",dataDir,files)
    for(file in files){
      #file <- files[1]
      ptm <- proc.time()
      cat("\n\nProcessing ", file,"\n")
      sample <- sub("\\.bam$","",file, ignore.case=T)
      
      if(inferExperiments){
        if( exists("substituteChromosomeNames") && !is.null(substituteChromosomeNames) ){
          inferredResults <- inferExperiment(file=file, gffFile=gffFile, substituteChromosomeNames=substituteChromosomeNames)
        } else {
          inferredResults <- inferExperiment(file=file, gffFile=gffFile)
        }
        if(inferredResults[["paired"]] != paired){
          cat(sep="","\nInferred paired result (",inferredResults[["paired"]],") overriding previous paired setting (",paired,").\n")
          paired <- inferredResults[["paired"]]
        } 
        if(inferredResults[["strandMode"]] != strandMode){
          cat(sep="","\nInferred strandedness result (",inferredResults[["strandMode"]],") overriding previous strandedness setting (",strandMode,").\n")
          strandMode <- inferredResults[["strandMode"]]
        }  
      }
      
      ## echo back some settings
      cat("\n  ### Settings ###\n")
      cat("  data folder :",dataDir, "\n")
      cat("  calculate bins per million :",TPM, "\n")
      cat("  strandMode :",strandMode, "\n")
      cat("  paired end :",paired, "\n")
      if((filter || countSinglets || onlyProperPairFragments) && !paired){
        cat("  require proper pairs : NA \n")
        cat("  include unpaired reads : NA \n")
        cat("  filter fragment sizes : NA \n")
        cat("  \t maximum size : NA \n")
      } else {
        cat("  require proper pairs :", onlyProperPairFragments, "\n")
        cat("  include unpaired reads :", countSinglets, "\n")
        cat("  filter fragment sizes :",filter, "\n")
        cat("  \t maximum size :",maxTLEN, "\n")
      }
      
      #generate genome index (only used in alternative export function)
      if( file.access(paste0(file,".bai"),4) == -1 ){
        cat("\n  No index file found... generating one (indexBam)\n")
        indexBam(file)
      }
      cat("\n  Reading chromosome information (scanBamHeader)\n")
      idx <- scanBamHeader(file)[[1]]$target
      idx <- data.frame(chr=names(idx), start=1,end=idx, stringsAsFactors = F)
      row.names(idx) <- NULL
      
      #read flagstats to find number of proper pairs or mapped reads
      cat("\n  Determining flag statistics (idxstatsBam / countBam)\n")
      unmappedReads <- sum(idxstatsBam(file)$unmapped)
      mappedReads <- sum(idxstatsBam(file)$mapped)
      propPairs <- countBam(file, param=ScanBamParam(flag=scanBamFlag(isProperPair = T) ))$records/2
      cat("    unmapped reads:  ", unmappedReads,"\n    mapped reads:    ", mappedReads,"\n    properly paired: ", propPairs*2,"\n")
      
      #import reads
      cat("\n  Importing read information (readGAlignments / readGAlignmentsPairs)\n")
      if(exists("last_time_per_mapped")){
        cat("    ETA ~", round(last_time_per_mapped * mappedReads / 60,1),"minutes.\n")
      } else {
        cat("    ETA unknown.\n")
      }
      
      #import unpaired reads
      flags <- scanBamFlag(isDuplicate = F,
                           isSecondaryAlignment = F,
                           isPaired = F)
      param <- ScanBamParam(flag=flags,
                            what=c("flag"))
      reads <- list("singles"=readGAlignments(file, param=param) )
      
      #reverse strand if called for
      if(strandMode==2){
        reads[["singles"]] <- invertStrand(reads[["singles"]] )
      }
      
      #make the pairs section blank
      reads[["pairs"]] <- reads[["singles"]][0]
      
      if(paired){
        #import paired reads
        flags <- scanBamFlag(isDuplicate = F,
                             isSecondaryAlignment = F)
        param <- ScanBamParam(flag=flags,
                              what=c("flag"))
        reads <- list( "pairs" = readGAlignmentPairs(file, param=param, strandMode = strandMode) )
        
        #make the singles section blank
        reads[["singles"]] <- reads[["pairs"]][0]
        
        if(onlyProperPairFragments){
          if(countSinglets){
            #reads[["pairs"]][!bamFlagTest(first(reads[["pairs"]])@elementMetadata$flag,"isProperPair")] 
            reads[["singles"]] <- append(reads[["singles"]], after=length(reads[["singles"]]), 
                                         first(reads[["pairs"]][!isProperPair(reads[["pairs"]])], real.strand=T))
            reads[["singles"]] <- append(reads[["singles"]], after=length(reads[["singles"]]), 
                                         last(reads[["pairs"]][!isProperPair(reads[["pairs"]])], real.strand=T))
          }
          reads[["pairs"]] <- reads[["pairs"]][isProperPair(reads[["pairs"]])]
        }
        
        if(filter){
          if(countSinglets){
            reads[["singles"]] <- append(reads[["singles"]], after=length(reads[["singles"]]), 
                                         first(reads[["pairs"]][width(GRanges(reads[["pairs"]])) >= maxTLEN], real.strand=T))
            reads[["singles"]] <- append(reads[["singles"]], after=length(reads[["singles"]]), 
                                         last(reads[["pairs"]][width(GRanges(reads[["pairs"]])) >= maxTLEN], real.strand=T))
          }
          reads[["pairs"]] <- reads[["pairs"]][width(GRanges(reads[["pairs"]])) < maxTLEN]
        }
        cat("\n  Determining final fragment total: ")
        cat("\n  Paired fragments: ",length(reads[["pairs"]]))
        cat("\n    + single reads: ",length(reads[["singles"]]))
        propPairs <- length(reads[["pairs"]]) + length(reads[["singles"]])
        cat("\n          in total: ",propPairs,"\n\n")
      }
      
      #strandMode is TRUE if it is 1 or 2, FALSE if 0
      if(strandMode){
        strands <- c("+","-")
      } else {
        strands <- "*"
      }
      
      #split up by strand
      reads <- lapply(reads, function(x) {
        split(x, strand(x))  
      })
      
      #calculate coverage of reads and pairs (simple start-end)
      bedCoords <- list()
      for(strand in strands){
        #strand <- "*"
        bed <- list()
        for(type in names(reads) ){
          #type <- "pairs"
          bed[[type]] <- coverage( GRanges(reads[[type]][[strand]]) )
        }
        bedCoords[[strand]] <- bed[["pairs"]] + bed[["singles"]]
      }
      rm(bed, reads)
      
      #calculate unstranded coverage if not already present
      if(!"*" %in% names(bedCoords)){
        bedCoords[["*"]] <- bedCoords[["+"]] + bedCoords[["-"]]
      }
      
      #replace unwanted chromosome names, again
      if( exists("substituteChromosomeNames") && !is.null(substituteChromosomeNames) ){
        if( !is.data.frame(substituteChromosomeNames) || ncol(substituteChromosomeNames) != 2 ){
          stop("Please supply 'substituteChromosomeNames' as a data.frame, where the first column is the query name and the second column the replacement name")
        }
        
        #substitute genome names
        names(substituteChromosomeNames) <- c("original","replacement")
        
        for(strand in names(bedCoords) ){
          #strand <- names(bedCoords)[1]
          strains <- names(bedCoords[[strand]])
          for(subID in 1:nrow(substituteChromosomeNames) ){
            strains[grepl(substituteChromosomeNames$original[subID],strains)] <- as.character(substituteChromosomeNames$replacement[subID])
            idx[grepl(substituteChromosomeNames$original[subID],idx$chr), "chr"] <- as.character(substituteChromosomeNames$replacement[subID])
          }
          #check for non-unique chromosome names now
          if( length(strains) != length(unique(strains)) ){
            stop("non-unique chromosome names!",
                 "\n\t names: ",toString(strains))
          }
          names(bedCoords[[strand]]) <- strains
        }
        
      }
      
      #export bigwig files
      for( strand in names(bedCoords) ){
        strandName <- switch(strand,
                             "+"="_plus",
                             "-"="_minus",
                             "*"="")
        
        #export counts
        cat("  Saving strand: ", strandName,"\n")
        
        ##easy RLE export
        export.bw(bedCoords[[strand]], paste0(outDir,"counts/",basename(sample),strandName,".bw"))
        
        # ##alternative export; this produces smaller files due to fixed step size
        # df <- as.data.frame(bedCoords[[strand]])
        # names(df) <- c("pos","chr","score")
        # for(chr in unique(df$chr) ){
        #   df[df$chr==chr,"pos"] <- 1:nrow(df[df$chr==chr,])
        # }
        # dfGR <- makeGRangesFromDataFrame(df,
        #                                  keep.extra.columns=T,
        #                                  ignore.strand=T,
        #                                  seqinfo=Seqinfo(idx$chr, seqlengths=idx$end, isCircular=NA),
        #                                  seqnames.field=c("chr"),
        #                                  start.field="pos",
        #                                  end.field="pos",
        #                                  starts.in.df.are.0based=F
        # )
        # export.bw(dfGR, paste0("coverage/counts/", basename(sample), "_", strandName, "_aE.bw"))
        
        
        #export TPM normalize counts if requested
        if(TPM){
          if(paired){
            tpmDenom <- propPairs
          } else{
            tpmDenom <- mappedReads    
          }
          
          if(tpmDenom != 0){
            cat("    Normalization factor: ", 1/tpmDenom * 1000000 ,"\n")
            bedCoords[[strand]] <- bedCoords[[strand]] / tpmDenom * 1000000
          } else {
            cat("  Error: skipping normalization... tpmDenom value is invalid : ", tpmDenom,"\n")
          }
          
          ##easy RLE export of normalized counts
          export.bw(bedCoords[[strand]], paste0(outDir,"TPM/",basename(sample),strandName,".bw"))
          
          # ##alternative export; this produces smaller files due to fixed step size
          # df <- as.data.frame(bedCoords[[strand]])
          # names(df) <- c("pos","chr","score")
          # for(chr in unique(df$chr) ){
          #   df[df$chr==chr,"pos"] <- 1:nrow(df[df$chr==chr,])
          # }
          # dfGR <- makeGRangesFromDataFrame(df,
          #                                  keep.extra.columns=T,
          #                                  ignore.strand=T,
          #                                  seqinfo=Seqinfo(idx$chr, seqlengths=idx$end, isCircular=NA),
          #                                  seqnames.field=c("chr"),
          #                                  start.field="pos",
          #                                  end.field="pos",
          #                                  starts.in.df.are.0based=F
          # )
          # export.bw(dfGR, paste0("coverage/TPM/", basename(sample), "_", strandName, "_aE.bw"))
          
        }
      }
      
      #enforce clean up
      rm(bedCoords)
      
      cat("\nSystem stats:\n")
      print.table(proc.time() - ptm)
      last_time <- (proc.time() - ptm)["elapsed"]
      last_time_per_mapped <- last_time / mappedReads
    }
#  })
}


### stacked coverage graphing - main
plotStacked <- function(target=NULL, strains=NULL, palette=NULL, 
                        labelGroups=T, height=3, width=3, binsPerInch=200, textScaling=1){
  
  #labelGroups=T; binsPerInch=200; textScaling=1
  #labelGroups=T; height=3; width=3; binsPerInch=200; textScaling=1; strands=NULL
  
  #check textScaling
  stopifnot(is.numeric(textScaling))
  stopifnot(textScaling > 0)
  
  #check if group names are unique
  if( length(names(strains)) != length(unique(names(strains))) ){
    stop("Please provide unique group names.")
  }
  
  #check if strain names are unique
  if( any(!sapply(strains,function(x) length(names(x))) == sapply(strains,function(x) length(unique(names(x))))) ){
    stop("Please provide unique strain names (within each group).")
  }
  
  #check if target is proper
  if( class(target) != "GRanges" ){
    stop("must supply GRanges object as target.")
  } else if( length(target) != 1){
    stop("must supply single GRanges sequence range.")
  } else if( length(seqnames(target)) != 1 ){
    stop("must supply a single GRanges seqname")
  }
  left <- start(target)
  right <- end(target)
  chromosome <- as.character(seqnames(target))
  
  #check if strains is semi-reasonable input
  if( any( unlist(lapply(strains, function(x) {is.null(x) | length(x) == 0})) ) ){
    stop("strains must be supplied.")
  } else if( any(!sapply(strains, function(x) {is.list(x)})) | any( !unlist(sapply(strains, function(x) sapply(x, function(y) {is.character(y)}))) ) ){
    stop("strains must be supplied as a list of groups, where each group is a list of character vectors specifying the files, \n eg: strains=list(
           \t\t\"group1\"=list(
           \t\t\t\"strain1\"=c(\"strain1_file1.bw\",\"strain1_file2.bw\"),
           \t\t\t\"strain1\"=c(\"strain1_file1.bw\",\"strain1_file2.bw\")
           \t\t\t)
           \t\t) ")
  } 
  
  ##make a data.frame of strain attributes
  strainInfo <- as.data.frame(matrix(nrow=sum(sapply(strains, length)),
                                     ncol=3,
                                     dimnames=list(NULL,c("id","group","strain"))) )
  row <- 1
  for( group in names(strains) ){
    for( strain in names(strains[[group]])){
      strainInfo[row,] <- c(paste0(group,"-",strain), group, strain)
      row <- row + 1
    }
  }
  #place default values
  strainInfo$strandRequest <- as.character(strand(target))
  strainInfo$strand <- NA
  strainInfo$strainFontFace <- "plain"
  strainInfo$covPal <- viridis_pal(end=0.8)( sum(sapply(strains, length)) )
  
  ##pull out any attributes set on the strains list, then apply them top-to-bottom (priority to bottom)
  RCovAttrVector <- c("strandRequest","strainFontFace")
  #apply top level attributes
  for( RCovAttr in RCovAttrVector ){
    if( !is.null(attr(strains, paste0("RCov.",RCovAttr) )) ) { 
      strainInfo[,RCovAttr] <- attr(strains,paste0("RCov.",RCovAttr)) 
    }
  }
  #next apply group and then strain attributes
  for( group in names(strains) ){
    for( RCovAttr in RCovAttrVector ){
      if( !is.null(attr(strains[[group]], paste0("RCov.",RCovAttr) )) ) { 
        strainInfo[strainInfo$group == group,RCovAttr] <- attr(strains[[group]],paste0("RCov.",RCovAttr)) 
      }
    }
    for( strain in names(strains[[group]])){
      for( RCovAttr in RCovAttrVector ){
        if( !is.null(attr(strains[[group]][[strain]], paste0("RCov.",RCovAttr) )) ) { 
          strainInfo[strainInfo$group == group & strainInfo$strain == strain,RCovAttr] <- attr(strains[[group]][[strain]],paste0("RCov.",RCovAttr)) 
        }
      }
    }
  }
  
  ##check if supplied palette covers all strains
  if(!is.null(palette)){
    if( length(palette) < nrow(strainInfo) ){
      warning("Palette was insufficient for requested strains (need ", sum(sapply(strains, length)) ," colors). Using default colors instead.")
    } else {
      strainInfo$covPal <- palette[1:nrow(strainInfo)]
    }
  }
  
  #pull in the data for this region from all available files
  strands <- strainInfo$strandRequest
  dataReps <- extractCov(left,right,strands,chromosome,binsPerInch)
  
  #convert list to a data frame
  # THIS ONLY WORKS PROPERLY IF THE COLUMNS ALL HAVE THE SAME NAME (which they should, as they are hard-coded above)
  slice <- do.call(rbind.data.frame, dataReps)
  row.names(slice) <- NULL
  
  #place trimmed data set in list of each group
  covData <- list()
  for( group in unique(slice$group) ){
    #group <- unique(slice$group)[1]
    covSlice <- slice[slice$group %in% group,]
    covSlice$strain <- factor(covSlice$strain, levels=names(strains[[group]]) )
    covData[group] <- list(covSlice)
  }
  #add list name as an attribute
  for (i in 1:length(covData)){
    if (is.null(attributes(covData[[i]]))){
      attributes(covData[[i]]) <- list(ref=names(covData)[i])
    } else {
      attributes(covData[[i]]) <- c(attributes(covData[[i]]), ref=names(covData)[i])
    }
  }
  
  #record actual strand retrieved in strainInfo data.frame    
  for(group in names(strains)){
    strainStrandPairs <- unique(covData[[group]][,c("strain","strand")])
    strainStrandPairs <- as.data.frame( lapply(strainStrandPairs, as.character), stringsAsFactors=F)
    for(strain in strainStrandPairs$strain){
      strainInfo[strainInfo$group == group & strainInfo$strain == strain,"strand"] <- strainStrandPairs[strainStrandPairs==strain,"strand"]
    }
  }
  
  #realign left/right limits to actual coordinate limits
  left <- min(slice$start)
  right <- max(slice$end)
  
  #set up gene track gene names
  exonsSlice <- exons[ (exons$left < right & exons$right > left),]
  if( nrow(exonsSlice) > 0 ){
    #save original exon boundaries
    exonsSlice$exonLeft <- exonsSlice$left
    exonsSlice$exonRight <- exonsSlice$right
    
    #truncate to plotting area
    exonsSlice[exonsSlice$left < left,"left"] <- left
    exonsSlice[exonsSlice$right > right,"right"] <- right
    
    #convert from factor
    exonsSlice$strand <- as.character(exonsSlice$strand)
    exonsSlice$locus <- as.character(exonsSlice$locus)
    
    #at the requested plot width, how many inches/nt?
    # 0.0787402 == sum of plot border (2 mm) in inches
    # X*textScaling; X == width of side bar (0.2 in) and spacer (0.01 in)
    inPerX <- (width-0.0787402-0.21*textScaling) / (right-left)
    
    #process each exon name
    for(row in 1:nrow(exonsSlice)){
      #row <- 1
      #determine length of gene box
      maxIn <- (exonsSlice[row,"right"]-exonsSlice[row,"left"])*inPerX
      
      #pad the name with a space on each side and determine the length of this gene name
      exonsSlice[row,"locus"] <- paste0(" ",exonsSlice[row,"locus"]," ")
      labelIn <- strwidth(exonsSlice[row,"locus"], units="in",  cex=(8*textScaling)/12 )
      
      #check if gene name is too long to fit in the gene box
      if( maxIn < labelIn ){
        #if so, remove it
        exonsSlice[row,"locus"] <- ""
      }
      
      #if gene area is too small, don't draw the arrow
      if(maxIn < 0.05){
        exonsSlice[row,"arrow"] <- F
      } else {
        exonsSlice[row,"arrow"] <- T
      }
      
      #if gene is truncated in the view, don't draw the arrow
      if(exonsSlice[row,"exonRight"] > right & exonsSlice[row,"strand"] == "+"){
        exonsSlice[row,"arrow"] <- F
      }
      if(exonsSlice[row,"exonLeft"] < left & exonsSlice[row,"strand"] == "-"){
        exonsSlice[row,"arrow"] <- F
      }
    }
  } else {
    exonsSlice <- data.frame(left=left, right=right, locus=NA)
  }
  
  #review data before plotting
  head(exonsSlice)
  lapply(covData, head)
  strainInfo
  
  verticalSpacing <- 0.02 #for null unit spacing; distance between each coverage graph
  
  ggList <- list()
  #ggList <- lapply(covData, function(groupData){
  for( group in names(covData)){
    #group <- names(covData)[[1]]
    groupData <- covData[[group]]
    tickTop <- round(max(groupData[,"covMeans"]), digits=1)
    
    #retrieve matching data from strainInfo
    ggPal <- strainInfo[strainInfo$group==group,"covPal"]
    ggStrains <- strainInfo[strainInfo$group==group,"strain"]
    ggStrands <- strainInfo[strainInfo$group==group,"strand"]
    ggStrainFontFace <- strainInfo[strainInfo$group==group,"strainFontFace"]
    
    #prepare list of ggplots, one per group
    if(labelGroups){
      plotMargin <- unit(c(verticalSpacing,0,verticalSpacing,0.5), c("null","null","null","mm"))
      strandLabel <- paste(ggStrains, paste0("(",ggStrands,")"))
      strainFontFace <- ggStrainFontFace
    } else {
      plotMargin <- unit(c(verticalSpacing/2,0,verticalSpacing/2,0.5), c("null","null","null","mm"))
      strandLabel <- paste0("(",ggStrands,")")
      strainFontFace <- "plain"
    }
    strandLabel <- sub(" \\(\\*\\)$","",strandLabel)
    
    #alternative data labeling approach
    ggStrainInfo <- strainInfo[strainInfo$group==group,c("strain","strainFontFace")]
    ggStrainInfo$strainFontFace <- strainFontFace
    ggStrainInfo$strandLabel <- strandLabel
    ggStrainInfo$strain <- factor(levels(groupData$strain),levels=levels(groupData$strain))
    
    levels(groupData$strain)
    levels(ggStrainInfo$strain)
    
    ggCov <- 
      ggplot(data=groupData) +
      geom_rect(aes(xmin=start, xmax=end, ymin=0,ymax=covMeans, fill=strain)) +
      scale_fill_manual(values = ggPal ) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(breaks=c(tickTop), limits=c(0,tickTop*1.01)) +
      coord_cartesian(xlim=c(left,right)) +
      facet_wrap(facets=~strain, dir="v", ncol=1) +
      theme_grey() +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "mm"),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        panel.spacing = unit(c(verticalSpacing),c("null")),
        plot.margin = plotMargin,
        plot.background = element_blank(),
        panel.background =  element_rect(fill = "grey95"),
        panel.border = element_blank(),
        panel.grid = element_blank()
      ) +
      labs(x=NULL, y=NULL) +
      annotate("label", # coverage minimum scale ('0') label
               x=left,
               y=0,
               label=0, size=2.2*textScaling, hjust=0, vjust=0, label.r=unit(0,"cm"),
               fill="grey95", alpha =0.7, label.size=0, label.padding=unit(c(0.05),"lines")) +
      annotate("label", # coverage maximum scale label
               x=left,
               y=tickTop,
               label=tickTop, size=2.2*textScaling, hjust=0, vjust=1, label.r=unit(0,"cm"),
               fill="grey95", alpha =0.7, label.size=0, label.padding=unit(c(0.05),"lines")) +
      geom_label( # strain / strand label
               data=ggStrainInfo,
               mapping=aes(x=right, y=tickTop, label=strandLabel, fontface=strainFontFace),
               size=2.2*textScaling, hjust=1, vjust=1, label.r=unit(0,"cm"),
               fill="grey95", alpha =0.7, label.size=0, label.padding=unit(c(0.05),"lines"))
    ggList[[group]] <- ggplotGrob(ggCov)
  }
  
  #prepare ggplot of gene track
  ggGenes <- ggplot(data=exonsSlice) + 
    theme_classic() + 
    scale_x_continuous(breaks=c(left, right), expand=c(0,0))+
    coord_cartesian(xlim=c(left,right)) +
    theme(      
      plot.margin = unit(c(verticalSpacing,0,0,0.5), c("null","null","null","mm")),
      panel.spacing = unit(verticalSpacing, "null"),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0, "cm"),
      axis.line = element_blank(), 
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=7*textScaling, hjust=c(0,1)),
      strip.background = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank()) +
    labs(x=NULL, y=NULL) +
    annotate("rect", 
             xmin=left, 
             xmax=right,
             ymin=1,
             ymax=2,
             fill=NA)
  
  for( row in 1:nrow(na.omit(exonsSlice)) ){
    exonName <- na.omit(exonsSlice)[row,"locus"]
    exonFill <- ifelse(na.omit(exonsSlice)[row,"strand"] == "+","darkgreen","darkblue")
    
    #set up coordinates for a rectangle or arrow
    if(na.omit(exonsSlice)[row,"arrow"] || !na.omit(exonsSlice)[row,"strand"] %in% c("-","+") ){
      if( na.omit(exonsSlice)[row,"strand"] == "+" ){
        #right arrow; note that inches must be substracted from npc unit to get the right-hand arrow
        exonPolyX <- unit(c(1,1,0,0,1),c("npc")) - unit(c(0,0.04,0,0,0.04),c("npc","in","npc","npc","in"))
        exonPolyY <- unit(c(0.5,0,0,1,1),c("npc"))
      } else {
        #left arrow
        exonPolyX <- unit(c(0,0.04,1,1,0.04),c("npc","in","npc","npc","in"))
        exonPolyY <- unit(c(0.5,1,1,0,0),c("npc"))
      }
    } else {
      #simple rectangle
      exonPolyX <- unit(c(0,1,1,0),c("npc"))
      exonPolyY <- unit(c(1,1,0,0),c("npc"))
    }
    
    exonGrob <- grobTree(
      polygonGrob( x=exonPolyX,
                   y=exonPolyY, 
                   gp=gpar(col="white", fill=exonFill)),
      textGrob(exonName, gp=gpar(fontsize=8*textScaling, col="white", fontface="italic"), vp=viewport(clip="on")) )
    
    ggGenes <- ggGenes + 
      annotation_custom(exonGrob, 
                        xmin=na.omit(exonsSlice)[row,"left"],
                        xmax=na.omit(exonsSlice)[row,"right"],
                        ymin=1,
                        ymax=2)
  }
  
  ##manual grid composition to add some labeling
  grid.newpage()
  #first make a single central viewport bordered by 1 mm on each side 
  pushViewport(viewport(layout = grid.layout( 3, 3,
                                              widths=unit(c(1,1,1),c("mm","null","mm")),
                                              heights=unit(c(0.5,1,0.5),c("mm","null","mm"))  )))
  
  
  #next layout the coverage plot structure within the central region
  pushViewport(viewport(
    layout = grid.layout(length(ggList)+1, 3, 
                         widths=unit(c(0.2*textScaling,0.01*textScaling,1),c("in","in","null")), 
                         heights=unit(c(sapply(strains, function(x) length(names(x)) ),0.33*textScaling), 
                                      c(rep("null",length(ggList)),"in")) ),
    layout.pos.row=2, layout.pos.col=2
  ))
  
  #to see the layout scheme use:
  #grid.show.layout( grid.layout(...) )
  
  #define convenience function for setting viewport column and row
  vplayout <- function(x,y) viewport(layout.pos.col=x,layout.pos.row=y,clip="on")
  
  #iterate over each group / individual ggplot graphs
  for(row in 1:length(ggList)){
    #determine the number of strains graphed in this single ggplot grob
    rows <- length(names(strains[[row]]))
    
    #print the ggplot graph
    pushViewport(vp=vplayout(3,row))
    grid.draw(ggList[[row]])
    popViewport()
    
    ##draw group or strain label background on the left
    if(labelGroups){
      groupHeights <- unit(c(verticalSpacing/rows, 1, verticalSpacing/rows), c("null") )
    } else {
      groupHeights <- unit(c(0, 1, 0), c("null") )
    }
    
    pushViewport(viewport(layout = grid.layout( 3, 1,
                                                heights=groupHeights ),
                          layout.pos.row=row, layout.pos.col=1))
    
    #either print group names or each strain name, depending on labelGroups settings
    if(labelGroups){
      grid.rect(vp=vplayout(1,2), gp=gpar(fill="grey80", col=NA))
      grid.text(names(ggList[row]), vp=vplayout(1,2),gp=gpar(fontsize=8), rot=90)  
    } else {
      ggStrainFontFace <- strainInfo[strainInfo$group == names(ggList)[row], "strainFontFace"]
      
      pushViewport(viewport(layout = grid.layout( length(strains[[row]]), 1),
                            layout.pos.row=2, layout.pos.col=1))
      for( strainRow in 1:rows ){
        pushViewport(viewport(layout = grid.layout( 3, 1,
                                                    heights=unit(c(verticalSpacing/2,1,verticalSpacing/2),c("null"))  ),
                              layout.pos.row=strainRow, layout.pos.col=1))
        grid.rect(vp=vplayout(1,2), gp=gpar(fill="grey80", col=NA))
        grid.text( names(strains[[row]][strainRow]), vp=vplayout(1,2),gp=gpar(fontsize=8*textScaling, fontface=ggStrainFontFace[strainRow]), rot=90)
        popViewport() 
      }
      popViewport() 
    }
    popViewport()
  }
  print(ggGenes,vp=vplayout(3, 1+length(ggList)))
}

### extract coverage data function
extractCov <- function(left,right,strands,chromosome,binsPerInch){
  
  #check coordinate input
  if(left != as.integer(left) | right != as.integer(right)){
    stop("left and right coordinates must be whole numbers (integers).")
  }
  if(left >= right){
    stop("left coordinate must be less than right coordinate.")
  }
  
  #check & standardize strand input
  if( ! all(strands %in% c("+","-","*")) ){
    stop("strands must be \"+\", \"-\", or \"*\"")
  }
  
  #check coordinate input
  stopifnot(is.character(chromosome))
  
  #check binsPerInch input
  stopifnot(is.numeric(binsPerInch))
  
  ###extract the desired region from bigwig files
  
  ##check for matching non-bigwig files and give a warning if found
  testFiles <- unique(unlist(strains))
  
  bwCheck <- vapply(testFiles, function(fileVAP) {
    grepl(".bw$",fileVAP, ignore.case = T)
  }, logical(1))
  if(!all(bwCheck)){
    warning("Supplied files are not all bigwig files (\".bw\")... please review inputs carefully for any errors.")
  }
  
  #also drop non .bw files from the list
  bigwigFiles <- testFiles[bwCheck]
  
  ##filter by strandedness... discards opposite strand but keeps all
  ## or only keeps all if neither strand
  #First make a fileInfo data frame with meta data for each file
  fileInfo <- data.frame(id=character(0),
                         group=character(0),
                         strain=character(0),
                         file=character(0)
  )
  #fill out fileInfo with group, strain, id (group-strain), and relative file path
  for( group in names(strains) ){
    #group <- names(strains)[1]
    for( strain in names(strains[[group]])){
      #strain <- names(strains[[group]])[1]
      reps <- strains[[group]][[strain]]
      fileInfo <- rbind(fileInfo, 
                        data.frame(
                          id=rep(paste0(group,"-",strain), length(reps)),
                          group=rep(group, length(reps)),
                          strain=rep(strain, length(reps)),
                          file=reps
                        ))
    }
  }
  
  #drop files that failed .bw file check earlier (absent from bigwigFiles)
  fileInfo <- fileInfo[fileInfo$file %in% bigwigFiles,]
  
  #categorize each file by strandeness (according to file name suffix)
  # "_minus.bw" == minus; "_plus.bw" == plus; neither == unstranded
  fileInfo$strand <- 
    vapply(fileInfo$file, function(fileNameVAP) {
      c("+","-","*")[match(
        sub(".*_(plus|minus)\\.bw","\\1",fileNameVAP, ignore.case = T),
        c("plus","minus"), nomatch=3)]
    }, character(1))
  
  #match up strand requests with available files, reverting to unstranded if requested unavailable
  for( id in unique(fileInfo$id) ){
    #id <- unique(fileInfo$id)[1]
    
    #retrieve requested strand for this group
    strand <- strands[match(id, unique(fileInfo$id))]
    
    #if no files in this group match the request, try again with unstranded
    if( !strand %in% fileInfo[which(fileInfo$id == id),"strand"] ){
      strand <- "*"
      warning("No files in group-strain combination \"",id,"\" matched requested strand (",strand,"). 
                Attempting to use unstranded files.")
      #check again, then fail out if nothing matches
      if( !strand %in% fileInfo[which(fileInfo$id == id),"strand"] ){
        stop("No files in group-strain combination \"",id,"\" matches unstranded file request.
               Please review files and file names.")
      }
    }
    #drop any files in this group that don't have the requested strand (or the unstranded backup request)
    if( any(fileInfo$id == id & fileInfo$strand != strand) ){
      fileInfo <- fileInfo[-which(fileInfo$id == id & fileInfo$strand != strand),]      
    }
  }
  #convert from factor to character
  fileInfo <- as.data.frame( apply(fileInfo, 2, as.character), stringsAsFactors=F )
  
  ### read in files at the requested coordinates
  data <- list()
  for( id in unique(fileInfo$id) ){
    # id <- unique(fileInfo$id)[1]
    for(file in fileInfo[fileInfo$id == id,"file"] ){
      # file <- fileInfo[fileInfo$id == id,"file"][1]
      
      tryCatch({
        fileExtract <- import.bw(file, 
                                 which = GRanges(seqnames = c(chromosome), 
                                                 ranges = IRanges::IRanges( start = left+1, end = right) # +1 to make the query a 1-based coordinate request
                                 ))
      }, 
      error= function(e) stop("Could not import requested data range.\n Error: ", e))
      
      if( length(fileExtract) == 0 ){
        message("Nothing was imported from BigWig file... attempting to retrieve chromosome names.")
        tryCatch({fileExtract <- seqinfo(import.bw(file)) },
                 error= function(e) stop("Could not import BigWig file.\n Error: ", e)
        )
        print(fileExtract)
        stop("Make sure requested chromosome names and ranges to match those above.")
      }
      #apply strand info to the coverage data (since bigwig files are imported unstranded)
      strand(fileExtract) <- fileInfo[fileInfo$id == id & fileInfo$file==file,"strand"]
      
      #drop any seqnames that are not being used (otherwise there will be errors later)
      seqlevels(fileExtract) <- as.character(droplevels(seqnames(fileExtract))[1])
      
      #shove this into a list
      data[[id]][file] <- list(fileExtract)
    }
  }
  
  #find average coverage over all replicates, then 
  # format output data frame of id, group, strain, strand, start, end, covMeans
  #also apply binning if requested 
  dataReps <- list()
  for(id in names(data) ){
    #id <- names(data)[1]
    
    strand <- unique(sapply(data[[id]], function(x) as.character(unique(strand(x))) ))
    if(length(strand) != 1){stop("error tracking strandedness properly (more than one strand after filtering files)")}
    
    covMeans <- 0
    for(rep in 1:length(data[[id]]) ){
      #rep <-1
      covMeans <- covMeans + coverage(data[[id]][[rep]], weight=data[[id]][[rep]]$score)
    }
    covMeans <- covMeans / length(data[[id]])
    
    #'trim' this Rle to just the region of interest with Views function
    covMeans <- Views(covMeans,(left+1):right)[[1]]
    
    #generate output dat.frame
    if( binsPerInch > 0 ){
      #construct the data table using binned averages (100 bins / inch)
      dataOut <- data.frame(
        id = id,
        group = fileInfo[fileInfo$id == id,"group"][1],
        strain = fileInfo[fileInfo$id == id,"strain"][1],
        strand = strand,
        start = aggregate(start(covMeans):end(covMeans), list(cut(start(covMeans):end(covMeans), width*binsPerInch, labels=F)), min)[,2] - 1,
        end   = aggregate(start(covMeans):end(covMeans), list(cut(start(covMeans):end(covMeans), width*binsPerInch, labels=F)), max)[,2],
        covMeans = aggregate(as.vector(covMeans[[1]]), list(cut(start(covMeans):end(covMeans),width*binsPerInch, labels=F)), mean)[,2]
      )
    } else {
      #construct the data table using the cumulative rle length to increment the start and end positions
      dataOut <- data.frame(
        id = id,
        group = fileInfo[fileInfo$id == id,"group"][1],
        strain = fileInfo[fileInfo$id == id,"strain"][1],
        strand = strand,
        start = c(start(covMeans), start(covMeans) + cumsum(runLength(covMeans)[[1]])[-length(runLength(covMeans)[[1]])] ) - 1,
        end = c(start(covMeans) + cumsum(runLength(covMeans)[[1]]) ) - 1,
        covMeans = runValue(covMeans)[[1]]
      )
    }
    
    dataReps[id] <- list(  dataOut  )
  }
  return(dataReps)
}

