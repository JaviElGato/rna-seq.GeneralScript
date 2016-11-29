#USAGE: bsub -q normal -G ichip -J allPC -R "select[mem>20000] rusage[mem=20000]" -M20000 -o /lustre/scratch115/projects/pneumo-inf/TranscriptCountFromfastqFiles.kallisto/allPCs.out -e /lustre/scratch115/projects/pneumo-inf/TranscriptCountFromfastqFiles.kallisto/allPCs.err /software/R-3.3.0/bin/R CMD BATCH --vanilla '--args Arg1="Value1" Arg2="Value2" Arg3="Value3" Arg4="Value4"' /FolderWithRscript/rnaseq.generalTools.R  /FolderForRout/file.Rout


#### Example: for i in {2..22};do bsub -q normal -G ichip -J chr$i -R "select[mem>2000] rusage[mem=2000]" -M2000 -o chr$i.out -e chr$i.err /software/R-3.1.2/bin/R CMD BATCH --vanilla '--args toModif="/lustre/scratch113/projects/crohns/iibdgc_meta/IIBDGC/ibd/'$i'.assoc" pattern="/lustre/scratch113/projects/crohns/iibdgc_meta/GWAS3/ibd/'$i'.assoc" output="'$i'" which="only1"' /nfs/team152/javi/ToolBoox/ScriptsR/equalizeAllelesFreq.Katie.R equalizeAlleles.Katie.Chr$i.Rout;done

# usage for txImportData: 
# /software/R-3.3.0/bin/R CMD BATCH --vanilla '--args '
# module="tximport" 
# pathToDir="path/to/dir/with/files" 
# pathTosampleFile="listOfSamples" 
# typeOfTranscript="kallisto/salmon" 
# t2g="transcript2geneFile" 
# output="outputFile"
#R CMD BATCH --vanilla '--args module="tximport" pathToDir="/Users/jga/001_Projects/Project3_LAIV_RNAseq/salmon.rel85.nasal.ExcludeExperimentalCarriers/" pathTosampleFile="../listSeqFiles.MasterFile.nasal.ExcludeExperimentalCarriers" typeOfTranscript="salmon" t2g="../transcript2gene.GRCh38.rel85.Modif.Transc_GeneEns.txt" output="txDataTest"' ~/001_Projects/ToolBox/rna-seq.GeneralScript/rnaseq.generalTools.WithArguments.R rnaseq.generalTools.WithArguments.Rout

# usage for dgeAnalysis:
# /software/R-3.3.0/bin/R CMD BATCH --vanilla '--args '
# module="dgeAnalysis"
# dataFromTxImportData(counts Matrix if this function is used alone)
# pathTosampleFile
# offset=c("txImport, calcNormFactors", "both", "none")
# cpmToFilter=5
# samplesToFilter=4
# package=c("edgeR, deseq2, limma")
# CountToRemove=0
# meanCounts=0
# outplots=NULL
# method=c("lm, robust")
# rounds = 0
# contrast=FALSE
# onlyVoom=FALSE
# plots=FALSE
# allSample="all"
# output="outputFile"

library(ggplot2)
library(batch)
#library(sleuth)
library(tximport)
library(readr)
library(biomaRt)
library(edgeR)
library(limma)
library(EDASeq)
library(statmod)
library(NOISeq)
library(data.table)
library(dtplyr)
library(DESeq2)
library(pcaExplorer)
library(pheatmap)
library(RColorBrewer)
library(calibrate)
# source("/Users/jga/001_Projects/ToolBox/R/qq_plot_javi.R")
# args=(commandArgs(TRUE))

# Determine default arguments
# fullAnalysis = "Yes"
transcript2gene = "/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/transcript2gene.txt" # rel 79
transcript2gene.85 = "/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/transcript2gene.GRCh38.rel85.txt" # rel 85
genesGC = "/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/genesGC.txt"
genesLength = "/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/genesLength.txt"
genesBiotype = "/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/genesBiotype.txt"
genesChrStartEnd = "/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/genesChrStartEnd.txt"
dispersionAnalysis = "Yes"
onlyNOISeq = "Yes"



args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
        print("No arguments supplied.")
        ##supply default values
        break()
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

print(module) 
# pathToDir="" 
# pathTosampleFile="" 
# typeOfTranscript="" 
# t2g="" 
# output=""


### This part was writen at the beginning, but now I have the files so it's easy/fast to load them. I'll keep it in the case of updated data needed
# print("Retreiving data from Ensembl")
# mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = 'ensembl.org')
# print("Retreiving t2g")
# t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",  "external_gene_name", "percentage_gc_content"), mart = mart) entrezgene
# t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, gcContent= percentage_gc_content )
# print("Retreiving biotypes")
# genesBiotype = biomaRt::getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=mart)

# print("Loading data info: Genes %GC, length, biotype")
# transcript2gene = fread(transcript2gene)
# print("transcript2gene loaded")
# genesGC = fread(genesGC)
# print("genesGC loaded")
# genesLength = fread(genesLength)
# print("genesLength loaded")
# genesBiotype = fread(genesBiotype)
# print("genesBiotype loaded")
# genesChrStartEnd = fread(genesChrStartEnd)
# print("genesChrStartEnd loaded")
# print("Done loading data")


## import data TXIMPORT
txImportData = function(pathToDir, pathTosampleFile, typeOfTranscript, t2g, toReturn=c("counts", "tpm", "rpkm")){
  print(paste("Importing", typeOfTranscript, "data using tximport from", pathToDir, sep=" "))

  if(is.data.frame(pathTosampleFile)){
    print("Great...  Sample file is already a dataframe")
  } else {
  pathTosampleFile = fread(pathTosampleFile)
  }

  if(is.data.frame(t2g)){
    print("Great...  T2G is already a dataframe")
  } else {
  t2g = fread(t2g)
  }

  # t2g = t2g[,!c("ext_gene", "gcContent"),with=FALSE]
  if (typeOfTranscript == "kallisto"){
  files =  file.path(pathToDir, pathTosampleFile$Sample_ID2, "abundance.tsv") ### Check this pathTosampleFile$sample_ID1 for the Column name
  } else if (typeOfTranscript == "salmon") {
  files =  file.path(pathToDir, pathTosampleFile$Sample_ID1, "quant.sf")
  }
  
  KallistoTximport = tximport(files, type = typeOfTranscript, tx2gene = t2g, reader = read_tsv)
  colnames(KallistoTximport$counts) = pathTosampleFile$Sample_ID1
  return(KallistoTximport)
}

## import data counts Matrix for DGE analysis
dgeImportData = function(dataFromTxImportData, pathTosampleFile, offset=c("txImport, calcNormFactors", "both", "none"), cpmToFilter=5, samplesToFilter=4, package=c("edgeR, deseq2, limma"), CountToRemove=0, meanCounts=0, outplots=NULL, method="", rounds = 0, contrast=FALSE, onlyVoom=FALSE, plots=FALSE, allSample="all"){
  print("Importing tximport into edgeR and calculate offset")
  #Calculate offset 1
  if(is.data.frame(pathTosampleFile)){
    print("Great... Sample file is already a dataframe")
  } else {
  pathTosampleFile = fread(pathTosampleFile)
  }
  print(head(pathTosampleFile))

  if(is.matrix(dataFromTxImportData)){
    print("Great... The counts are already a matrix")
    counts = dataFromTxImportData
    } else {
    counts = dataFromTxImportData$counts
  }
  print(head(counts))
  dataFromTxImportData.bk = dataFromTxImportData

if (CountToRemove > 0 ){
  print("Before filtering")
  print(dim(counts))
  print(paste("Extracting genes with more than", CountToRemove, "Counts", sep=" "))
  print("After filtering")
  print(dim(counts))
} else if (meanCounts > 0){
  print("Before filtering")
  print(dim(dataFromTxImportData$counts))
  filter = apply(dataFromTxImportData$counts,1,function(x) mean(x)>10)
  dataFromTxImportData$counts = as.matrix(dataFromTxImportData$counts[filter,])
  dataFromTxImportData$abundance = as.matrix(dataFromTxImportData$abundance[filter,])
  dataFromTxImportData$length = as.matrix(dataFromTxImportData$length[filter,])
  print("After filtering")
  print(dim(dataFromTxImportData$counts))
  } else{
  print(paste("Extracting genes with more than", cpmToFilter, "in at least", samplesToFilter, "samples", sep=" "))
  print("Before filtering")
  print(dim(counts))
  counts = counts[rowSums(cpm(counts)>cpmToFilter) >= samplesToFilter, ]
  print("After filtering")
  print(dim(counts))
  dataFromTxImportData$counts = counts[rowSums(cpm(counts)>cpmToFilter) >= samplesToFilter, ]
  GenesAfterfilter = rownames(counts)
  counts = counts[counts = GenesAfterfilter, ]
  print(class(counts))
}
  #### Commented for Pickrell Analysis
  # lengthGenes = dataFromTxImportData$length
  # lengthGenes = lengthGenes[lengthGenes = GenesAfterfilter, ]
  # dataFromTxImportData$length = lengthGenes
  #### Commented for Pickrell Analysis

  if ( package == "edgeR" && offset == "calcNormFactors"){
    print("converting to edgeR with calcNormFactors")
    predesignMatrix = factor(paste(pathTosampleFile$Tissue, pathTosampleFile$Condition, sep='.'))
    designMatrix = model.matrix(~0 + Sample_ID5 + Tissue + Condition + Tissue:Condition, pathTosampleFile)
    print(head(designMatrix))
    print(rowSums(designMatrix))

    forDgeData = DGEList(counts)
    forDgeData = calcNormFactors(forDgeData)
    forDgeData = estimateDisp(forDgeData, designMatrix, robust=TRUE)
    print(class(forDgeData))

  } else if ( package == "edgeR" && offset == "both"){
    normMat = lengthGenes/exp(rowMeans(log(lengthGenes)))
    offset = log(calcNormFactors(counts/normMat)) + log(colSums(counts/normMat))
    forDgeData = DGEList(counts)
    forDgeData$offset = t(t(log(normMat)) + offset)
    forDgeData = calcNormFactors(forDgeData, method="upperquartile")

  } else if (package == "deseq2" ){
    print("converting to DESeq2")
    # coldata = data.frame(Gender = pathTosampleFile$Gender)
    coldata = data.frame(Condition = pathTosampleFile$Condition)
    # coldata = data.frame(Tissue = pathTosampleFile$Tissue, Condition = pathTosampleFile$Condition )
    # coldata = data.frame(Tissue = pathTosampleFile$Tissue,  Condition = pathTosampleFile$Condition, InterGroup = pathTosampleFile$sample_ID2)
    # rownames(coldata) = pathTosampleFile$IDs
    # print(head(coldata))
    rownames(coldata) = pathTosampleFile$Sample_ID3
    # print(head(dataFromTxImportData$counts))
    # print(head(dataFromTxImportData$length))
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~ 0+InterGroup+Condition)
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~ Tissue)
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~ Gender)
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~  Tissue + Tissue:Condition)
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~ Tissue+Condition)
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~ Tissue:Condition)   => Not possible to construct:
          #   Error in checkFullRank(modelMatrix) :
          # the model matrix is not full rank, so the model cannot be fit as specified.
          # One or more variables or interaction terms in the design formula are linear
          # combinations of the others and must be removed.
          # Please read the vignette section 'Model matrix not full rank':
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~ 0+Tissue+Condition+Tissue:Condition) # So called "Full Model"
    # forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~ Tissue + Tissue:InterGroup + Tissue:Condition)
    forDgeData = DESeqDataSetFromTximport(dataFromTxImportData, colData=coldata, design= ~Condition)
    # forDgeData = DESeqDataSetFromMatrix(counts, colData=coldata, design= ~ Gender) ### Used for Pickrell_Montgomery_Data analysis

  } else if (package == "limma" && offset == "none") {
    if (plots){
      pdf(paste(outplots,".pdf",sep=""))
    } else {
      print("Not Plotting")
    }
    print("starting Limma")

    if(allSample == "all") {
      #### Matrix for full Data ########################################
      predesignMatrix = factor(paste(pathTosampleFile$Tissue, pathTosampleFile$Condition, sep='.'))
      # predesignMatrix = factor(paste(pathTosampleFile$Tissue, pathTosampleFile$CarrierStatusAll, sep='.'))
      # designMatrix = model.matrix(~ 0 + predesignMatrix)
      # designMatrix = model.matrix(~ Tissue + Condition, pathTosampleFile)
      # print(head(designMatrix))
      # colnames(designMatrix) = c("Blood.Minus5","Blood.Plus2","Nasal.Minus5","Nasal.Plus2")
      # designMatrix = model.matrix(~ CarrierStatusAll, pathTosampleFile)
      # designMatrix = model.matrix(~ Tissue + CarrierStatusAll, pathTosampleFile)
      # designMatrix = model.matrix(~ Condition + CarrierStatusAll, pathTosampleFile)
      designMatrix = model.matrix(~ Tissue + Condition + CarrierStatusAll, pathTosampleFile)
      # designMatrix = model.matrix(~ Tissue + Condition + CarrierStatusNatural, pathTosampleFile)
      # designMatrix = model.matrix(~ Tissue + Condition + CarrierStatusExperimental, pathTosampleFile)
      print(head(designMatrix))
      forDgeData = DGEList(dataFromTxImportData$counts, group=predesignMatrix)
      print(head(forDgeData$samples))
      ### close Matrix for full Data ########################################
      plotMDS(forDgeData, labels=1:56, col=as.numeric(as.factor(forDgeData$sample$group)), main="MDS plot vaccines pilot - Blood+Nasal (Exclude Exp. carriers)")
      # plotMDS(forDgeData, labels=1:22, col=as.numeric(as.factor(forDgeData$sample$group)), main="MDS plot vaccines pilot - Blood+Nasal Only Exp. Carriers")
      # plotMDS(forDgeData, labels=1:6, col=as.numeric(as.factor(forDgeData$sample$group)), main="MDS plot vaccines pilot - Blood+Nasal Only Nat. Carriers")
      # legend("topleft", legend=c("Blood Minus5", "Blood Plus2", "Nasal Minus2", "Nasal Plus2"), pch=15, col=1:4)
      legend("topright", legend=c("Blood Minus5", "Blood Plus2", "Nasal Minus2", "Nasal Plus2"), pch=15, col=1:4)

    } else {
      ### Matrix for Nasal or Blood subset ########################################
      predesignMatrix = factor(paste(pathTosampleFile$Condition, pathTosampleFile$CarrierStatusAll, sep='.'))
      # designMatrix = model.matrix(~ CarrierStatusAll, pathTosampleFile)
      designMatrix = model.matrix(~ CarrierStatusAll, pathTosampleFile)
      # designMatrix = model.matrix(~ Condition + CarrierStatusAll, pathTosampleFile)
      # designMatrix = model.matrix(~ Condition + CarrierStatusNatural, pathTosampleFile)
      # designMatrix = model.matrix(~ Condition + CarrierStatusExperimental, pathTosampleFile)
      # colnames(designMatrix) = c("Minus5","Plus2")
      print(head(designMatrix))
      forDgeData = DGEList(dataFromTxImportData$counts, group=predesignMatrix)
      # forDgeData$sample$group = pathTosampleFile$CarrierStatusAll
      # forDgeData$sample$group = pathTosampleFile$Condition
      ### Close Matrix for Nasal or Blood subset ########################################
      if(allSample == "nasal") {
      plotMDS(forDgeData, labels=1:24, col=as.numeric(as.factor(forDgeData$sample$group)), main="MDS plot vaccines pilot Nasal (Excluding Exp. Carriers)") #Excluding Nat. Carriers #Only Exp. Indvs
      # legend("topleft", legend=c("M5-NonCarrier", "P2-NonCarrier", "P2-Carrier"), pch=15, col=c("black", "green", "red"))
      legend("bottomright", legend=c("M5-NonCarrier", "P2-NonCarrier", "M5-Carrier", "P2-Carrier"), pch=15, col=c("red", "blue", "black", "green"))
      # legend("bottomright", legend=c("M5-NonCarrier", "P2-Carrier"), pch=15, col=c("black","red"))
    } else if(allSample == "blood"){
      plotMDS(forDgeData, labels=1:32, col=as.numeric(as.factor(forDgeData$sample$group)), main="MDS plot vaccines pilot Blood (Excluding Exp. Carriers)") # Only Experimental Indivs
      legend("topright", legend=c("M5-NonCarrier", "P2-NonCarrier", "M5-Carrier", "P2-Carrier"), pch=15, col=c("red", "blue", "black", "green"))
      # legend("bottomright", legend=c("M5-NonCarrier", "P2-NonCarrier", "P2-Carrier"), pch=15, col=c("black", "green", "red"))
      # legend("bottomright", legend=c("M5-NonCarrier", "P2-Carrier"), pch=15, col=c("black","red"))
      }
    }

    print("calcNormFactors")
    # forDgeData.dgelist = calcNormFactors(forDgeData, method="TMM")
    forDgeData.dgelist = calcNormFactors(forDgeData)
    print(head(forDgeData.dgelist$samples))
    # nf = calcNormFactors(forDgeData)


    if ( rounds == 1 ){
      print("voom 1")

      if(onlyVoom){
      print("Only voom normalization (without lmFit and eBayes)")
      # forDgeData <- voom(forDgeData, designMatrix, plot=TRUE)
      forDgeData.voom <- voomWithQualityWeights(forDgeData.dgelist, designMatrix, plot=TRUE,  normalization="none")
      dev.off()
      return(forDgeData.voom)

      } else {
      # forDgeData <- voom(forDgeData, designMatrix, plot=TRUE)
      forDgeData.voom <- voomWithQualityWeights(forDgeData.dgelist, designMatrix, plot=TRUE,  normalization="none")
      # forDgeData <- voom(dataFromTxImportData$counts, designMatrix, plot=TRUE,  lib.size=colSums(dataFromTxImportData$counts)*nf)
      print("Calc Dups")
      forDgeData.dups = duplicateCorrelation(forDgeData.voom, designMatrix, block=pathTosampleFile$Sample_ID5)
      print(forDgeData.dups$consensus)

      forDgeData2.voom = forDgeData.voom
      forDgeData.dups2 = forDgeData.dups
    }

  } else if (rounds == 2) {
    print("voom 1")
    # forDgeData <- voom(forDgeData, designMatrix, plot=TRUE)
    forDgeData.voom <- voomWithQualityWeights(forDgeData.dgelist, designMatrix, plot=TRUE,  normalize.method="none")
    # forDgeData <- voom(dataFromTxImportData$counts, designMatrix, plot=TRUE,  lib.size=colSums(dataFromTxImportData$counts)*nf)
    print("Calc Dups")
    forDgeData.dups = duplicateCorrelation(forDgeData.voom, designMatrix, block=pathTosampleFile$Sample_ID5)
    print(forDgeData.dups$consensus)

    # Second round voom and duplicateCorrelation
    print("voom 2")
    forDgeData2.voom <- voomWithQualityWeights(forDgeData.dgelist, designMatrix, plot=FALSE, block=pathTosampleFile$Sample_ID5, correlation=forDgeData.dups$consensus)
    print("Calc Dups No. 2")
    forDgeData.dups2 = duplicateCorrelation(forDgeData2.voom, designMatrix, block=pathTosampleFile$Sample_ID5)
    print(forDgeData.dups2$consensus)
  }

    if(method == "ls") {
      print("Fitting the Model ls")
    forDgeData.fit = lmFit(forDgeData2.voom, designMatrix, block=pathTosampleFile$Sample_ID5, correlation=forDgeData.dups2$consensus, method="ls")
  } else if (method == "robust") {
      print("Fitting the Model robust")
    forDgeData.fit = lmFit(forDgeData2.voom, designMatrix, block=pathTosampleFile$Sample_ID5, correlation=forDgeData.dups2$consensus, method="robust")
  }
    # forDgeData.fit = lmFit(forDgeData, designMatrix, method="robust")

    ########################################
    if (contrast){
    contrastToTest =  makeContrasts(diffPlus2Minus5 = (Blood.Plus2-Blood.Minus5)-(Nasal.Plus2-Nasal.Minus5), levels=designMatrix)
    # contrastToTest =  makeContrasts(diffPlus2Minus5 = (Plus2-Minus5), levels=designMatrix)
    ########################################

    print("Fitting the Model No. 2")
    forDgeData.fit2 = contrasts.fit(forDgeData.fit, contrastToTest)
    forDgeData.ebayes = eBayes(forDgeData.fit2, robust=T)
    # forDgeData.ebayes = treat(forDgeData.fit2, robust=T)

    # return(forDgeData.dup)
  } else {
    forDgeData.ebayes = eBayes(forDgeData.fit, robust=T)
  }

  print(summary(decideTests(forDgeData.ebayes)))
    # qqt(forDgeData.ebayes$t, df=forDgeData.ebayes$df.prior + forDgeData.ebayes$df.residual, pch=19,cex=0.2)
  # abline(0,1)

  # forDgeData.ebayes.df = toptable(forDgeData.ebayes, coef="ConditionPlus2", n=Inf )
  forDgeData.ebayes.df = toptable(forDgeData.ebayes, coef="CarrierStatusAllNonCarrier", n=Inf )
  # forDgeData.ebayes.df = toptable(forDgeData.ebayes, coef="CarrierStatusNaturalNonCarrier", n=Inf )
  # forDgeData.ebayes.df = toptable(forDgeData.ebayes, coef="CarrierStatusExperimentalNonCarrier", n=Inf )
  forDgeData.ebayes.df$genes = rownames(forDgeData.ebayes.df)
  forDgeData.ebayes.AverExp = data.frame(AverLogExp=forDgeData.ebayes$Amean, genes=names(forDgeData.ebayes$Amean))
  forDgeData.ebayes.df.tmp = merge(forDgeData.ebayes.AverExp, forDgeData.ebayes.df, by="genes")
  p = ggplot(forDgeData.ebayes.df.tmp, aes(AverLogExp, logFC, colour=(adj.P.Val <= 0.05))) + geom_point() + theme(legend.position="none")
  print(p)
  # qqplotFunction(toptable(forDgeData.ebayes, n=Inf, coef="ConditionPlus2")$P.Value, outplots)
  qqplotFunction(toptable(forDgeData.ebayes, n=Inf, coef="CarrierStatusAllNonCarrier")$P.Value, outplots)
  # qqplotFunction(toptable(forDgeData.ebayes, n=Inf, coef="CarrierStatusNaturalNonCarrier")$P.Value, outplots)
  # qqplotFunction(toptable(forDgeData.ebayes, n=Inf, coef="CarrierStatusExperimentalNonCarrier")$P.Value, outplots)
  dev.off()
  return(forDgeData.ebayes.df)

    # forDgeData <- voom(forDgeData, designMatrix, plot=TRUE,  lib.size=colSums(counts)*nf$samples$norm.factors)
    # forDgeData = calcNormFactors(forDgeData)
    # forDgeData <- voom(forDgeData, designMatrix, plot=TRUE)


  } else if(offset == "txImport"){
    print("converting to edgeR with offset and calcNormFactors")
    normMat = lengthGenes/exp(rowMeans(log(lengthGenes)))
    offset = log(calcNormFactors(counts/normMat)) + log(colSums(counts/normMat))
    forDgeData = DGEList(counts)
    forDgeData$offset = t(t(log(normMat)) + offset)

  }
  # else {
  #   print("You might like to normalize your data somehow")
  # }

  return(forDgeData)

}


###################
###################
###################
# DGE analysis with edgeR - UNDER CONSTRUCTION - ADD THE SECTION SHOULD BE RE-CHECKED ###################

testCPM = function(output, dataFromTxImportData, samplesToFilter, filterRange){
  # Quick evaluation to determine the threshold of genes to remove because of low counts
  pdf(paste("DensityLogCPM.filters.", output ,".pdf", sep=""))
  for(filter in 1:filterRange){
  edgeRData.all = dgeImportData(dataFromTxImportData, cpmToFilter=filter , offset="both",package="edgeR", samplesToFilter=samplesToFilter)
  counts = edgeRData.all$counts
  counts.cpm = cpm(counts)
  # counts = counts[rowSums(cpm(counts)>filter) >= 4, ]
  toKeep = rowSums(counts.cpm > filter) >= 4
  edgeRData.all = edgeRData.all[toKeep, , keep.lib.sizes=FALSE]
  lcounts = cpm(edgeRData.all$counts, log=TRUE)
  genesleft = dim(lcounts)[1]
  plot(density(lcounts[,1]), col=col[1], lwd=2, las=2, ylim=c(0, 0.3), xlab="log(CPM)", main=paste("filter CPM >", filter, "Genes left =", genesleft, sep=" "))
  for (i in 2:67){den = density(lcounts[,i]); lines(den$x, den$y, col=col[i], lwd=2)}
  }
  dev.off()
}

makePlots = function(DgeData, TopGenes=NULL, output){
  DgeData.res = results(DgeData, alpha=0.05)
  # DgeData.res = results(DgeData, alpha=0.05, contrast=list("TissueBlood.ConditionPlus2", "TissueNasal.ConditionPlus2"))
  print("res Done")
  DgeData.resOrdered = DgeData.res[order(DgeData.res$padj),]
  print("resOrdered Done")
  DgeData.resOrdered.rownanes = rownames(DgeData.resOrdered)
  DgeData.resOrdered$Genes = DgeData.resOrdered.rownanes
  print("resOrdered Genes Done")
  print(head(DgeData.resOrdered))
  DgeData.rownames = rownames(DgeData.resOrdered[DgeData.resOrdered$padj < 0.05,])
  # DgeData.select = order(rowMeans(counts(DgeData, normalized=T)), decreasing=T)[1:5] # In case of the most expressed Genes
  DgeData.select = c()
  for(i in 1:length(DgeData.rownames)){DgeData.select = c(DgeData.select,match(DgeData.rownames[i], rownames(counts(DgeData))))}
  print(DgeData.select)
  if(length(deseqData.Nasal.Condition.select) > 20){
    deseqData.Nasal.Condition.select = deseqData.Nasal.Condition.select[1:20]
  }
  DgeData.nt = normTransform(DgeData)
  DgeData.log2_norm_counts = assay(DgeData.nt)[DgeData.select,]
  DgeData.df = as.data.frame(colData(DgeData)[,c("Condition", "Tissue")])
  print("VST Done")
  DgeData.vcd = varianceStabilizingTransformation(DgeData, blind=F)
  print("Dist Done")
  DgeData.sampleDist = dist(t(assay(DgeData.vcd)))
  DgeData.sampleDistMatrix = as.matrix(DgeData.sampleDist)
  rownames(DgeData.sampleDistMatrix) = paste(DgeData.vcd$Tissue, DgeData.vcd$InterGroup, sep=".")
  colnames(DgeData.sampleDistMatrix) = NULL
  colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  write.table(DgeData.resOrdered, file=paste(output,".txt", sep=""), quote=F, row.names=T, col.names=T)
  #
  #
  pdf(paste(output,".plots.pdf", sep=""))
  DESeq2::plotMA(DgeData.res)
  plotDispEsts(DgeData)
  # qqplotFunction(DgeData.res$pvalue)
  # qqplotFunction(DgeData.res$padj)

  # plotCounts(DgeData, gene="ENSG00000211959.2", intgroup=c("Tissue","Condition"), pch=19, main="ENSG00000211959.2 - IGHV4-39 - pval=7.21e-11 - padj=1.13e-06")
  # plotCounts(DgeData, gene="ENSG00000282960.1", intgroup=c("Tissue","Condition"), pch=19, main="ENSG00000282960.1 - AL513412.1 (KDM4C iso)- pval=5.72e-06 - padj=0.044")
  # plotCounts(DgeData, gene="ENSG00000131778.17", intgroup=c("Tissue","Condition"), pch=19, main="ENSG00000131778.17 - CHD1L - pval=1.64e-05 - padj=0.086")
  # plotCounts(DgeData, gene="ENSG00000107485.15", intgroup=c("Tissue","Condition"), pch=19, main="ENSG00000282657.2 - GATA3 - pval=3.06e-05 -  padj=0.099")
  # plotCounts(DgeData, gene="ENSG00000282657.2", intgroup=c("Tissue","Condition"), pch=19, main="ENSG00000282657.2 - IGHM - pval=3.18e-05 - padj=0.099")

  pheatmap(DgeData.log2_norm_counts, cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col=DgeData.df)
  pheatmap(DgeData.log2_norm_counts, cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col=DgeData.df)
  pheatmap(DgeData.log2_norm_counts, cluster_rows=F, show_rownames=T, cluster_cols=T, annotation_col=DgeData.df)
  pheatmap(DgeData.log2_norm_counts, cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col=DgeData.df)
  pheatmap(DgeData.sampleDistMatrix, clustering_distance_rows=DgeData.sampleDist, clustering_distance_cols=DgeData.sampleDist, col=colors)
  # Make a basic volcano plot
  with(DgeData.resOrdered, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot")) #xlim=c(-2.5,2.5)
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(DgeData.resOrdered, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
  with(subset(DgeData.resOrdered, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  # with(subset(DgeData.resOrdered, padj<.1), textxy(log2FoldChange, -log10(pvalue), labs=Genes, cex=.8))

  dev.off()
}

getDEnamesEdgeR = function(edgeRobject, type=c("lrt, qlf")) {
  if(type == "lrt") {
  decideTest = decideTestsDGE(object=paste(edgeRobject,".lrt", sep=""))
} else if (type == "qlf"){
    decideTest = decideTestsDGE(object=paste(edgeRobject,".qlf", sep=""))
}
  isDE = as.logical(decideTest)
  DEnames = rownames(edgeRobject)[isDE]
  return(DEnames)
}


rpkm2tpm = function(data){
  if(class(data[,1]) == "character"){
    geneNames = data[,1]
    colnames = data
    dataToWork = data[,c(-1)]
  } else {
    dataToWork = data
  }
  numberRows = dim(dataToWork)[1]
  numberCols = dim(dataToWork)[2]
  matrixToFill = matrix(c(0), nrow=numberRows , ncol=numberCols)
  vectorOfSums = colSums(dataToWork)
  for(cols in 1:numberCols){
    matrixToFill[,cols] = dataToWork[,cols]/vectorOfSums[cols]*10^6
  }

  if(class(data[,1]) == "character"){
    rownames(matrixToFill) = data[,1]
    colnames(matrixToFill) = data[,2:numberCols]
  }

  return(matrixToFill)
}

# QQ plot

qqplotFunction = function(pvalues, title){
  PVAL=-log10(pvalues)
  N <- length(PVAL) ## number of p-values

  ## create the null distribution
  ## (-log10 of the uniform)
  null <- -log(1:N/N,10)
  MAX <- max(c(PVAL,null),na.rm=T)


  ## create the confidence intervals
  c95 <- rep(0,N)
  c05 <- rep(0,N)

  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)

  for(i in 1:N){
    c95[i] <- qbeta(0.95,i,N-i+1)
    c05[i] <- qbeta(0.05,i,N-i+1)
  }

  ## plot the two confidence lines
  plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",lty=2,
  axes=FALSE, xlab="", ylab="")
  par(new=T)
  plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", lty=2,
  axes=FALSE, xlab="", ylab="")
  ## add the diagonal
  abline(0,1,col="red",lwd=2)
  par(new=T)

  ## add the qqplot
  qqplot(null,PVAL, ylim=c(0,MAX),xlim=c(0,MAX), pch=19,main=title,
  xlab=expression(paste("Theoretical ", -log[10](italic(p)))), ylab= expression(paste("Observed ",-log[10](italic(p)))))

  #END
}

#############################################################################################################################

if (module == "tximport"){

  txDatatemp = txImportData(pathToDir=pathToDir, 
               pathTosampleFile=pathTosampleFile, 
               typeOfTranscript=typeOfTranscript, 
               t2g=t2g)

  print(head(txDatatemp$counts))

  write.table(txDatatemp$counts, file=output, col.names=TRUE, row.names=TRUE, quote=FALSE, sep='\t')

} else if (module == "dgeAnalysis"){

  dgeAnalysisTemp = dgeAnalysis(dataFromTxImportData=dataFromTxImportData, 
              pathTosampleFile=pathTosampleFile, 
              offset=offset, 
              cpmToFilter=cpmToFilter, 
              samplesToFilter=samplesToFilter, 
              package=package, 
              CountToRemove=0, 
              meanCounts=0, 
              outplots=NULL, 
              method=method, 
              rounds = 1, 
              contrast=FALSE, 
              onlyVoom=FALSE, 
              plots=FALSE, 
              allSample="all")

 write.table(dgeAnalysisTemp, file=output, col.names=TRUE, row.names=TRUE, quote=FALSE, sep='\t')

} else if (module == "importAndAnalysis"){
  txDataTemp = txImportData(pathToDir=pathToDir, 
               pathTosampleFile=pathTosampleFile, 
               typeOfTranscript=typeOfTranscript, 
               t2g=t2g)

  dgeAnalysis(dataFromTxImportData=txDataTemp, 
              pathTosampleFile=pathTosampleFile, 
              offset=offset, 
              cpmToFilter=cpmToFilter, 
              samplesToFilter=samplesToFilter, 
              package=package, 
              CountToRemove=0, 
              meanCounts=0, 
              outplots=NULL, 
              method=method, 
              rounds = 1, 
              contrast=FALSE, 
              onlyVoom=FALSE, 
              plots=FALSE, 
              allSample="all")

  write.table(dgeAnalysisTemp, file=output, col.names=TRUE, row.names=TRUE, quote=FALSE, sep='\t')
}

