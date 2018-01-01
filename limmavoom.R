#!/usr/local/apps/R/gcc_4.9.1/3.2.3/bin/Rscript

#if you don't have limma installed, install it.
if(!require('limma')) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
}
if(!require('edgeR')) {
  install.packages('edgeR')
}

# Set threshold FDR and FC values!
MaxFDR = .1
MinAbsFC = 2

# Get Args # RScript --vanilla WD GroupAFiles GroupBFiles Associations GroupAName GroupBName
args = commandArgs(trailingOnly=TRUE)

if(length(args)!=6) stop("Faulty Arguments Supplied\nPlease use the following syntax:\n
RScript --vanilla [YourScript.R] WD GroupAFiles GroupBFiles Associations GroupAName GroupBName
\n Files are comma delimited.")

#set wd
WD = args[1]
setwd(WD)

writeLines('Executing command with args:')
writeLines(paste('\tWorking Directory:',args[1]))
writeLines(paste('\tGroup ',args[5],':',args[2]))
writeLines(paste('\tGroup ',args[6],':',args[3]))
writeLines(paste('\tGene Associations File:',args[4]))


getCountsMatrix = function(groupFiles,groupName){
  #this function reads the bash input for file paths and puts counts data into one counts matrix
  sepGroupFiles = unlist(strsplit(groupFiles, ","))
  group = list()
  for(i in 1:length(sepGroupFiles)) {
    group[[i]] = read.table(sepGroupFiles[i], header=TRUE, sep="\t")
    # Rename counts column to clearly reflect sample name
    names(group[[i]])[7] = paste(groupName, i, sep=" ")
  }
  # Merge all count data into a counts matrix
  for(i in c(1:length(group))[-2]){
    if(i == 1) {
      CountMatrix = merge(group[[i]][,c("Geneid", paste(groupName, i, sep=" "))], group[[i+1]][,c("Geneid", paste(groupName, i+1, sep=" "))], by ="Geneid")
    }
    else {CountMatrix = merge(CountMatrix, group[[i]][,c("Geneid", paste(groupName, i, sep=" "))], by ="Geneid")}
  }
  CountMatrix = cbind(CountMatrix, group[[1]][6])
  return(CountMatrix)
}

#runs the getCountsMatrix function for both groups, using args
groupAName = args[5]
groupBName = args[6]
filePrefix = paste(gsub("[[:punct:] ]+", "",groupAName),'-',gsub("[[:punct:] ]+", "",groupBName))
groupA = getCountsMatrix(args[2], groupAName)
groupB = getCountsMatrix(args[3], groupBName)
counts = merge(groupA, groupB, by=c("Geneid","Length"))
groupASampleCount = ncol(groupA) - 2
groupBSampleCount = ncol(groupB) - 2
#makes the gene ids the rownames of the countsmatrix
row.names(counts) = counts[,1]
counts = counts[,-1]
#writes it to the working directory
write.csv(counts, paste(filePrefix,'_','counts.csv', sep=""))

#split up counts and gene length data
length = counts[,1]
counts = counts[,2:length(counts)]

#reads in the Gene Association file
associations = read.table(args[4], sep='\t', quote="", comment.char="", header=F, skip=24, fill=TRUE)[,c(5,11)]
associations$V11 = as.character(associations$V11)
associations$V11 = sapply(strsplit(associations$V11, "\\|"), function(v) v[1])

##------ LimmaVoom Analysis Starts ------
#the magic happens down here...

#counts <- read.csv("counts.csv",header=T,stringsAsFactors=F,row.names=1) #not necessary since created earlier
filtered_counts <- counts[apply(counts, 1, function(x) any(x>0)),]
group = factor(c(rep("A",groupASampleCount),rep("B",groupBSampleCount)))
design = model.matrix(~0+group)
colnames(design) <- levels(group)

#put MV trend graph in a pdf
pdf(paste(filePrefix,'_','voom.pdf', sep=""))
v = voom(as.matrix(filtered_counts), design, plot=TRUE, normalize="quantile")
dev.off()

fit <- lmFit(v,design)
contrast.matrix <- makeContrasts(contrasts="A-B",levels=design)
fit <- contrasts.fit(fit, contrast.matrix)
ebayes.fit=eBayes(fit)
results = topTable(ebayes.fit, coef="A-B", number=nrow(ebayes.fit))
results$FC <- ifelse(results$logFC<0, -1/(2^results$logFC), 2^results$logFC)
write.csv(results, paste(filePrefix,'_','voom.csv', sep=""))


##----- Get GO Associations ----

GO_all = c()
for(i in 1:nrow(results)) {
  gene_name = row.names(results)[i]
  GO = associations[(associations$V11==gene_name & !is.na(associations$V11)),1]
  GO_all[i] = paste(unique(GO),collapse="; ") #for one particular gene...
}
results=cbind(results, GO_Enrichment = GO_all)

##----- Get RPKM ------------

RPKM = rpkm(counts, length, log=FALSE, normalized.lib.sizes=TRUE)
write.csv(RPKM, paste(filePrefix,'_','rpkm.csv', sep=""))
results = cbind(RPKM[match(rownames(results), rownames(RPKM)),1:ncol(RPKM)],results)

##----- Print all_res ----
write.csv(results, paste(filePrefix,'_','all.csv', sep=""))

##----- Keep only results with FDR above threshold ------
sig = all[which(all$adj.P.Val<MaxFDR),]

##----- Keep only results with fold change above threshold
up = sig[which(sig$FC< (-1*MinAbsFC)),]
dn = sig[which(sig$FC> MinAbsFC),]

##----- Print significantly upregulated and downregulated results -------
write.csv(up, paste(filePrefix,'_','up.csv', sep=""))
write.csv(dn, paste(filePrefix,'_','dn.csv', sep=""))

##----- Save the workspace image in case we need to troubleshoot later -------
save.image(paste(filePrefix,'_',"image.RData", sep=""))
