library("bpcp")
library("exact2x2")


#makes a pre-allocated (faster) output dataframe
numRows=nrow(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice)
numCols=9

#specify column types
outData<-data.frame(geneA=character(numRows),
                       geneB=character(numRows),
                       AandB=numeric(numRows),
                       AnotB=numeric(numRows),
                       BnotA=numeric(numRows),
                       notAnotB=numeric(numRows),
                       lowerInt=numeric(numRows),
                       upperInt=numeric(numRows),
                       pvalue=numeric(numRows),stringsAsFactors=FALSE)
                       

#creates 2x2 matrix; 1st=AandB, 2nd=AandNotB, 3rd=notAandB, 4th=notAandNotB
#just to illustrate format of genePairs so it is legible:
#genePairs<-matrix(c(dataset$AandB[row], dataset$AnotB[row], dataset$BnotA[row], dataset$notAnotB[row]),
#                 nrow=2, dimnames=list(c(paste(dataset$GeneA[row], "+"),
#                                         paste(dataset$GeneA[row], "-")),
#                                       c(paste(dataset$GeneB[row], "+"),
#                                         paste(dataset$GeneB[row], "-"))))
#the first parameters in the matrix (ie dataset$AandB[row], etc...) will extract an integer values of the number of
#patients that have one of the four possible gene mutation combinations
#nrow=2 will format it into a 2x2 matrix for proper fisher exact test format
#paste will properly concantenate the name of the two genes and whether there is a mutation in one (+) or not (-)



#will iterate through all rows of dataset and put into a 2x2 matrix

for (i in 1:nrow(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice)){
  genePairs<-matrix(c(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$AandB[i],
                      co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$AnotB[i],
                      co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$BnotA[i],
                      co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$notAnotB[i]),
                    nrow=2, dimnames=list(c(paste(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$GeneA[i],"+"),
                                            paste(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$GeneA[i],"-")),
                                          c(paste(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$GeneB[i],"+"),
                                            paste(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$GeneB[i],"-"))))
  
  #calculates the fisher exact test to determine if gene pairs are signicant for mutual exclusion and co-occurence
  calculation<-fisher.exact(genePairs, conf.int=TRUE, conf.level=0.95)
  
  #unlists the confidence interval values of the Fisher's exact test so sub-class attributes can be extracted
  confidenceInterval<-unlist(attributes(calculation$conf.int))
  
  #extracts lower and upper confidence intervals
  lowerInt<-confidenceInterval[2]
  upperInt<-confidenceInterval[3]
  #only extracts p-value from Fisher exact
  pvalue<-fisher.exact(genePairs, conf.int=TRUE, conf.level=0.95)$p.value
  
  #since the gene names are factors, need to convert to characters
  outData[i,]<-c(as.character(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$GeneA[i]),
                 as.character(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$GeneB[i]),
                 as.numeric(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$AandB[i]),
                 as.numeric(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$AnotB[i]),
                 as.numeric(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$BnotA[i]),
                 as.numeric(co.occurence.mutual.exlusion.counts.PDAC.only.minCADD15.starts.stops.gains.losses.nonsyn.splice$notAnotB[i]),
                 lowerInt, upperInt, pvalue)
  }
#counts number of gene pairs that have p-value less than 0.05
with(outData, sum(pvalue < 0.05))

#counts number of gene pairs that have p-value less than 0.01
with(outData, sum(pvalue < 0.01))

#counts number of gene pairs that have p-value less than designated cut-off and have more than 1 patient that has a mutation
pVal=0.01
count=0
for (i in 1:numRows){
  if ( (sum(as.numeric(outData$AandB[i]), as.numeric(outData$AnotB[i]), as.numeric(outData$BnotA[i]))>1 ) & (outData$pvalue[i]<pVal)){
    count<-count+1
    print(paste(outData$geneA[i], outData$geneB[i]))
  }
}
print(count)

