import sys
from collections import Counter
import pandas
import seaborn as sns

#sys.argv[1] input data format

def makeCluster():
	#key=patID, value=list of predicted protein function loss
	patientMutations={}
	#stores list of all genes (including redundant), predicted for loss of protein function
	numberDiffentGenes=[]
	with open(sys.argv[1]) as input:
		for line in input:
			line=line.split('\t')
			#line[3] is the name of the mutated gene
			numberDiffentGenes.append(line[3])
			#line[0] is the patientID
			if line[0] in patientMutations:
				patientMutations[line[0]]=patientMutations[line[0]]+[line[3]]
			else:
				patientMutations[line[0]]=[line[3]]

	#prints the number of unique genes seen			
	print len(set(numberDiffentGenes))
	#counts the number of times a gene is predicted to be deleterious in predicted protein function
	totalGenesMuts=dict(Counter(numberDiffentGenes))

	#removes all genes that only appear once, to help reduce noise in data set
	removeSingletons={}
	for key in totalGenesMuts:
		if totalGenesMuts[key]<=2:
			continue;
		else:
			removeSingletons[key]=totalGenesMuts[key]

	print len(removeSingletons) 
	genesToConsider=[key for key in removeSingletons]
	#print patientMutations
	#sortByFreq=Counter(removeSingletons).most_common(101)
	#orderToSort=[sortByFreq[x][0] for x in range(0, len(sortByFreq))]

	listOfpatsWithMutCounts=[]
	for key in patientMutations:
		for i in range(0, len(genesToConsider)):
			if genesToConsider[i] in patientMutations[key]:
				entry=[str(key), genesToConsider[i], patientMutations[key].count(genesToConsider[i])]
				listOfpatsWithMutCounts.append(entry)
			else:
				entry=[str(key), genesToConsider[i], 0]
				listOfpatsWithMutCounts.append(entry)

	dataframe=pandas.DataFrame(listOfpatsWithMutCounts, columns=['patientID', 'gene_mutated', 'frequency'])
	#print dataframe
	dataframe=dataframe.pivot("patientID", "gene_mutated", "frequency")
	clusterData=sns.clustermap(dataframe, method='single', metric='euclidean')
	sns.plt.show()

def input_ROCKclustering_in_R():

	print "unfinished"


if __name__=='__main__':
	makeCluster();