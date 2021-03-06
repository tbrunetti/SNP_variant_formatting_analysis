import sys
from collections import Counter
from operator import itemgetter
import pandas
import seaborn as sns

#sys.argv[1] input data format is output from match_and_format_variants_with_CADD
#make sure 0=patientID, 7=gene, 16=CADD score, if all other fields are padded, code should work
def makeCluster():
	#key=patID, value=list of predicted protein function loss
	patientMutations={}
	#stores list of all genes (including redundant), predicted for loss of protein function
	patientMutationsCADD={}
	numberDiffentGenes=[]
	with open(sys.argv[1]) as input:
		headers=next(input)
		for line in input:
			line=line.split('\t')
			#line[7] is the name of the mutated gene
			numberDiffentGenes.append(line[7])
			#line[0] is the patientID
			if line[0] in patientMutations:
				patientMutations[line[0]]=patientMutations[line[0]]+[line[7]]
				patientMutationsCADD[line[0]]=patientMutationsCADD[line[0]]+[(str(line[7]),float(line[16].rstrip('\n')))]
			else:
				patientMutations[line[0]]=[line[7]]
				patientMutationsCADD[line[0]]=[(str(line[7]),float(line[16].rstrip('\n')))]

	for key in patientMutationsCADD:
		print patientMutationsCADD[key]
		#print Counter(patientMutationsCADD[key][0] for key in patientMutationsCADD)


	#prints the number of unique genes seen			
	print len(set(numberDiffentGenes))
	#counts the number of times a gene is predicted to be deleterious in predicted protein function
	totalGenesMuts=dict(Counter(numberDiffentGenes))

	#removes all genes that only appear once, to help reduce noise in data set
	removeSingletons={}
	for key in totalGenesMuts:
		if totalGenesMuts[key]<=1:
			continue;
		else:
			removeSingletons[key]=totalGenesMuts[key]

	print len(removeSingletons) 
	genesToConsider=[key for key in removeSingletons]
	for i in range(len(genesToConsider)):
		print genesToConsider[i]
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
	clusterData=sns.clustermap(dataframe, method='weighted', metric='jaccard')
	sns.plt.show()

	#cluster data similar to above, but hue is now CADD score
	listOfpatsCADD=[]
	for key in patientMutationsCADD:
		for i in range(0, len(genesToConsider)):
			#checks if patient had gene mutation, if not, the bottom else statment will give it a CADD of zero
			if genesToConsider[i] in [patientMutationsCADD[key][tuples][0] for tuples in range(len(patientMutationsCADD[key]))]:
				caddScore=[patientMutationsCADD[key][tuples][1] for tuples in range(len(patientMutationsCADD[key])) if patientMutationsCADD[key][tuples][0]==genesToConsider[i]]
				#if there are multiple mutations in the same gene for a patient, take the highest CADD score
				if len(caddScore)>1:
					entryCADD=[str(key), genesToConsider[i], max(caddScore)]
					print entryCADD
					listOfpatsCADD.append(entryCADD)
				else:
					entryCADD=[str(key), genesToConsider[i], caddScore[0]]
					listOfpatsCADD.append(entryCADD)
			else:
				entryCADD=[str(key), genesToConsider[i], 0]
				listOfpatsCADD.append(entryCADD)
	#print listOfpatsCADD
	dataframeCADD=pandas.DataFrame(listOfpatsCADD, columns=['patientID', 'gene_mutated', 'CADD_score'])
	dataframeCADD=dataframeCADD.pivot('patientID', 'gene_mutated', 'CADD_score')
	clusterCADDdata=sns.clustermap(dataframeCADD, method='single', metric='euclidean', cmap='hot_r')
	sns.plt.show()

#this function properly formats data for use in PCA, MCA, and some clustering algorithms implemented in R
def input_clustering_PCA_MCA_in_R():
	#key=patID, value=list of predicted protein function loss
	patientMutations={}
	#stores list of all genes (including redundant), predicted for loss of protein function
	numberDiffentGenes=[]
	with open(sys.argv[1]) as input:
		for line in input:
			line=line.split('\t')
			#line[7] is the name of the mutated gene
			numberDiffentGenes.append(line[7])
			#line[0] is the patientID
			if line[0] in patientMutations:
				patientMutations[line[0]]=patientMutations[line[0]]+[line[7]]
			else:
				patientMutations[line[0]]=[line[7]]

	#prints the number of unique genes seen			
	print len(set(numberDiffentGenes))
	#counts the number of times a gene is predicted to be deleterious in predicted protein function
	totalGenesMuts=dict(Counter(numberDiffentGenes))

	#removes all genes that only appear once or set int, to help reduce noise in data set
	removeSingletons={}
	for key in totalGenesMuts:
		if totalGenesMuts[key]<=1:
			continue;
		else:
			removeSingletons[key]=totalGenesMuts[key]

	#a list of genes that are being considered for downstream analysis, dependend upon above loop filtering
	genesToConsider=[key for key in removeSingletons]

	f1=open('PDAC-patients-mutations-profile-minCADD_15_nonsyn-splice-gain-loss-stops-starts.txt', 'w')
	f1.write('patientID'+'\t')
	for x in range(0, len(genesToConsider)-1):
		f1.write(genesToConsider[x]+'\t')
	f1.write(genesToConsider[len(genesToConsider)-1]+'\n')

	
	#writes tab delimited file with one row per patient and columns are all possible gene mutations
	#1 for mutation present 0 for no mutation called
	for key in patientMutations:
		f1.write(str(key)+'\t')
		for i in range(0, len(genesToConsider)-1):
			if genesToConsider[i] in patientMutations[key]:
				f1.write('1'+'\t')
			else:
				f1.write('0'+'\t')
		if genesToConsider[len(genesToConsider)-1] in patientMutations[key]:
			f1.write('1'+'\n')
		else:
			f1.write('0'+'\n')


if __name__=='__main__':
	#makeCluster();
	input_clustering_PCA_MCA_in_R();
