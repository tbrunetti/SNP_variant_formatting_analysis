import sys
import itertools

#sys.argv[1]=tab-delimited txt file where each line represents a single snp from patient
#tab[0]=patientID, tab[7]=gene mutated; same file will suffice for all methods/functions listed

def bothGeneAandGeneB(patientMuts, geneComparisons):
	#output is 2 files, one that records the fequency (not a percent) of co-occurence between two genes,
	#the 2nd file lists all the co-occuring genes per patient
	f=open('co-occuring_gene_fequency.txt', 'w')
	f2=open('patientID_with_co-occurence.txt', 'w')
	
	#iterates through all possible gene comparisons and the number of patients that have mutations in both genes
	for i in range(0, len(geneComparisons)):
		coOccuring=0
		for key in patientMuts:
			if (geneComparisons[i][0] in patientMuts[key]) and (geneComparisons[i][1] in patientMuts[key]):
				coOccuring+=1
				#records the patientID and the co-occurence genes
				f2.write(str(key)+'\t'+str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\n')
		#records the co-occurence genes and the total number of patients where this co-occurence exists
		f.write(str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\t'+str(coOccuring)+'\n')


def geneAandNotGeneB(patientMuts, geneComparisons):
	
	f=open('geneA_and_not_geneB_frequency.txt', 'w')
	f2=open('patientID_with_geneA_but_Not_geneB.txt', 'w')
	
	for i in range(0, len(geneComparisons)):
		#counts occurrences of having geneA and not having geneB
		geneAonly=0
		for key in patientMuts:
			if (geneComparisons[i][0] in patientMuts[key]) and (geneComparisons[i][1] not in patientMuts[key]):
				geneAonly+=1
				#records the patientID and the gene tuple in which patient only has 1st gene in tuple and
				#does not have mutation in 2nd gene in tuple
				f2.write(str(key)+'\t'+ str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\n')

		#records the gene tuple under investigation and the frequency of patients that have a called muatation
		#in just the first gene in the tuple, but not mutation in the 2nd gene listed in the tuple
		f.write(str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\t'+str(geneAonly)+'\n')

def notGeneAandGeneB(patientMuts, geneComparisons):
	
	f=open('not_geneA_but_with_geneB_frequency.txt', 'w')
	f2=open('patientID_Not_geneA_but_with_geneB.txt', 'w')
	
	for i in range(0, len(geneComparisons)):
		#counts only patients that have a mutation in geneB but no mutation in geneA
		geneBonly=0
		for key in patientMuts:
			if (geneComparisons[i][0] not in patientMuts[key]) and (geneComparisons[i][1] in patientMuts[key]):
				geneBonly+=1
				#records the patientID and the gene tuple in which patient only had 2nd gene in tuble and
				#does not have mutation in 1st gene in tuple
				f2.write(str(key)+'\t'+str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\n')

		#records the gene tuple under investigation, and the frequency (number of patients) that have a
		#called mutation in the 2nd gene in the tuple but not in the first gene listed in the tuple
		f.write(str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\t'+str(geneBonly)+'\n')

def notGeneAandNotGeneB(patientMuts, geneComparisons):
	
	f=open('not_geneA_and_not_geneB.txt', 'w')
	f2=open('patientID_Not_geneA_and_Not_geneB.txt', 'w')

	for i in range(0, len(geneComparisons)):
		#counts the number of patients that do not have a mutation in either A nor B
		neitherAorB=0
		for key in patientMuts:
			if (geneComparisons[i][0] not in patientMuts[key]) and (geneComparisons[i][1] not in patientMuts[key]):
				neitherAorB+=1
				#records patientID and the tuple under investigation in which the patient does not have
				#either mutation listed in the tuple
				f2.write(str(key)+'\t'+str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\n')

		#records the gene tuple under investigation and the number of patients that do not have a mutation
		#in either gene A or geneB
		f.write(str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\t'+str(neitherAorB)+'\n')


if __name__=='__main__':
	#a list of all gene mutations of all patients, there will be several duplicates
	genesRepresented=[]
	#key=patientID
	#value=list of genes mutated
	patientMuts={}
	with open(sys.argv[1]) as input:
		#if no headers, this can be commented out
		headers=next(input)
		for line in input:
			line=line.split('\t')
			genesRepresented.append(line[7])
			if line[0] in patientMuts:
				patientMuts[line[0]]=patientMuts[line[0]]+[line[7]]
			else:
				patientMuts[line[0]]=[line[7]]

	#This will generate a list of all the gene comparisons that need to be made		
	geneComparisons=list(itertools.combinations(set(genesRepresented), 2))	
	print "The number of comparisons to make are "+str(len(geneComparisons))
	print "The number of patient genomes analyzed are "+str(len(patientMuts))
	
	#calls made to methods/funnctions
	bothGeneAandGeneB(patientMuts, geneComparisons);
	geneAandNotGeneB(patientMuts, geneComparisons);
	notGeneAandGeneB(patientMuts, geneComparisons);
	notGeneAandNotGeneB(patientMuts, geneComparisons);