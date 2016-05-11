import sys
import itertools
#sys.argv[1]=tab-delimited txt file where each line represents a single snp from patient
#tab[0]=patientID, tab[7]=gene mutated
def bothGeneAandGeneB():
	#a list of all gene mutations of all patients, there will be several duplicates
	genesRepresented=[]
	#key=patientID
	#value=list of genes mutated
	patientMuts={}
	with open(sys.argv[1]) as input:
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

	f=open('co-occuring_gene_fequency.txt', 'w')
	f2=open('patientID_with_co-occurence.txt', 'w')
	
	#iterates through all possible gene comparisons and the number of patients that have mutations in both genes
	for i in range(0, len(geneComparisons)):
		coOccuring=0
		for key in patientMuts:
			if geneComparisons[i][0] in patientMuts[key] and geneComparisons[i][1] in patientMuts[key]:
				coOccuring+=1
				#records the patientID and the co-occurence genes
				f2.write(str(key)+'\t'+str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\n')
		#records the co-occurence genes and the total number of patients where this co-occurence exists
		f.write(str(geneComparisons[i][0])+'\t'+str(geneComparisons[i][1])+'\t'+str(coOccuring)+'\n')




def geneAandNotGeneB():
	print "non-functional"

def notGeneAandGeneB():
	print "non-functional"

def notGeneAandNotGeneB():
	print "non-functional"

if __name__=='__main__':
	bothGeneAandGeneB();
