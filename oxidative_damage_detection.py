import sys
import numpy as np
from collections import Counter

#FFPE refers to oxidative damage induced by formalin fixation of samples
#FFPE footprint C->T or G->A
def detectFFPE():
	print "unfinished"

#OxoG refers to oxidative damage induced by shearing DNA
#oxoG footprint typically C->A or G->T
#sys.argv[1]=tab delimited list of SNPs in VCF format
def detectOxoG():
	#matrix setup for each nt identifer: (each pos in matrix corresponds to context on both sides of mut)
	#[T_T	T_C		T_A		T_G]
	#[C_T	C_C		C_A		C_G]
	#[A_T	A_C		A_A		A_G]
	#[G_T	G_C		G_A		G_G]
	
	#the two nt identifiers for zero arrays correspond to reference nt change to alt nt
	ct=np.zeros((4,4), dtype=np.int)
	ca=np.zeros((4,4), dtype=np.int)
	cg=np.zeros((4,4), dtype=np.int)
	ag=np.zeros((4,4), dtype=np.int)
	ac=np.zeros((4,4), dtype=np.int)
	at=np.zeros((4,4), dtype=np.int)
	#ct[0][0]=ct[0][0]+5
	#print ct
	#key=mutation, values=5' and 3' nt surrounding mutation
	mutWithContext={'C->T':[], 'C->A':[], 'C->G':[], 'A->G':[], 'A->C':[], 'A->T':[]}
	with open(sys.argv[1]) as input:
		totalMutations=0
		for line in input:
			totalMutations+=1
			line=line.split('\t')
			#for ct matrix
			if (line[2]=='C' and line[3]=='T') or (line[2]=='G' and line[3]=='A'):
				#5' end pre-SNP + 3' end post-SNP
				mutWithContext['C->T']=mutWithContext['C->T']+[str(line[12][2])+str(line[12][4])]			
			
			#for ca matrix
			elif (line[2]=='C' and line[3]=='A') or (line[2]=='G' and line[3]=='T'):
				mutWithContext['C->A']=mutWithContext['C->A']+[str(line[12][2])+str(line[12][4])]
			
			#for cg matrix
			elif (line[2]=='C' and line[3]=='G') or (line[2]=='G' and line[3]=='C'):
				mutWithContext['C->G']=mutWithContext['C->G']+[str(line[12][2])+str(line[12][4])]

			#for ag matrix
			elif (line[2]=='A' and line[3]=='G') or (line[2]=='T' and line[3]=='C'):
				mutWithContext['A->G']=mutWithContext['A->G']+[str(line[12][2])+str(line[12][4])]

			#for ac matrix
			elif (line[2]=='A' and line[3]=='C') or (line[2]=='T' and line[3]=='G'):
				mutWithContext['A->C']=mutWithContext['A->C']+[str(line[12][2])+str(line[12][4])]

			#for at matrix
			elif (line[2]=='A' and line[3]=='T') or (line[2]=='T' and line[3]=='A'):
				mutWithContext['A->T']=mutWithContext['A->T']+[str(line[12][2])+str(line[12][4])]




		print "unfinished"

if __name__=='__main__':
	#detectFFPE();
	detectOxoG();