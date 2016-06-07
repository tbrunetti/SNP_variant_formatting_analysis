import sys
import os

#sys.argv[1]=exact names of files to be analyzed (SIFT, PROVEAN, CADD must be present in columns) (files listed, should be tab-delimited files)
def format_data(pathToData):
	f=open('UC-retrospective-patients-nonSyn-stopStart-codon-changes-spliceSites-SIFT-and-or-PROVEAN-deleterious-CADD.txt', 'w')
	f.write('patient_ID'+'\t'+'chr'+'\t'+'pos'+'\t'+'ref'+'\t'+'alt'+'\t'+'pct_alt_normal'+'\t'+'pct_alt_tumor'+'\t'+'gene'+'\t'+'mutation_type'+'\t'+'amino_acid_change'+'\t'+'aa_length'+'\t'+'SIFT_protein_call'+'\t'+'SIFT_score'+'\t'+'PROVEAN_protein_call'+'\t'+'PROVEAN_score'+'\t'+'CADD-raw'+'\t'+'CADD-scaled'+'\n')
	with open(sys.argv[1]) as input:
		os.chdir(pathToData)
		#note: here "line" will be the file name
		for line in input:
			patientID=''
			line=line.rstrip('\n')
			for x in range(0, len(line)):
				if line[x]=='-':
					break;
				patientID=patientID+line[x]
		
			with open(line) as input:
				os.chdir('/home/tonya/pan_can_project/')
				snpCount=0
				for snp in input:
					snpCount+=1
					if snpCount==1:
						continue;
					snp=snp.split('\t')
					#only selects non-synonymous mutations and stop gained/lost mutations
					if snp[14]=='NON_SYNONYMOUS_CODING' or snp[14]=='STOP_GAINED' or snp[14]=='STOP_LOST' or snp[14]=='SPLICE_SITE_REGION' or snp[14]=='START_GAINED':
						#snp[31] is second confidence flag, if 31=='NA' then the confidence value has been resolved
						#all combinations below indicate amino acid substitution could not be accurately assessed
						if snp[22]=='NA' and snp[31]!='NA' and (snp[32]=='NA' or snp[32]=='NOT_ENOUGH_SEQUENCES' or snp[32]=='ERROR_WRONG_AA_CALL' or snp[32]=='ISOFORM_NOT_FOUND'):
							continue;
						#select mutations that have a confident call produced by SIFT and PROVEAN
						print patientID
						#the 'NA' is in f.write because that column is amino acid change which is empty if snp is one of the following effects below
						if snp[14]=='STOP_GAINED' or snp[14]=='STOP_LOST' or snp[14]=='START_GAINED' or snp[14]=='SPLICE_SITE_REGION':
								f.write(str(patientID)+'\t'+str(snp[0])+'\t'+str(snp[1])+'\t'+str(snp[2])+'\t'+str(snp[3])+'\t'+str(snp[6])+'\t'+str(snp[9])+'\t'+str(snp[11])+'\t'+str(snp[14])+'\t'+'NA'+'\t'+str(snp[18])+'\t'+str(snp[22])+'\t'+str(snp[23].rstrip('\n'))+'\t'+str(snp[33])+'\t'+str(snp[32])+'\t'+str(snp[35])+'\t'+str(snp[36]))
						
						if (snp[31]=='NA' or snp[33]!='NA') and (snp[22]=='AFFECT_PROTEIN_FUNCTION' or snp[27]=='AFFECT_PROTEIN_FUNCTION' or snp[33]=='DELETERIOUS'):	
							if snp[26]=='NA' and snp[22]=='AFFECT_PROTEIN_FUNCTION':
								f.write(str(patientID)+'\t'+str(snp[0])+'\t'+str(snp[1])+'\t'+str(snp[2])+'\t'+str(snp[3])+'\t'+str(snp[6])+'\t'+str(snp[9])+'\t'+str(snp[11])+'\t'+str(snp[14])+'\t'+str(snp[17])+'\t'+str(snp[18])+'\t'+str(snp[22])+'\t'+str(snp[23].rstrip('\n'))+'\t'+str(snp[33])+'\t'+str(snp[32])+'\t'+str(snp[35])+'\t'+str(snp[36]))
							elif snp[27]=='TOLERATED':
								continue;
							else:
								f.write(str(patientID)+'\t'+str(snp[0])+'\t'+str(snp[1])+'\t'+str(snp[2])+'\t'+str(snp[3])+'\t'+str(snp[6])+'\t'+str(snp[9])+'\t'+str(snp[11])+'\t'+str(snp[14])+'\t'+str(snp[17])+'\t'+str(snp[18])+'\t'+str(snp[27])+'\t'+str(snp[28].rstrip('\n'))+'\t'+str(snp[33])+'\t'+str(snp[32])+'\t'+str(snp[35])+'\t'+str(snp[36]))
			os.chdir(pathToData)

if __name__=='__main__':
	pathToData='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/match-SIFT-CADD'
	format_data(pathToData);