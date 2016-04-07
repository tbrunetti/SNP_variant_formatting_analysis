import sys
import os

#sys.argv[1]=exact names of files to be analyzed
def format_data(pathToData):
	f=open('PDAC-patients-nonSyn-stop-codon-changes-SIFT-score-bw-0-and-point-1.txt', 'w')
	f.write('patient_ID'+'\t'+'pct_alt_normal'+'\t'+'pct_alt_tumor'+'\t'+'gene'+'\t'+'mutation_type'+'\t'+'amino_acid_change'+'\t'+'SIFT_protein_call'+'\t'+'SIFT_score'+'\n')
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
					#print snp
					#only selects non-synonymous mutations and stop gained/lost mutations
					if snp[14]=='NON_SYNONYMOUS_CODING' or snp[14]=='STOP_GAINED' or snp[14]=='STOP_LOST':
						if snp[22]=='NA' or snp[22]=='NOT_ENOUGH_SEQUENCES':
							continue;
						#select mutations that have a particular SIFT score range
						print patientID
						if snp[23]=='' or float(snp[23])<=0.10:
							f.write(str(patientID)+'\t'+str(snp[6])+'\t'+str(snp[9])+'\t'+str(snp[11])+'\t'+str(snp[14])+'\t'+str(snp[17])+'\t'+str(snp[22])+'\t'+str(snp[23].rstrip('\n'))+'\n')
			os.chdir(pathToData)

if __name__=='__main__':
	pathToData='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/pathology_dept_formatting'
	format_data(pathToData);