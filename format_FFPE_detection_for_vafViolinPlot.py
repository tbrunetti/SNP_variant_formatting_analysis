#output of this script is used as input into vafViolinPlot.R
import sys
import os


#sys.argv[1]=a single column list of exact VCF names to analyze
def FFPE_mutation_detection_format(pathToVCF):
	#creates the final output file
	#can change name/directory as needed
	f=open('/home/tonya/pancan/vaf_by_substitution.txt', 'w')
	os.chdir(pathToVCF)
	#open the list of VCFs to analyze
	with open(sys.argv[1]) as input:
		for line in input:
			line=line.rstrip('\n')
			#to associate each patient with their SNPs
			for i in range(0, len(line))
				patientID=''
				if line[i]=='-':
					break;
				patientID=patientID+line[i]
			#analyzes every SNV called in the patient and formats in the following columns:
			#patientID, vaf, refAllele>altAllele
			with open(line) as input:
				#directory where file 'vaf_by_substitution.txt' is located
				#change as needed
				os.chdir('/home/tonya/pancan/')
				for variant in input:
					variant=variant.split('\t')
					vaf=float(variant[8])/(float(variant[7])+float(variant[8]))
					f.write(str(patientID)+'\t'+str(vaf)+'\t'+str(variant[2])+'>'+str(variant[3])+'\n')
			os.chdir(pathToVCF)

if __name__=='__main__':
	pathToVCF='/home/tonya/pancan/panData_UCretro...'
	FFPE_mutation_detection_format(pathToVCF);