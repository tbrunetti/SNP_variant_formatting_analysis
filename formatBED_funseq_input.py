import sys
import os

#sys.argv[1] is a list of exact file names to reformat, one file name per line
def formatBed_funseq(pathToVCF):

	with open(sys.argv[1]) as input:
		for line in input:
			os.chdir(pathToVCF)
			with open(line.rstrip('\n')) as input:
				os.chdir(pathToOutput)
				f=open(str(line[:-5].rstrip('\n'))+'-'+'input-funseq.txt', 'w')
				for variant in input:
					variant=variant.split('\t')
					f.write(str(variant[0])+'\t'+str((int(variant[1])-1))+'\t'+str(int(variant[1]))+'\t'+str(variant[2])+'\t'+str(variant[3])+'\n')


if __name__=='__main__':
	#this is path to 'patient_files_offtargets_removed' folder with no headers
	pathToVCF='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/patient-files-offtargets-removed'
	pathToOutput='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/bedFormatted_files_funseq'
	formatBed_funseq(pathToVCF);