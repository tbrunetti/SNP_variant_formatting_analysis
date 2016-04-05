import sys
import os
from subprocess import call
import pandas
#this script should be used to format a BAF file for input into ADTEx CNV algorithm
#input files are VCF format for each matched tumor-normal pair per patient

#storeInfo reads in VCFs and matches them to patient IDs and stores them in an easy to access format
def storeInfo(pathToVCFs):
	#change to directory where VCF files are located
	os.chdir(pathToVCFs)
	#store patient IDs with corresponding VCFs
	storeFiles={}
	#extract file names from current working directory
	for files in os.listdir(os.getcwd()):
		temp=[]
		fileName=''
		with open(files) as input:
			#establishes patient study number and BID
			for i in range(0, len(files)):
				fileName=fileName+str(files[i])
				if fileName.count('-')==3:
					fileName=fileName[:-1]
					break;
			
			#reads in patient VCF file and removes EOL characters	
			for line in input:
				line=line.split('\t')
				line[len(line)-1]=line[len(line)-1].rstrip('\n')
				temp.append(line)
		
		#extracts headers from files
		headers=[temp[0][x] for x in range(len(temp[0]))]
		temp.pop(0)
		dataframe=pandas.DataFrame(temp, columns=headers)
		#key=studyID-BID
		#value=dataframe of VCF
		storeFiles[fileName]=dataframe

	print 'There are '+str(len(storeFiles)) +' patients with matched VCF files'
	return storeFiles

def formatBAF(pathToOutput, storeFiles):
	os.chdir(pathToOutput)
	for key in storeFiles:
		f=open(str(key)+'-BAF-ADTEx-Input.txt', 'w')
		#header names are chosen by ADTEx and must be kept as is
		f.write('chrom'+'\t'+'SNP_loc'+'\t'+'control_BAF'+'\t'+'tumor_BAF'+'\t'+'control_doc'+'\t'+'tumor_doc'+'\n')
		#.shape returns a tuple (row, columns), only interested in rows
		for i in range (0, int(storeFiles[key].shape[0])):
			f.write(str(storeFiles[key].loc[i, 'chr'])+'\t')
			f.write(str(storeFiles[key].loc[i, 'pos'])+'\t')
			f.write(str(float(storeFiles[key].loc[i, 'pct_normal_alt'].rstrip('%'))/100)+'\t')
			f.write(str(float(storeFiles[key].loc[i, 'pct_tumor_alt'].rstrip('%'))/100)+'\t')
			f.write(str(int(storeFiles[key].loc[i, 'normal_ref_count'])+int(storeFiles[key].loc[i, 'normal_alt_count']))+'\t')
			
			#if/else is here so an additional EOL character is not added to the end of the file
			if i==(int(storeFiles[key].shape[0])-1):
				f.write(str(int(storeFiles[key].loc[i, 'tumor_ref_count'])+int(storeFiles[key].loc[i, 'tumor_alt_count'])))
			else:
				f.write(str(int(storeFiles[key].loc[i, 'tumor_ref_count'])+int(storeFiles[key].loc[i, 'tumor_alt_count']))+'\n')



if __name__=='__main__':
	#change this to path of VCFs and to where output files should be stored
	pathToVCFs='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/patient-files-offtargets-removed-headers-added'
	pathToOutput='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/BAF-ADTEx-files'
	storeFiles=storeInfo(pathToVCFs);
	formatBAF(pathToOutput, storeFiles)
