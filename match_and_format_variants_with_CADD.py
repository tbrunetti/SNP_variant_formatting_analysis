import sys
import pandas
import os
from operator import itemgetter
#sys.argv[1] is a txt file with a list of names of VCFs to format, headers included; one name per line
#designate location of VCFs in pathToVCF variable

#format VCFs for input into CADD database
def formatForCADDinput(pathToVCF, pathToOutput):
	#sets default output path to pathToVCF if path is not specified
	if pathToOutput=='':
		pathToOutput=pathToVCF
	#changes directory to VCF location
	os.chdir(pathToVCF)
	with open(sys.argv[1]) as input:
		for fileVCF in input:
			if fileVCF.rstrip('\n') in os.listdir(pathToVCF):
				with open(fileVCF.rstrip('\n')) as input:
					temp=[]
					#checks to see if header has been recorded
					headerExtracted=0
					for line in input:
						line=line.split('\t')
						headerExtracted+=1
						if headerExtracted==1:
							header=[line[i] for i in range(len(line))]
						else:
							temp.append(line)
					dataframe=pandas.DataFrame(temp, columns=header)
					#only outputs the chr, pos, dnp id, ref, and alt columns in that order for CADD
					#does not write out headers per CADD interface request
					dataframe.to_csv(str(pathToOutput)+'input-CADD-'+str(fileVCF[:-5])+'.txt', sep='\t', cols=['chr', 'pos', 'dbsnp_id', 'ref', 'alt'], header=False, index=False)
			else:
				print str(line)+' does not exist in directory'

#sys.argv[1]=txt file of list of VCFs to add CADD scores and matching CADD output file name
#format -> exact VCF file name,exact CADD output file name
#NOTE: file from CADD database is required!!!!

def findCADD(pathToVCF, pathToCADDfile, pathToMergedFiles):
	with open(sys.argv[1]) as input:
		for line in input:
			line=line.split(',')
			os.chdir(pathToVCF)
			tempVCF=[]
			tempCADD=[]
			with open(line[0]) as input:
				#keep in mind, first row is header information
				for row in input:
					row=row.split('\t')
					#removes EOL
					row[len(row)-1]=row[len(row)-1].rstrip('\n')
					tempVCF.append(row)
				#changes directory to final file output so new file can be created
				os.chdir(pathToMergedFiles)
				f=open(str(line[0][:-4])+'-matched-CADD-scores.txt', 'w')
				#adds the header to the file
				for i in range(0, len(tempVCF[0])-1):
					f.write(str(tempVCF[0][i])+'\t')
				f.write(str(tempVCF[0][len(tempVCF[0])-1])+'\t'+'CADD_rawScore'+'\t'+'CADD_PHRED'+'\n')
				#removes the header index from VCF, so not included in CADD match
				tempVCF.pop(0)
			#sort tempVCF by chromosome, to make comparisons run faster in code
			sortedTempVCF=sorted(tempVCF, key=itemgetter(0))
			

			os.chdir(pathToCADDfile)
			with open(line[1].rstrip('\n')) as input:
				for row in input:
					#skips the header lines in the CADD output
					if '#' in row:
						continue;
					row=row.split('\t')
					#to keep consistent formatting with VCF files, making search easier
					row[0]='chr'+str(row[0])
					row[len(row)-1]=row[len(row)-1].rstrip('\n')
					tempCADD.append(row)
			#sorted, just to make comparisons run faster in code
			sortedTempCADD=sorted(tempCADD, key=itemgetter(0))
			

			os.chdir(pathToMergedFiles)
			for i in range(0, len(sortedTempVCF)):
				#checkpoint to determine if a CADD score has been identified
				checkpoint=0
				for x in range(0, len(sortedTempCADD)):
					if sortedTempVCF[i][0]==sortedTempCADD[x][0] and sortedTempVCF[i][1]==sortedTempCADD[x][1]:
						for z in range(0, len(sortedTempVCF[i])):
							f.write(str(sortedTempVCF[i][z])+'\t')
						f.write(str(sortedTempCADD[x][4])+'\t'+str(sortedTempCADD[x][5])+'\n')
						#critical to have checkpoint to make sure those SNPs that did not have an
						#associated CADD score are still listed in the final output
						checkpoint+=1
						break;
				#means CADD score not available for particular SNP
				if checkpoint==0:
					for j in range(0, len(sortedTempVCF)-1):
						f.write(str(sortedTempVCF[i][j])+'\t')
					f.write(str(sortedTempVCF[i][len(sortedTempVCF)-1])+'\n')


if __name__=='__main__':
	#change to VCF location
	pathToVCF='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/patient-files-offtargets-removed-headers-added/'
	pathToOutput='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/CADD-formatted-off-targets-removed/'
	pathToCADDfile='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/CADD-output/'
	pathToMergedFiles='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/CADD-matched-with-VCF/'
	#formatForCADDinput(pathToVCF, pathToOutput);
	findCADD(pathToVCF, pathToCADDfile, pathToMergedFiles);