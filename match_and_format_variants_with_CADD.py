import sys
import pandas
import os
#sys.argv[1] is a txt file with a list of names of VCFs to format, headers included; one name per line
#designate location of VCFs in pathToVCF variable

def formatForCADDinput(pathToVCF):
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
					dataframe.to_csv('input-CADD-'+str(fileVCF[:-4])+'.txt', sep='\t', cols=['chr', 'pos', 'dbsnp_id', 'ref', 'alt'], header=False)
			else:
				print str(line)+' does not exist in directory'

def findCADD():
	print "unfinished"


if __name__=='__main__':
	#change to VCF location
	pathToVCF='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/patient-files-offtargets-removed-headers-added/'
	formatForCADDinput(pathToVCF);
	#findCADD();