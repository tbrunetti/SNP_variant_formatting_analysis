#make sure KAVIAR_trim is downloaded to get list of named variants
#check genome build and make sure same as build used to create VCF

import sys
import os

#sys.argv[1] is a txt file that lists the names of VCFs in pathToVCF to anlayze (VCF files should have headers)
#should be formatted as one file name per line
#sys.argv[2] is KAVIAR vcf library
def findVariant(pathToVCF):
	os.chdir(pathToVCF)
	for vcfFile in sys.argv[1]:
		if vcfFile in os.listdir(pathToVCF):
			


if __name__=='__main__':
	pathToVCF='/mnt/'
	findVariant(pathToVCF);
