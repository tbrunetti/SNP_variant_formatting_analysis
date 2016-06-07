import os
import sys
import pandas
import seaborn as sns

def isoformCheck(path):
	gene='TP53'
	patientIso={}
	noPatIDonlyTups=[]
	os.chdir(path)
	for files in os.listdir(path):
		with open(files) as input:
			#extracts patient ID from file name
			patientID=''
			for i in range(len(files)):
				if files[i]=='-':
					break;
				else:
					patientID=patientID+files[i]
			#comment out header line if no header is included in input
			header=next(input)
			for line in input:
				line=line.split('\t')
				if len(line)<35:
					print line
					break;
				if (line[11]==gene) and (float(line[36].rstrip('\n'))>=15):
					noPatIDonlyTups.append([str(line[14]), str(line[17]), str(line[18])])
					if patientID in patientIso:
						patientIso[patientID]=patientIso[patientID]+[(str(line[14]), str(line[17]), str(line[18]))]
					else:
						patientIso[patientID]=[(str(line[14]), str(line[17]), str(line[18]))]


	print patientIso
	print len(patientIso)
	print noPatIDonlyTups
	for x in range(0, len(noPatIDonlyTups)):
		if noPatIDonlyTups[x][0]=='NON_SYNONYMOUS_CODING':
			noPatIDonlyTups[x][0]='NON_SYN'
	#print sorted(noPatIDonlyTups, key=lambda x: x[0])
	dataframe=pandas.DataFrame(noPatIDonlyTups, columns=['mutation type', 'mutation', 'isoform'])
	print dataframe

	with sns.plotting_context("notebook", font_scale=1.5):
		sns.countplot(y="mutation type", hue="isoform", data=dataframe, palette="Set2")
		sns.plt.show()

if __name__=='__main__':
	#path to VCF with added CADD columns
	path='/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/match-SIFT-CADD'
	pathOutput=''
	isoformCheck(path);