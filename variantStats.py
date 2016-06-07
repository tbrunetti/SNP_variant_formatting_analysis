import os
import sys
import numpy as np
import plotly
import plotly.plotly as plotly
import plotly.graph_objs as go
import plotly.offline
import seaborn as sns
import pandas
from collections import Counter

#plotly.sign_in('tbrunetti', '72r07ya8dr')


def variantStats(directoryPatientVariants):
	totalVarsPerPatient=[]
	os.chdir(directoryPatientVariants)
	for allFiles in os.listdir(directoryPatientVariants):
		with open(allFiles) as input:
			totalVariants=sum(1 for line in input)
		totalVarsPerPatient.append(totalVariants)

	os.chdir('/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/patient-files-offtargets-removed')
	for allFiles in os.listdir('/home/tonya/pan_can_project/panData_UCretro_DNA_mutations/patient-files-offtargets-removed'):
		with open(allFiles) as input:
			totalVariants=sum(1 for line in input)
		totalVarsPerPatient.append(totalVariants)

	print len(totalVarsPerPatient)
	print len(np.log10(totalVarsPerPatient))
	log10values=np.log10(totalVarsPerPatient)
	log10=[log10values[i] for i in range(len(log10values))]
	print log10

	print "The min number of variants called in a patient is " + str(min(totalVarsPerPatient))
	print "The max number of variants called in a patient is " + str(max(totalVarsPerPatient))
	print "The mean number of variants in a patient is "+ str(np.average(totalVarsPerPatient))
	print "The median number of variants in a patient is "+ str(np.median(totalVarsPerPatient))

	return totalVarsPerPatient

def makeBoxplot(totalVarsPerPatient):
	
	plotly.offline.plot({
		"data":[go.Box(
			y=np.log10(totalVarsPerPatient),
			name='Retrospective UC pancreatic cancer patients',
			boxpoints='Outliers')],
		"layout":go.Layout(
			title='Retrosepctive Pancreatic Cancer Variants from University of Chicago',
			yaxis=dict(title='log10(variant number)'))
		})

#sys.argv[1]=tab-delimited file, each line is one variant of a patient (PDAC)
#sys.argv[2]=tab-delimited file, each line is one variant of a patient (nonPDAC)
#first column is patientID, 8th column is gene mutation, rest of columns need to have space hold
def deleteriousVariants():
	patientVariantsPDAC={}
	with open(sys.argv[1]) as input:
		headersPDAC=next(input)
		for line in input:
			line=line.split('\t')
			if line[0] in patientVariantsPDAC:
				patientVariantsPDAC[line[0]]=patientVariantsPDAC[line[0]]+[line[8]]
			else:
				patientVariantsPDAC[line[0]]=[line[8]]

	patientVariantsNonPDAC={}
	with open(sys.argv[2]) as input:
		headersNonPDAC=next(input)
		for line in input:
			line=line.split('\t')
			if line[0] in patientVariantsNonPDAC:
				patientVariantsNonPDAC[line[0]]=patientVariantsNonPDAC[line[0]]+[line[8]]
			else:
				patientVariantsNonPDAC[line[0]]=[line[8]]

	numberOfDelMutationsPDAC=[len(patientVariantsPDAC[key]) for key in patientVariantsPDAC]
	numberOfDelMutationsNonPDAC=[len(patientVariantsNonPDAC[key]) for key in patientVariantsNonPDAC]

	trace0=go.Box(
		y=numberOfDelMutationsPDAC,
		boxpoints='Outliers',
		name='PDAC')

	trace1=go.Box(
		y=numberOfDelMutationsNonPDAC,
		boxpoints='Outliers',
		name='non-PDAC')

	plotly.offline.plot({
		"data":[trace0, trace1],
		"layout":go.Layout(
			title='Deleterious variant calls of UC retrosepctive patients',
			xaxis=dict(tickfont=dict(
				size=18)),
			yaxis=dict(title='total variants', 
				titlefont=dict(size=18)),
			margin=go.Margin(
				l=450,
				r=450))
		})

	print "The mean number of variants in a PDAC patient is "+ str(np.average(numberOfDelMutationsPDAC))
	print "The median number of variants in a PDAC patient is "+ str(np.median(numberOfDelMutationsPDAC))
	print "The mean number of variants in a non-PDAC patient is "+ str(np.average(numberOfDelMutationsNonPDAC))
	print "The median number of variants in a non-PDAC patient is "+ str(np.median(numberOfDelMutationsNonPDAC))
	

#sys.argv[1]=tab-delimited file, each line is one variant of a patient
#first column is patientID, 8th column is gene mutation and 9th is mutation type, rest of columns need to have space hold
#return a graph of the number of patients that have a particular mutation type in a given gene,
#also prints out total number of patients and the NON-SYNONYMOUS mutation change is available
def breakdown_of_specific_deleterious_vars():
	geneInterest='SMAD4'
	mutationTypes=[]
	patientMutation={}
	patientMutNonSyn={}
	with open(sys.argv[1]) as input:
		headers=next(input)
		for line in input:
			line=line.split('\t')

			if line[7]==geneInterest:
				mutationTypes.append(line[8])
				if line[0] in patientMutation:
					patientMutation[line[0]]=patientMutation[line[0]]+[line[8]]
				else:
					patientMutation[line[0]]=[line[8]]

			if line[7]==geneInterest and line[8]=='NON_SYNONYMOUS_CODING':
				if line[0] in patientMutNonSyn:
					patientMutNonSyn[line[0]]=patientMutNonSyn[line[0]]+[line[9]]
				else:
					patientMutNonSyn[line[0]]=[line[9]]	

	print patientMutation
	print patientMutNonSyn
	print len(patientMutation)

	dataframe=pandas.DataFrame(mutationTypes, columns=['mutation type'])
	
	with sns.plotting_context("notebook", font_scale=2):
		sns.countplot(x='mutation type', data=dataframe)
		sns.plt.show()

def breakdown_of_deleterious_vars_between_groups():
	pdacMuts=[]
	nonPDACmuts=[]
	patientCompleteProfilePDAC={}
	patientCompleteProfileNonPDAC={}
	with open(sys.argv[1]) as input:
		headers=next(input)
		for line in input:
			line=line.split('\t')
			if line[0] in patientCompleteProfilePDAC:
				patientCompleteProfilePDAC[line[0]]=patientCompleteProfilePDAC[line[0]]+[line[7]]
			if line[0] not in patientCompleteProfilePDAC:
				patientCompleteProfilePDAC[line[0]]=[line[7]]

	with open(sys.argv[2]) as input:
		headers=next(input)
		for line in input:
			line=line.split('\t')
			if line[0] in patientCompleteProfileNonPDAC:
				patientCompleteProfileNonPDAC[line[0]]=patientCompleteProfileNonPDAC[line[0]]+[line[7]]
			if line[0] not in patientCompleteProfileNonPDAC:
				patientCompleteProfileNonPDAC[line[0]]=[line[7]]

	#print patientCompleteProfilePDAC
	#print patientCompleteProfileNonPDAC
	temp1=[]
	temp2=[]
	for key in patientCompleteProfilePDAC:
		temp1=temp1+list(set(patientCompleteProfilePDAC[key]))
	for key in patientCompleteProfileNonPDAC:
		temp2=temp2+list(set(patientCompleteProfileNonPDAC[key]))

	#print len(patientCompleteProfilePDAC)
	#print len(patientCompleteProfileNonPDAC)
	#print Counter(temp1).most_common(10)[i][j]
	#print Counter(temp2).most_common(10)

	forDataFrame=[]
	for i in range(0, 10):
		forDataFrame.append([Counter(temp1).most_common(10)[i][0], Counter(temp1).most_common(10)[i][1], 'PDAC'])
		forDataFrame.append([Counter(temp2).most_common(10)[i][0], Counter(temp2).most_common(10)[i][1], 'NON-PDAC'])

	dataframe=pandas.DataFrame(forDataFrame, columns=['gene', 'frequency', 'diagnosis'])
	print dataframe

	geneOrder=['KRAS', 'TP53', 'GNAQ', 'CHEK2', 'SMAD4', 'KMT2C', 'PRDM16', 'LTBP2', 'NOTCH3', 'OBSCN', 'PTGS1', 'TBX22', 'PRSS1', 'RICTOR', 'TSC2', 'OPTN']
	with sns.plotting_context("notebook", font_scale=1.5):
		sns.barplot(x='gene', y='frequency', hue='diagnosis', data=dataframe, order=geneOrder)
	#sns.set_context("poster", font_scale=5)
		sns.set(font_scale=2.0)
		sns.plt.legend(loc='upper right')
		sns.plt.show()

if __name__=='__main__':
	directoryPatientVariants='/home/tonya/pan_can_project/panData_UCprospective_DNA_mutations/variant-viewer-files-offtargets-removed'
	#totalVarsPerPatient=variantStats(directoryPatientVariants);
	#makeBoxplot(totalVarsPerPatient);
	#deleteriousVariants();
	breakdown_of_specific_deleterious_vars();
	#breakdown_of_deleterious_vars_between_groups();