import os
import sys
import numpy as np
import plotly
import plotly.plotly as plotly
import plotly.graph_objs as go
import plotly.offline
 
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
#first column is patientID, 8th column is mutation type, rest of columns need to have space hold
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
	

	
if __name__=='__main__':
	directoryPatientVariants='/home/tonya/pan_can_project/panData_UCprospective_DNA_mutations/variant-viewer-files-offtargets-removed'
	#totalVarsPerPatient=variantStats(directoryPatientVariants);
	#makeBoxplot(totalVarsPerPatient);
	deleteriousVariants();