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
		

	
if __name__=='__main__':
	directoryPatientVariants='/home/tonya/pan_can_project/panData_UCprospective_DNA_mutations/variant-viewer-files-offtargets-removed'
	totalVarsPerPatient=variantStats(directoryPatientVariants);
	#makeBoxplot(totalVarsPerPatient);