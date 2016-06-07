import sys

#sys.argv[1]:tab-delimited file of lists of int(chromosome number) with int(location of SNP)
#sys.argv[2]:bed file to check against, i.e. your capture panel=int(chr), int(start), int(end)
def capturePanel_To_SNParray():
	snpArray={}
	capture={}
	with open(sys.argv[1]) as input:
		for line in input:
			line=line.split('\t')
			if line[0] in snpArray:
				snpArray[line[0]]=snpArray[line[0]]+[line[1].rstrip('\n')]
			else:
				snpArray[line[0]]=[line[1].rstrip('\n')]

	with open(sys.argv[2]) as input:
		for line in input:
			line=line.split('\t')
			if line[0] in capture:
				capture[line[0]].append([line[1], line[2].rstrip('\n')])
			else:
				capture[line[0]]=[[line[1], line[2].rstrip('\n')]]

	
	f=open('/home/tonya/Desktop/illumina-snps-oncoArray500-matched-CPCI-capture-panel.txt', 'w')
	matchedRegions=[]
	for key in snpArray:
		try:
			for x in range(0, len(snpArray[key])):
				for i in range(0, len(capture[key])):
					try:
						if int(capture[key][i][0])<=int(snpArray[key][x])<=int(capture[key][i][1]):
							matchedRegions.append(str(key)+'-'+str(capture[key][i][0])+'-'+str(capture[key][i][0]))
							f.write(str(key)+'\t'+str(capture[key][i][0])+'\t'+str(capture[key][i][1])+'\t'+snpArray[key][x]+'\n')
					except KeyError:
						print "Error, key does not exist"
		except KeyError:
			print "Error, key "+ str(key)+ " does not exist (2)"

	print 'The total number capture regions covered by the SNP array is '+ str(len(set(matchedRegions)))

if __name__=='__main__':
	capturePanel_To_SNParray();