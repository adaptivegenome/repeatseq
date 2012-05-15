#!/usr/bin/python
#

import sys
import math

MAX_REF = 70
REF_BIN_SIZE = 15
MAX_FLANK = 5
ERROR = [[[0 for bq in range(5)] for unitsize in range(5)] for ref in range(MAX_REF/REF_BIN_SIZE + 1)]
TRUE = [[[0 for bq in range(5)] for unitsize in range(5)] for ref in range(MAX_REF/REF_BIN_SIZE + 1)]
TOTAL = [[[0 for bq in range(5)] for unitsize in range(5)] for ref in range(MAX_REF/REF_BIN_SIZE + 1)]
InvestigateRepeat = False
majority = numReads = refLength = unitSize = 0

filelist = open(sys.argv[1], 'r')
for filename in filelist:
	file = open(filename.rstrip(), 'r')
	for line in file:
		if line[0] == '~': # it's a header line
			InvestigateRepeat = False
			line = line.split(':')
			numReads = int(line[6].split(' ')[0])
			if numReads >= 25 and float(line[4].split(' ')[0]) != 1: InvestigateRepeat = True
			else: continue
			refLength = int(line[2].split(' ')[0])
			if refLength > MAX_REF: refLength = MAX_REF
			unitSize = int(line[1].split(' ')[1].split('_')[0])
			majority = line[3].split(' ')[0].split('r')[0].split('[')[0] #majority allele
			if majority == "NA": 
				InvestigateRepeat = False
				continue
			else: majority = int(majority)
		
		elif InvestigateRepeat == True and len(line.split(' ')) > 3: 
			# we care about the reads & it's not the reference sequence
			line = line.split(' ')
	
			if unitSize > 5: continue
			BQ = -30*math.log10(float(line[7].split(':')[-1]))
			if BQ < 0: BQ = 0
			elif BQ > 4: BQ = 4
			BQ = int(BQ)

			allele = len(line[1]) - line[1].count('-')
		
			# increment counters
			# print refLength/25, unitSize-1, int(5*BQ), minFlank
			TOTAL[refLength/REF_BIN_SIZE][unitSize - 1][BQ] += 1
			if allele != majority: ERROR[refLength/REF_BIN_SIZE][unitSize - 1][BQ] += 1
			else: TRUE[refLength/REF_BIN_SIZE][unitSize - 1][BQ] += 1
		
print "{",
for unit in range(5):
	print "{",
	for ref in range(MAX_REF/REF_BIN_SIZE + 1):
		print "{",
		for bq in range(5):
			while TRUE[ref][unit][bq] >= 1000 or ERROR[ref][unit][bq] >= 1000: 
				TRUE[ref][unit][bq] /= 10
				ERROR[ref][unit][bq] /= 10
			print "{", TRUE[ref][unit][bq], ",", ERROR[ref][unit][bq], "}",
			
			#if TOTAL[ref][unit][bq] >= 1000: print float(ERROR[ref][unit][bq])/TOTAL[ref][unit][bq],
			#else: print "NA",
			if bq != 4: print ",",
		print "}",
		if ref != MAX_REF/REF_BIN_SIZE: print ",",
		print #newline
	print "}",
	if unit != 4: print ",",
	print #newline
print "}"
		
