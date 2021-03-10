#!/usr/bin/python
import sys,os
tab = '/home/zdh/Projects/04.CF_bats/00.toHuman/OCG2protein.tab'
infile = sys.argv[1]
back = sys.argv[2]
forground = open('forground.txt','w')
background = open('background.txt','w')
dict1 = {}
with open(tab) as TAB:
	for line in TAB:
		line = line.strip()
		lines = line.split()
		ocg = lines[0].split('.')[0]
		dict1[ocg] = lines[-1]
		

with open(infile) as IN:
	for line in IN:
		lines = line.split()
		ocg = lines[0].split('.')[0]
		if ocg in dict1:
			forground.write('%s\n' %dict1[ocg])
		else:
			print(ocg)

for i in os.listdir(back):
	ocg = i.split('.')[0]
	if ocg in dict1:
		background.write('%s\n' %dict1[ocg])

forground.close()
background.close()

os.system('clusterprofiler.pl -a=forground.txt -b=background.txt')
