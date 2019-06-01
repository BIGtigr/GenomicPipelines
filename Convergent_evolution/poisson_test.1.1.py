# Zhang, J., and S. Kumar (1997) Detection of convergent and parallel evolution at the amino 
# acid sequence level. Mol Biol Evol 14: 527-536

import sys
from decimal import Decimal
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--i', help=('result file of probCal.py'))
parser.add_option('--rate', help=('path to 03_rate file'))

(options, args) = parser.parse_args()
usage = u'''`Usage:python %prog [options...] [result] [path]'''
# from scipy.misc import factorial
# m = 124#seq length
# Nc = 1 #observed number of convergent-change sites 
# mfc = 0.039 #The probability that a site is a convergent-change site (fc) is the sum of 
# probability of occurrence of all site configurations satisfying the condition
########################################################
numOpts = len(sys.argv)
if numOpts < 2:
    parser.print_help()
    sys.exit()

def factorial(N):
    a=1
    for i in range(2,N+1):
        a*=i
    return a


########################################################
# orginal equation
########################################################
def binomial(m,Nc,mfc):
	x0 = factorial(m)
	#print (a)
	x9 = 0
	for i in range(Nc):
        	x1 = m-i
	        x2 = factorial(x1)
	        x3 = factorial(i)
	        x4 = x2*x3
	        x5 = x0/x4
	        x6 = (mfc/m)**i
	        x7 = (1-(mfc/m))**(m-i)
	        x8 = x5*x6*x7
	        x9 += x8
	theta = 1 - x9
	return(theta)


#######################################################
# poisson test
#######################################################
def poisson(m,Nc,mfc):
	y5 = 0
	for i in range(Nc):	
		y1 = np.e**(-mfc)
		y2 = mfc**i
		y3 = y1*y2
		y4 = y3/factorial(i)
		y5 += y4
	theta = 1-y5	
	return(theta)


######################################################
#a1=fun1(m,Nc,mfc)
#a2=fun2(m,Nc,mfc)
#print (str(a1)+'\n'+str(a2))

with open(options.i) as IN:
	for line in IN:
		line = line.strip()
		(GeneName,branch,Nc,mfc) = line.split()
		if float(Nc) > 10:
			continue
		with open('%s/%s.rat' %(options.rate,GeneName)) as RAT:
			sites = RAT.readlines()
			m = len(sites)
		
		a1=binomial(m,int(float(Nc)),float(mfc))
		a2=poisson(m,int(float(Nc)),float(mfc))
		if a1 < 0.05 or a2 < 0.05:
			a3 = '%s\t%s\t%s\t%s\t%s\t%s' %(GeneName,branch,Nc,mfc,a1,a2)
			print (a3)
		
