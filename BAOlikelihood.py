'''
Add new results to BOSS DR12 BAO likelihood
'''
from math import *
import numpy as np
from Cosmo import distance #this is other python in LSSanalysis
#directory with DR12 DM/HZ likelihood
dirDR12 = '/Users/ashleyross/Downloads/ALAM_ET_AL_2016_consensus_and_individual_Gaussian_constraints/COMBINEDDR12_BAO_consensus_dM_Hz/'

zl = [0.38,0.38,0.51,0.51,0.61,0.61]
ml = [1512.39,81.2087,1975.22,90.9029,2306.68,98.9647]
covt = np.loadtxt(dirDR12+'BAO_consensus_covtot_dM_Hz.txt')
d = distance(.31,.69,h=0.676)

def testcosmo(alphat,alphar,z):
	'''
	The purpose of this is to confirm that Cosmo.py reproduces the alpha_perp/|| -> DM/H conversion
	alphat is the transverse alpha (alpha_perp), alphar is the radial alpha (alpha_||) 
	'''
	
	print(d.dc(z)*alphat)
	print(d.Hz(z)/alphar*100.)
	'''
	comparison with file with ROSS in name and Ross et al. 2017 gives matches all to within 0.05%; should be good enough
	'''
	
#def converttoalph():
fidl  = np.zeros(len(ml))
mtl = np.array(ml)
for i in range(0,len(zl),2):
	z = zl[i]
	fidl[i] = d.dc(z)
	fidl[i+1] = .01/d.Hz(z)
	mtl[i+1] = 1./ml[i+1]
fidl = np.array(fidl)	
alphl = mtl/fidl
print(alphl) #check that matches expectation
cova = covt
corr = np.zeros((6,6))
for i in range(0,len(covt)):
	for j in range(0,len(covt)):
		cova[i][j] /= ml[i]*ml[j]
		#corr[i][j] = cova[i][j]/sqrt(cova[i][i]*cova[j][j])
		#print(cova[i][j])
		if i == j:
			print(sqrt(cova[i][j]))
for i in range(0,len(covt)):
	for j in range(0,len(covt)):
		corr[i][j] = cova[i][j]/sqrt(cova[i][i]*cova[j][j])
		
print(corr)
#note positive correlation implies use 1/alpha_||
#correlation to be preserved is that between high and low redshift z bins
cap = cova[0][4]/sqrt(cova[0][0]*cova[4][4]) #transverse correlation
car = cova[1][5]/sqrt(cova[1][1]*cova[5][5]) #radial correlation
carp = cova[0][5]/sqrt(cova[0][0]*cova[5][5]) #cross corr
carp2 = cova[1][4]/sqrt(cova[1][1]*cova[4][4])
print(cap,car,carp,carp2)
#take mean to get D_A/H cross corr
carpm = (carp+carp2)/2.
canew = np.zeros((6,6))
#first four elements are same as old one
for i in range(0,4):
	for j in range(0,4):
		canew[i][j] = cova[i][j]
'''
NEW BAO measurement goes HERE
'''
newmeas = np.zeros((2,2))
newmeas[0][0] = 0.021*0.021 #Gaussian approximation to alpha_perp error
newmeas[1][1] = 0.027*0.027 #Gaussian approximation to 1/alpha_|| error
newmeas[1][0] = 0.39*0.021*0.027 #Gaussian approximation to cross term
newmeas[0][1] = 0.39*0.021*0.027
diralph = '/Users/ashleyross/Dropbox/eboss/BAO_dM_Hz_BOSSLRG/' #wherever you put the alphas
alphas = np.loadtxt(diralph+'LRGgaussalphas.txt') #do it this way so nothing is public
ap = alphas[0] #alpha_perp
ar = alphas[1] #alpha_||
zlrg = 0.72
zl[4] = zlrg
zl[5] = zlrg
ml[4] = d.dc(zlrg)*ap
print(ml[4])
ml[5] = d.Hz(zlrg)/ar*100.
print(ml[5])
#assign these to cov matrix
canew[4][4] = newmeas[0][0]
canew[5][5] = newmeas[1][1]
canew[4][5] = newmeas[0][1]
canew[5][4] = newmeas[1][0]
#assign new cross terms, first z bin 1 and zbin 3
canew[0][4] = cap*sqrt(canew[0][0]*canew[4][4])
canew[1][5] = car*sqrt(canew[1][1]*canew[5][5])
canew[1][4] = carpm*sqrt(canew[1][1]*canew[4][4])
canew[0][5] = carpm*sqrt(canew[0][0]*canew[5][5])
#now zbin 2 and zbin 3
canew[2][4] = cap*sqrt(canew[2][2]*canew[4][4])
canew[3][5] = car*sqrt(canew[3][3]*canew[5][5])
canew[3][4] = carpm*sqrt(canew[3][3]*canew[4][4])
canew[2][5] = carpm*sqrt(canew[2][2]*canew[5][5])
#apply symmetry
for i in range(0,6):
	for j in range(0,6):
		if canew[i][j] == 0:
			canew[i][j] = canew[j][i]

corrnew = np.zeros((6,6))
for i in range(0,len(canew)):
	for j in range(0,len(canew)):
		corrnew[i][j] = canew[i][j]/sqrt(canew[i][i]*canew[j][j])
		
print(corrnew)

#print(canew)
#check to see if we gained any information over BOSS DR12
print(sum(sum(np.linalg.inv(canew))),sum(sum(np.linalg.inv(cova))))

fo = open(diralph+'BAO_vals_DR16_dM_Hz.txt','w')
for i in range(0,len(zl),2):
	fo.write(str(zl[i])+' dM(rsfid/rs) '+str(ml[i])+'\n')
	fo.write(str(zl[i+1])+' Hz(rs/rsfid) '+str(ml[i+1])+'\n')
fo.close()
fo = open(diralph+'BAO_covtot_dM_Hz.txt','w')
for i in range(0,6):
	for j in range(0,6):
		c = canew[i][j]*ml[i]*ml[j]
		fo.write(str(c)+' ')
	fo.write('\n')
fo.close()
			
	
	