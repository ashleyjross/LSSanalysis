from math import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import rc
from matplotlib import rcParams


rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)


ebossdir = '/Users/ashleyross/Dropbox/eboss/' #where AJR puts correlation functions, writes out results

datxi = np.loadtxt(ebossdir +'dr16_qso_v7_1_NS_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()
BAOmod = np.loadtxt(ebossdir + 'ximodQSOdataJH7_15meBp0.40.44.03.08.015.00NScomb.dat').transpose()
BAOmodep0 = np.loadtxt(ebossdir + 'ximodQSOdataJH7_15meBp0.40.44.03.08.015.00NScombep0.dat').transpose()

def P2(mu):
	return .5*(3.*mu**2.-1.)


def plotQSOBAO_xi0():
	d = datxi[1][10:30]
	err = datxi[2][10:30]
	rl = datxi[0][10:30]
	mod = BAOmod[1][0:20]
	modsm = BAOmod[2][0:20]
	dat = d - modsm
	mod = mod-modsm
	plt.errorbar(rl,dat*1.e3,err*1.e3,fmt='o',mfc='w',ms=8,elinewidth=1.75,capsize=4,mew=1.5,color='k')
	plt.plot(rl,mod*1.e3,'k--')
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$10^3(\xi_0 - \xi_{\rm 0, smooth})$',labelpad=-0.5)
	plt.savefig(ebossdir+'QSOxiBAO_mono.png')
	
	plt.show()

def plotQSOBAO_xi2():
	d = datxi[3][10:30]
	err = datxi[4][10:30]
	rl = datxi[0][10:30]
	mod = BAOmod[1][20:40]
	modsm = BAOmodep0[1][20:40]
	dat = d - modsm
	mod = mod-modsm
	plt.errorbar(rl,dat*1.e3,err*1.e3,fmt='o',mfc='w',ms=8,elinewidth=1.75,capsize=4,mew=1.5,color='k')
	plt.plot(rl,mod*1.e3,'k--')
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$10^3(\xi_2 - \xi_2(\epsilon=0))$',labelpad=-0.5)
	plt.savefig(ebossdir+'QSOxiBAO_quad.png')
	
	plt.show()
	
def plotQSOBAO_xiring():
	d2 = datxi[3][10:30]
	d0 = datxi[1][10:30]
	kl = datxi[0][10:30]
	modsm0 = BAOmod[2][0:20]
	modsm2 = BAOmod[2][20:40]
	dat0 = d0-modsm0
	dat2 = d2-modsm2
	datm = np.zeros((2*(len(kl)+10),2*(len(kl)+10)))
	for i in range(0,2*(len(kl)+10)):
		#if i < len(kl):
		#	kp = kl[-1-i]
		#else:
		#	kp = kl[i-len(kl)]
		kp = -150+2.5+5.*i
		for j in range(0,2*(len(kl)+10)):
			#if j < len(kl):
			#	kr = kl[-1-j]
			#else:	
			#	kr = kl[j-len(kl)]
			kr = -150+2.5+5.*j
			kt = sqrt(kp**2.+kr**2.)
			if kt < 150 and kt > 50:
				mu = kr/kt
				kind = int((kt-50.)/5.)
				pt = dat0[kind]+dat2[kind]*P2(mu)
				datm[i][j] = pt
	#err = d[2]/d0[4]
	#mod = (d[3]-d[4])/d0[4]
	#plt.errorbar(kl,dat,err,fmt='ko')
	#plt.plot(kl,mod,'k--')
	plt.imshow(datm,origin='lower',interpolation='bicubic',cmap='gray')
	locs,labels = plt.xticks()
	locs = np.array([-.5 ,9.5,  19.5,  29.5,  39.5,49.5,59.5  ])
	labs = []
	for loc in locs:
		kval = abs(-150+2.5+5.*loc)
		labs.append(str(kval.round(3)))
		
	plt.xticks(locs,labs)
	plt.yticks(locs,labs)
	print(locs)
	plt.xlim(-1,len(datm))
	plt.ylim(-1,len(datm))	
	plt.xlabel(r'$s_{\perp}$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s_{||}$ ($h^{-1}$Mpc)')
	plt.savefig(ebossdir+'QSOxiBAOring.png')
	plt.show()
