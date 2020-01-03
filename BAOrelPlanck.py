import numpy as np
from math import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from matplotlib import rcParams
from pylab import *
from numpy import loadtxt as load
import matplotlib.cm as cm
from Cosmo import *
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)

key_dict = {'DR7Per':['Darkgoldenrod', 'goldenrod', '>', 'SDSS-DR7',0.22,1.02],
			  '6dFGS'  :['Brown', 'BlanchedAlmond', 'o', '6dFGS', 0.06,.943],
			  'Wznorec'  :['tomato', 'salmon', 'D', 'WiggleZ', 0.65,1.05],
			  'DR7rec':['Crimson', 'salmon', 'p', 'DR7-recon', 0.3,1.001],
			  'Wzrec1'   :['thistle', 'Linen', 'h', 'Wz-recon', 0.45,1.082],
			  'Wzrec2'   :['thistle', 'Linen', 'h', 'Wz-recon', 0.45,1.082],
			  'Wzrec3'   :['thistle', 'Linen', 'h', 'Wz-recon', 0.45,1.082],
			  'BOSSCMDR11':['tomato', 'MintCream', '<', 'DR11-CMASS', 0.34,.975],
			  'BOSSLZDR11':['tomato', 'MintCream', '<', 'DR11-LowZ', 0.34,.975],
			  'LyaXDR11'   :['Dodgerblue', 'lightblue', 'o', r'DR11-Ly$\alpha$', 2.,1.02],
			  'MGS'        :['hotpink', 'salmon', '>', 'SDSS-MGS', 0.093,1.084],
			  'BOSSDR12_1' :['k', 'grey', 's', 'BOSS-LRG', 0.34,.975],
			  'BOSSDR12_2' :['k', 'grey', 's', 'BOSS-LRG', 0.34,.975],
			  'BOSSDR12_3' :['k', 'grey', 's', 'BOSS-LRG', 0.34,.975],
			  'LyADR12':['DarkTurquoise', 'skyblue', 'd', r'DR12-Ly$\alpha$', 2.2,1.06],
			  'LyADR12C':['DarkTurquoise', 'skyblue', 'd', r'DR12-Ly$\alpha$', 2.2,1.06],
			  'QSODR14'     :['tomato', 'salmon', 'H', 'DR14-QSO', 1.4,1.06],
			  'LRGDR14':['tomato', 'salmon', 'H', 'DR14-LRG',0.6,1.025],
			  'DESY1'     :['purple', 'plum', '>', 'DES-y1',0.59,.955],
			  'LyADR14':['Dodgerblue', 'powderblue', 'o', r'DR14-Ly$\alpha-{\rm auto}$',2.,1.02],
			  'LyADR14X':['DarkTurquoise', 'powderblue', 'd', r'DR14-Ly$\alpha-{\rm cross}$',2.,1.02],
			  'LyADR14C':['Dodgerblue', 'powderblue', 'o', r'DR14-Ly$\alpha-{\rm comb}$',2.,1.02],
			  'QSODR16'   :['orange', 'gold', '*', 'eBOSS-QSO', 1.4,1.06],
			  'ELGDR16'   :['forestgreen', 'mediumseagreen', '*', 'eBOSS-ELG', 0.86,.97],
			  'LRGDR16'   :['firebrick', 'indianred', '*', 'DR16-LRG', 0.6,1.025],
			  'LyADR16C'   :['DarkTurquoise', 'powderblue', '*', r'eBOSS-Ly$\alpha$',2.26,1.01]}


def plotBAOrelBOSSfid_dv():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		dv = float(ln[4])
		if dv != -1:
			zl.append(z)
			sdv = float(ln[5])
			EH = int(ln[-1].strip('\n'))
			dvrsfidc = dcosfid.dV(z)/rsfidcamb
			if EH:
				dvrsfid = dcosfid.dV(z)/rsfidEH
			else:
				dvrsfid = dvrsfidc
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='o')
	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		pdvfidl.append(dpcos.dV(z)/prs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-')
	plt.show()	

def plotBAOrelBOSSfid_dv_preboss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	#pom = 0.3147
	#pobh = 0.02233
	#ph = 0.6737
	#prs = 147.18
	#use WMAP7
	pom = 0.271
	pobh = 0.02227
	ph = 0.703
	prs = 152.76
	dpcos = distance(pom,1-pom,ph,obhh=pobh)
	print(dpcos.rs)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[4])
			
		sdv = float(ln[5])
		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		if EH:
			dvrsfid = dcosfid.dV(z)/rsfidEH
		else:
			dvrsfid = dvrsfidc
		if lab == 'DR7Per':
			plt.text(0.2,.945,'SDSS DR7',fontsize=14,color='steelblue')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='<',color='steelblue',markeredgecolor='k',markersize=7)
		if lab == '6dFGS':
			plt.text(0.06,.935,'6dFGS',fontsize=14,color='green')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='gs',markeredgecolor='k',markersize=7)
		if lab == 'Wznorec':
			plt.text(0.65,1.05,'WiggleZ',fontsize=14,color='.5')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='>',color='.5',markeredgecolor='k',markersize=7)
		if lab == 'DR7rec':
			plt.text(0.4,.99,'SDSS DR7',fontsize=14,color='firebrick')	
			plt.text(0.4,.98,'LRGs rec',fontsize=14,color='firebrick')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='h',color='firebrick',markeredgecolor='k',markersize=7)
		

	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.dV(z)/rsfidEH
		pdvfidl.append(dpcos.dV(z)/dpcos.rs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-',label='WMAP7')
	plt.legend(loc='lower right')
	plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel('BAO distance/eBOSS fiducial')
	plt.title('Spherically averaged, before BOSS')
		
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAObeforeBOSS.png')
	plt.show()

def plotBAOrelBOSSfid_dv_boss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2015 best LCDM
	pom = 0.315
	pobh = 0.02205
	ph = 0.673
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[4])
			
		sdv = float(ln[5])
		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		if EH:
			dvrsfid = dcosfid.dV(z)/rsfidEH
		else:
			dvrsfid = dvrsfidc
		if lab == '6dFGS':
			plt.text(0.06,.935,'6dFGS',fontsize=14,color='green')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='gs',markeredgecolor='k',markersize=7)
		if lab == 'Wzrec1':
			plt.text(0.45,1.082,'WiggleZ (rec)',fontsize=14,color='.8')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='v',markersize=6,elinewidth=.5,color='.8',markeredgecolor='k')
		if lab == 'Wzrec3' or lab == 'Wzrec2':	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='v',markersize=6,elinewidth=.5,color='.8',markeredgecolor='k')
		if lab == 'MGS':
			plt.text(0.1,1.095,'SDSS',fontsize=14,color='purple')
			plt.text(0.1,1.084,'MGS',fontsize=14,color='purple')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='^',markeredgecolor='k',markersize=7,elinewidth=1.75,color='purple')
		if lab == 'BOSSDR12_1':
			plt.text(0.34,.975,'BOSS',fontsize=14,color='r')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'BOSSDR12_2' or lab == 'BOSSDR12_3':
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR12C':
			plt.text(2.,1.02,'BOSS',fontsize=14,color='k')	
			plt.text(1.95,1.01,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='H',markersize=7,elinewidth=1.75,color='magenta',markeredgecolor='k')
		

	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		pdvfidl.append(dpcos.dV(z)/prs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-',label='Planck 2015')
	plt.legend(loc='lower right')
	plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel('BAO distance/eBOSS fiducial')
	plt.title('Spherically averaged or best constrained, after BOSS')
		
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAOafterBOSS.png')
	plt.show()


def plotBAOrelBOSSfid_dv_eboss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[4])
			
		sdv = float(ln[5])
		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		if EH:
			dvrsfid = dcosfid.dV(z)/rsfidEH
		else:
			dvrsfid = dvrsfidc
		if lab == '6dFGS':
			plt.text(0.06,.935,'6dFGS',fontsize=14,color='green')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='gs',markeredgecolor='k',markersize=7)
		if lab == 'Wzrec1':
			plt.text(0.45,1.082,'WiggleZ (rec)',fontsize=14,color='.8')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='v',markersize=6,elinewidth=.5,color='.8',markeredgecolor='k')
		if lab == 'Wzrec3' or lab == 'Wzrec2':	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='v',markersize=6,elinewidth=.5,color='.8',markeredgecolor='k')
		if lab == 'MGS':
			plt.text(0.1,1.095,'SDSS',fontsize=14,color='purple')
			plt.text(0.1,1.084,'MGS',fontsize=14,color='purple')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='^',markeredgecolor='k',markersize=7,elinewidth=1.75,color='purple')
		if lab == 'BOSSDR12_1':
			plt.text(0.34,.975,'BOSS',fontsize=14,color='r')	
			plt.text(0.34,.965,'galaxies',fontsize=14,color='r')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'BOSSDR12_2':
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'LRGDR16':
			plt.text(0.6,1.025,'eBOSS',fontsize=14,color='orange')	
			plt.text(0.61,1.015,'LRGs',fontsize=14,color='orange')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='orange',markeredgecolor='k')
		if lab == 'DESY1':
			plt.text(0.5,.955,'DESY1',fontsize=14,color='lightblue')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='o',markersize=7,elinewidth=1.,color='lightblue',markeredgecolor='k')

		if lab == 'ELGDR16':
			plt.text(0.86,.97,'eBOSS',fontsize=14,color='b')	
			plt.text(0.87,.96,'ELGs',fontsize=14,color='b')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='b',markeredgecolor='k')
		if lab == 'QSODR16':
			plt.text(1.4,1.06,'eBOSS',fontsize=14,color='k')	
			plt.text(1.4,1.05,'quasars',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='cyan',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR16C':
			plt.text(2.,1.02,'eBOSS',fontsize=14,color='k')	
			plt.text(1.95,1.01,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		

	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		pdvfidl.append(dpcos.dV(z)/prs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-',label='Planck 2018')
	plt.legend(loc='lower right')
	plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel('BAO distance/eBOSS fiducial')
	plt.title('Spherically averaged or best constrained, after eBOSS')
	
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAOaftereBOSS.png')	
	plt.show()

def plotBAOrelBOSSfid_dvmh_eboss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[4])
			
		sdv = float(ln[5])
		dm = float(ln[6])			
		sdm = float(ln[7])
		dh = float(ln[8])			
		sdh = float(ln[9])

		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		mrs = 1.
		if EH:
			dvrsfid = dcosfid.dV(z)/rsfidEH
			mrs = rsfidEH/rsfidcamb
		else:
			dvrsfid = dvrsfidc
		if lab == '6dFGS':
			#plt.text(0.06,.935,'6dFGS',fontsize=14,color='green')	
			plt.errorbar(z,dv/np.sqrt(z)*mrs,sdv/np.sqrt(z)*mrs,fmt='k>',markeredgecolor='k',markersize=7)
		if lab == 'MGS':
			#plt.text(0.1,1.095,'SDSS',fontsize=14,color='purple')
			#plt.text(0.1,1.084,'MGS',fontsize=14,color='purple')	
			plt.errorbar(z,dv/np.sqrt(z)*mrs,sdv/np.sqrt(z)*mrs,fmt='k^',markeredgecolor='k',markersize=7)
		if lab == 'BOSSDR12_1':
			#plt.text(0.34,.975,'BOSS',fontsize=14,color='r')	
			#plt.text(0.34,.965,'galaxies',fontsize=14,color='r')
			plt.errorbar(z,dm/np.sqrt(z)*mrs,sdm/np.sqrt(z)*mrs,fmt='D',markersize=7,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*np.sqrt(z)*mrs,sdh*np.sqrt(z)*mrs,fmt='D',markersize=7,color='steelblue',markeredgecolor='k')
		if lab == 'BOSSDR12_2':
			plt.errorbar(z,dm/np.sqrt(z)*mrs,sdm/np.sqrt(z)*mrs,fmt='D',markersize=7,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*np.sqrt(z)*mrs,sdh*np.sqrt(z)*mrs,fmt='D',markersize=7,color='steelblue',markeredgecolor='k')
		if lab == 'LRGDR16':
			#plt.text(0.6,1.025,'eBOSS',fontsize=14,color='orange')	
			#plt.text(0.61,1.015,'LRGs',fontsize=14,color='orange')
			plt.errorbar(z,dm/np.sqrt(z)*mrs,sdm/np.sqrt(z)*mrs,fmt='*',markersize=14,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*np.sqrt(z)*mrs,sdh*np.sqrt(z)*mrs,fmt='*',markersize=14,color='steelblue',markeredgecolor='k')
		if lab == 'DESY1':
			#plt.text(0.5,.955,'DESY1',fontsize=14,color='lightblue')	
			plt.errorbar(z,dm/np.sqrt(z)*mrs,sdm/np.sqrt(z)*mrs,fmt='o',markersize=7,color='firebrick',markeredgecolor='k')

		if lab == 'ELGDR16':
			#plt.text(0.86,.97,'eBOSS',fontsize=14,color='b')	
			#plt.text(0.87,.96,'ELGs',fontsize=14,color='b')
			plt.errorbar(z,dv/np.sqrt(z)*mrs,sdv/np.sqrt(z)*mrs,fmt='*',markersize=14,color='k',markeredgecolor='k')
		if lab == 'QSODR16':
			#plt.text(1.4,1.06,'eBOSS',fontsize=14,color='k')	
			#plt.text(1.4,1.05,'quasars',fontsize=14,color='k')
			plt.errorbar(z,dm/np.sqrt(z)*mrs,sdm/np.sqrt(z)*mrs,fmt='*',markersize=14,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*np.sqrt(z)*mrs,sdh*np.sqrt(z)*mrs,fmt='*',markersize=14,color='steelblue',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR16C':
			#plt.text(2.,1.02,'eBOSS',fontsize=14,color='k')	
			#plt.text(1.95,1.01,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dm/np.sqrt(z)*mrs,sdm/np.sqrt(z)*mrs,fmt='*',markersize=14,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*np.sqrt(z)*mrs,sdh*np.sqrt(z)*mrs,fmt='*',markersize=14,color='steelblue',markeredgecolor='k')
		

	zls = zl.sort()
	pdvfidl = []
	pdmfidl = []
	pdhfidl = []
	for z in zl:
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		dmrsfidc = dcosfid.dc(z)/rsfidcamb
		dhrsfidc = dcosfid.cHz(z)/rsfidcamb
		pdvfidl.append(dvrsfidc/np.sqrt(z))
		pdmfidl.append(dmrsfidc/np.sqrt(z))
		pdhfidl.append(dhrsfidc*np.sqrt(z))
	plt.plot(zl,pdvfidl,'k-',label=r' $D_V$')
	plt.plot(zl,pdmfidl,'-',color='firebrick',label=r' $D_M$')
	plt.plot(zl,pdhfidl,'-',color='steelblue',label=r' $zD_H$')
	plt.legend(loc='lower right')
	#plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel(r'BAO distance/($\sqrt{z}r_s(z_{\rm drag})$')
	plt.ylabel(r'BAO distance/$D_V(\rm{eBOSS fid})$')
	#plt.title('Spherically averaged or best constrained, after eBOSS')
	
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAOVMHaftereBOSS.png')	
	plt.show()

def plotBAOrelBOSSfid_dvmhdv_eboss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[4])
			
		sdv = float(ln[5])
		dm = float(ln[6])			
		sdm = float(ln[7])
		dh = float(ln[8])			
		sdh = float(ln[9])

		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		mrs = 1.
		if EH:
			dvrsfid = dcosfid.dV(z)/rsfidEH
			mrs = rsfidEH/rsfidcamb
		else:
			dvrsfid = dvrsfidc
		if lab == '6dFGS':
			#plt.text(0.06,.935,'6dFGS',fontsize=14,color='green')	
			df = plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='k>',markeredgecolor='k',markersize=4,label='6dFGS')
		if lab == 'MGS':
			#plt.text(0.1,1.095,'SDSS',fontsize=14,color='purple')
			#plt.text(0.1,1.084,'MGS',fontsize=14,color='purple')	
			mg = plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='k^',markeredgecolor='k',markersize=4,label='SDSS MGS')
		if lab == 'BOSSDR12_1':
			#plt.text(0.34,.975,'BOSS',fontsize=14,color='r')	
			#plt.text(0.34,.965,'galaxies',fontsize=14,color='r')
			bs = plt.errorbar(2.5,1,.01,fmt='Dk',markersize=4,label='BOSS')
			plt.errorbar(z,dm/dvrsfid,sdm/dvrsfid,fmt='D',markersize=4,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*z/dvrsfid,sdh*z/dvrsfid,fmt='D',markersize=4,color='steelblue',markeredgecolor='k')
		if lab == 'BOSSDR12_2':
			plt.errorbar(z,dm/dvrsfid,sdm/dvrsfid,fmt='D',markersize=4,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*z/dvrsfid,sdh*z/dvrsfid,fmt='D',markersize=4,color='steelblue',markeredgecolor='k')
		if lab == 'LRGDR16':
			#plt.text(0.6,1.025,'eBOSS',fontsize=14,color='orange')	
			#plt.text(0.61,1.015,'LRGs',fontsize=14,color='orange')
			plt.errorbar(z,dm/dvrsfid,sdm/dvrsfid,fmt='*',markersize=10,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*z/dvrsfid,sdh*z/dvrsfid,fmt='*',markersize=10,color='steelblue',markeredgecolor='k')
		if lab == 'DESY1':
			#plt.text(0.5,.955,'DESY1',fontsize=14,color='lightblue')	
			ds = plt.errorbar(2.5,1,.01,fmt='ko',markersize=4,label='DES')
			plt.errorbar(z,dm/dvrsfid,sdm/dvrsfid,fmt='o',markersize=4,color='firebrick',markeredgecolor='k')

		if lab == 'ELGDR16':
			#plt.text(0.86,.97,'eBOSS',fontsize=14,color='b')	
			#plt.text(0.87,.96,'ELGs',fontsize=14,color='b')
			eb = plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=10,color='k',markeredgecolor='k',label='eBOSS')
		if lab == 'QSODR16':
			#plt.text(1.4,1.06,'eBOSS',fontsize=14,color='k')	
			#plt.text(1.4,1.05,'quasars',fontsize=14,color='k')
			plt.errorbar(z,dm/dvrsfid,sdm/dvrsfid,fmt='*',markersize=10,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*z/dvrsfid,sdh*z/dvrsfid,fmt='*',markersize=10,color='steelblue',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR16C':
			#plt.text(2.,1.02,'eBOSS',fontsize=14,color='k')	
			#plt.text(1.95,1.01,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dm/dvrsfid,sdm/dvrsfid,fmt='*',markersize=10,color='firebrick',markeredgecolor='k')
			plt.errorbar(z,dh*z/dvrsfid,sdh*z/dvrsfid,fmt='*',markersize=5,color='steelblue',markeredgecolor='k')
		
	first_legend = plt.legend(handles=[df,mg,eb,bs,ds],loc='center right',facecolor='w',framealpha=1)
	ax = plt.gca().add_artist(first_legend)
	zls = zl.sort()
	pdvfidl = []
	pdmfidl = []
	pdhfidl = []
	for z in zl:
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		dvrspc = dpcos.dV(z)/prs
		dmrspc = dpcos.dc(z)/prs
		dhrspc = dpcos.cHz(z)/prs
		pdvfidl.append(dvrspc/dvrsfidc)
		pdmfidl.append(dmrspc/dvrsfidc)
		pdhfidl.append(dhrspc*z/dvrsfidc)
	dv, = plt.plot(zl,pdvfidl,'k-',label=r' $D_V$')
	dm, = plt.plot(zl,pdmfidl,'-',color='firebrick',label=r' $D_M$')
	dh, = plt.plot(zl,pdhfidl,'-',color='steelblue',label=r' $zD_H$')
	#plt.legend(loc='lower left')
	plt.legend(handles=[dv,dm,dh],loc='lower left')
	#plt.ylim(.901,1.109)
	plt.xlim(0,2.4)
	plt.xlabel('redshift')
	plt.ylabel(r'BAO distance/eBOSS fiducial $D_V$')
	#plt.title('Spherically averaged or best constrained, after eBOSS')
	
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAOVMHdVftereBOSS.png')	
	plt.show()


def plotBAOrelBOSSfid_dm_eboss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[6])
			
		sdv = float(ln[7])
		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.dc(z)/rsfidcamb
		if EH:
			dvrsfid = dcosfid.dc(z)/rsfidEH
		else:
			dvrsfid = dvrsfidc
		if lab == 'BOSSDR12_1':
			plt.text(0.12,.975,'BOSS',fontsize=14,color='r')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'BOSSDR12_2':
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'LRGDR16':
			plt.text(0.35,1.025,'eBOSS',fontsize=14,color='orange')	
			plt.text(0.36,1.015,'LRGs',fontsize=14,color='orange')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='orange',markeredgecolor='k')
		if lab == 'DESY1':
			plt.text(0.85,.955,'DESY1',fontsize=14,color='lightblue')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='o',markersize=7,elinewidth=1.,color='lightblue',markeredgecolor='k')

		if lab == 'QSODR16':
			plt.text(1.4,1.07,'eBOSS',fontsize=14,color='k')	
			plt.text(1.4,1.06,'quasars',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='cyan',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR16C':
			plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
			plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		

	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.dc(z)/rsfidcamb
		pdvfidl.append(dpcos.dc(z)/prs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-',label='Planck 2018')
	plt.legend(loc='lower left')
	plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel('BAO distance/eBOSS fiducial')
	plt.title('Co-moving angular-diameter, after DR16')
		
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAO_DMaftereBOSS.png')
	plt.show()

def plotBAOrelBOSSfid_dm_boss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[6])
			
		sdv = float(ln[7])
		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.dc(z)/rsfidcamb
		if EH:
			dvrsfid = dcosfid.dc(z)/rsfidEH
		else:
			dvrsfid = dvrsfidc
		if lab == 'BOSSDR12_1':
			plt.text(0.34,.965,'BOSS',fontsize=14,color='r')
			plt.text(0.34,.955,'galaxies',fontsize=14,color='r')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'BOSSDR12_2' or lab == 'BOSSDR12_3':
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR12C':
			plt.text(2.,.985,'BOSS',fontsize=14,color='k')	
			plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='H',markersize=7,elinewidth=1.75,color='magenta',markeredgecolor='k')
		

	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.dc(z)/rsfidcamb
		pdvfidl.append(dpcos.dc(z)/prs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-',label='Planck 2018')
	plt.legend(loc='lower left')
	plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel('BAO distance/eBOSS fiducial')
	plt.title('Co-moving angular-diameter, after DR12')
		
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAO_DMafterBOSS.png')
	plt.show()


def plotBAOrelBOSSfid_Hz_eboss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[8])
			
		sdv = float(ln[9])
		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.cHz(z)/rsfidcamb
		if EH:
			dvrsfid = dcosfid.cHz(z)/rsfidEH
		else:
			dvrsfid = dvrsfidc
		if lab == 'BOSSDR12_1':
			plt.text(0.15,.975,'BOSS',fontsize=14,color='r')	
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'BOSSDR12_2':
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'LRGDR16':
			plt.text(0.8,.955,'eBOSS',fontsize=14,color='orange')	
			plt.text(0.81,.945,'LRGs',fontsize=14,color='orange')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='orange',markeredgecolor='k')

		if lab == 'QSODR16':
			plt.text(1.4,1.07,'eBOSS',fontsize=14,color='k')	
			plt.text(1.4,1.06,'quasars',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='cyan',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR16C':
			plt.text(2.,1.02,'eBOSS',fontsize=14,color='k')	
			plt.text(1.95,1.01,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		

	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.cHz(z)/rsfidcamb
		pdvfidl.append(dpcos.cHz(z)/prs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-',label='Planck 2018')
	plt.legend(loc='lower right')
	plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel('BAO distance/eBOSS fiducial')
	plt.title('Expansion rate, after DR16')
	
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAOHzaftereBOSS.png')	
	plt.show()

def plotBAOrelBOSSfid_Hz_boss():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	pdvfidl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dv = float(ln[8])
			
		sdv = float(ln[9])
		EH = int(ln[-1].strip('\n'))
		dvrsfidc = dcosfid.cHz(z)/rsfidcamb
		if EH:
			dvrsfid = dcosfid.cHz(z)/rsfidEH
		else:
			dvrsfid = dvrsfidc
		if lab == 'BOSSDR12_1':
			plt.text(0.12,.975,'BOSS',fontsize=14,color='r')	
			plt.text(0.12,.965,'galaxies',fontsize=14,color='r')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		if lab == 'BOSSDR12_2' or lab == 'BOSSDR12_3':
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='D',markersize=7,elinewidth=1.75,color='r',markeredgecolor='k')
		#if lab == 'LyADR14':
		#	plt.text(2.,.985,'eBOSS',fontsize=14,color='k')	
		#	plt.text(1.95,.975,r'Lyman-$\alpha$',fontsize=14,color='k')
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		#if lab == 'LyADR14X':
		#	plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='*',markersize=20,elinewidth=1.75,color='magenta',markeredgecolor='k')
		if lab == 'LyADR12C':
			plt.text(2.,1.02,'BOSS',fontsize=14,color='k')	
			plt.text(1.95,1.01,r'Lyman-$\alpha$',fontsize=14,color='k')
			plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='H',markersize=7,elinewidth=1.75,color='magenta',markeredgecolor='k')
		

	zls = zl.sort()
	for z in zl:
		dvrsfidc = dcosfid.cHz(z)/rsfidcamb
		pdvfidl.append(dpcos.cHz(z)/prs/dvrsfidc)
	plt.plot(zl,pdvfidl,'k-',label='Planck 2018')
	plt.legend(loc='lower right')
	plt.ylim(.901,1.109)
	plt.xlabel('redshift')
	plt.ylabel('BAO distance/eBOSS fiducial')
	plt.title('Expansion rate, after DR12')
	
	plt.savefig('/Users/ashleyross/Dropbox/BAOwebossDR16/BAOHzafterBOSS.png')	
	plt.show()


def plotBAOrelBOSSfid_dv_aniYear():
	#reads in dv/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)

	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	szl = []
	pdvfidl = []
	dvdvl = []
	sdvdvl = []
	yearl = []
	keyl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		dv = float(ln[4])
		if dv != -1:
			zl.append(z)
			szl.append(z) 
			sdv = float(ln[5])
			EH = int(ln[-1].strip('\n'))
			dvrsfidc = dcosfid.dV(z)/rsfidcamb
			if EH:
				dvrsfid = dcosfid.dV(z)/rsfidEH
			else:
				dvrsfid = dvrsfidc
			#plt.errorbar(z,dv/dvrsfid,sdv/dvrsfid,fmt='o')
			dvdvl.append(dv/dvrsfid)
			sdvdvl.append(sdv/dvrsfid)
			yearl.append(int(ln[0]))
			keyl.append(ln[2])
	zls = szl.sort()
	keyl = np.array(keyl)
	for z in szl:
		dvrsfidc = dcosfid.dV(z)/rsfidcamb
		pdvfidl.append(dpcos.dV(z)/prs/dvrsfidc)
	fig, ax = plt.subplots()
	fig.set_tight_layout(True)
	ax.plot(szl,pdvfidl,'k-',label='Plack 2018 best-fit prediction')
	ax.legend()
	zl = np.array(zl)
	dvdvl = np.array(dvdvl)
	sdvdvl = np.array(sdvdvl)
	yearl = np.array(yearl)
	def plotyear(year):
		w = yearl == year
		#for key in keyl[w]:
		for i in range(0,len(keyl[w])):
			key = keyl[w][i]
			ax.errorbar(zl[w][i],dvdvl[w][i],sdvdvl[w][i],fmt=key_dict[key][2],ms=8,elinewidth=1.75,
									 color=key_dict[key][0], mec=key_dict[key][0],
									 mfc=key_dict[key][1],capsize=4,mew=1.5)

			plt.text(key_dict[key][4], key_dict[key][5], key_dict[key][3], color=key_dict[key][0])
		#ax.errorbar(zl[w],dvdvl[w],sdvdvl[w],fmt='o',label=str(year))
		#ax.legend()
		#ax.title(str(year))
	#writer = ImageMagickFileWriter()	
	anim = FuncAnimation(fig,plotyear, frames=np.unique(yearl), interval=500)	
	#anim.save('baoyear.gif',writer='imagemagick')#, fps=10,writer='imagemagick',bitrate=-1)
	#plt.savefig('baoyear.gif')
	#plt.show()
	#plt.plot(zl,pdvfidl,'k-')


	plt.show()	


def plotBAOrelBOSSfid_H():
	#reads in c/H/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	zl = []
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		zl.append(z)
		dH = float(ln[-6])
		if dH != -1:
			sdH = float(ln[-5])
			EH = int(ln[-1].strip('\n'))
			if EH:
				dHrsfid = dcosfid.cHz(z)/rsfidEH
			else:
				dHrsfid = dcosfid.cHz(z)/rsfidcamb
			plt.errorbar(z,dH/dHrsfid,sdH/dHrsfid,fmt='o')
	zl = np.array(zl)
	#planck 2018 best LCDM
	pom = 0.3147
	pobh = 0.02233
	ph = 0.6737
	prs = 147.18
	dpcos = distance(pom,1-pom,ph,obhh=pobh)
	pdvrs = dpcos.cHz(zl)/prs
	plt.plot(zl,pdvrs,'k-')
	plt.show()	

def plotBAOrelBOSSfid_DM():
	#reads in c/H/rs measurements, plots them relative to the expectation in the BOSS fiducial cosmology
	dcosfid = distance(.31,.69,.676,obhh=0.022)
	rsfidcamb = 147.77 #r_s(z_drag) from camb and fiducial cosmology
	rsfidEH =  dcosfid.rs
	dir = '/Users/ashleyross/Dropbox/BAOwebossDR16/'	
	meas = open(dir+'BAOcomp.txt').readlines()
	#columns in meas are
	#year, ref, label, zeff, dvrs, sigdv, dmrs, sigdm, hrs, sigh, omfid, hfid, omb2fid, rsEH
	for i in range(4,len(meas)):
		ln = meas[i].split(',')
		lab = ln[2]
		#if np.isin(lab,labels):
		z = float(ln[3])
		dM = float(ln[-8])
		if dM != -1:
			sdM = float(ln[-7])
			EH = int(ln[-1].strip('\n'))
			if EH:
				dMrsfid = dcosfid.dc(z)/rsfidEH
			else:
				dMrsfid = dcosfid.dc(z)/rsfidcamb
			plt.errorbar(z,dM/dMrsfid,sdM/dMrsfid,fmt='o')
	plt.show()	


def BAOrelPlanckCurrent(wo='current',xmax=2.5,Lya=True,BOSS=False,BOSSDR12=True,MGS=True,wz=True,sdss=False,df6=True,QSODR14=True,LRGDR14=True,des=False,desy1=True,eboss=False,desi=False):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.axes as ax
	import numpy as np
	
	from numpy import loadtxt as load

	pp = PdfPages('BAOrelPlanck'+wo+'.pdf')
	plt.clf()
	fig=plt.figure()
	ax0=fig.add_subplot(1,1,1)
	ax0.yaxis.set_ticks(np.arange(0.95,1.1,.05))
	plt.minorticks_on()
	ax.Axes.axis(ax0,'auto')
	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.ylim ( 0.911, 1.09 )
	#pe = np.loadtxt('/Users/ashleyross/DR7VAC/DVordPlanck.txt').transpose()
	pe = []
	f = open('Dvrserr.dat').readlines()
	fo = open('DVout.dat','w')
	zlf = []
	dvlf = []
	rmd = 147.78
	for i in range(0,len(f)):
		ln= f[i].split()
		pe.append(float(ln[2])/float(ln[1]))
		zlf.append(.05*(1.+i))
		dvlf.append(float(ln[1])*rmd)
	#x = pe[0]
	x=np.arange(0.05,3.,0.05)
	y=np.ones(len(x))
	yu = np.zeros(len(x))
	yd = np.zeros(len(x))
	yu2 = np.zeros(len(x))
	yd2 = np.zeros(len(x))
	omwig = .27
	hwig = .71
	obp = .02205
	obw = .0225
	omp = 0.315
	hp = 0.673
	hpu = 0.679
	delh = 0.012
	delm = .017
	hpd = 0.6665
	df = distance(omp,1.-omp,hp,obhh=.022)
	dh = distance(omp,1.-omp,hp+.01,obhh=.022)
	dm = distance(omp+.01,1.-omp-.01,hp,obhh=.022)
	wigf = (alph(.44,omp,hp,obp,omwig,hwig,obw)[0],alph(.6,omp,hp,obp,omwig,hwig,obw)[0],alph(.73,omp,hp,obp,omwig,hwig,obw)[0])
	print(wigf)
	#ru = du.rs/df.rs
	#rd = dd.rs/df.rs
	#re = (ru-rd)/2.
	#fo = open('planckdvrserr_om.dat','w')
	for i in range(0,len(x)):
		#dvrsf = df.dV(x[i])/df.rs
		#dd = delm*(dm.dV(x[i])/dm.rs/dvrsf-1.)/.01
		#fo.write(str(x[i])+' '+str(dd)+'\n')
		yu[i] = y[i]+pe[i]#pe[1][i]#+dd#sqrt(pe[1][i]**2.)#+re**2.)
		yd[i] = y[i]-pe[i]#pe[1][i]#dd#sqrt(pe[1][i]**2.)#+re**2.)
	plt.fill_between(x,yd,yu,color='0.75')
	#plt.fill_between(x,yd2,yu2,color='0.75')
	plt.plot(x,1.0*y,'k-')
	#plt.plot(x,1.0*yu,'k-')
	#plt.plot(x,yd,'k-')
	if df6:
		plt.text(0.06,.925,'6dFGS',fontsize=18,color='green')
		x6 = [0.1]
		y6 = [0.987]
		e6 = [0.045]
		plt.errorbar(x6, y6,e6,fmt='s',markersize=7,elinewidth=1.75,color='green',markeredgecolor='k')
		dv = np.interp(.1,zlf,dvlf)
		fo.write(str(0.1)+' '+str(dv*.987)+' '+str(dv*0.045)+'\n')
	if BOSS:
		plt.text(0.3,.945,'BOSS DR11',fontsize=18)
		xl = [0.32,0.57]
		yl = [0.974,0.983]
		el = [0.020,0.01]
		plt.errorbar(xl, yl,el,fmt='o',markersize=7,elinewidth=1.75,color='k')
	if BOSSDR12:
		plt.text(0.38,.96,'BOSS',fontsize=18,color='purple')
		xl = [0.35,0.61]
		yl = [0.999*.995,0.985*.996]
		el = [0.01,0.01]
		plt.errorbar(xl, yl,el,fmt='D',markersize=7,elinewidth=1.75,color='purple',markeredgecolor='k')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	#plt.text(0.26,.9,'LOWZ',fontsize=18)
	#plt.text(0.51,.92,'BOSS',fontsize=18)
	#plt.text(0.5,.9,'CMASS',fontsize=18)
	if MGS:
		plt.text(0.1,1.085,'SDSS',fontsize=18,color='r')
		plt.text(0.1,1.074,'MGS',fontsize=18,color='r')
		xl = [0.15]
		yl = [1.04*.9957]
		el = [0.037]
		plt.errorbar(xl,yl,el,fmt='^',markeredgecolor='k',markersize=7,elinewidth=1.75,color='r')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if wz:
		plt.text(0.45,1.072,'WiggleZ',fontsize=18,color='.5')
		ywl = [1.061*wigf[0],1.065*wigf[1],1.039*wigf[2]]
		xwl = [.44,.6,.73] 
		ewl = [0.048,0.045,0.034]
		plt.errorbar(xwl, ywl,ewl,fmt='v',markersize=6,elinewidth=1.,color='.5',markeredgecolor='k')
		for i in range(0,len(xwl)):
			dv = np.interp(xwl[i],zlf,dvlf)
			fo.write(str(xwl[i])+' '+str(dv*ywl[i])+' '+str(dv*ewl[i])+'\n')

	if sdss:		
		xsl = [0.275,0.35]
		ysl = [0.972,0.969]
		esl = [0.027,0.018]
		plt.errorbar(xsl, ysl,esl,fmt='s',markersize=6,elinewidth=1.,mfc='w',color='k')
	#plt.text(0.95,1.07,r'$\bar{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_avg)+'\pm'+'{:.3f}'.format(delta_alpha_std)+'$')
	#plt.text(0.95,1.04,r'$\tilde{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_med)+'^{'+'{:+.3f}'.format(delta_alpha_upp)+'}_{'+'{:+.3f}'.format(delta_alpha_low)+'}$')
	if eboss:
		plt.text(1.,.97,'eBOSS projected',fontsize=18,color='r')
		xl = [0.75,0.8,1.5]
		yl = [.99,.99,.99]
		el = [0.013,0.026,0.021]
		plt.errorbar(xl,yl,el,fmt='D',markeredgecolor='r',markersize=7,elinewidth=1.75,color='r')
	if QSODR14:
		plt.text(1.2,1.045,'eBOSS quasars',fontsize=18,color='orange')
		#plt.text(1.2,.94,'(predicted)',fontsize=18,color='b')
		xl = [1.5]
		yl = [0.996]
		el = [0.044]
		plt.errorbar(xl,yl,el,color='orange',markeredgecolor='k',markersize=8,elinewidth=2.25,markeredgewidth=0,linewidth=0)
		plt.plot(xl,yl,'*',color='orange',markeredgecolor='k',markersize=18,markeredgewidth=1.5)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if LRGDR14:
		plt.text(.68,.93,'eBOSS LRGs',fontsize=18,color='firebrick')
		xl = [.71]
		yl = [0.968]
		el = [0.026]
		plt.errorbar(xl,yl,el,fmt='x',markeredgecolor='k',markersize=7,elinewidth=1.75,color='firebrick')

	if des:
		plt.text(.85,.98,'DES Y3',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.8]
		yld = [1.0]
		eld = [0.025]
		#yu = [1.025,1.025]
		#yd = [0.975,0.975]
		plt.errorbar(xld,yld,eld,fmt='h',markeredgecolor='k',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
	if desy1:
		plt.text(.85,1.015,'DES',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.8]
		yld = [0.991]
		eld = [0.041]
		yu = [1.019,1.019]
		yd = [0.981,0.981]
		plt.errorbar(xld,yld,eld,fmt='-s',markeredgecolor='steelblue',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
		plt.ylim(.9,1.1)
	if desi:
		plt.text(1.4,1.03,'DESI',fontsize=18,color='purple')		
		d = load('/Users/ashleyross/BAOpredictionscode/BAO_BigBOSS.dat').transpose()
		xl = d[0]
		yl = np.ones((len(d[0])))
		el = d[-1]
		plt.errorbar(xl,yl,el,fmt='-*',markeredgecolor='purple',markersize=7,elinewidth=1.75,color='purple')
	if Lya:
		lyf = [alph(2.33,omp,hp,obp,0.1426/.6731/.6731,.6731,0.02222)[0],alph(2.35,omp,hp,obp,.27,.7,.0227)[0]]
		print(lyf)
		lyb = [1.026/lyf[0],0.988/lyf[1]]
		lyt = (lyb[0]/.026**2.+lyb[1]/.022**2.)/(1/.026**2.+1/.022**2.)
		lyst = sqrt(1./(1/.026**2.+1/.022**2.))
		print(lyt,lyst)
		xl = [2.34]
		yl = [lyt]
		el = [lyst]
		plt.errorbar(xl,yl,el,fmt='o',color='b',markersize=8,markeredgecolor='k')
		plt.text(1.85,.97,r'(e)BOSS Ly$\alpha$',color='b',fontsize=18)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	plt.xlim ( 0.0, xmax )

	plt.xlabel ('Redshift', fontsize=18)
	#plt.ylabel (r'$(D_{\rm V}/r_{\rm d})/(D_{\rm V}/r_{\rm d})_{\rm Planck}$', fontsize=18)
	plt.ylabel (r'Distance/Distance(Planck$\Lambda$CDM)', fontsize=18)
	pp.savefig()
	pp.close()
	fo.close()
	return True

def BAOrelPlanckNear(wo='Near',xmax=2.5,Lya=True,BOSS=False,BOSSDR12e=True,MGS=True,wz=True,sdss=False,df6=True,QSO=True,QSODR14=False,ELG=True,LRG=True,LRGDR14=False,des=True,desy1=False,eboss=False,desi=False):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.axes as ax
	import numpy as np
	#from Cosmo import distance,a
	from numpy import loadtxt as load

	pp = PdfPages('BAOrelPlanck'+wo+'.pdf')
	plt.clf()
	fig=plt.figure()
	ax0=fig.add_subplot(1,1,1)
	ax0.yaxis.set_ticks(np.arange(0.95,1.1,.05))
	plt.minorticks_on()
	ax.Axes.axis(ax0,'auto')
	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.ylim ( 0.911, 1.09 )
	#pe = np.loadtxt('/Users/ashleyross/DR7VAC/DVordPlanck.txt').transpose()
	pe = []
	f = open('Dvrserr.dat').readlines()
	fo = open('DVout.dat','w')
	zlf = []
	dvlf = []
	rmd = 147.78
	for i in range(0,len(f)):
		ln= f[i].split()
		pe.append(float(ln[2])/float(ln[1]))
		zlf.append(.05*(1.+i))
		dvlf.append(float(ln[1])*rmd)
	#x = pe[0]
	x=np.arange(0.05,3.,0.05)
	y=np.ones(len(x))
	yu = np.zeros(len(x))
	yd = np.zeros(len(x))
	yu2 = np.zeros(len(x))
	yd2 = np.zeros(len(x))
	omwig = .27
	hwig = .71
	obp = .02205
	obw = .0225
	omp = 0.315
	hp = 0.673
	hpu = 0.679
	delh = 0.012
	delm = .017
	hpd = 0.6665
	df = distance(omp,1.-omp,hp,obhh=.022)
	dh = distance(omp,1.-omp,hp+.01,obhh=.022)
	dm = distance(omp+.01,1.-omp-.01,hp,obhh=.022)
	wigf = (alph(.44,omp,hp,obp,omwig,hwig,obw)[0],alph(.6,omp,hp,obp,omwig,hwig,obw)[0],alph(.73,omp,hp,obp,omwig,hwig,obw)[0])
	print(wigf)
	#ru = du.rs/df.rs
	#rd = dd.rs/df.rs
	#re = (ru-rd)/2.
	#fo = open('planckdvrserr_om.dat','w')
	for i in range(0,len(x)):
		#dvrsf = df.dV(x[i])/df.rs
		#dd = delm*(dm.dV(x[i])/dm.rs/dvrsf-1.)/.01
		#fo.write(str(x[i])+' '+str(dd)+'\n')
		yu[i] = y[i]+pe[i]#pe[1][i]#+dd#sqrt(pe[1][i]**2.)#+re**2.)
		yd[i] = y[i]-pe[i]#pe[1][i]#dd#sqrt(pe[1][i]**2.)#+re**2.)
	plt.fill_between(x,yd,yu,color='0.75')
	#plt.fill_between(x,yd2,yu2,color='0.75')
	plt.plot(x,1.0*y,'k-')
	#plt.plot(x,1.0*yu,'k-')
	#plt.plot(x,yd,'k-')
	if df6:
		plt.text(0.06,.925,'6dFGS',fontsize=18,color='green')
		x6 = [0.1]
		y6 = [0.987]
		e6 = [0.045]
		plt.errorbar(x6, y6,e6,fmt='s',markersize=7,elinewidth=1.75,color='green',markeredgecolor='k')
		dv = np.interp(.1,zlf,dvlf)
		fo.write(str(0.1)+' '+str(dv*.987)+' '+str(dv*0.045)+'\n')
	if BOSS:
		plt.text(0.3,.945,'BOSS DR11',fontsize=18)
		xl = [0.32,0.57]
		yl = [0.974,0.983]
		el = [0.020,0.01]
		plt.errorbar(xl, yl,el,fmt='o',markersize=7,elinewidth=1.75,color='k')
	if BOSSDR12e:
		plt.text(0.35,.967,'BOSS',fontsize=18,color='purple')
		xl = [0.35,0.52]
		yl = [0.999*.995,0.99*.996]
		el = [0.01,0.01]
		plt.errorbar(xl, yl,el,fmt='D',markersize=7,elinewidth=1.75,color='purple',markeredgecolor='k')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	#plt.text(0.26,.9,'LOWZ',fontsize=18)
	#plt.text(0.51,.92,'BOSS',fontsize=18)
	#plt.text(0.5,.9,'CMASS',fontsize=18)
	if MGS:
		plt.text(0.13,1.08,'SDSS',fontsize=18,color='r')
		plt.text(0.13,1.069,'MGS',fontsize=18,color='r')
		xl = [0.15]
		yl = [1.04*.9957]
		el = [0.037]
		plt.errorbar(xl,yl,el,fmt='^',markeredgecolor='k',markersize=7,elinewidth=1.75,color='r')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if wz:
		plt.text(0.45,1.072,'WiggleZ',fontsize=18,color='.5')
		ywl = [1.061*wigf[0],1.065*wigf[1],1.039*wigf[2]]
		xwl = [.44,.6,.73] 
		ewl = [0.048,0.045,0.034]
		plt.errorbar(xwl, ywl,ewl,fmt='v',markersize=6,elinewidth=1.,color='.5',markeredgecolor='k')
		for i in range(0,len(xwl)):
			dv = np.interp(xwl[i],zlf,dvlf)
			fo.write(str(xwl[i])+' '+str(dv*ywl[i])+' '+str(dv*ewl[i])+'\n')

	if sdss:		
		xsl = [0.275,0.35]
		ysl = [0.972,0.969]
		esl = [0.027,0.018]
		plt.errorbar(xsl, ysl,esl,fmt='s',markersize=6,elinewidth=1.,mfc='w',color='k')
	#plt.text(0.95,1.07,r'$\bar{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_avg)+'\pm'+'{:.3f}'.format(delta_alpha_std)+'$')
	#plt.text(0.95,1.04,r'$\tilde{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_med)+'^{'+'{:+.3f}'.format(delta_alpha_upp)+'}_{'+'{:+.3f}'.format(delta_alpha_low)+'}$')
	if eboss:
		plt.text(1.,.97,'eBOSS projected',fontsize=18,color='r')
		xl = [0.75,0.8,1.5]
		yl = [.99,.99,.99]
		el = [0.013,0.026,0.021]
		plt.errorbar(xl,yl,el,fmt='D',markeredgecolor='r',markersize=7,elinewidth=1.75,color='r')
	if QSODR14:
		plt.text(1.2,1.045,'eBOSS quasars',fontsize=18,color='orange')
		#plt.text(1.2,.94,'(predicted)',fontsize=18,color='b')
		xl = [1.5]
		yl = [0.996]
		el = [0.044]
		plt.errorbar(xl,yl,el,color='orange',markeredgecolor='k',markersize=8,elinewidth=2.25,markeredgewidth=0,linewidth=0)
		plt.plot(xl,yl,'*',color='orange',markeredgecolor='k',markersize=18,markeredgewidth=1.5)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')
	if QSO:
		plt.text(1.2,1.02,'eBOSS quasars',fontsize=18,color='orange')
		#plt.text(1.2,.94,'(predicted)',fontsize=18,color='b')
		xl = [1.5]
		yl = [1.005]
		el = [0.02]
		plt.errorbar(xl,yl,el,color='orange',markeredgecolor='k',markersize=8,elinewidth=2.25,markeredgewidth=0,linewidth=0)
		plt.plot(xl,yl,'*',color='orange',markeredgecolor='k',markersize=18,markeredgewidth=1.5)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if LRGDR14:
		plt.text(.68,.93,'eBOSS LRGs',fontsize=18,color='firebrick')
		xl = [.71]
		yl = [0.968]
		el = [0.026]
		plt.errorbar(xl,yl,el,fmt='x',markeredgecolor='k',markersize=7,elinewidth=1.75,color='firebrick')

	if LRG:
		plt.text(.68,.95,'eBOSS+BOSS LRGs',fontsize=18,color='firebrick')
		xl = [.71]
		yl = [0.975]
		el = [0.015]
		plt.errorbar(xl,yl,el,fmt='x',markeredgecolor='k',markersize=7,elinewidth=1.75,color='firebrick')

	if ELG:
		plt.text(.68,1.03,'eBOSS ELGs',fontsize=18,color='k')
		xl = [.81]
		yl = [.99]
		el = [0.03]
		plt.errorbar(xl,yl,el,fmt='o',markeredgecolor='k',markersize=7,elinewidth=1.75,color='k')
		plt.plot(xl,yl,'wo')

	if des:
		plt.text(.85,.98,'DES',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.85]
		yld = [0.991]
		eld = [0.02]
		#yu = [1.025,1.025]
		#yd = [0.975,0.975]
		plt.errorbar(xld,yld,eld,fmt='h',markeredgecolor='k',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
	if desy1:
		plt.text(.85,1.015,'DES',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.8]
		yld = [0.991]
		eld = [0.041]
		yu = [1.019,1.019]
		yd = [0.981,0.981]
		plt.errorbar(xld,yld,eld,fmt='-s',markeredgecolor='steelblue',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
		plt.ylim(.9,1.1)
	if desi:
		plt.text(1.4,1.03,'DESI',fontsize=18,color='purple')		
		d = load('/Users/ashleyross/BAOpredictionscode/BAO_BigBOSS.dat').transpose()
		xl = d[0]
		yl = np.ones((len(d[0])))
		el = d[-1]
		plt.errorbar(xl,yl,el,fmt='-*',markeredgecolor='purple',markersize=7,elinewidth=1.75,color='purple')
	if Lya:
		lyf = [alph(2.33,omp,hp,obp,0.1426/.6731/.6731,.6731,0.02222)[0],alph(2.35,omp,hp,obp,.27,.7,.0227)[0]]
		print(lyf)
		lyb = [1.026/lyf[0],0.988/lyf[1]]
		lyt = (lyb[0]/.026**2.+lyb[1]/.022**2.)/(1/.026**2.+1/.022**2.)
		lyst = sqrt(1./(1/.026**2.+1/.022**2.))
		print(lyt,lyst)
		xl = [2.34]
		yl = [lyt]
		el = [lyst*.8]
		plt.errorbar(xl,yl,el,fmt='o',color='b',markersize=8,markeredgecolor='k')
		plt.text(1.85,.97,r'(e)BOSS Ly$\alpha$',color='b',fontsize=18)
		#for i in range(0,len(xl)):
		#	dv = np.interp(xl[i],zlf,dvlf)
		#	fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	plt.xlim ( 0.0, xmax )

	plt.xlabel ('Redshift', fontsize=18)
	#plt.ylabel (r'$(D_{\rm V}/r_{\rm d})/(D_{\rm V}/r_{\rm d})_{\rm Planck}$', fontsize=18)
	plt.ylabel (r'Distance/Distance(Planck$\Lambda$CDM)', fontsize=18)
	pp.savefig()
	pp.close()
	fo.close()
	return True

	
if __name__ == "__main__":
	 BAOrelPlanckNear()		
