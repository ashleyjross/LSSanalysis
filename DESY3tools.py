dir = '/Users/ashleyross/Dropbox/DESY3/'
dirs = '/Users/ashleyross/Dropbox/BAO-DES-Y3/Data/ACF-Santi/'
import fitsio
from healpix import thphi2radec,radec2thphi#,healpix,ang2pix_ring,pix2ang_ring# #comes from AJR's healpix routines
try:
    import healpy as hp
    hpm = True
except:
    print 'no healpy, this will cause problems '
    hpm = False
from math import *
from numpy import loadtxt as load
from numpy import array,zeros
import pylab as plb
from random import random
from math import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.cm as cm

rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)


mask = '_footprint_xcorr_4096_gt22_nimgriz_'

def mkgalmapY3ac(res,zr,md,fore,wm='',syscut=''):
	gl = []
	for i in range(0,12*res*res):
		gl.append(0)
	f = fitsio.read(dir+'dr1_lss_red_'+zr+'_v0_redux.fits.gz',ext=1)
	ngt = 0
	w = 1.
	zem = 0
	fw = ''
	if fore == 'fore':
		fw = '_fore'
	#if fore == 'auto':
	#	fw = '_auto'
	if md == 'nodepth':
		md = '_none'
	else:
		md = '_'+md	
	for i in range(0,len(f)):
		ra,dec = f[i]['RA'],f[i]['DEC']
		
		if f[i]['v0'+md+fw] == 1.:
		
		#if wm != '':
		#	w = float(ln[4])
		#if z > zmin and z < zmax:
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi,nest=True)
			gl[p] += w
			ngt += w
	print len(gl),ngt
	return gl

def maskd(res,md='mode',fore='fore'):
    #degrade mask
	f = fitsio.read(dir+md+mask+fore+'.fits.gz')
	mo = []
	for i in range(0,12*res*res):
		mo.append(0)
	frac = (res/4096.)**2.
	for i in range(0,len(f)):
		p = f[i]['PIXEL']
		#mv = float(ln[1])
		#if mv > 0:
		th,phi = hp.pix2ang(4096,p,nest=True)
		a = frac
		po = hp.ang2pix(res,th,phi,nest=True)
		mo[po] += a
	fo = open(dir+'Y3'+md+fore+str(res)+'nest.dat','w')
	for i in range(0,len(mo)):
		if mo[i] > 0:
			fo.write(str(i)+' '+str(mo[i])+'\n')
	fo.close()
	return True


def mkgalodensY3ac4w(zr,res=256,tol=.2,md='mode',fore='fore'):
	from healpix import pix2ang_ring,thphi2radec
	np = 12*res*res
	ml = []
	gl = []
	for i in range(0,np):
		ml.append(0)
		gl.append(0)
		
	mf = open(dir+'Y3'+md+fore+str(res)+'nest.dat').readlines()
	print 'mask read'
	for i in range(0,len(mf)):
		p = int(mf[i].split()[0])
		ml[p] = float(mf[i].split()[1])
	gl = mkgalmapY3ac(res,zr,md,fore)
		

	ng = 0
	np = 0 
	for i in range(0,len(gl)):
		if ml[i] > tol:
			m = 0
			#if maskst == 'y':
			#th,phi = pix2ang_ring(res,i)
			if m == 0:
				ng += gl[i]
				np += ml[i]
	print ng,np
	ave = ng/np
	print ave
	cw = ''
	fo = open(dir+'galY3ac'+md+fore+zr+str(res)+'odenspczw.dat','w')
	ft = open(dir+'galY3ac'+md+fore+zr+str(res)+'rdodens.dat','w')
	no = 0		
	for i in range(0,len(ml)):
		if ml[i] > tol:
			th,phi = hp.pix2ang(res,i,nest=True)
			m = 0
			#dec = -180./pi*th+90.
			ra,dec = thphi2radec(th,phi)
			if m == 0:
				sra = sin(phi)
				cra = cos(phi)
				sdec = sin(-1.*(th-pi/2.))
				cdec = cos(-1.*(th-pi/2.))
				od = gl[i]/(ave*ml[i]) -1.
				fo.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(od)+' '+str(ml[i])+'\n')
				ft.write(str(th)+' '+str(phi)+' '+str(od)+' '+str(ml[i])+'\n')
	print no
	#ft.close()
	fo.close()
	return True
	
def plotY3compY1(zr,md='autofore',baor=(0,0)):
	if zr == '0607':
		zr1 = '0.60to0.70'
	if zr == '0708':
		zr1 = '0.70to0.80'
	if zr == '0809':
		zr1 = '0.80to0.90'
	if zr == '0910':
		zr1 = '0.90to1.00'
	pp = PdfPages(dir+'DR1compY1'+zr+'.pdf')	
	#d = load(dir+'galY3ac'+md+'fore'+zr+'5122ptPixclb6.dat').transpose()
	d = load(dirs+'dr1_lss_red_'+zr+'_v0_redux.masked.2PC_exact_fine').transpose()
	d1 = load('/Users/ashleyross/Dropbox/BAO-DES/Y1-Data/ACF_Santi/Y1redLSS_Y1_DNFv1.0_masked.dat.z'+zr1+'.w2PCF_exact_th0.30').transpose()
	thl = []
	thl1 = []
	wl = []
	wl1 = []
	for i in range(0,len(d[0])):
		if d[0][i]< baor[0] or d[0][i] > baor[1]:
			thl.append(d[0][i])
			wl.append(d[1][i])
	for i in range(0,len(d1[0])):
		if d1[0][i]< baor[0] or d1[0][i] > baor[1]:
			thl1.append(d1[0][i])
			wl1.append(d1[1][i])
	thl = np.array(thl)
	thl1 = np.array(thl1)
	wl = np.array(wl)
	wl1 = np.array(wl1)		
	plt.plot(thl,thl*wl*1.e3,thl1,thl1*wl1*1.e3)
	plt.xlabel(r'$\theta$ (degrees)')
	plt.ylabel(r'$10^3\theta w(\theta)$')
	plt.title(zr+', comparing to Y1 (orange)')
	pp.savefig()
	pp.close()
	plt.clf()
	return True

def plotY3compmd(zr,baor=(2.1,4.4)):
	pp = PdfPages(dir+'DR1compmask'+zr+'.pdf')
	d = load(dir+'galY3acmodefore'+zr+'5122ptPixclb6.dat').transpose()
	d1 = load(dir+'galY3acautofore'+zr+'5122ptPixclb6.dat').transpose()
	d2 = load(dir+'galY3acnodepthfore'+zr+'5122ptPixclb6.dat').transpose()
	thl = []
	thl1 = []
	wl = []
	wl1 = []
	thl2 = []
	wl2 = []
	for i in range(0,len(d[0])):
		if d[0][i]< baor[0] or d[0][i] > baor[1]:
			thl.append(d[0][i])
			wl.append(d[1][i])
	for i in range(0,len(d1[0])):
		if d1[0][i]< baor[0] or d1[0][i] > baor[1]:
			thl1.append(d1[0][i])
			wl1.append(d1[1][i])
	for i in range(0,len(d2[0])):
		if d2[0][i]< baor[0] or d2[0][i] > baor[1]:
			thl2.append(d2[0][i])
			wl2.append(d2[1][i])
	thl = np.array(thl)
	thl1 = np.array(thl1)
	thl2 = np.array(thl2)
	wl = np.array(wl)
	wl1 = np.array(wl1)	
	wl2 = np.array(wl2)	
	plt.plot(thl,thl*wl*1.e3,thl1,thl1*wl1*1.e3,thl2,thl2*wl2*1.e3)
	plt.xlabel(r'$\theta$ (degrees)')
	plt.ylabel(r'$10^3\theta w(\theta)$')
	plt.title(zr+', comparing three masks')
	#plt.show()
	pp.savefig()
	pp.close()
	plt.clf()
	return True
		
