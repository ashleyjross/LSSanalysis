dir = '/Users/ashleyross/Dropbox/DESY3/'
dirs = '/Users/ashleyross/Dropbox/BAO-DES-Y3/Data/ACF-Santi/'
import fitsio
from healpix import thphi2radec,radec2thphi,healpix,ang2pix_ring,pix2ang_ring# #comes from AJR's healpix routines
try:
    import healpy as hp
    hpm = True
except:
    print 'no healpy, this will cause problems '
    hpm = False
from math import *
import numpy as np
from numpy import loadtxt as load
from numpy import array,zeros
import pylab as plb
from random import random
from math import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.cm as cm

rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)


#mask = '_footprint_xcorr_4096_gt22_nimgriz_' #used for DR1
#mask = 'mask_v0_lssred'
#mask = '_no2massfaint'
mask = '_c08'

def compY1Y3():
	f1 = fitsio.read('/Users/ashleyross/Dropbox/DESY3/Y1mag22.fits')
	fs_1 = fitsio.read('/Users/ashleyross/Dropbox/DESY3/Y3red_1.fits')
	fs_2 = fitsio.read('/Users/ashleyross/Dropbox/DESY3/Y3red_2.fits')
	maskf = open('/Users/ashleyross/DESY1/mask_Y1redBAO_mean_z_bpz_VFF_4096ring.dat')
	npix = 12*4096*4096
	mask = []
	pixl1 = np.zeros(f1.size)
	pixl3_1 = np.zeros(fs_1.size)
	pixl3_2 = np.zeros(fs_2.size)
	for i in range(0,npix):
		mask.append(0)
	for line in maskf:
		pix = int(float(line.split()[0]))
		mask[pix] = 1
	print 'mask done'	
	ng = 0
	ng3 = 0
	for i in range(0,f1.size):
		ra,dec = f1[i]['RA'],f1[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = hp.ang2pix(4096,th,phi)
		if mask[int(p)] == 1:
			pixl1[i] = 1
	print 'f1 done'
	for i in range(0,fs_1.size):
		ra,dec = fs_1[i]['RA'],fs_1[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = hp.ang2pix(4096,th,phi)
		if mask[int(p)] == 1:
			pixl3_1[i] = 1
	print 'f3_1 done'		
	for i in range(0,fs_2.size):
		ra,dec = fs_2[i]['RA'],fs_2[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = hp.ang2pix(4096,th,phi)
		if mask[int(p)] == 1:
			pixl3_2[i] = 1
	print 'f3_2 done'		
	w1 = (pixl1 == 1) & (f1['MAG_AUTO_I']-f1['MAG_AUTO_Z'] +2.*(f1['MAG_AUTO_R']-f1['MAG_AUTO_I'])>1.7) & (f1['MAG_AUTO_I'] >17.5)
	f1m = f1[w1]
	w3_1 = (pixl3_1 == 1)
	f3_1m = fs_1[w3_1]
	w3_2 = (pixl3_2 == 1)
	f3_2m = fs_2[w3_2]
	print f1m.size, f3_1m.size+f3_2m.size 

	

def compSOFMA():
	fma = fitsio.read('/Users/ashleyross/Dropbox/DESY3/allgal22magauto_1.fits')
	fs = fitsio.read('/Users/ashleyross/Dropbox/DESY3/allgal22_1.fits')
	f1 = fitsio.read('/Users/ashleyross/Dropbox/DESY3/Y1mag22.fits')
	facs = 5.3
	facm = 6.
	fac1 = 3.
	print np.mean(fma['SOF_CM_MAG_CORRECTED_I']),np.mean(fma['MAG_AUTO_I']-1.569*fma['EBV_SFD98'])
	print fma.size*facm,fs.size*facs,f1.size*fac1,f1.size
	wa = (fma['SOF_CM_MAG_CORRECTED_I']-fma['SOF_CM_MAG_CORRECTED_Z'] +2.*(fma['SOF_CM_MAG_CORRECTED_R']-fma['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fma['SOF_CM_MAG_CORRECTED_I']>17.5)
	ws =  (fs['SOF_CM_MAG_CORRECTED_I']-fs['SOF_CM_MAG_CORRECTED_Z'] +2.*(fs['SOF_CM_MAG_CORRECTED_R']-fs['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fs['SOF_CM_MAG_CORRECTED_I']>17.5)
	w1 = (f1['MAG_AUTO_I']-f1['MAG_AUTO_Z'] +2.*(f1['MAG_AUTO_R']-f1['MAG_AUTO_I'])>1.7) & (f1['MAG_AUTO_I'] >17.5)
	print fma[wa].size*6.,fs[ws].size*facs,f1[w1].size*fac1,f1[w1].size
	wa = (fma['SOF_CM_MAG_CORRECTED_I']-fma['SOF_CM_MAG_CORRECTED_Z'] +2.*(fma['SOF_CM_MAG_CORRECTED_R']-fma['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fma['SOF_CM_MAG_CORRECTED_I']>17.5)\
	& ((fma['SOF_CM_MAG_CORRECTED_G'] - fma['SOF_CM_MAG_CORRECTED_R'])>-1.) & ((fma['SOF_CM_MAG_CORRECTED_G'] - fma['SOF_CM_MAG_CORRECTED_R'])<3.)
	ws =  (fs['SOF_CM_MAG_CORRECTED_I']-fs['SOF_CM_MAG_CORRECTED_Z'] +2.*(fs['SOF_CM_MAG_CORRECTED_R']-fs['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fs['SOF_CM_MAG_CORRECTED_I']>17.5)\
	& ((fs['SOF_CM_MAG_CORRECTED_G'] - fs['SOF_CM_MAG_CORRECTED_R'])>-1.) & ((fs['SOF_CM_MAG_CORRECTED_G'] - fs['SOF_CM_MAG_CORRECTED_R'])<3.)

	w1 = (f1['MAG_AUTO_I']-f1['MAG_AUTO_Z'] +2.*(f1['MAG_AUTO_R']-f1['MAG_AUTO_I'])>1.7) & (f1['MAG_AUTO_I'] >17.5)\
	& ((f1['MAG_AUTO_G'] - f1['MAG_AUTO_R'])>-1.) & ((f1['MAG_AUTO_G'] - f1['MAG_AUTO_R'])<3.)

	print fma[wa].size*6.,fs[ws].size*facs,f1[w1].size*fac1,f1[w1].size

	wa = (fma['SOF_CM_MAG_CORRECTED_I']-fma['SOF_CM_MAG_CORRECTED_Z'] +2.*(fma['SOF_CM_MAG_CORRECTED_R']-fma['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fma['SOF_CM_MAG_CORRECTED_I']>17.5)\
	& ((fma['SOF_CM_MAG_CORRECTED_G'] - fma['SOF_CM_MAG_CORRECTED_R'])>-1.) & ((fma['SOF_CM_MAG_CORRECTED_G'] - fma['SOF_CM_MAG_CORRECTED_R'])<3.)\
	& ((fma['SOF_CM_MAG_CORRECTED_R'] - fma['SOF_CM_MAG_CORRECTED_I'])>-1.) & ((fma['SOF_CM_MAG_CORRECTED_R'] - fma['SOF_CM_MAG_CORRECTED_I'])<2.5)

	ws =  (fs['SOF_CM_MAG_CORRECTED_I']-fs['SOF_CM_MAG_CORRECTED_Z'] +2.*(fs['SOF_CM_MAG_CORRECTED_R']-fs['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fs['SOF_CM_MAG_CORRECTED_I']>17.5)\
	& ((fs['SOF_CM_MAG_CORRECTED_G'] - fs['SOF_CM_MAG_CORRECTED_R'])>-1.) & ((fs['SOF_CM_MAG_CORRECTED_G'] - fs['SOF_CM_MAG_CORRECTED_R'])<3.)\
	& ((fs['SOF_CM_MAG_CORRECTED_R'] - fs['SOF_CM_MAG_CORRECTED_I'])>-1.) & ((fs['SOF_CM_MAG_CORRECTED_R'] - fs['SOF_CM_MAG_CORRECTED_I'])<2.5)

	w1 = (f1['MAG_AUTO_I']-f1['MAG_AUTO_Z'] +2.*(f1['MAG_AUTO_R']-f1['MAG_AUTO_I'])>1.7) & (f1['MAG_AUTO_I'] >17.5)\
	& ((f1['MAG_AUTO_G'] - f1['MAG_AUTO_R'])>-1.) & ((f1['MAG_AUTO_G'] - f1['MAG_AUTO_R'])<3.)\
	& ((f1['MAG_AUTO_R'] - f1['MAG_AUTO_I'])>-1.) & ((f1['MAG_AUTO_R'] - f1['MAG_AUTO_I'])<2.5)
	print fma[wa].size*6.,fs[ws].size*facs,f1[w1].size*fac1,f1[w1].size

	wa = (fma['SOF_CM_MAG_CORRECTED_I']-fma['SOF_CM_MAG_CORRECTED_Z'] +2.*(fma['SOF_CM_MAG_CORRECTED_R']-fma['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fma['SOF_CM_MAG_CORRECTED_I']>17.5)\
	& ((fma['SOF_CM_MAG_CORRECTED_G'] - fma['SOF_CM_MAG_CORRECTED_R'])>-1.) & ((fma['SOF_CM_MAG_CORRECTED_G'] - fma['SOF_CM_MAG_CORRECTED_R'])<3.)\
	& ((fma['SOF_CM_MAG_CORRECTED_R'] - fma['SOF_CM_MAG_CORRECTED_I'])>-1.) & ((fma['SOF_CM_MAG_CORRECTED_R'] - fma['SOF_CM_MAG_CORRECTED_I'])<2.5)\
	& ((fma['SOF_CM_MAG_CORRECTED_I'] - fma['SOF_CM_MAG_CORRECTED_Z'])>-1.) & ((fma['SOF_CM_MAG_CORRECTED_I'] - fma['SOF_CM_MAG_CORRECTED_Z'])<2.)
	
	ws =  (fs['SOF_CM_MAG_CORRECTED_I']-fs['SOF_CM_MAG_CORRECTED_Z'] +2.*(fs['SOF_CM_MAG_CORRECTED_R']-fs['SOF_CM_MAG_CORRECTED_I'])>1.7) & (fs['SOF_CM_MAG_CORRECTED_I']>17.5)\
	& ((fs['SOF_CM_MAG_CORRECTED_G'] - fs['SOF_CM_MAG_CORRECTED_R'])>-1.) & ((fs['SOF_CM_MAG_CORRECTED_G'] - fs['SOF_CM_MAG_CORRECTED_R'])<3.)\
	& ((fs['SOF_CM_MAG_CORRECTED_R'] - fs['SOF_CM_MAG_CORRECTED_I'])>-1.) & ((fs['SOF_CM_MAG_CORRECTED_R'] - fs['SOF_CM_MAG_CORRECTED_I'])<2.5)\
	& ((fs['SOF_CM_MAG_CORRECTED_I'] - fs['SOF_CM_MAG_CORRECTED_Z'])>-1.) & ((fs['SOF_CM_MAG_CORRECTED_I'] - fs['SOF_CM_MAG_CORRECTED_Z'])<2.)

	f1 = fitsio.read('/Users/ashleyross/Dropbox/DESY3/Y1mag22_chbpz.fits')

	w1 = (f1['mag_auto_i']-f1['mag_auto_z'] +2.*(f1['mag_auto_r']-f1['mag_auto_i'])>1.7) & (f1['mag_auto_i'] >17.5)\
	& ((f1['mag_auto_g'] - f1['mag_auto_r'])>-1.) & ((f1['mag_auto_g'] - f1['mag_auto_r'])<3.)\
	& ((f1['mag_auto_r'] - f1['mag_auto_i'])>-1.) & ((f1['mag_auto_r'] - f1['mag_auto_i'])<2.5)\
	& ((f1['mag_auto_i'] - f1['mag_auto_z'])>-1.) & ((f1['mag_auto_i'] - f1['mag_auto_z'])<2.)
	print fma[wa].size*6.,fs[ws].size*facs,f1[w1].size*fac1,f1[w1].size
	fmaw = fma[wa]
	fsw = fs[ws]
	f1w = f1[w1]
	maskf = open('/Users/ashleyross/DESY1/mask_Y1redBAO_mean_z_bpz_VFF_4096ring.dat')
	npix = 12*4096*4096
	mask = []
	for i in range(0,npix):
		mask.append(0)
	for line in maskf:
		pix = int(float(line.split()[0]))
		mask[pix] = 1
	ng = 0
	ng3 = 0
	for i in range(0,f1w.size):
		ra,dec = f1w[i]['ra'],f1w[i]['dec']
		th,phi = radec2thphi(ra,dec)
		p = hp.ang2pix(4096,th,phi)
		if mask[int(p)] == 1:
			ng += 1.
	for i in range(0,fsw.size):
		ra,dec = fsw[i]['RA'],fsw[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = hp.ang2pix(4096,th,phi)
		if mask[int(p)] == 1:
			ng3 += 1.

	print ng,ng3

	



	ws = (fsw['SOF_CM_MAG_CORRECTED_I']<19.+(3*fsw['DNF_ZMEAN_SOF']))
	w1 = (f1w['mag_auto_i'] < 19.+(3.*f1w['mean_z_bpz']))
	print fsw[ws].size*facs,f1w[w1].size*fac1,f1w[w1].size

	ws = (fsw['SOF_CM_MAG_CORRECTED_I']<19.+(3*fsw['DNF_ZMEAN_SOF'])) & (fsw['DNF_ZMEAN_SOF'] > 0.6) & (fsw['DNF_ZMEAN_SOF'] < 1.)
	w1 = (f1w['mag_auto_i'] < 19.+(3.*f1w['mean_z_bpz']))  & (f1w['mean_z_bpz']>0.6) & (f1w['mean_z_bpz']<1.)
	print fsw[ws].size*facs,f1w[w1].size*fac1,f1w[w1].size

	wa = (fma['DNF_ZMEAN_SOF'] > 0.6) & (fma['DNF_ZMEAN_SOF'] < 1.0)
	ws = (fs['DNF_ZMEAN_SOF'] > 0.6) & (fs['DNF_ZMEAN_SOF'] < 1.0)

	print fma[wa].size*6.,fs[ws].size*facs
	wa = (fma['DNF_ZMEAN_SOF'] > 0.6) & (fma['DNF_ZMEAN_SOF'] < 1.0) & (fma['SOF_CM_MAG_CORRECTED_I']-fma['SOF_CM_MAG_CORRECTED_Z'] +2.*(fma['SOF_CM_MAG_CORRECTED_R']-fma['SOF_CM_MAG_CORRECTED_I'])>1.7)
	ws = (fs['DNF_ZMEAN_SOF'] > 0.6) & (fs['DNF_ZMEAN_SOF'] < 1.0) & (fs['SOF_CM_MAG_CORRECTED_I']-fs['SOF_CM_MAG_CORRECTED_Z'] +2.*(fs['SOF_CM_MAG_CORRECTED_R']-fs['SOF_CM_MAG_CORRECTED_I'])>1.7)
	print fma[wa].size*6.,fs[ws].size*facs
	wa = (fma['DNF_ZMEAN_SOF'] > 0.6) & (fma['DNF_ZMEAN_SOF'] < 1.0) & (fma['SOF_CM_MAG_CORRECTED_I']-fma['SOF_CM_MAG_CORRECTED_Z'] +2.*(fma['SOF_CM_MAG_CORRECTED_R']-fma['SOF_CM_MAG_CORRECTED_I'])>1.7)\
	& (fma['MAG_AUTO_I']-1.569*fma['EBV_SFD98'] <19.+(3*fma['DNF_ZMEAN_SOF']))
	ws = (fs['DNF_ZMEAN_SOF'] > 0.6) & (fs['DNF_ZMEAN_SOF'] < 1.0) & (fs['SOF_CM_MAG_CORRECTED_I']-fs['SOF_CM_MAG_CORRECTED_Z'] +2.*(fs['SOF_CM_MAG_CORRECTED_R']-fs['SOF_CM_MAG_CORRECTED_I'])>1.7)\
	& (fs['SOF_CM_MAG_CORRECTED_I']<19.+(3*fs['DNF_ZMEAN_SOF']))
	print fma[wa].size*6.,fs[ws].size*facs
	print np.mean(fma['SOF_CM_MAG_CORRECTED_I'][wa]),np.mean(fma['MAG_AUTO_I'][wa]-1.569*fma['EBV_SFD98'][wa])


def mkgalmapY3ac(res,zr,gz='.gz',md='',fore='',wm='',syscut=''):
	gl = []
	for i in range(0,12*res*res):
		gl.append(0)
	#f = fitsio.read(dir+'dr1_lss_red_'+zr+'_v0_redux.fits.gz',ext=1)
	f = fitsio.read(dir+'test'+zr+mask+'.fits'+gz,ext=1)
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
	#else:
	#	md = '_'+md	
	for i in range(0,len(f)):
		ra,dec = f[i]['RA'],f[i]['DEC']
		
		#if f[i]['v0'+md+fw] == 1.:
		
		#if wm != '':
		#	w = float(ln[4])
		#if z > zmin and z < zmax:
		th,phi = radec2thphi(ra,dec)
		p = hp.ang2pix(res,th,phi,nest=True)
		gl[p] += w
		ngt += w
	print len(gl),ngt
	return gl

def calczerr(zr,md='',fore='',wm='',syscut=''):
	
	f = fitsio.read(dir+'test'+zr+mask+'.fits',ext=1)
	dtot = sum((f['DNF_ZMEAN_SOF']-f['DNF_ZMC_SOF'])**2.)
	detot = sum(f['DNF_ZSIGMA_SOF'])
	n = float(len(f))
	print sqrt(dtot/n),detot/n
	dnoout = 0
	nno = 0
	for i in range(0,len(f)):
		zmc = f[i]['DNF_ZMC_SOF']
		if zmc > .2 and zmc < 1.5:
			dnoout += (f[i]['DNF_ZMEAN_SOF']-zmc)**2.
			nno += 1.
	print nno,n		
	print sqrt(dnoout/nno)		
	return True

def maskd(res,gz='.gz'):
    #degrade mask
	f = fitsio.read(dir+'mask'+mask+'_lssred.fits'+gz)
	mo = []
	for i in range(0,12*res*res):
		mo.append(0)
	frac = (res/4096.)**2.
	for i in range(0,len(f)):
		p = f[i]['PIXEL']
		#mv = float(ln[1])
		#if mv > 0:
		th,phi = hp.pix2ang(4096,p,nest=True)
		#a = frac
		po = hp.ang2pix(res,th,phi,nest=True)
		mo[po] += frac
	fo = open(dir+'Y3mask'+mask+'lssred'+str(res)+'nest.dat','w')
	for i in range(0,len(mo)):
		if mo[i] > 0:
			fo.write(str(i)+' '+str(mo[i])+'\n')
	fo.close()
	return True


def mkgalodensY3ac4w(zr,res=256,tol=.2,md='',fore=''):
	from healpix import pix2ang_ring,thphi2radec
	np = 12*res*res
	ml = []
	gl = []
	for i in range(0,np):
		ml.append(0)
		gl.append(0)
		
	mf = open(dir+'Y3mask'+mask+'lssred'+str(res)+'nest.dat').readlines()
	print 'mask read'
	for i in range(0,len(mf)):
		p = int(mf[i].split()[0])
		ml[p] = float(mf[i].split()[1])
	gl = mkgalmapY3ac(res,zr)#,md,fore)
		

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
	fo = open(dir+'galY3'+mask+'lssred'+zr+str(res)+'odenspczw.dat','w')
	ft = open(dir+'galY3'+mask+'lssred'+zr+str(res)+'rdodens.dat','w')
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

def mkgalodensY3ac4w_split(zr,res=256,tol=.2,splt=1.,md='',fore=''):
	from healpix import pix2ang_ring,thphi2radec
	np = 12*res*res
	ml = []
	gl = []
	for i in range(0,np):
		ml.append(0)
		gl.append(0)
		
	mf = open(dir+'Y3mask'+mask+'lssred'+str(res)+'nest.dat').readlines()
	print 'mask read'
	for i in range(0,len(mf)):
		p = int(mf[i].split()[0])
		ml[p] = float(mf[i].split()[1])
	gl = mkgalmapY3ac(res,zr)#,md,fore)
		

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
	gll = zeros((len(gl)))
	glh = zeros((len(gl)))
	ngl = 0
	ngh = 0
	for i in range(0,len(gl)):
		if ml[i] > tol:
			if gl[i]/ml[i] > ave*splt:
				glh[i] = gl[i]
				ngh += gl[i]
			else:
				gll[i] = gl[i]
				ngl += gl[i]
	avel = ngl/np
	aveh = ngh/np			
	cw = ''
	fol = open(dir+'galY3'+mask+'lssred'+zr+str(res)+'lowodenspczw.dat','w')
	ftl = open(dir+'galY3'+mask+'lssred'+zr+str(res)+'lowrdodens.dat','w')
	foh = open(dir+'galY3'+mask+'lssred'+zr+str(res)+'highodenspczw.dat','w')
	fth = open(dir+'galY3'+mask+'lssred'+zr+str(res)+'highrdodens.dat','w')
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
				odl = gll[i]/(avel*ml[i]) -1.
				odh = glh[i]/(aveh*ml[i]) -1.
				fol.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(odl)+' '+str(ml[i])+'\n')
				ftl.write(str(th)+' '+str(phi)+' '+str(odl)+' '+str(ml[i])+'\n')
				foh.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(odh)+' '+str(ml[i])+'\n')
				fth.write(str(th)+' '+str(phi)+' '+str(odh)+' '+str(ml[i])+'\n')
	print no
	#ft.close()
	fol.close()
	foh.close()
	return True

	
def plotY3compY1(zr,b=1.8,res=512,baor=(0,0),f=.82,offset=0):
	if zr == '0607':
		zr1 = '0.60to0.70'
	if zr == '0708':
		zr1 = '0.70to0.80'
	if zr == '0809':
		zr1 = '0.80to0.90'
	if zr == '0910':
		zr1 = '0.90to1.00'
	pp = PdfPages(dir+'Y3compY1'+zr+str(res)+mask+'.pdf')	
	d = load(dir+'galY3'+mask+'lssred'+zr+str(res)+'2ptPixclb6.dat').transpose()
	#d = load(dirs+'dr1_lss_red_'+zr+'_v0_redux.masked.2PC_exact_fine').transpose()
	d1 = load('/Users/ashleyross/Dropbox/BAO-DES/Y1-Data/ACF_Santi/Y1redLSS_Y1_DNFv1.0_masked.dat.z'+zr1+'.w2PCF_exact_th0.30').transpose()
	dr1 = load(dir+'galY3acautofore'+zr+'5122ptPixclb6.dat').transpose()
	beta = f/b
	#		w1 = self.wthl[ind][1][indd]+self.beta*(2/3.*self.wthl[ind][1][indd]+4/3.*self.wthl[ind][2][indd])+(self.beta)**2.*(.2*self.wthl[ind][1][indd]+4/7.*self.wthl[ind][2][indd]+8/35.*self.wthl[ind][3][indd])
	
	if zr == '0708':
		zb = '2'
	if zr == '0607':
		zb = '1'
	if zr == '0809':
		zb = '3'
	if zr == '0910':
		zb = '4'
	dth = load('/Users/ashleyross/Dropbox/BAO-DES-Y3/Templates/w_template/wtheta_model_D5.2_TD1_MockData5_cosmo1_dtheta0.01_phi1_nzbins4_Ver4bin_'+zb+'.txt').transpose() #blinded template from DR1
	thl = []
	thl1 = []
	wl = []
	wl1 = []
	wl2 = []
	for i in range(0,len(d[0])):
		if d[0][i]< baor[0] or d[0][i] > baor[1]:
			thl.append(d[0][i])
			wl.append(d[1][i])
			wl2.append(dr1[1][i])
	for i in range(0,len(d1[0])):
		if d1[0][i]< baor[0] or d1[0][i] > baor[1]:
			thl1.append(d1[0][i])
			wl1.append(d1[1][i])
			
	thl = np.array(thl)
	thl1 = np.array(thl1)
	wl = np.array(wl)
	wl1 = np.array(wl1)		
	plt.plot(thl,thl*(wl-offset)*1.e3,thl1,thl1*wl1*1.e3)#,thl,thl*wl2*1.e3)
	thv = b*b*(dth[1]+beta*(2/3.*dth[1]+4/3.*dth[2])+beta**2.*(.2*dth[1]+4/7.*dth[2]+8/35.*dth[3]))
	#plt.plot(dth[0],dth[0]*thv*1.e3,'k--')
	plt.xlabel(r'$\theta$ (degrees)')
	plt.ylabel(r'$10^3\theta w(\theta)$')
	plt.title(zr+', comparing to Y1 (orange)')
	plt.xlim(0,6)
	pp.savefig()
	pp.close()
	plt.clf()
	return True

def plotY3compmd(zr,baor=(0,0)):
	pp = PdfPages(dir+'DR1compmask'+zr+'.pdf')
	d = load(dir+'galY3_no2massfaintlssred'+zr+'5122ptPixclb6.dat').transpose()
	d1 = load(dir+'galY3acautofore'+zr+'5122ptPixclb6.dat').transpose()
	d2 = load(dir+'galY3mask_v0_lssred'+zr+'5122ptPixclb6.dat').transpose()
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
	#plt.plot(thl,thl*wl*1.e3,thl1,thl1*wl1*1.e3,thl2,thl2*wl2*1.e3)
	plt.plot(thl,thl*wl*1.e3,thl2,thl2*wl2*1.e3)
	plt.xlabel(r'$\theta$ (degrees)')
	plt.ylabel(r'$10^3\theta w(\theta)$')
	plt.title(zr+', comparing two masks')
	#plt.show()
	pp.savefig()
	pp.close()
	plt.clf()
	return True
		
def ngvext(file,mask,res=256,mc=.8,extmax=.15,nbin=10):
	extmap = np.loadtxt('maps/healSFD_r_256_fullsky.dat')
	extc = 2.751
	h = healpix()
	ml = zeros((4096*4096*12))
	mask = fitsio.read(dir+mask+'.fits.gz')
	for i in range(0,len(mask)):
		ml[mask[i]['PIXEL']] = mask[i]['SIGNAL']
	data = fitsio.read(dir+file+'.fits.gz')
	ngl = np.zeros((12.*res*res))
	ngt = 0
	for i in range(0,len(data)):
		if data[i]['HPIX_4096'] > mc:
			th,phi = radec2thphi(data[i]['RA'],data[i]['DEC'])
			pix = h.ang2pix_nest(res,th,phi)
			ngl[pix] += 1.
			ngt += 1.
	print ngt, len(data)
	mlr = zeros((12.*res*res))
	for i in range(0,len(mask)):
		if 	mask[i]['SIGNAL'] > mc:
			th,phi = h.pix2ang_nest(4096,mask[i]['PIXEL'])
			pix = h.ang2pix_nest(res,th,phi)
			mlr[pix] += (res/4096.)**2.
	ave = sum(ngl)/sum(mlr)
	print ave
	bing = zeros((nbin))
	binr = zeros((nbin))
	
	for i in range(0,12*res*res):
		extv = extmap[i]/extc
		if extv < extmax:
			bin = int(extv/extmax*nbin) 		
			bing[bin] += ngl[i]
			binr[bin] += mlr[i]		
	print bing, binr,bing/(binr*ave)
	fo = open(dir+file+'vsext.dat','w')
	for i in range(0,nbin):
		extv = i*extmax/float(nbin)+extmax/float(2.*nbin)
		fo.write(str(extv)+' '+str(bing[i]/(binr[i]*ave))+' '+str(sqrt(bing[i])/(binr[i]*ave))+'\n')
	fo.close()
	return True
			