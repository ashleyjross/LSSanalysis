import fitsio #needed to read data
import numpy as np
from math import *
from optimize import fmin
import os
import pyfits

#try: mkesample_ver = os.environ['MKESAMPLE_VER']
#except: mkesample_ver = None
#mkesample_ver='1.0'
#dirsys = os.environ['MKESAMPLE_DIR']+'/inputFiles/' #change to local directory where maps are
#dir = os.environ['EBOSS_LSS_CATALOGS']+'/'+mkesample_ver+'/' #change to where your catalog files are
#dir = os.environ['MKESAMPLE_SCRATCH']+'/lsscats/'+mkesample_ver+'/' #change to where your catalog files are
dir = '/Users/ashleyross/fitsfiles/'
outdir = '/Users/ashleyross/eboss/'
#outdir=os.environ['EBOSS_LSS_CATALOGS']+'/'+mkesample_ver+'/' #to where files will be written

def ngvsys_ran(sampl,NS,ver,sys,sysmin,sysmax,zmin,zmax,band=-1,wm=''):
	#sample is the sample being used, e.g. 'lrg'
	#uses files with imaging properties filled
	#NS is either 'N' or 'S'
	#ver is the version, e.g., 'v1.0_IRt'
	#sys is a string containing the name of the systematic to be tested
	#sysmin is the minimum value of the systematic to be tested, ~25 is a good value to use for stars
	#sysmax is the maximum value of the systematic to be tested, ~200 is a good value to use for stars
	#res is Nside for the healpix map, default should be 512 for sky,seeing,air and 256 for extinction and stars
	#zmin is the minimum redshift to use, 0.6 is minimum used for lrgs, 0.9 for qsos
	#zmax is the maximum redshift to used, 1.0 for lrgs and 2.2 for qsos
	
	stl = []
	wstl = []
	errl = []
	if wm == 'wgdepthext': #this is used to set the weights for extinction
		b,m = findlinmb(sampl,NS,ver,'IMAGE_DEPTH_EXT1',.8,2.2,wm='nosys')
	if wm == 'wgdepthextext': #this can be used to test that the weights work for everything
		b,m = findlinmb(sampl,NS,ver,'IMAGE_DEPTH_EXT1',.8,2.2)
		be,me = findlinmb(sampl,NS,ver,'EB_MINUS_V-1',.8,2.2,wm='wgdepthext')
		print 'fit parameters for g-band depth'
		print b,m
		print 'fit parameters for EB_MINUS_V (after weighting for depth)'
		print be,me
	binng = []
	binnr = []
	nsysbin = 10 #10 bins are being used
	for i in range(0,nsysbin):
		binng.append(0)
		binnr.append(0)

	
	f = fitsio.read(dir+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.ran.fits') #read random file
	print 'reading randoms from ' + dir+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.ran.fits'
	nr = 0
	sysm = float(nsysbin)/(sysmax-sysmin)
	bs = 0
	bsr = 0
	extc = [4.239,3.303,2.285,1.698,1.263]
	for i in range (0,len(f)):
		if band == -1:
			sysv = f[i][sys]
		else:
			if sys == 'IMAGE_DEPTH_EXT':
				sysv = f[i]['IMAGE_DEPTH'][band]
				sysv = luptm(sysv,band)-extc[band]*f[i]['EB_MINUS_V']
			else:
				sysv = f[i][sys][band]	 
				if sys == 'IMAGE_DEPTH':
					sysv = luptm(sysv,band)

		bins = int((sysv-sysmin)*sysm)
		if bins >= 0 and bins < nsysbin:
			binnr[bins] += 1.
		else:
			bsr += 1.

		nr += 1.
	print nr

	f = fitsio.read(dir+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.dat.fits') #read quasar file
	print 'reading quasars from: ' + dir+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.dat.fits'
	no = 0
	zm = 0
	nt = 0
	for i in range (0,len(f)):
		z = f[i]['Z']
		gc = True
		um = f[i]['MODELMAG'][0]-f[i]['EXTINCTION'][0]
		gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
		rm = f[i]['MODELMAG'][2]-f[i]['EXTINCTION'][2]
		im = f[i]['MODELMAG'][3]-f[i]['EXTINCTION'][3]
		if z > zmin and z < zmax:
			no += 1
			#w = 1.
			#if wm == '':
			
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT'] #weight in current catalog
			if wm == 'nosys':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP'] #compare to result with no systematic weights
			if sys == 'SKYFLUX':
				sys = 'SKY_FLUX'
			if band == -1:
				sysv = f[i][sys]
			else:
				if sys == 'IMAGE_DEPTH_EXT':
					sysv = f[i]['IMAGE_DEPTH'][band]
					sysv = luptm(sysv,band)-extc[band]*f[i]['EB_MINUS_V']
				else:
					sysv = f[i][sys][band]	 
					if sys == 'IMAGE_DEPTH':
						sysv = luptm(sysv,band)

			bins = int((sysv-sysmin)*sysm)
			if wm == 'wgdepthext':
				sysw = f[i]['IMAGE_DEPTH'][1]
				sysw = luptm(sysw,1)-extc[1]*f[i]['EB_MINUS_V']
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*(1./(b+m*sysw))
			if wm == 'wgdepthextext':
				sysw = f[i]['IMAGE_DEPTH'][1]
				sysw = luptm(sysw,1)-extc[1]*f[i]['EB_MINUS_V']
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*(1./(b+m*sysw))
				ext = f[i]['EB_MINUS_V']
				w = w*(1./(be+me*ext))
					
			if bins >= 0 and bins < nsysbin:
				binng[bins] += 1.*w
			else:
				bs += w #count numbers outside of sysmin/sysmax

			zm += w*z
			nt += w
	print 'total number, weighted number'
	print no,nt
	print 'mean redshift'
	print zm/nt

	print 'total number of randoms/objects '+str(nr)+'/'+str(nt)
	print 'number of randoms/objects outside tested range '+str(bsr)+'/'+str(bs)			
	ave = nt/nr
	print 'average number of objects per random is '+ str(ave)
	fs = open(dir+'n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'v'+sys+str(band)+'.dat','w')
	xl = []
	yl = []
	el = []
	for i in range(0,nsysbin):
		sysv = sysmin + 1./(2.*sysm) + i/sysm
		if binnr[i] > 0:
			ns = binng[i]/binnr[i]/ave
			nse = sqrt(binng[i]/(binnr[i])**2./(ave)**2.+(binng[i]/ave)**2./(binnr[i])**3.) #calculate poisson error
		else:
			ns = 1. #write out 1.0 1.0 if no pixels at given value of sys
			nse = 1.		
		fs.write(str(sysv)+' '+str(ns)+' '+str(nse)+'\n')
		xl.append(sysv)
		yl.append(ns)
		el.append(nse)
	fs.close()
	chin = sum((np.array(yl)-1.)**2./np.array(el)**2.)
	print 'chi2 for null test'
	print chin
	xl = np.array(xl)
	yl = np.array(yl)
	el = np.array(el)
	#plotvssys_simp(xl,yl,el,sys)
	return xl,yl,el



def findlinmb(sampl,NS,ver,sys,zmin,zmax,wm='',res=''):
	#finds linear fit parameters (depth or stellar density relationships hard-coded to expect given resolutions)
	if sys == 'star':
		res = '256'
	if sys == 'depth':
		res = '512'
	d = np.loadtxt(dir+'n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+str(res)+'v'+sys+'.dat').transpose()	
	lf = linfit(d[0],d[1],d[2])
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	return b0,m0

class linfit:
	def __init__(self,xl,yl,el):
		self.xl = xl
		self.yl = yl
		self.el = el
		
	def chilin(self,bml):
		chi = 0
		b = bml[0]
		m = bml[1]
		for i in range(0,len(self.xl)):
			y = b+m*self.xl[i]
			chi += (self.yl[i]-y)**2./self.el[i]**2.
		return chi	


def luptm(nmag,bnd):
	#calculates SDSS magnitudes from fluxes
	b = []
	b.append(1.4*10**-10.)
	b.append(.9*10**-10.)
	b.append(1.2*10**-10.)
	b.append(1.8*10**-10.)
	b.append(7.4*10**-10.)
	return -2.5/log(10.)*(asinh((nmag/10.**9.)/(2.*b[bnd]))+log(b[bnd]))

def calcweights(sample,NS,version,zmin=.8,zmax=2.2,app='.fits'):
	#note, zmin/zmax assume LRG sample these need to be change for QSO files
	outfile = outdir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat'+app
	b,m = findlinmb(sample,NS,version,'IMAGE_DEPTH_EXT1',.8,2.2,wm='nosys')
	be,me = findlinmb(sample,NS,version,'EB_MINUS_V-1',.8,2.2,wm='wgdepthext')
	print 'fit parameters for g-band depth'
	print b,m
	print 'fit parameters for EB_MINUS_V (after weighting for depth)'
	print be,me
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat'+app) #read galaxy/quasar file
	#no = 0
	#fo = open(os.environ['HOME']+'/mkEsample_work/depthextweights'+sample+'_'+NS+version+'_mz'+str(zmin)+'xz'+str(zmax)+'.dat','w')
# 	for i in range(0,len(f)):
# 		sysw = f[i]['IMAGE_DEPTH'][1]
# 		sysw = luptm(sysw,1)-3.303*f[i]['EB_MINUS_V']
# 		wd = (1./(b+m*sysw))
# 		ext = f[i]['EB_MINUS_V']
# 		we = (1./(be+me*ext))
# 		wtot = wd*we
# 		#fo.write(str(ws)+'\n')
# 		f[i]['WEIGHT_SYSTOT'] = wtot
# 	#fo.close()
# 	pyfits.writeto(outfile,f,clobber=True)
# 	print 'Wrote file with systematic weights to: ' + outfile
	#print no
	return True

# if __name__ == '__main__':
# 	zmin = .6
# 	zmax = 1.
# 	sample = 'lrg'
# 	NS = 'N'
# 	version = 'v0.8_IRc'
# 	ngvsys(sample,NS,version,zmin,zmax)
# 	calcstarweight(sample,NS,version,zmin=.6,zmax=1.)
