import fitsio
import numpy as np
from healpix import thphi2radec,radec2thphi#,healpix,ang2pix_ring,pix2ang_ring# #comes from AJR's healpix routines
try:
    import healpy as hp
    hpm = True
except:
    print 'no healpy, this will cause problems '
    hpm = False
from math import *

inputdir = '/Users/ashleyross/DESY1/'
samp = '5bins_hidens_hilum_higherlum_jointmask_0.15-0.9'
zlim = [0.15,0.3,0.45,0.6,0.75,0.9]

def maskd(res):
    #degrade mask
	f = fitsio.read(inputdir+samp+'_mask.fits.gz')
	mo = []
	for i in range(0,12*res*res):
		mo.append(0)
	frac = (res/4096.)**2.
	for i in range(0,len(f)):
		p = f[i]['HPIX']
		fr = f[i]['FRACGOOD']
		if fr > 0:
			th,phi = hp.pix2ang(4096,p)
			a = frac*fr
			po = hp.ang2pix(res,th,phi)
			mo[po] += a
	fo = open(inputdir+'Y1RMmask'+str(res)+'ring.dat','w')
	for i in range(0,len(mo)):
		if mo[i] > 0:
			fo.write(str(i)+' '+str(mo[i])+'\n')
	fo.close()
	return True


def mkRMmap(zmin,zmax,res=4096,pixmin=0,pixmax=False,wm='',wo=''):
	if pixmax == False:
		npix = 12*res*res- pixmin
	else:
		npix = 1+pixmax-pixmin	
	pixl = np.zeros((npix))
	f = fitsio.read(inputdir+samp+'_sample.fits.gz')
	if wo != '':
		fo = open('Y1redmagic'+str(zmin)+str(zmax)+wm+'.dat','w')
	if wm != '':		
		sysb = wm.split('_')
		if sysb[0] == 'seerw':
			seemap = mkmap('r','FWHM_MEAN',pixmin=pixmin,pixmax=pixmax)
			sys = sysb[1]
			band = sysb[2]
			b,m = findlinmb(band,sys,zmin,zmax)
			if sys == 'maglimit3':
				tp = '_'
				pr = 1.
			if sys == 'EXPTIME':
				tp = '_total'
				pr = 0.5
			sysmap = mkmap(band,sys,tp,pixmin=pixmin,pixmax=pixmax,pr=pr)
			
		else:
			sys = sysb[0]
			band = sysb[1]
			b,m = findlinmb(band,sys,zmin,zmax)
			if sys == 'maglimit3':
				tp = '_'
				pr = 1.
			if sys == 'EXPTIME':
				tp = '_total'
				pr = 0.5
			sysmap = mkmap(band,sys,tp,pixmin=pixmin,pixmax=pixmax,pr=pr)
	A = 0.505
	bs = 5.2
	sig = .7
	for i in range(0,len(f)):
		z = f[i]['ZREDMAGIC']
		w = 1.
		if z > zmin and z <= zmax:
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi)-pixmin
			
			if wm != '':
				if sysb[0] == 'seerw':
					see = seemap[p]
					sys = sysmap[p]
					rsee = A*(1.-erf((see-bs)/sig))
					rsys = m*sys+b
					w = 1./(rsee*rsys)
				else:
					sys = sysmap[p]
					w = 1./(m*sys+b)
			pixl[p] += w
		if wo != '':
			fo.write(str(f[i]['RA'])+' '+str(f[i]['DEC'])+' '+str(z)+' '+str(w)+'\n')
	fo.close()	
	return pixl

def findlinmb(band,sys,zmin,zmax,wm=''):
	#finds linear fit parameters 
	from optimize import fmin
	d = np.loadtxt(inputdir+'Y1RM42x3pt'+str(zmin)+str(zmax)+'v'+sys+band+wm+'jackerr.dat').transpose()
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
	


def mkmap(band,sys,type='coaddweights3_mean',res=4096,pixmin=0,pixmax=-999,pr=1.):
	if sys == 'nstar':
		starl = hp.read_map(inputdir+'y1a1_gold_1.0.2_stars_nside0512.fits')
		return starl
	if sys == 'ZMAX_hilum' or sys == 'ZMAX_hidens' or sys == 'ZMAX_higherlum':
		f = fitsio.read(inputdir+samp+'_mask.fits.gz')
		starl = []
		if pixmax == -999:
			npix = 12*res*res-pixmin
		else:
			npix = 1+pixmax-pixmin	
		for i in range(0,npix):
			starl.append(0)
		for i in range(0,len(f)):
			p = f[i]['HPIX']-pixmin
			if p >= 0 and p < npix:
				starl[p] = f[i][sys]**pr
		return starl
	f = fitsio.read(inputdir+'Y1A1NEW_COADD_SPT_band_'+band+'/Y1A1NEW_COADD_SPT_band_'+band+'_nside4096_oversamp4_'+sys+'_'+type+'.fits.gz')
	starl = []
	if pixmax == -999:
		npix = 12*res*res
	else:
		npix = 1+pixmax-pixmin	
	for i in range(0,npix):
		starl.append(0)
	for i in range(0,len(f)):
		p = f[i]['PIXEL']-pixmin
		if p >= 0 and p < npix:
			starl[p] = f[i]['SIGNAL']**pr
	f = fitsio.read(inputdir+'Y1A1NEW_COADD_STRIPE82/nside4096_oversamp4/Y1A1NEW_COADD_STRIPE82_band_'+band+'_nside4096_oversamp4_'+sys+'_'+type+'.fits.gz')
	for i in range(0,len(f)):
		p = f[i]['PIXEL']-pixmin
		if p >= 0 and p < npix:
			starl[p] = f[i]['SIGNAL']**pr

	return starl

def ngalvnstarsysl(gl,sl,mask,res=4096,t=.2,nbin=10,smin=50,smax=600,pixmin=0):
	from healpix import thphi2radec,pix2ang_ring
	sysl = []
	for i in range(0,nbin*2):
		sysl.append([])
	minv = 1000
	maxv = -100
	print smin,smax
	nb = 0
# 	for i in range(0,len(sl)):
# 		m = ml[i]
# 		if m > 0:
# 			if m > t:
# 				if sl[i] > maxv and sl[i] < smax:
# 					maxv = sl[i]
# 				if sl[i] < minv and sl[i] > smin:
# 					minv = sl[i]
	for i in range(0,len(mask)):
		if res == 4096:
			p = mask[i]['HPIX']-pixmin
			fr = mask[i]['FRACGOOD']
		else:
			p = int(mask[i][0])
			fr = mask[i][1]
		if fr > t:
			slv = sl[p]
			if slv > maxv and slv < smax:
				maxv = sl[p]
			if slv < minv and slv > smin:
				minv = slv
			
	print 'min,max,nb'
	print minv,maxv,nb
	gtl = []
	pl = []
	for i in range(0,nbin):
		gtl.append(0)
		pl.append(0)
	gt = 0
	pt = 0
	slt = 0
	print len(gl)
# 	for i in range(0,len(gl)):
# 		m = ml[i]
# 		if m > 0:
# 			if m > t:
# 				slv = sl[i]
	for i in range(0,len(mask)):
		if res == 4096:
			p = mask[i]['HPIX']-pixmin
			fr = mask[i]['FRACGOOD']
		else:
			p = mask[i][0]
			fr = mask[i][1]
		if fr > t:
			slv = sl[p]

			if slv >= minv and slv <= maxv:
				pt += fr
				gt += gl[p]

				slt += sl[p]
				bin = int((slv-minv)/(maxv*1.00001-minv)*nbin)
				if bin < nbin:
					sysl[bin*2].append(gl[p])
					sysl[bin*2+1].append(fr)
				gtl[bin] += slv*fr
				pl[bin] += fr
	bcl = []
	for i in range(0,nbin):
		bcl.append(gtl[i]/pl[i])
	print gt,pt,t,slt
	print slt/pt
	return sysl,gt/pt,minv,maxv,bcl
    
def putngalvsysjack(band,sys,type='coaddweights3_mean',zmin=0,zmax=2,smin=-999,smax=999,njack=20,res=4096,t=.2,nbin=10,wm='',pr=1.):
	if res == 4096:
		mask = fitsio.read(inputdir+samp+'_mask.fits.gz')
		pixmin = min(mask['HPIX'])
		pixmax = max(mask['HPIX'])
		npix = 1+pixmax-pixmin
		ml = np.zeros((npix))
		npixv = 0
		for i in range(0,len(mask)):
			p = mask[i]['HPIX']-pixmin
				#if res == 4096:
				#fr = 1.
				#else:
			fr = mask[i]['FRACGOOD']
			ml[p] = fr
			npixv += 1
	else:
		pixmin = 0
		pixmax = 12*res*res-1
		npix = 1+pixmax-pixmin
		mask = np.loadtxt(inputdir+'Y1RMmask'+str(res)+'ring.dat')

	gall = mkRMmap(zmin,zmax,res,pixmin=pixmin,pixmax=pixmax,wm=wm)
	sysmap = mkmap(band,sys,type,pixmin=pixmin,pixmax=pixmax,pr=pr)
	if smin == -999:
		smin = min(sysmap)
	if smax == 999:
		smax = max(sysmap)
	print smin,smax	
	sysl,mt,smin,smax,bcl = ngalvnstarsysl(gall,sysmap,mask,res=res,t=t,nbin=nbin,smin=smin,smax=smax,pixmin=pixmin)
	print mt
	sl = []
	ml = []
	print len(sysl)
	for i in range(0,2*nbin,2):
		print len(sysl[i])
	for i in range(0,2*nbin,2):
		std = 0
		ng = 0
		npix = 0
		print len(sysl[i]),len(sysl[i+1])
		for j in range(0,len(sysl[i])):
			ng += sysl[i][j]
			npix += sysl[i+1][j]
		mean = ng/npix/mt
		print i,mean
		ml.append(mean)
		if len(sysl[i]) < njack:
			for k in range(0,len(sysl[i])):
				std += (sysl[i][k]/sysl[i+1][k]/mt-mean)**2.
			std = sqrt(std/(len(sysl[i])-1.))/sqrt(len(sysl[i])-1.)   # I think it should be outside the loop?
		else:
			jkf = len(sysl[i])/njack
			for k in range(0,njack):
				ng = 0
				npix = 0
				minj = jkf*k
				maxj = jkf*(k+1)
				for j in range(0,len(sysl[i])):
					if j < minj or j >= maxj:
						ng += sysl[i][j]
						npix += sysl[i+1][j]
				mj = ng/npix/mt
				std += (mj-mean)**2.
			std = sqrt((njack-1.)/float(njack)*std)
		print i,mean,std
		sl.append(std)
	rw = ''
	fo = open(inputdir+'Y1RM42x3pt'+str(zmin)+str(zmax)+rw+'v'+sys+band+wm+'jackerr.dat','w')
	for i in range(0,nbin):
		slb = sl[i]
		st = i*(smax-smin)/float(nbin)+(smax-smin)/nbin/2.+smin
		#fo.write(str(st)+' '+str(ml[i])+' '+str(slb)+'\n')
		fo.write(str(bcl[i])+' '+str(ml[i])+' '+str(slb)+'\n')
	fo.close()
	return True
    
def plotsys(sys,band,zmin,zmax,wm='',pr=1.):
	from matplotlib import pyplot as plt
	import pylab
	plt.clf()
	pylab.ion()
	d = np.loadtxt(inputdir+'Y1RM42x3pt'+str(zmin)+str(zmax)+'v'+sys+band+wm+'jackerr.dat').transpose()
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(inputdir+'RMsysplots/nY1RM42x3pt'+str(zmin)+str(zmax)+'vsys'+sys+band+wm+'.pdf')
	xl = d[0]**pr
	null = sum((d[1]-1.)**2./d[2]**2.)
	ol = np.ones((len(d[0])))
	plt.plot(xl,ol,'k--')
	plt.errorbar(xl,d[1],d[2],fmt='ko')
	x = max(xl)-.3*(max(xl)-min(xl))
	y = 1.15
	plt.text(x,y,r'$\chi^2$/dof = '+str(null)[:4]+'/'+str(len(ol)))
	plt.ylim(.8,1.2)
	if pr != 1.:
		sys = '('+sys+')^'+str(pr)
	plt.xlabel(band+'-band '+sys)
	plt.ylabel('number density/average number density')
	plt.title('Redmagic galaxies for 3x2pt, '+str(zmin)+'< z <'+str(zmax))
	plt.show()
	pp.savefig()
	pp.close()
	return True

def testseeplot(sys,band,zmin,zmax,wm='',A=.505,sig=.7,bs=5.2):
	from matplotlib import pyplot as plt
	import pylab
	plt.clf()
	pylab.ion()
	d = np.loadtxt(inputdir+'Y1RM42x3pt'+str(zmin)+str(zmax)+'v'+sys+band+wm+'jackerr.dat').transpose()
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(inputdir+'RMsysplots/nY1RM42x3pt'+str(zmin)+str(zmax)+'vsys'+sys+band+wm+'.pdf')
	xl = d[0]
	null = sum((d[1]-1.)**2./d[2]**2.)
	ol = np.ones((len(d[0])))
	tl = []
	for i in range(0,len(d[0])):
		tl.append(A*(1.-erf((d[0][i]-bs)/sig)))
	print tl[0]
	expchi = sum((d[1]-np.array(tl))**2./d[2]**2.)
	print expchi
	plt.plot(xl,ol,'k--')
	plt.plot(xl,tl,'k-')
	plt.errorbar(xl,d[1],d[2],fmt='ko')
	x = max(xl)-.2*(max(xl)-min(xl))
	y = 1.15
	plt.text(x,y,r'$\chi^2$/dof = '+str(null)[:4])
	plt.ylim(.6,1.2)
	#if pr != 1.:
	#	sys = '('+sys+')^'+str(pr)
	plt.xlabel(band+'-band '+sys)
	plt.ylabel('number density/average number density')
	plt.title('Redmagic galaxies for 3x2pt, '+str(zmin)+'< z <'+str(zmax))
	plt.show()
	pp.savefig()
	pp.close()
	return True
