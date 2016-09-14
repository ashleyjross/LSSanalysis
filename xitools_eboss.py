dirsci = '/mnt/lustre/ashleyr/eboss/' #where AJR puts eboss catalogs, change this to wherever you have put catalogs
dirsys = '/Users/ashleyross/ngalvsys_wiki/' #change to local directory where ngalvsys from wiki was put, note star map and depth map included
dirfits = '/Users/ashleyross/fitsfiles/' #change to where your catalog files are
ebossdir = '/Users/ashleyross/eboss/' #where AJR puts correlation functions, writes out results

import fitsio #needed to read data
import numpy as np
from math import *



def P2(mu):
	return .5*(3.*mu**2.-1.)
	
def P4(mu):
	return .125*(35.*mu**4.-30.*mu**2.+3.)

def P6(mu):
	return 1./16.*(231.*mu**6.-315.*mu**4.+105.*mu**2.-5.)

def P8(mu):
	return 1./128.*(6435.*mu**8.-12012.*mu**6.+6930.*mu**4.-1260.*mu**2.+35.)


def mksubfile_Dmufbfjack(file1,file2,jack,off=0,a=''):
	fo = open('sub'+str(jack+off)+a+'.sh','w')
	fo.write('#!/bin/bash\n')
	fo.write('#$ -V -cwd\n')
	for i in range(0,20):
		fo.write('~/pp2pt_Dmufb '+file1+str(i) +' '+file2+str(jack) +' -o output.'+file1+' \n')
	fo.close()
	return True

def mksuball_nran_Dmufbfjack(ranf,galf,nran,wr,njack=20):
	fo = open('suball.sh','w')
	fo.write('#!/bin/bash\n')
	for i in range(0,njack):
		mksubfile_Dmufbfjack(galf,galf,i)
		fo.write('qsub sub'+str(i)+'.sh \n')
		fo.write('sleep 1 \n')
	for j in range(0,nran):	
		off = 20+40*j
		for i in range(0,njack):
			mksubfile_Dmufbfjack(ranf+str(j)+wr,galf,i,off)
			fo.write('qsub sub'+str(i+off)+'.sh \n')
			fo.write('sleep 1 \n')
		off = 40+40*j
		for i in range(0,njack):
			mksubfile_Dmufbfjack(ranf+str(j)+wr,ranf+str(j)+wr,i,off)
			fo.write('qsub sub'+str(i+off)+'.sh \n')
			fo.write('sleep 1 \n')

	fo.close()
	return True

def mkran1mil(sample,NS,version,N=0,c='sci',app='.fits',compmin=0):
	dirout = ''
	dir = ''
	if c == 'sci':
		dir = dirsci
		dirout = dir
	minc = N*10**6
	maxc = (N+1)*10**6 #will become relevant once files are big enough
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.ran'+app)#[minc:maxc]
	if maxc > len(f):
		maxc = len(f)
	#f = fitsio.FITS(dir+'eboss_'+version+'-'+sample+'-'+NS+'eboss_'+version+'.ran'+app)[1]
	no = 0
	fo = open(dirout+'reboss'+sample+'_'+NS+version+'_'+str(N)+'.dat','w')
	#for i in range(0,len(f)):		
	for i in range(minc,maxc):
		comp = f[i]['COMP']
		if comp > compmin:
			fo.write(str(f[i]['RA'])+' '+str(f[i]['DEC'])+'\n') #rest of the information comes from matching to galaxy catalog, this gives set of randoms matching angular footprint
		#fo.write(str(float(f[i]['RA']))+' '+str(float(f[i]['DEC']))+'\n') #to be used if opened and read from fitsio.FITS
	fo.close()
	return True

def mkjackf(sample,NS,version,Njack=20):
	#defines jack-knifes
	mf = open('ranHeal_pix256eboss'+sample+'_'+NS+version+'.dat').readlines()
	fo = open('jackhpixeboss'+sample+'_'+NS+version+str(Njack)+'.dat','w')
	for i in range(0,Njack-1):
		fo.write(str(mf[(len(mf)/Njack)*(i+1)].split()[0])+'\n')
	fo.close()
	return True

def ranHealp(sample,NS,version,res=256,rad=''):
	#pixelizes random file to create jack-knifes
	import gzip
	dir = dirsci
	angm = 1.
	if rad == 'rad':
		angm = 180./pi
	from healpix import healpix,radec2thphi,thphi2radec
	pixl = []
	h = healpix()
	np = 12*res**2
	for i in range(0,np):
		pixl.append(0)
	f = open(dirsci+'reboss'+sample+'_'+NS+version+'_0'+'.dat')	
	for line in f:
		if line[0] != '#' and line[1] != '#' and line.split()[0] != 'az(d)':
			ln = line.split()
			ra,dec = float(ln[0])*angm,float(ln[1])*angm
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			pixl[p] += 1.
	fo = open('ranHeal_pix'+str(res)+'eboss'+sample+'_'+NS+version+'.dat','w')
	for i in range(0,len(pixl)):
		if pixl[i] > 0:
			fo.write(str(i)+' '+str(pixl[i])+'\n')
	fo.close()
	return True


def mkgal4xi(sample,NS,version,zmin=.6,zmax=1.,c='sci',app='.fits',wm='',compmin=0,gmin=0,gmax=30):
	#note, zmin/zmax assume LRG sample these need to be change for QSO files
	from healpix import healpix, radec2thphi
	from optimize import fmin
	if wm == 'wstar' or wm == 'cpstar':
		wsys = np.loadtxt('allstars17.519.9Healpixall256.dat')
		b,m = findlinmb(sample,NS,version,'star',zmin,zmax)
		h = healpix()
	if wm == 'wdepth' or wm == 'cpdepth':
		wsys = np.loadtxt('healdepthinm512.dat').transpose() #note, file is zipped in directory	
		#b,m = findlinmb(sample,NS,version,'depth',zmin,zmax)
		ds = np.loadtxt('ngebossQSO_S'+version+'_mz0.9xz2.2gm'+str(gmin)+str(gmax)+'512vdepth.dat').transpose()
		dn = np.loadtxt('ngebossQSO_N'+version+'_mz0.9xz2.2gm'+str(gmin)+str(gmax)+'512vdepth.dat').transpose()
		dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
		e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
		lf = linfit(ds[0],dt,e)
		inl = np.array([1.,0])
		b,m = fmin(lf.chilin,inl)

		h = healpix()
	if c == 'sci': #AJR uses this define directory for machine he uses
		dir = dirsci
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat'+app) #read galaxy/quasar file
	no = 0
	gw = ''
	if gmax != 30:
		gw = 'gx'+str(gmax)
	fo = open(dir+'geboss'+sample+'_'+NS+version+'_mz'+str(zmin)+'xz'+str(zmax)+gw+wm+'4xi.dat','w')
	for i in range(0,len(f)):
		z = f[i]['Z']
		comp = f[i]['COMP']
		gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
		if z > zmin and z < zmax and comp > compmin and gm < gmax:
			no += 1
			#w = 1.
			#if wm == '':
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT'] #standard weight to use 
			if wm == 'nfkp':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.) #do this if you want to see difference FKP weights make
			if wm == 'cp':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
			if wm == 'st':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_STAR']*f[i]['WEIGHT_FKP']
			if wm == 'see':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']*f[i]['WEIGHT_FKP']
			if wm == 'seenfkp':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']
			if wm == 'stsee':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']*f[i]['WEIGHT_STAR']*f[i]['WEIGHT_FKP']
			if wm == 'cpstar' or wm == 'cpdepth':
				w = (f[i]['WEIGHT_CP'])*f[i]['WEIGHT_FKP']
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			if wm == 'wstar' or wm == 'cpstar':
				pix2 = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws
			if wm == 'wdepth' or wm == 'cpdepth':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				

			fo.write(str(f[i]['RA'])+' '+str(f[i]['DEC'])+' '+str(z)+' '+str(w)+'\n')
	fo.close()
	print no
	return True

def mkran4xi(sample,NS,version,N=0,wm='',zmin=.6,zmax=.1,comp = 'sci',gmax=30):
	from random import random
	if comp == 'sci':
		dir = dirsci 
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gw = ''
	if gmax != 30:
		gw = 'gx'+str(gmax)

	gf = np.loadtxt(dir+'geboss'+sample+'_'+NS+version+'_'+wz+gw+wm+'4xi.dat').transpose()
	fr = np.loadtxt(dir+'reboss'+sample+'_'+NS+version+'_'+str(N)+'.dat').transpose()
	fo = open(dir+'reboss'+sample+'_'+NS+version+'_'+str(N)+wz+gw+wm+'4xi.dat','w')
	n = 0
	for i in range(0,len(fr[0])):
		indz = int(random()*len(gf[0]))
		fo.write(str(fr[0][i])+' '+str(fr[1][i])+' '+str(gf[2][indz])+' '+str(gf[3][indz])+'\n')
		#redshifts, and galaxy weight at that redshift, are randomly sampled and added to random file
		n += 1.
	print n #just helps to know things worked properly
	fo.close()
	return True

def mkran4xifit(sample,NS,version,N=0,wm='',zmin=.6,zmax=.1,comp = 'sci'):
	from random import random
	if comp == 'sci':
		dir = dirsci 
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.ran.fits')
	fo = open(dir+'reboss'+sample+'_'+NS+version+'_'+str(N)+wz+wm+'4xi.dat','w')
	n = 0
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z > zmin and z < zmax:
			fo.write(str(f[i]['RA'])+' '+str(f[i]['DEC'])+' '+str(f[i]['Z'])+' '+str(f[i]['WEIGHT_FKP'])+'\n')
			n += 1.
	print n #just helps to know things worked properly
	fo.close()
	return True


def createSourcesrd_ad(file,md='deg'):
	#takes the 4xi files and put them into the format AJR's pair-counting code uses
	dir = dirsci
	from healpix import healpix,radec2thphi
	from random import random
	from Cosmo import distance
	d = distance(.31,.69) #cosmology assumed in final BOSS analyses, make sure this conforms to current
	h = healpix()
	fo = open(file+'pcadw.dat','w')
	f = open(dir+file+'4xi.dat')
	count = 0
	for line in f:
		if line[0] != '#':
			cols = line.split()
			ra = float(cols[0])         # The longitude coordinate
			dec = float(cols[1])       # The latitude coordinate
			if md == 'deg':
				th,phi = radec2thphi(ra,dec)
			if md == 'rad':
				th,phi = radec2thphi(ra*180./pi,dec*180./pi)
			pix = int(h.ang2pix_nest(256,th,phi))
			z = float(cols[2])
			cd = d.dc(z)
			w = float(cols[3])
			if md == 'deg':
				sra = sin(radians(ra))
				cra = cos(radians(ra))
				sdec = sin(radians(dec))
				cdec = cos(radians(dec))
			if md == 'rad':
				sra = sin(ra)
				cra = cos(ra)
				sdec = sin(dec)
				cdec = cos(dec)
			
			fo.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(cd)+' '+str(w)+' '+str(pix)+'\n')
			count += 1.*w

	fo.close()
	return True

def createSourcesrd_adJack(file,jack,NS='N2',Njack=20):
	#divides the file created by previous module into 20 jack-knife samples
	from healpix import healpix,radec2thphi
	from random import random
	import gzip
	fo = open(file+str(jack)+'pcadw.dat','w')
	f = open(file+'pcadw.dat')
	fjack = open('jackhpix'+NS+str(Njack)+'.dat').readlines()
	count = 0
	jackmin = int(fjack[jack-1])
	if jack == Njack-1:
		jackmax = 10000000000
	else:
		jackmax = int(fjack[jack])
	if jack == 0:
		jackmin= 0
	h = healpix()
	for line in f:
		if line[0] != '#':
			cols = line.split()
			pix = int(cols[-1])
			if pix > jackmin and pix <= jackmax:
				#if zm == '':
				fo.write(line)

	fo.close()
	return True

def createalladfilesfb(sample,NS,version,nran=1,wm='',zmin=.6,zmax=1.,gmax=30):
	#after defining jack-knifes, this makes all of the divided files and the job submission scripts
	#./suball.sh sends all of the jobs to the queue on the system I use
	mkgal4xi(sample,NS,version,wm=wm,zmin=zmin,zmax=zmax,gmax=gmax)
	gw = ''
	if gmax != 30:
		gw = 'gx'+str(gmax)

	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gf = 'geboss'+sample+'_'+NS+version+'_'+wz+gw+wm
	createSourcesrd_ad(gf)
	for i in range(0,20):
		createSourcesrd_adJack(gf,i,'eboss'+sample+'_'+NS+version)

	rf = 'reboss'+sample+'_'+NS+version+'_'
	for rann in range(0,nran):	
		print rann
		mkran4xi(sample,NS,version,N=rann,wm=wm,zmin=zmin,zmax=zmax,gmax=gmax)
		#mkran4xifit(sample,NS,version,N=rann,wm=wm,zmin=zmin,zmax=zmax)
		rfi = rf+str(rann)+wz+gw+wm
		createSourcesrd_ad(rfi)
		for i in range(0,20):
			createSourcesrd_adJack(rfi,i,'eboss'+sample+'_'+NS+version)
	mksuball_nran_Dmufbfjack(rf,gf,nran,wr=wz+gw+wm)
	return True

def ppxilcalc_LSDfjack_bs(sample,NS,version,jack,mom,zmin=.6,zmax=1.,wm='',bs=5,start=0,rmax=250,mumin=0,mumax=1.,nranf=1,njack=20,wf='n',wmu = ''):
	#finds xi, for no jack-knife, set jack = -1, otherwise can be used to calculate jack-knife xi for 0 <= jack < Njack
	from numpy import zeros
	from time import time
	#DDnl = zeros((nbin,njack),'f')
	rf = 'reboss'+sample+'_'+NS+version+'_'
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gf = 'geboss'+sample+'_'+NS+version+'_'+wz+wm

	DDnl = []	
	DDnorml = 0
	DDnormt = 0
	DRnl = []
	DRnorml = 0
	DRnormt = 0
	RRnl = []
	RRnl0 = []
	nmubin = 100
	nbin = rmax/bs
	for i in range(0,rmax*nmubin):
		DDnl.append(0)
		DRnl.append(0)
		RRnl.append(0)
		RRnl0.append(0)
	RRnorml = 0
	RRnormt = 0
	pl = []
#	pl.append((0,0,0,0))
	if wf == 'y':
		fD = open('DDcounts'+file+'.dat','w')
		fDR = open('DRcounts'+file+'.dat','w')
		fR = open('RRcounts'+file+'.dat','w')
	nmut = 0
	for i in range(0,nmubin):
		mu = i/float(nmubin)+.5/float(nmubin)
		mub = int(mu*nmubin)
		#print mu
		if mu > mumin and mu < mumax:
			pl.append((1.,P2(mu),P4(mu),P6(mu),P8(mu),mub))
			nmut += 1.
		else:
			pl.append((0,0,0,0,0,0))	
	#print len(pl)
	for i in range(0,njack):
		#print i
		for j in range(0,njack):
			if jack != i and jack != j:
				fdp = open(gf+str(j)+gf+str(i)+'2ptdmu.dat').readlines()
				DDnormt += float(fdp[0])
				fdnp = open(rf+'0'+wz+wm+str(j)+gf+str(i)+'2ptdmu.dat').readlines()
				fr = open(rf+'0'+wz+wm+str(j)+rf+'0'+wz+wm+str(i)+'2ptdmu.dat').readlines()
				DRnormt += float(fdnp[0])
				RRnormt += float(fr[0])
				for k in range(1,len(fdp)):
					dp = float(fdp[k])
					dr = float(fdnp[k])
					rp = float(fr[k])
					DDnl[k-1] += dp
					DRnl[k-1] += dr
					RRnl[k-1] += rp
					
					
	for nr in range(1,nranf):
		for i in range(0,njack):
			#print i
			for j in range(0,njack):
				if jack != i and jack != j:
					fdnp = open(rf+str(nr)+wz+wm+str(j)+gf+str(i)+'2ptdmu.dat').readlines()
					fr = open(rf+str(nr)+wz+wm+str(j)+rf+str(nr)+wz+wm+str(i)+'2ptdmu.dat').readlines()
					DRnormt += float(fdnp[0])
					RRnormt += float(fr[0])
					for k in range(1,len(fdp)):						
						dr = float(fdnp[k])
						rp = float(fr[k])
						DRnl[k-1] += dr
						RRnl[k-1] += rp
		
	#print RRnl
	
	if wf == 'y':
		fD.write(str(DDnormt)+'\n')
		fDR.write(str(DRnormt)+'\n')
		fR.write(str(RRnormt)+'\n')
		for j in range(0,100):
			for i in range(0,rmax):
				fD.write(str(DDnl[j+100*i])+' ')
				fDR.write(str(DRnl[j+100*i])+' ')
				fR.write(str(RRnl[j+100*i])+' ')
			fD.write('\n')
			fDR.write('\n')
			fR.write('\n')
		fD.close()
		fDR.close()
		fR.close()
	xil = zeros((nbin),'f')
	for i in range(start,rmax,bs):
		xi = 0
		dd = 0
		dr = 0
		rr = 0
		for j in range(0,nmubin):
			if wmu != 'counts':
				dd = 0
				dr = 0
				rr = 0
			for k in range(0,bs):
				bin = nmubin*(i+k)+j			
				if bin < len(RRnl):
					if RRnl[bin] == 0:
						pass
				
					else:
						dd += DDnl[bin]
						rr += RRnl[bin]
						dr += DRnl[bin]
			#if rr != 0 and wm == 'muw':			
			if wmu != 'counts':
				xi += pl[j][mom]/float(nmut)*(dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr
		if wmu == 'counts':
			xi = (dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr		
		if i/bs < nbin:
			xil[i/bs] = xi
	return xil


def ppxilfile_bs(sample,NS,version,mom,zmin=.6,zmax=1.,wm='',bs=5,start=0,rmax=250,nranf=1,njack=20,mumin=0,mumax=1.,wmu=''):
	#write out xi to a file, no jack-knife errors
	#for quadrupole, set mom = 1, for hexadecapole set mom = 2, ect. (odd moments not supported)
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gf = 'geboss'+sample+'_'+NS+version+'_'+wz+wm
	ans = ppxilcalc_LSDfjack_bs(sample,NS,version,-1,mom,zmin=zmin,zmax=zmax,wm=wm,mumin=mumin,mumax=mumax,bs=bs,start=start,rmax=rmax,nranf=nranf,wmu=wmu)
	print ans
	if mumin != 0:
		gf += 'mum'+str(mumin)
	if mumax != 1.:
		gf += 'mux'+str(mumax)
	fo = open('xi'+str(2*mom)+gf+wmu+str(bs)+'st'+str(start)+'.dat','w')
	print gf
	for i in range(0,rmax/bs):
		r = float(bs)/2.+float(bs)*i+float(start)
		fo.write(str(r)+' '+str(ans[i]*(4.*mom+1.))+'\n')
	fo.close()
	return True

def mkODmap(sample='lrg',NS='N',version='v1.0_IRt',res=256,zmin=.6,zmax=1.):
	#make an 2D over-density map to be used for w(theta) calculation
	from healpix import healpix,radec2thphi
	h = healpix()
	ff = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.ran.fits')
	fdf = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat.fits')
	nran = len(ff)
	np = 12*res*res
	ml = [] #pixel list for randoms/mask
	gl = [] #pixel list for galaxies
	for i in range(0,np):
		ml.append(0)
		gl.append(0)
	for i in range(0,len(ff)):
		ra,dec = ff[i]['RA'],ff[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = h.ang2pix_nest(res,th,phi)
		ml[p] += 1.
	for i in range(0,len(fdf)):
		z = fdf[i]['Z']
		if z > zmin and z < zmax:
			ra,dec = fdf[i]['RA'],fdf[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			p = h.ang2pix_nest(res,th,phi)
			w = (fdf[i]['WEIGHT_NOZ']+fdf[i]['WEIGHT_CP']-1.)*fdf[i]['WEIGHT_FKP']*fdf[i]['WEIGHT_SYSTOT'] #standard weight to use 
			#multiply by your own weight here
			gl[p] += w
	ave = sum(gl)/sum(ml)
	print ave
	fo = open('geboss'+sample+'_'+NS+version+'_mz'+str(zmin)+'xz'+str(zmax)+str(res)+'odenspczw.dat','w') #this is the file for the code
	ft = open('geboss'+sample+'_'+NS+version+'_mz'+str(zmin)+'xz'+str(zmax)+str(res)+'rdodens.dat','w') #this is with theta,phi coordinates, in case you want to plot it
	no = 0		
	for i in range(0,len(ml)):
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(res,i)
			sra = sin(phi)
			cra = cos(phi)
			sdec = sin(-1.*(th-pi/2.))
			cdec = cos(-1.*(th-pi/2.))
			od = gl[i]/(ave*ml[i]) -1.
			fo.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(od)+' '+str(ml[i])+'\n')
			ft.write(str(th)+' '+str(phi)+' '+str(od)+' '+str(ml[i])+'\n')
	ft.close()
	fo.close()
	
###Routines below are for systematic analysis, etc.

def ngvsys(sampl,NS,ver,sys,sysmin,sysmax,res,zmin,zmax,wm='',gmag=False):
	#sample is the sample being used, e.g. 'lrg'
	#NS is either 'N' or 'S'
	#ver is the version, e.g., 'v1.0_IRt'
	#sys is a string containing the name of the systematic to be tested
	#sysmin is the minimum value of the systematic to be tested, ~25 is a good value to use for stars
	#sysmax is the maximum value of the systematic to be tested, ~200 is a good value to use for stars
	#res is Nside for the healpix map, default should be 512 for sky,seeing,air and 256 for extinction and stars
	#zmin is the minimum redshift to use, 0.6 is minimum used for lrgs, 0.9 for qsos
	#zmax is the maximum redshift to used, 1.0 for lrgs and 2.2 for qsos
	from healpix import healpix, radec2thphi
	h = healpix()
	
	stl = []
	wstl = []
	errl = []
	if wm == 'wstar':
		wsys = np.loadtxt(dirsys+'allstars17.519.9Healpixall256.dat')
		b,m = findlinmb(sampl,NS,ver,'star',zmin,zmax)
	if wm == 'wdepth':
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()	
		b,m = findlinmb(sampl,NS,ver,'depth',zmin,zmax)
	if sys != 'star' and sys != 'ext' and sys != 'depth':
		fsys = np.loadtxt(dirsys+'heal'+sys+'nm'+str(res)+'.dat')
	if sys == 'ext':
		fsys = np.loadtxt(dirsys+'healSFD_r_'+str(res)+'_fullsky.dat')*2.285/2.751 #r-band extinction
	if sys == 'star':
		fsys = np.loadtxt(dirsys+'allstars17.519.9Healpixall256.dat')
	if sys == 'depth':
		try:
			fsys = np.loadtxt(dirsys+'healdepthinm'+str(res)+'.dat').transpose()
		except:
			#from DJS email Oct 7th to eboss targeting:
			# Flux(5-sigma) = anorm * PSF_FWHM * sqrt(SKYFLUX) * 10^(0.4*kterm*AIRMASS)
			# anorm = 0.387, 0.218, 0.241, 0.297, 0.665
			# kterm =  [0.49, 0.17, 0.10, 0.06, 0.06] 
			# current files are for i band so
			anorm = 0.297
			kterm = 0.06
			fsee = np.loadtxt(dirsys+'healseenm'+str(res)+'.dat')
			fsky = np.loadtxt(dirsys+'healskynm'+str(res)+'.dat')
			fair = np.loadtxt(dirsys+'healairnm'+str(res)+'.dat')
			fsys = []
			febv = np.loadtxt(dirsys+'healSFD_r_256_fullsky.dat')/2.751 #ebv map at Nside 256
			ebv5 = []
			for i in range(0,12*512*512):
				th,phi = h.pix2ang_nest(512,i)
				if fsee[i] > 0 and fsee[i] < 2.5:
					pix2 = h.ang2pix_nest(256,th,phi)
					ebv5.append(febv[pix2])
				else:
					ebv5.append(0)	
			for i in range(0,len(fsee)):
				see = float(fsee[i])
				sky = float(fsky[i])
				air = float(fair[i])
			
				dp = luptm(anorm*see*sqrt(sky)*10.**(.4*kterm*air),3)-ebv5[i]*1.698
				fsys.append(dp)
			fo = open(dirsys+'healdepthinm'+str(res)+'.dat','w')
			for i in range(0,len(fsys)):
				fo.write(str(fsys[i])+'\n')
			fo.close()			

	ml = []
	npo = 12*res**2

	print min(fsys),max(fsys)
	h = healpix()
	pixlg = []
	pixlr = []
	for i in range(0,npo):
		pixlg.append(0)
		pixlr.append(0)

	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.ran.fits') #read galaxy/quasar file
	nr = 0
	for i in range (0,len(f)):
		ra,dec = f[i]['RA'],f[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(res,th,phi))
		pixlr[p] += 1.
		nr += 1.
	print nr

	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.dat.fits') #read galaxy/quasar file
	no = 0
	zm = 0
	nt = 0
	for i in range (0,len(f)):
		z = f[i]['Z']
		gc = True
		if gmag != False:
			gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if gm < gmag[0] or gm > gmag[1]:
				gc = False
		if z > zmin and z < zmax and gc:
			no += 1
			#w = 1.
			#if wm == '':
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP'] #standard weight to use if no systematic weights have been defined
			if wm == 'fid':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT']
			if wm == 'nfkp':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.) #do this if you want to see difference FKP weights make
			if wm == 'cp':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
			if wm == 'st':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_STAR']*f[i]['WEIGHT_FKP']
			if wm == 'see':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']*f[i]['WEIGHT_FKP']
			if wm == 'seenfkp':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']
			if wm == 'stsee':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']*f[i]['WEIGHT_STAR']*f[i]['WEIGHT_FKP']
			if wm == 'wstar':
				pix2 = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws
			if wm == 'wdepth':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws				
			pixlg[p] += 1.*w
			zm += w*z
			nt += w
	print 'total number, weighted number'
	print no,nt
	print 'mean redshift'
	print zm/nt
	binnbs = []
	binns = []
	nsysbin = 10 #10 bins are being used
	for i in range(0,nsysbin):
		binnbs.append(0)
		binns.append(0)

	nbt = 0
	nt = 0
	bs = 0
	bsr = 0
	n0 = 0
	sysm = float(nsysbin)/(sysmax-sysmin)
	ng0 = 0
	for i in range(0,npo):
		sysv = float(fsys[i])
		if sysv != 0: #the maps are not perfect, entries with 0s shouldn't be used
			nt += pixlr[i]
			nbt += pixlg[i]
			bins = int((sysv-sysmin)*sysm)
			if bins >= 0 and bins < nsysbin:
				binnbs[bins] += pixlg[i]
				binns[bins] += pixlr[i]
			else:
				bs += pixlg[i] #count numbers outside of sysmin/sysmax
				bsr += pixlr[i]
		else:
			n0 += pixlr[i] #count numbers inside bad pixels in sys map
			ng0 += pixlg[i]
					
	print 'total number of randoms/objects '+str(nt)+'/'+str(nbt)
	print 'number of randoms/objects where sys = 0 '+str(n0)+'/'+str(ng0)
	print 'number of randoms/objects outside tested range '+str(bsr)+'/'+str(bs)			
	ave = nbt/nt
	print 'average number of objects per random is '+ str(ave)
	if gmag != False:
		wm += 'gm'+str(gmag[0])+str(gmag[1])
	fs = open(ebossdir+'n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+str(res)+'v'+sys+'.dat','w')
	xl = []
	yl = []
	el = []
	for i in range(0,nsysbin):
		sysv = sysmin + 1./(2.*sysm) + i/sysm
		if binns[i] > 0:
			ns = binnbs[i]/binns[i]/ave
			nse = sqrt(binnbs[i]/(binns[i])**2./(ave)**2.+(binnbs[i]/ave)**2./(binns[i])**3.) #calculate poisson error
		else:
			ns = 1. #write out 1.0 1.0 if no pixels at given value of sys
			nse = 1.		
		fs.write(str(sysv)+' '+str(ns)+' '+str(nse)+'\n')
		xl.append(sysv)
		yl.append(ns)
		el.append(nse)

	return xl,yl,el

def findlinmb(sampl,NS,ver,sys,zmin,zmax,wm='',res='512'):
	#finds linear fit parameters (depth or stellar density relationships hard-coded to expect given resolutions)
	from optimize import fmin
	if sys == 'star':
		res = '256'
	if sys == 'depth':
		res = '512'
	d = np.loadtxt('n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+str(res)+'v'+sys+'.dat').transpose()	
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

### Below is for theoretical covariance matrix
###!Do NOT use for actual science!

def covcalcxi0(r1,r2,version,sample,amp,zmin,zmax,dr=10.,sp=1.,nzbs=20,a='',v='y',gam=-1.7,file='Challenge_matterpower',dir='/Users/ashleyross/DR12/',mun=0,beta=0.4,sfog=3.0,amc=0.0,sigt=6.,sigr=10.,mult=1.,sigs=15.):
	#calculates very approximate covariance, terms added in rather ad-hoc way, DO NOT USE for actual science
	from Cosmo import distance
	d = distance(.31,.69)
	#vol = veff
	fk = open(dir+file+'.dat').readlines()
	k0 = float(fk[0].split()[0])
	k1 = float(fk[1].split()[0])
	ldk = log(k1)-log(k0)
	sumxi = 0
	nbN = np.loadtxt(ebossdir+'nbar-eboss_'+version+'-'+sample+'-N-eboss_'+version+'.dat').transpose()
	nbS = np.loadtxt(ebossdir+'nbar-eboss_'+version+'-'+sample+'-S-eboss_'+version+'.dat').transpose()
	veff = 0
	nl = []
	vl = []
	for i in range(0,len(nbN[0])):
		z = nbN[0][i]
		if z > zmin and z < zmax:
			nl.append((nbN[3][i]+nbS[3][i])/2.)
			vl.append(nbN[5][i]+nbS[5][i])
	sigi = sqrt((sigr**2.+2.*sigt**2.)/3.)
		
	for kb in range(0,len(fk)-1):
		sigk = 0
		k = float(fk[kb].split()[0])
		pk1 = float(fk[kb].split()[1])*amp
		sigv21 = sigi**2./4.
		damp1 = exp(-1.*k*k*sigv21)	
		C1 = damp1
		pk1 = C1**2.*pk1
		dk = ldk*k
		isigk = 0
		for i in range(0,len(nl)):
			isigk += vl[i]/(pk1**2.+2.*pk1/nl[i]) #sum inverse variance for each redshift shell
		sigk = dk/(r1*r2)/isigk*(sin(k*r1)*sin(k*r2))
		sumxi += sigk*2. #factor of 4 makes BOSS results ~correct
	sumxi = sumxi/(2.*pi**2.)#/vol
	if r1 == r2: #pure shot-noise piece has been separated out since this depends on the r bin size
		npp = 0
		for i in range(0,len(nl),nzbs):
			vt = 0
			nave = 0
			for j in range(0,nzbs):
				vt += vl[i+j]
				nave += nl[i+j]*vl[i+j]
			nave = nave/vt
			#npp += 2.*pi*r1**2.*dr*vl[i]*nl[i]**2.
			#npp += 2.*pi*r1**2.*dr*vt*nave**2.
			npp += 1./3.*pi*((r1+dr/2.)**3.-(r1-dr/2.)**3.)*vt*nave**2.
		sumxi += 1./npp
		print r1,sqrt(1./npp)*r1
	return sumxi

def covcalcxi0boss(r1,r2,amp,zmin,zmax,dr=10.,sp=1.,a='',v='y',gam=-1.7,file='Challenge_matterpower',dir='/Users/ashleyross/DR12/',mun=0,beta=0.4,sfog=3.0,amc=0.0,sigt=6.,sigr=10.,mult=1.,sigs=15.):
	#same as covcalcxi0, but for BOSS catalogs
	from Cosmo import distance
	d = distance(.31,.69)
	#vol = veff
	fk = open(dir+file+'.dat').readlines()
	k0 = float(fk[0].split()[0])
	k1 = float(fk[1].split()[0])
	ldk = log(k1)-log(k0)
	sumxi = 0
	nbN = np.loadtxt(dir+'nbar-cmasslowz-dr12v4-N-Reid-om0p31.dat').transpose()
	nbS = np.loadtxt(dir+'nbar-cmasslowz-dr12v4-S-Reid-om0p31.dat').transpose()
	veff = 0
	nl = []
	vl = []
	for i in range(0,len(nbN[0])):
		z = nbN[0][i]
		if z > zmin and z < zmax:
			nl.append((nbN[3][i]+nbS[3][i])/2.)
			vl.append(nbN[5][i]+nbS[5][i])
	sigi = sqrt((sigr**2.+2.*sigt**2.)/3.)
		
	for kb in range(0,len(fk)-1):
		sigk = 0
		k = float(fk[kb].split()[0])
		pk = float(fk[kb].split()[1])*amp
		sigv21 = sigi**2./4.
		damp1 = exp(-1.*k*k*sigv21)	
		C1 = damp1
		pk1 = C1**2.*pk
		dk = ldk*k
		isigk = 0
		for i in range(0,len(nl)):
			isigk += vl[i]/(pk1**2.+2.*pk1/nl[i])#*exp(-1.*k**2.)
			#isigk += vl[i]/(pk1+1./nl[i]/sqrt(dr))**2.
		sigk = dk/(r1*r2)/isigk*(sin(k*r1)*sin(k*r2))
		sumxi += sigk*4.
	sumxi = sumxi/(2.*pi**2.)#/vol
	if r1 == r2: #pure shot-noise piece has been separated out since this depends on the r bin size
		npp = 0
		for i in range(0,len(nl)):
			npp += 2.*pi*r1**2.*dr*vl[i]*nl[i]**2.
		sumxi += 1./npp
	return sumxi


def mkcovxiNS(sample,version,amp,zmin,zmax,bs=10.,nzbs=20):
	#writes our analytic covariance matrix
	#Uses very approximate covariance, with terms added in rather ad-hoc way, DO NOT USE for actual science
	wm = 'mz'+str(zmin)+'xz'+str(zmax)
	fo = open(ebossdir+'covxiNS'+sample+version+wm+str(bs)+'.dat','w')
	r0 = bs/2.
	r1 = r0
	while r1 < 200:
		print r1
		r2 = r0
		while r2 < 200:
			if sample == 'BOSS':
				a = covcalcxi0boss(r1,r2,amp,zmin,zmax,dr=bs)
			else:
				a = covcalcxi0(r1,r2,version,sample,amp,zmin,zmax,dr=bs,nzbs=nzbs)
			fo.write(str(a)+' ')
			r2 += bs
		fo.write('\n')
		r1 += bs
	fo.close()
	return True

def docovsandbao():
	#make both of the covariance matrices for the current data, then does BAO
	mkcovxiNS('QSO','v0.7',1.3,.9,2.2) #QSO amplitude is ~1.3 compared to z = 0 clustering
	mkcovxiNS('lrg','v0.8_IRc',2.,.6,1.) #LRG amplitude is 2.0 compared to z = 0
	xibao('QSO','v0.7',0.9,2.2,'wdepth',md=1.,m=1.)
	xibao('lrg','v0.8_IRc',0.6,1.,'wstar')

### Find the BAO scale

def doallBAOqsomocks(sig=1,covmd='mock',bs=10):
	ma = 0
	sa = 0
	siga = 0
	chia = 0
	n = 0
	fo = open('BAOmockQPM_QSOv1.0'+covmd+str(bs)+'.dat','w')
	for i in range(1,101):
		zer = ''
		if i < 1000:
			zer += '0'
		if i < 100:
			zer += '0'
		if i < 10:
			zer += '0'
		mn = zer+str(i)	
		a = xibao('QPM_QSO',.9,2.2,mockn=mn,covmd=covmd,bs=bs)
		for j in range(0,len(a)):
			fo.write(str(a[j])+' ')
		fo.write('\n')
		if sig == 1:
			s1b = float(a[1]),float(a[2])
		if sig == 2:
			s1b = float(a[3]),float(a[4])	
		if s1b[0] > .8 and s1b[1] < 1.2:

			ma += a[0]
			sa += a[0]**2.
			siga += (float(a[2])-float(a[1]))/2.
			chia += a[-1]
			n += 1.
	ma = ma/n
	sa = sqrt(sa/n-ma**2.)
	siga = siga/n
	chia = chia/n
	return ma,sa,siga,chia,n

def putallBAOqsomocks(sig=1,covmd='mock',bs=10):
	d = np.loadtxt('BAOmockQPM_QSOv1.0'+covmd+str(bs)+'.dat')
	ma = 0
	sa = 0
	siga = 0
	chia = 0
	n = 0
	print len(d)
	for i in range(0,len(d)):
		a = d[i]
		if sig == 1:
			s1b = float(a[1]),float(a[2])
		if sig == 2:
			s1b = float(a[3]),float(a[4])	
		if s1b[0] > .8 and s1b[1] < 1.2:

			ma += a[0]
			sa += a[0]**2.
			siga += (float(a[2])-float(a[1]))/2.
			chia += a[-1]
			n += 1.
	ma = ma/n
	sa = sqrt(sa/n-ma**2.)
	siga = siga/n
	chia = chia/n
	return ma,sa,siga,chia,n
	

def xibao(sample,zmin,zmax,version='v1.0',wm='',bs=10,start=0,rmin=30,rmax=150.,md=1.,m=1.,mb='',Bp=.4,v='n',mockn='',covmd='mock'):
	#does baofits, set mb='nobao' to do no BAO fit
	from baofit_pubtest import doxi_isolike
	from Cosmo import distance
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	bsst = str(bs)+'st'+str(start)
	if sample == 'lrg' or sample == 'QSO':
		dn = np.loadtxt(ebossdir+'xi0geboss'+sample+'_N'+version+'_'+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0geboss'+sample+'_S'+version+'_'+wz+wm+bsst+'.dat').transpose()
	if sample == 'QPM_QSO':	
		dn = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'NGC'+version+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'SGC'+version+wz+wm+bsst+'.dat').transpose()
	if sample == 'aveQPM_QSO':	
		dn = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'NGC'+version+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'SGC'+version+wz+wm+bsst+'.dat').transpose()
	if sample == 'lrg':
		wt = (ds[1]*1.8+1.3*dn[1])/3.1
	if sample == 'QSO' or sample == 'QPM_QSO' or sample == 'aveQPM_QSO':
		if version == 'v1.2' or version == 'v1.5':
			wt = (ds[1]*.52+dn[1]*.66)/(.52+.66)
		else:
			wt = (ds[1]+dn[1])/2.
	if sample == 'LRGmod' or sample == 'QSOmod':
		dn = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.06.010.015.00_10st0.dat').transpose()	
	if sample == 'LRGmod':
		wt = dn[1]*2.
		sample = 'lrg'
	if sample == 'QSOmod':
		wt = dn[1]*1.3
		sample = 'QSO'
	if sample == 'BOSS':
		if zmin == 0.5:
			zbin = '3'
		dn = np.loadtxt(ebossdir+'dr12v5pre_bin'+zbin+'.xi').transpose()
		wt = dn[1]	
	dd = wt
	rl = dn[0]
	#print rl
	if mb == 'nobao':
		#mod = np.loadtxt('/Users/ashleyross/DR12/xi0smChallenge_matterpower0.43.06.010.015.00.dat').transpose()[1]
		mod = np.loadtxt('BAOtemplates/xi0smChallenge_matterpower0.43.02.55.015.00.dat').transpose()[1]
	else:
		#mod = np.loadtxt('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()[1]
		if sample == 'QPM_QSO':
			mod = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.02.55.015.00.dat').transpose()[1]
		else:
			mod = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.02.55.015.00.dat').transpose()[1]


	csample = sample
	if sample == 'aveQPM_QSO' or sample == 'QPM_QSO':
		if covmd == 'an':
			csample = 'QSO'
		if covmd == 'mock':
			csample = 'QPM_QSO'	
	if sample == 'QSO' and covmd == 'mock':
		csample = 'QPM_QSO'	
	if covmd == 'an':	
		cov = np.loadtxt(ebossdir+'covxiNS'+csample+version+wz+str(float(bs))+'.dat')
	if covmd == 'mock':
		cov = np.loadtxt(ebossdir+'cov0'+csample+version+'NScomb'+str(bs)+'st'+str(start)+'.dat')
	cov = cov*m	
	if md != 1:
		for i in range(0,len(cov)):
			cov[i][i] = cov[i][i]*md
				
	
	chil = doxi_isolike(dd,cov,mod,rl,rmin=rmin,rmax=rmax,v=v,wo=sample+version+mb,Bp=Bp)
	fo = open(ebossdir+'BAOxichil'+sample+mockn+version+mb+str(Bp)+'.dat','w')
	for i in range(0,len(chil)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chil[i])+'\n')
	fo.close()
	a = sigreg_c12(ebossdir+'BAOxichil'+sample+mockn+version+mb+str(Bp))
	#print a
	return a

def sigreg_c12(file,fac=1.):
	#report the confidence region +/-1 for chi2
	dir = ''
	f = open(file+'.dat').readlines()
	chil = []
	chim = 1000
	
	fl = []
	for i in range(0,len(f)):
		a = float(f[i].split()[0])
		#if a > min and a < max:
		chiv = float(f[i].split()[-1])*fac
		chil.append((chiv,a))
		if chiv < chim:
			#better to fit a parabola to get these values
			chim = chiv	
			im = i
			am = a
	a1u = 2.
	a1d = 0
	a2u = 2.
	a2d = 0
	oa = 0
	ocd = 0
	s0 = 0
	s1 = 0
	for i in range(im+1,len(chil)):
		chid = chil[i][0] - chim
		if chid > 1. and s0 == 0:
			a1u = (chil[i][1]/abs(chid-1.)+oa/abs(ocd-1.))/(1./abs(chid-1.)+1./abs(ocd-1.))
			s0 = 1
		if chid > 4. and s1 == 0:
			a2u = (chil[i][1]/abs(chid-4.)+oa/abs(ocd-4.))/(1./abs(chid-4.)+1./abs(ocd-4.))
			s1 = 1
		ocd = chid	
		oa = chil[i][1]
	oa = 0
	ocd = 0
	s0 = 0
	s1 = 0
	for i in range(1,im):
		chid = chil[im-i][0] - chim
		if chid > 1. and s0 == 0:
			a1d = (chil[im-i][1]/abs(chid-1.)+oa/abs(ocd-1.))/(1./abs(chid-1.)+1./abs(ocd-1.))
			s0 = 1
		if chid > 4. and s1 == 0:
			a2d = (chil[im-i][1]/abs(chid-4.)+oa/abs(ocd-4.))/(1./abs(chid-4.)+1./abs(ocd-4.))
			s1 = 1
		ocd = chid	
		oa = chil[im-i][1]
	if a1u < a1d:
		a1u = 2.
		a1d = 0
	if a2u < a2d:
		a2u = 2.
		a2d = 0
			
	return am,a1d,a1u,a2d,a2u,chim	


def compmat(NS='NScomb',m=1.):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	from numpy import loadtxt as load
	
	ma = load(ebossdir+'covxiNSQSOv1.0mz0.9xz2.210.0.dat')	
	mm = load(ebossdir+'cov0QPM_QSOv1.0'+NS+'10st0.dat')
	da = []
	dm = []
	th = []
	for i in range(0,len(ma)):
		da.append(sqrt(ma[i][i]))
		dm.append(sqrt(mm[i][i]))
		th.append(5.+10.*i)
	da = np.array(da)
	dm = np.array(dm)
	th = np.array(th)	
	plt.plot(th,da*th*m,'k-',th,dm*th,'r-')
	plt.show()
	md = np.zeros((15,15))
	for i in range(0,15):
		ri = 5.+i*10.
		for j in range(0,15):
			rj = 5.+j*10.
			md[i][j] = (ma[i][j]-mm[i][j])*ri*rj
	fig = plt.figure()
	ax = fig.add_subplot(111)
	res = ax.matshow(md)
	fig.colorbar(res)
	plt.show()
	return True

def testzweights(p=1.):
	from Cosmo import distance
	from random import gauss
	#test difference in recovered signal to noise with and without weighting
	d = distance(.3,.7)
	url = []
	wrl = []
	for j in range(0,1000):
		z = .905
		l = []
		wl = []
		while z < 2.2:
			w = 1./d.D(z)**p
			sigma = sqrt(1./w)
			for i in range(0,10):
				v = gauss(0,sigma)
				l.append(v)
				wl.append(w)
			z += .01
		u = np.mean(l)
		wr = sum(np.array(l)*np.array(wl))/sum(wl)
		url.append(u)
		wrl.append(wr)
	print len(url),len(wrl),np.mean(url),np.mean(wrl)	
	vu = np.mean(np.array(url)**2.)-np.mean(url)**2.
	vw = np.mean(np.array(wrl)**2.)-np.mean(wrl)**2.
	return vu,vw

def nzQSOgal(sample='QSO',version='v1.5',zmin=.9,zmax=2.2,zb=.05):
	from matplotlib import pyplot as plt
	fn = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-N'+'-eboss_'+version+'-full.dat.fits')
	fs = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-S'+'-eboss_'+version+'-full.dat.fits')
	nzls = []
	nzln = []
	for i in range(0,int(1/float(zb)*zmax)):
		nzls.append(0)
		nzln.append(0)
	sumzw = 0
	sumw = 0
	for i in range(0,len(fn)):
		z = fn[i]['Z']
		im = fn[i]['IMATCH']
		if z > zmin and z < zmax and im == 9:
			zind = int(1/float(zb)*z)
			nzln[zind] += 1.#*ff[i]['WEIGHT_FKP']
			sumzw += z#*ff[i]['WEIGHT_FKP']
			sumw += 1.#*ff[i]['WEIGHT_FKP']
	print sumzw/sumw
	sumzws = 0
	sumws = 0

	for i in range(0,len(fs)):
		z = fs[i]['Z']
		im = fs[i]['IMATCH']
		if z > zmin and z < zmax and im == 9:
			zind = int(1/float(zb)*z)
			w = 1.
			nzls[zind] += w
			sumzws += z*w
			sumws += w
	print sumzws/sumws
# 	for i in range(0,len(fd[0])):
# 		z = fd[2][i]
# 		if z > zmin and z < zmax:
# 			zind = int(1/float(zb)*z)
# 			nzld[zind] += fd[3][i]
	nzls = np.array(nzls)/sum(nzls)/zb		
	#nzld = np.array(nzld)/sum(nzld)
	nzln = np.array(nzln)/sum(nzln)/zb
	#print nzld
	nzl = np.arange(0,zmax,zb)
	plt.plot(nzl,nzln,'r-',nzl,nzls,'b-')
	plt.xlim(zmin,zmax)
	plt.xlabel('redshift')
	plt.ylabel('normalized QSO-galaxy density')
	plt.show()
	return True	


def compnz(sample='QSO',NS='N',version='v1.2',zmin=.9,zmax=2.2,zb=.05):
	from matplotlib import pyplot as plt
	ff = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.ran.fits')
	fdf = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat.fits')
	#fd = np.loadtxt(ebossdir+'reboss'+sample+'_'+NS+version+'_0mz'+str(zmin)+'xz'+str(zmax)+'4xi.dat').transpose()
	nzlf = []
	nzld = []
	nzldd = []
	for i in range(0,int(1/float(zb)*zmax)):
		nzlf.append(0)
		nzld.append(0)
		nzldd.append(0)
	sumzw = 0
	sumw = 0
	for i in range(0,len(ff)):
		z = ff[i]['Z']
		if z > zmin and z < zmax:
			zind = int(1/float(zb)*z)
			nzlf[zind] += 1.#*ff[i]['WEIGHT_FKP']
			sumzw += z#*ff[i]['WEIGHT_FKP']
			sumw += 1.#*ff[i]['WEIGHT_FKP']
	print sumzw/sumw
	sumzwd = 0
	sumwd = 0

	for i in range(0,len(fdf)):
		z = fdf[i]['Z']
		if z > zmin and z < zmax:
			zind = int(1/float(zb)*z)
			w = (fdf[i]['WEIGHT_NOZ']+fdf[i]['WEIGHT_CP']-1.)*fdf[i]['WEIGHT_SYSTOT']#*fdf[i]['WEIGHT_FKP']
			nzldd[zind] += w
			sumzwd += z*w
			sumwd += w
	print sumzwd/sumwd
# 	for i in range(0,len(fd[0])):
# 		z = fd[2][i]
# 		if z > zmin and z < zmax:
# 			zind = int(1/float(zb)*z)
# 			nzld[zind] += fd[3][i]
	nzlf = np.array(nzlf)/sum(nzlf)		
	#nzld = np.array(nzld)/sum(nzld)
	nzldd = np.array(nzldd)/sum(nzldd)
	#print nzld
	nzl = np.arange(0,zmax,zb)
	plt.plot(nzl,nzlf/nzldd,'r-')#,nzl,nzld/nzldd,'b-')
	plt.xlim(zmin,zmax)
	plt.xlabel('redshift')
	plt.ylabel('normalized random density / QSO density')
	plt.show()
	return True	

def compnzfib(sample='QSO',version='v1.3',zmin=.9,zmax=2.2,zb=.05,rz=.1):
	from matplotlib import pyplot as plt
	from random import random,gauss
	fn = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+'N-eboss_'+version+'.dat.fits')
	fs = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+'S-eboss_'+version+'.dat.fits')
	nzlf = []
	nzld = []
	nzldd = []
	for i in range(0,int(1/float(zb)*zmax)):
		nzlf.append(0)
		nzld.append(0)
		nzldd.append(0)
	n = 0
	for i in range(0,len(fn)):
		z = fn[i]['Z']
		if z > zmin and z < zmax:
			zind = int(1/float(zb)*z)
			ranz = random()*1.3+.9
			zindr = int(1/float(zb)*ranz)
			rf = random()
			zindf = zind
			if rf < rz:
				zindf = zindr
			#w =1.
			w = (fn[i]['WEIGHT_CP'])*fn[i]['WEIGHT_SYSTOT']#*fdf[i]['WEIGHT_FKP']
			fid = fn[i]['FIBERID']
			s = 0
			if fid < 0:
				n += 1
				s = 1
			if fid <= 50 and fid > 0:
				nzlf[zind] += w
				s =1
			if fid >450 and fid <= 550:
				nzlf[zind] += w
				s =1
			if fid > 950 and fid <=1000:
				nzlf[zind] += w	
				s = 1
			if s == 0:
				nzldd[zindf] += w
	print n
	for i in range(0,len(fs)):
		z = fs[i]['Z']
		if z > zmin and z < zmax:
			zind = int(1/float(zb)*z)
			ranz = random()*1.3+.9
			zindr = int(1/float(zb)*ranz)
			rf = random()
			zindf = zind
			if rf < rz:
				zindf = zindr
			#w =1.
			w = (fs[i]['WEIGHT_CP'])*fs[i]['WEIGHT_SYSTOT']#*fdf[i]['WEIGHT_FKP']
			fid = fs[i]['FIBERID']
			s = 0
			if fid < 0:
				s = 1
			if fid <= 50 and fid > 0:
				nzlf[zind] += w
				s =1
			if fid >450 and fid <= 550:
				nzlf[zind] += w
				s =1
			if fid > 950 and fid <=1000:
				nzlf[zind] += w	
				s = 1
			if s == 0:
				nzldd[zindf] += w
# 	for i in range(0,len(fd[0])):
# 		z = fd[2][i]
# 		if z > zmin and z < zmax:
# 			zind = int(1/float(zb)*z)
# 			nzld[zind] += fd[3][i]
	nzfe = (np.array(nzlf))**.5/sum(nzlf)
	nzlf = np.array(nzlf)/sum(nzlf)		
	#nzld = np.array(nzld)/sum(nzld)
	nzldd = np.array(nzldd)/sum(nzldd)
		
	#print nzld
	nzl = np.arange(zb/2.,zmax,zb)
	c = np.ones((len(nzl)))
	c = c*.98
	plt.plot(nzl,nzlf,'r-',nzl,nzldd,'b-')
	plt.show()
	plt.plot(nzl,c,'k:')
	plt.errorbar(nzl,nzlf/nzldd,nzfe/nzldd,fmt='r-')#,nzl,nzld/nzldd,'b-')
	plt.xlim(zmin-zb/2.,zmax+zb/2.)
	plt.xlabel('redshift')
	plt.ylabel('normalized edge fiber density / other density')
	plt.show()
	return True	

def compnzknown(sample='QSO',version='v1.3',zmin=.9,zmax=2.2,zb=.05,rz=.1):
	from matplotlib import pyplot as plt
	from random import random,gauss
	fn = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+'N-eboss_'+version+'.dat.fits')
	fs = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+'S-eboss_'+version+'.dat.fits')
	nzlf = []
	nzld = []
	nzldd = []
	for i in range(0,int(1/float(zb)*zmax)):
		nzlf.append(0)
		nzld.append(0)
		nzldd.append(0)
	n = 0
	for i in range(0,len(fn)):
		z = fn[i]['Z']
		if z > zmin and z < zmax:
			zind = int(1/float(zb)*z)
			#w =1.
			w = (fn[i]['WEIGHT_CP']+fn[i]['WEIGHT_NOZ']-1.)*fn[i]['WEIGHT_SYSTOT']#*fdf[i]['WEIGHT_FKP']
			fid = fn[i]['FIBERID']
			nzldd[zind] += w
			if fid < 0:
				nzlf[zind] += w
			if fid > 0:
				nzld[zind] += w
	print n
	for i in range(0,len(fs)):
		z = fs[i]['Z']
		if z > zmin and z < zmax:
			zind = int(1/float(zb)*z)
			#w =1.
			w = (fs[i]['WEIGHT_CP']+fs[i]['WEIGHT_NOZ']-1.)*fs[i]['WEIGHT_SYSTOT']#*fdf[i]['WEIGHT_FKP']
			fid = fs[i]['FIBERID']
			nzldd[zind] += w
			if fid < 0:
				nzlf[zind] += w
			if fid > 0:
				nzld[zind] += w
	nzfe = (np.array(nzlf))**.5/sum(nzlf)
	nzlf = np.array(nzlf)/sum(nzlf)		
	nzld = np.array(nzld)/sum(nzld)
	nzldd = np.array(nzldd)/sum(nzldd)
		
	#print nzld
	nzl = np.arange(zb/2.,zmax,zb)
	c = np.ones((len(nzl)))
	c = c*.98
	plt.plot(nzl,nzlf,'r-',nzl,nzldd,'k-',nzl,nzld,'b-')
	plt.show()
# 	plt.plot(nzl,c,'k:')
# 	plt.errorbar(nzl,nzlf/nzldd,nzfe/nzldd,fmt='r-')#,nzl,nzld/nzldd,'b-')
# 	plt.xlim(zmin-zb/2.,zmax+zb/2.)
# 	plt.xlabel('redshift')
# 	plt.ylabel('normalized edge fiber density / other density')
# 	plt.show()
	return True	


def compextfib2(sample='lrg',NS='N',version='v1.0_IRt'):
	fdf = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat.fits')
	mgl = []
	el = []
	for i in range(0,len(fdf)):
		magf = luptm(fdf[i]['FIBER2FLUX'][-1],4)
		mgl.append(magf)
		el.append(fdf[i]['EXTINCTION'][4])
	mgl = np.array(mgl)
	el = np.array(el)
	print sum(mgl)/float(len(mgl))
	print (sum(mgl)-sum(el))/float(len(mgl))
	print max(mgl)
	print max(mgl-el)
	print max(el)
	print sum(el)/float(len(mgl))
	from matplotlib import pyplot as plt
	plt.hist((mgl-el),range=(21,max(mgl-el)))
	plt.show()
	return True

def compextfib2NS(sample='lrg',version='v1.0_IRt'):
	fdfn = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-N-eboss_'+version+'.dat.fits')
	fdfs = fitsio.read(dirfits+'eboss_'+version+'-'+sample+'-S-eboss_'+version+'.dat.fits')
	mgln = []
	eln = []
	for i in range(0,len(fdfn)):
		magf = luptm(fdfn[i]['FIBER2FLUX'][-1],4)
		mgln.append(magf)
		eln.append(fdfn[i]['EXTINCTION'][4])
	mgln = np.array(mgln)
	eln = np.array(eln)
	mgls = []
	els = []
	for i in range(0,len(fdfs)):
		magf = luptm(fdfs[i]['FIBER2FLUX'][-1],4)
		mgls.append(magf)
		els.append(fdfs[i]['EXTINCTION'][4])
	mgls = np.array(mgls)
	els = np.array(els)
	#from matplotlib import pyplot as plt
	import pylab as P
	#plt.hist((mgln-eln),range=(21,21.6),normed=1)
	#plt.hist((mgls-els),range=(21,21.6),normed=1)
	P.figure()
	P.hist([(mgln-eln),(mgls-els)],range=(21,21.6),normed=1,label=['NGC','SGC'])
	P.title('extinction correction to magnitude from flux')
	P.legend()
	P.show()
	P.figure()
	P.hist([(mgln),(mgls)],range=(21,21.7),normed=1,label=['NGC','SGC'])
	P.title('No extinction correction to magnitude from flux')
	P.legend()
	P.show()

	#plt.clf()
	#plt.hist((mgln),range=(21,21.7),normed=1)
	#plt.hist((mgls),range=(21,21.7),normed=1)
	#plt.show()
	#plt.clf()
	return True
	
	
### Some modules to count redshift failures, etc.

def countzfcp(sample,NS,version,zmin=.6,zmax=1.,c='',app='.fits',wm=''):
	#note, zmin/zmax assume LRG sample these need to be change for QSO files
	from healpix import healpix, radec2thphi
	if wm == 'wstar' or wm == 'cpstar':
		wsys = np.loadtxt('allstars17.519.9Healpixall256.dat')
		b,m = findlinmb(sample,NS,version,'star',zmin,zmax)
		h = healpix()
	if wm == 'wdepth' or wm == 'cpdepth':
		wsys = np.loadtxt('healdepthinm512.dat').transpose() #note, file is zipped in directory	
		b,m = findlinmb(sample,NS,version,'depth',zmin,zmax)
		h = healpix()
	if c == 'sci': #AJR uses this define directory for machine he uses
		dir = dirsci
	else:
		dir = dirfits	
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat'+app) #read galaxy/quasar file
	nzf = 0
	ncp = 0
	nt = 0
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z > zmin and z < zmax:
			nt += (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)
			nzf += (f[i]['WEIGHT_NOZ']-1.)
			ncp += (f[i]['WEIGHT_CP']-1.)
	print nt,nzf/nt,ncp/nt
	return True


### Below is for plots


def plotxiLRGNS(bs='10st0',v='v0.8_IRc',wm=''):
	#Plots comparison between NGC and SGC clustering for and to theory LRGs, no stellar density correction
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiLRGNS'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	cov = np.loadtxt(ebossdir+'covxiNSlrgv0.8_IRcmz0.6xz1.010.dat')
	et = []
	for i in range(0,len(cov)):
		et.append(sqrt(cov[i][i]))
	#plt.fill_between(dt[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	ds = np.loadtxt(ebossdir+'xi0gebosslrg_S'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi0gebosslrg_N'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	
	t = (ds[1]*.8+.75*dn[1])/1.55
	plt.plot(ds[0],ds[0]**2.*ds[1],'b-s')
	plt.plot(dn[0],dn[0]**2.*dn[1],'r-d')
	plt.plot(dn[0],dn[0]**2.*t,'k-o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	plt.xlim(20,165)
	plt.ylim(-49,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,'SGC',color='b')
	plt.text(30,170,'NGC',color='r')
	plt.text(30,160,'Combined',color='k')
	plt.title(r'Correlation function of v0.8_IRc eboss only LRGs, 0.6 < z < 1.0')
	pp.savefig()
	pp.close()
	return True


def plotxiQSONS(bs='10st0',v='v0.7'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONS'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.02.55.015.00.dat').transpose()
	t = (ds[1]*.8+dn[1])/1.8
	plt.plot(ds[0],ds[0]**2.*ds[1],'b-s')
	plt.plot(dn[0],dn[0]**2.*dn[1],'r-d')
	plt.plot(dn[0],dn[0]**2.*t,'k-o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*1.3,'k:')
	plt.xlim(20,165)
	plt.ylim(-49,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,'SGC',color='b')
	plt.text(30,170,'NGC',color='r')
	plt.text(30,160,'Combined',color='k')
	plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True

def plotxiQSOvcomp(bs='10st0',v1='v0.7',v2='v1.0'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONS'+v1+'comp'+v2+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	ds1 = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v1+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dn1 = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v1+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	t1 = (ds1[1]*1.+dn1[1])/2.
	ds2 = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v2+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dn2 = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v2+'_mz0.9xz2.2'+bs+'.dat').transpose()
	t2 = (ds2[1]*1.+dn2[1])/2.
	plt.plot(dn1[0],dn1[0]**2.*t1,'b-s')
	plt.plot(dn1[0],dn1[0]**2.*t2,'r-d')
	plt.plot(dt[0],dt[0]**2.*dt[1]*1.3,'k:')
	plt.xlim(20,200)
	plt.ylim(-49,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,v1,color='b')
	plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True

def plotxiQSOcompw(bs='5st0',v='v1.3'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONScompdepthw'+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	ds1 = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dn1 = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	t1 = (ds1[1]*.5+dn1[1]*.66)/1.16
	ds2 = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dn2 = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	t2 = (ds2[1]*.5+dn2[1]*.66)/1.16
	ds3 = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v+'_mz0.9xz2.2gx22wdepth'+bs+'.dat').transpose()
	dn3 = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v+'_mz0.9xz2.2gx22wdepth'+bs+'.dat').transpose()
	t3 = (ds3[1]*.5+dn3[1]*.66)/1.16

	plt.plot(dn1[0],dn1[0]**2.*t1,'b-s')
	plt.plot(dn1[0],dn1[0]**2.*t2,'r-d')
	plt.plot(dn1[0],dn1[0]**2.*t3,'g-o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*1.3,'k:')
	plt.xlim(20,250)
	plt.ylim(-49,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,180,v1,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True


def plotxiQSOvcomps(bs='10st0',v1='v0.7',v2='v1.0'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xiQSONS'+v1+'comps'+v2+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	ds1 = np.loadtxt('xi0gebossQSO_S'+v1+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dn1 = np.loadtxt('xi0gebossQSO_N'+v1+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dt = np.loadtxt('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	t1 = (ds1[1]*1.+dn1[1])/2.
	ds2 = np.loadtxt('xi0gebossQSO_S'+v2+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dn2 = np.loadtxt('xi0gebossQSO_N'+v2+'_mz0.9xz2.2'+bs+'.dat').transpose()
	t2 = (ds2[1]*1.+dn2[1])/2.
	plt.plot(dn1[0],dn1[0]*t1,'b-s')
	plt.plot(dn1[0],dn1[0]*t2,'r-d')
	plt.plot(dt[0],dt[0]*dt[1]*1.3,'k:')
	plt.xlim(0,150)
	plt.ylim(-1,8)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s\xi(s)$ ($h^{-1}$Mpc)',size=16)
	plt.text(30,7.5,v1,color='b')
	plt.text(30,7,v2,color='r')
	plt.title(r'Correlation function of quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True


def plotxiQSOv(bs='10st0',v='v1.2'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONS'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	ds1 = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dn1 = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dna = np.loadtxt(ebossdir+'xi0gaveQPM_QSONGCv1.0mz0.9xz2.2'+bs+'.dat').transpose()
	dsa = np.loadtxt(ebossdir+'xi0gaveQPM_QSOSGCv1.0mz0.9xz2.2'+bs+'.dat').transpose()

	dt = np.loadtxt('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.02.55.015.00.dat').transpose()
	#if v == 'v1.2':
	wn = 0.66
	ws = .53
		
	if v == 'v1.0':
		wn = .5
		ws = .5
	t1 = (ds1[1]*ws+dn1[1]*wn)/(wn+ws)
	ta = (dsa[1]*1.+dna[1])/2.
	te = dna[2]/sqrt(2.)
	plt.errorbar(dn1[0],dn1[0]**2.*t1,dn1[0]**2.*te,fmt='bs')
	plt.plot(dt[0],dt[0]**2.*dt[1]*1.3,'k:')
	plt.errorbar(dn1[0],dn1[0]**2.*ta,dn1[0]**2.*te/10.,fmt='k-')
	plt.xlim(20,250)
	plt.ylim(-49,100)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True


def plotxiLRGNSwstar(bs='10st0',v='v1.0_IRt'):
	#Plots comparison between NGC and SGC clustering and to theory for LRGs, with stellar density correction
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xiLRGNSwstaro'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dt = np.loadtxt('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	cov = np.loadtxt('covxiNSlrgv0.8_IRcmz0.6xz1.010.dat')
	et = []
	for i in range(0,len(cov)):
		et.append(sqrt(cov[i][i]))
	#plt.fill_between(dt[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	#dswz = np.loadtxt('xi0gebosslrg_S'+v+'_mz0.6xz1.0cpstar'+bs+'.dat').transpose()
	#dnwz = np.loadtxt('xi0gebosslrg_N'+v+'_mz0.6xz1.0cpstar'+bs+'.dat').transpose()
	#wtz = (dswz[1]*1.8+dnwz[1])/2.8
	if v == 'v0.9_IRc':
		dsw = np.loadtxt('xi0gebosslrg_S'+v+'_mz0.6xz1.0wstar'+bs+'.dat').transpose()
		dnw = np.loadtxt('xi0gebosslrg_N'+v+'_mz0.6xz1.0wstar'+bs+'.dat').transpose()
		wt = (dsw[1]*2.2+1.7*dnw[1])/3.9
	else:			
		dsw = np.loadtxt('xi0gebosslrg_S'+v+'_mz0.6xz1.0st'+bs+'.dat').transpose()
		dnw = np.loadtxt('xi0gebosslrg_N'+v+'_mz0.6xz1.0st'+bs+'.dat').transpose()
		wt = (dsw[1]+dnw[1])/2.


	#ds = np.loadtxt('xi0gebosslrg_S'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	#dn = np.loadtxt('xi0gebosslrg_N'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	
	#t = (ds[1]*1.8+1.3*dn[1])/3.1
	#plt.plot(ds[0],ds[0]**2.*ds[1],'b--')
	#plt.plot(dn[0],dn[0]**2.*dn[1],'r--')
	#plt.plot(dnw[0],dnw[0]**2.*t,'k--')
	plt.plot(dsw[0],dsw[0]**2.*dsw[1],'b-s')
	plt.plot(dnw[0],dnw[0]**2.*dnw[1],'r-d')
	plt.errorbar(dnw[0][:20],dnw[0][:20]**2.*wt[:20],dnw[0][:20]**2.*et,fmt='k-o')
	#plt.plot(dsw[0],dsw[0]**2.*dswz[1],'b--s')
	#plt.plot(dnw[0],dnw[0]**2.*dnwz[1],'r--d')
	#plt.errorbar(dnw[0][:20],dnw[0][:20]**2.*wtz[:20],dnw[0][:20]**2.*et,fmt='k--o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	plt.xlim(20,165)
	plt.ylim(-50,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,'SGC',color='b')
	plt.text(30,170,'NGC',color='r')
	plt.text(30,160,'Combined',color='k')
	plt.text(90,175,'Corrected for stellar density')
	plt.title(r'Correlation function of v1.0_IRt eboss only LRGs, 0.6 < z < 1.0')
	pp.savefig()
	pp.close()
	return True

def plotxiQSONSwdepth(bs='10st0',v='v0.9'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, with depth correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xiQSONSwdeptho'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	cov = np.loadtxt('covxiNSQSOv0.7mz0.9xz2.210.dat')
	#dswz = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2cpdepth'+bs+'.dat').transpose()
	#dnwz = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2cpdepth'+bs+'.dat').transpose()
	#wtz = (dswz[1]*1.+dnwz[1])/2.
	
	dsw = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dnw = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()

	et = []
	for i in range(0,len(cov)):
		et.append(sqrt(cov[i][i]))
	#plt.fill_between(dsw[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	wt = (dsw[1]*1.+dnw[1])/2.
	#ds = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	#dn = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dt = np.loadtxt('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	#t = (ds[1]*1.+dn[1])/2.
	#plt.plot(ds[0],ds[0]**2.*ds[1],'b--')
	#plt.plot(dn[0],dn[0]**2.*dn[1],'r--')
	#plt.plot(dnw[0],dnw[0]**2.*t,'k--')
	plt.plot(dsw[0],dsw[0]**2.*dsw[1],'b-s')
	plt.plot(dnw[0],dnw[0]**2.*dnw[1],'r-d')
	plt.errorbar(dnw[0][:20],dnw[0][:20]**2.*wt[:20],dnw[0][:20]**2.*et,fmt='k-o')
	#plt.plot(dsw[0],dsw[0]**2.*dswz[1],'b--s')
	#plt.plot(dnw[0],dnw[0]**2.*dnwz[1],'r--d')
	#plt.errorbar(dnw[0][:20],dnw[0][:20]**2.*wtz[:20],dnw[0][:20]**2.*et,fmt='k--o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*1.3,'k:')
	plt.xlim(20,165)
	plt.ylim(-50,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,'SGC',color='b')
	plt.text(30,170,'NGC',color='r')
	plt.text(30,160,'Combined',color='k')
	plt.text(90,175,'Corrected for depth')
	plt.title(r'Correlation function of v0.9 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True

def plotxiLRGNSbaofit(bs='10st0',v='v0.8_IRc'):
	#Plots comparison between LRG clustering and best-fit BAO theory
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xiLRGNSbaofit'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dt = np.loadtxt('ximodlrgv0.8_IRc.dat').transpose()
	dtn = np.loadtxt('ximodlrgv0.8_IRcnobao.dat').transpose()
	cov = np.loadtxt('covxiNSlrgv0.8_IRcmz0.6xz1.010.dat')
	et = []
	for i in range(0,len(cov)):
		et.append(sqrt(cov[i][i]))
	#plt.fill_between(dt[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	dsw = np.loadtxt('xi0gebosslrg_S'+v+'_mz0.6xz1.0wstar'+bs+'.dat').transpose()
	dnw = np.loadtxt('xi0gebosslrg_N'+v+'_mz0.6xz1.0wstar'+bs+'.dat').transpose()
	wt = (dsw[1]*1.8+dnw[1])/2.8
	ds = np.loadtxt('xi0gebosslrg_S'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	dn = np.loadtxt('xi0gebosslrg_N'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	
	t = (ds[1]*1.8+1.3*dn[1])/3.1
	plt.errorbar(dnw[0][:20],dnw[0][:20]**2.*wt[:20],dnw[0][:20]**2.*et,fmt='ko')
	plt.plot(dt[0],dt[0]**2.*dt[1],'k-')
	plt.plot(dt[0],dt[0]**2.*dtn[1],'k--')
	plt.xlim(20,165)
	plt.ylim(-49,150)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,120,r'$\alpha=0.937\pm0.027$',color='k',size=16)
	plt.text(35,105,r'$\chi^2$/dof = 6.9/7',color='k',size=16)
	plt.title(r'BAO best-fit for v0.8_IRc eboss only LRGs, 0.6 < z < 1.0')
	pp.savefig()
	pp.close()
	return True

def plotxiQSONSbaofit(bs='10st0',v='v1.0'):
	#Plots comparison between QSO clustering and best-fit BAO theory
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xiQSONSbaofit'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	if bs == '10st0':
		bsc = '10.0'
	if bs == '5st0':
		bsc = '5.0'
	dt = np.loadtxt('ximodQSO'+v+'.dat').transpose()
	dtn = np.loadtxt('ximodQSO'+v+'nobao.dat').transpose()
	cov = np.loadtxt('covxiNSQSO'+v+'mz0.9xz2.2'+bsc+'.dat')
	et = []
	for i in range(0,len(cov)):
		et.append(sqrt(cov[i][i]))
	#plt.fill_between(dt[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	if v != 'v1.0':
		dsw = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
		dnw = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	else:
		dsw = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
		dnw = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	fac = int(10/float(bsc))	
	wt = (dsw[1]*1.+dnw[1])/2.
	plt.errorbar(dnw[0][:20*fac],dnw[0][:20*fac]**2.*wt[:20*fac],dnw[0][:20*fac]**2.*et,fmt='ko')
	plt.plot(dt[0],dt[0]**2.*dt[1],'k-')
	plt.plot(dt[0],dt[0]**2.*dtn[1],'k--')
	plt.xlim(20,170)
	plt.ylim(-35,80)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,71,r'$\alpha=1.044\pm0.042$',color='k',size=16)
	#plt.text(35,64,r'$\chi^2$/dof = 5.8/7',color='k',size=16)
	plt.text(30,71,r'$\alpha=1.038\pm0.045$',color='k',size=16)
	plt.text(35,64,r'$\chi^2$/dof = 19/19',color='k',size=16)
	plt.title(r'BAO best-fit for v1.0 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True


def plotxiQSONSwdeptha(bs='10st0',v='v0.7'):
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xiQSONSwdeptho'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	cov = np.loadtxt('covxiNSQSOv0.7mz0.9xz2.210.dat')
	dsw = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	dnw = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()

	et = []
	for i in range(0,len(cov)):
		et.append(sqrt(cov[i][i]))
	#plt.fill_between(dsw[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	wt = (dsw[1]*1.+dnw[1])/2.
	ds = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dn = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2'+bs+'.dat').transpose()
	dt = np.loadtxt('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	t = (ds[1]*1.+dn[1])/2.
	#plt.plot(ds[0],ds[0]**2.*ds[1],'b--')
	#plt.plot(dn[0],dn[0]**2.*dn[1],'r--')
	#plt.plot(dnw[0],dnw[0]**2.*t,'k--')
	plt.plot(dsw[0],dsw[0]**2.*dsw[1],'b-s')
	plt.plot(dnw[0],dnw[0]**2.*dnw[1],'r-d')
	plt.errorbar(dnw[0][:20],dnw[0][:20]**2.*wt[:20],dnw[0][:20]**2.*et,fmt='k-o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*1.3,'k:')
	plt.xlim(20,165)
	plt.ylim(-50,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,'SGC',color='b')
	plt.text(30,170,'NGC',color='r')
	plt.text(30,160,'Combined',color='k')
	plt.text(90,175,'Corrected for depth')
	plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True


def plotLRGNSvsstar(v='v0.8_IRc'):
	#plots N_LRG vs. N_star
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('nLRGNS'+v+'vstar.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt('ngebosslrg_S'+v+'_mz0.6xz1.0256vstar.dat').transpose()
	dn = np.loadtxt('ngebosslrg_N'+v+'_mz0.6xz1.0256vstar.dat').transpose()
	
	plt.errorbar(ds[0],ds[1],ds[2],fmt='bs')
	plt.errorbar(ds[0],dn[1],dn[2],fmt='rd')
	plt.xlabel(r'$N_{\rm star}$ (per Nside=256 pixel)',size=16)
	plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
	plt.ylim(.76,1.19)
	plt.text(60,.8,'SGC',color='b')
	plt.text(60,.78,'NGC',color='r')
	plt.title(r'galaxy density vs. stellar density for v0.8_IRc eboss only LRGs, 0.6 < z < 1.0')
	pp.savefig()
	pp.close()
	return True


def plotQSONSvsdepth(v='v0.7'):
	#plots N_QSO vs. i-band depth
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('nQSONS'+v+'vdepth.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt('ngebossQSO_S'+v+'_mz0.9xz2.2512vdepth.dat').transpose()
	dn = np.loadtxt('ngebossQSO_N'+v+'_mz0.9xz2.2512vdepth.dat').transpose()
	
	plt.errorbar(ds[0],ds[1],ds[2],fmt='bs')
	plt.errorbar(ds[0],dn[1],dn[2],fmt='rd')
	plt.xlabel(r'5$\sigma$ $i$-band depth, averaged in Nside=512 pixels (magnitudes)',size=16)
	plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
	plt.ylim(.76,1.19)
	plt.text(22,1.1,'SGC',color='b')
	plt.text(22,1.08,'NGC',color='r')
	plt.title(r'galaxy density vs. $i$-band depth for v0.7 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True

def plotQSOgmagNSvsdepth(v='v1.3'):
	#plots N_QSO vs. i-band depth
	from optimize import fmin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'nQSONS'+v+'vdepth.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2gm020.5512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2gm020.5512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='ko')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl = []
	ml = []
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'k--')
	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2gm20.521.0512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2gm20.521.0512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='rd')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'r--')

	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2gm21.021.5512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2gm21.021.5512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='bs')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'b--')

	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2gm21.522.0512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2gm21.522.0512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='g^')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'g--')
	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2gm22.030512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2gm22.030512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='^',color='purple')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'--',color='purple')
	plt.text(22.6,.81,r'$g < 20.5$',fontsize=18,color='k')
	plt.text(22.6,.77,r'$20.5 < g < 21$',fontsize=18,color='r')
	plt.text(22.6,.73,r'$21 < g < 21.5$',fontsize=18,color='b')
	plt.text(22.6,.69,r'$21. 5 < g < 22$',fontsize=18,color='g')
	plt.text(22.6,.65,r'$g > 22$',fontsize=18,color='purple')
	plt.xlabel(r'5$\sigma$ $i$-band depth, averaged in Nside=512 pixels (magnitudes)',size=16)
	plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
	plt.ylim(.6,1.4)
	plt.title(r'galaxy density vs. $i$-band depth for v1.3 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return bl,ml

	
def plotLRGNSbaolike(v='v0.8_IRc'):
	#plot bao likelihood for LRGs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xiLRGNSbaolik'+v+'.pdf')
	plt.clf()
	plt.minorticks_on()
	db = np.loadtxt('BAOxichillrgv0.8_IRc0.4.dat').transpose()
	dnb = np.loadtxt('BAOxichillrgv0.8_IRcnobao0.4.dat').transpose()
	chim = min(db[1])
	plt.plot(db[0],db[1]-chim,'k-')
	plt.plot(db[0],dnb[1]-chim,'k--')
	#plt.xlim(20,165)
	#plt.ylim(-49,150)
	plt.xlabel(r'$\alpha_{\rm BAO}$',size=16)
	plt.ylabel(r'$\Delta\chi^2$',size=16)
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	plt.title(r'BAO likelihood for v0.8_IRc eboss only LRGs, 0.6 < z < 1.0')
	pp.savefig()
	pp.close()
	return True

def plotQSONSbaolike(v='v1.0'):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONSbaolik'+v+'.pdf')
	plt.clf()
	plt.minorticks_on()
	db = np.loadtxt(ebossdir+'BAOxichilQSO'+v+'0.4.dat').transpose()
	dnb = np.loadtxt(ebossdir+'BAOxichilQSO'+v+'nobao0.4.dat').transpose()
	chim = min(db[1])
	plt.plot(db[0],db[1]-chim,'k-')
	plt.plot(db[0],dnb[1]-chim,'k--')
	#plt.xlim(20,165)
	#plt.ylim(-49,150)
	plt.xlabel(r'$\alpha_{\rm BAO}$',size=16)
	plt.ylabel(r'$\Delta\chi^2$',size=16)
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	plt.title(r'BAO likelihood for '+v+' eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True
		