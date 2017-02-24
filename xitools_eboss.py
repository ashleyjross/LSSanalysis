dirsci = '/mnt/lustre/ashleyr/eboss/' #where AJR puts eboss catalogs, change this to wherever you have put catalogs
dirsys = '/Users/ashleyross/ngalvsys_wiki/' #change to local directory where ngalvsys from wiki was put, note star map and depth map included
dirfits = '/Users/ashleyross/fitsfiles/' #change to where your catalog files are
ebossdir = '/Users/ashleyross/eboss/' #where AJR puts correlation functions, writes out results
dirscio = '/mnt/lustre/ashleyr/eboss/mocks/'

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

def mkran1mil(sample,NS,version,cm='',N=0,c='sci',app='.ran.fits',compmin=0):
	dirout = ''
	dir = ''
	if c == 'sci':
		dir = dirsci
		dirout = dir
	minc = N*10**6
	maxc = (N+1)*10**6 #will become relevant once files are big enough
	f = fitsio.read(dir+cm+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+app)#[minc:maxc]
	if maxc > len(f):
		maxc = len(f)
	#f = fitsio.FITS(dir+'eboss_'+version+'-'+sample+'-'+NS+'eboss_'+version+'.ran'+app)[1]
	no = 0
	#if app == '.2.ran.fits':
	#	N = N+2
	fo = open(dirout+'reboss'+cm+sample+'_'+NS+version+'_'+str(N)+'.dat','w')
	#for i in range(0,len(f)):		
	for i in range(minc,maxc):
		#comp = f[i]['COMP']
		#if comp > compmin:
		fo.write(str(f[i]['RA'])+' '+str(f[i]['DEC'])+'\n') #rest of the information comes from matching to galaxy catalog, this gives set of randoms matching angular footprint
		#else:
		#	no += 1
		#fo.write(str(float(f[i]['RA']))+' '+str(float(f[i]['DEC']))+'\n') #to be used if opened and read from fitsio.FITS
	print no
	fo.close()
	return True

def mkjackf(sample,NS,version,cm='',Njack=20):
	#defines jack-knifes
	mf = open('ranHeal_pix256eboss'+cm+sample+'_'+NS+version+'.dat').readlines()
	fo = open('jackhpixeboss'+cm+sample+'_'+NS+version+str(Njack)+'.dat','w')
	for i in range(0,Njack-1):
		fo.write(str(mf[(len(mf)/Njack)*(i+1)].split()[0])+'\n')
	fo.close()
	return True

def ranHealp(sample,NS,version,cm='',res=256,rad=''):
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
	f = open(dirsci+'reboss'+cm+sample+'_'+NS+version+'_0'+'.dat')	
	for line in f:
		if line[0] != '#' and line[1] != '#' and line.split()[0] != 'az(d)':
			ln = line.split()
			ra,dec = float(ln[0])*angm,float(ln[1])*angm
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			pixl[p] += 1.
	fo = open('ranHeal_pix'+str(res)+'eboss'+cm+sample+'_'+NS+version+'.dat','w')
	for i in range(0,len(pixl)):
		if pixl[i] > 0:
			fo.write(str(i)+' '+str(pixl[i])+'\n')
	fo.close()
	return True

def mkQSO4xi_allweights(NS,version,sample='QSO',cm='',zmin=.6,zmax=1.,c='sci',app='.fits',wm='',compmin=0,gmin=0,gmax=30):
	#note, zmin/zmax assume LRG sample these need to be change for QSO files
	from healpix import healpix, radec2thphi
	from optimize import fmin
	depthm = np.loadtxt('healdepthinm512.dat').transpose() 	
	#b,m = findlinmb(sample,NS,version,'depth',zmin,zmax)
	ds = np.loadtxt('ngebossQSO_S'+version+'_mz0.9xz2.2gm'+str(gmin)+str(gmax)+'512vdepth.dat').transpose()
	dn = np.loadtxt('ngebossQSO_N'+version+'_mz0.9xz2.2gm'+str(gmin)+str(gmax)+'512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	bd,md = fmin(lf.chilin,inl)
	h = healpix()
	slpl = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
	node = 22.45
	be,me = findlinmb(sampl,NS,ver,'ext',zmin,zmax,wm='wdepthgmag')
	wext = np.loadtxt('healSFD_r_'+str(256)+'_fullsky.dat')/2.751
	d = np.loadtxt('ngebossQSO_'+NS+version+'_mz0.9xz2.2fid256vext.dat').transpose()
	dt = d[1]
	e = d[2]	
	lf = linfit(d[0],dt,e)
	inl = np.array([1.,0])
	beo,meo = fmin(lf.chilin,inl)

	if c == 'sci': #AJR uses this define directory for machine he uses
		dir = dirsci
	f = fitsio.read(dir+cm+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat'+app) #read galaxy/quasar file
	no = 0
	gw = ''
	if gmax != 30:
		gw = 'gx'+str(gmax)
	fo = open(dir+'geboss'+cm+sample+'_'+NS+version+'_mz'+str(zmin)+'xz'+str(zmax)+gw+wm+'4xiweights.dat','w')

	for i in range(0,len(f)):
		z = f[i]['Z']
		gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
		if z > zmin and z < zmax and gm < gmax:
			no += 1
			wfkp = f[i]['WEIGHT_FKP']
			wfcp = wfkp*f[i]['WEIGHT_CP']
			wfzf = wfkp*f[i]['WEIGHT_NOZ']
			wfcpzf = wfkp*(f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)
			wfczst = wfkp*(f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SYS']
			wfczsee = wfkp*(f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']
			wfczss = wfkp*(f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_STAR']*f[i]['WEIGHT_SEEING']
			fo.write(str(f[i]['RA'])+' '+str(f[i]['DEC'])+' '+str(z)+' '+str(wfkp)+' '+str(wfcp)+' '+str(wfzf)+' '+str(wfcpzf)+' '+str(wfczst)+' '+str(wfczsee)+' '+str(wfczss)+'\n')
	fo.close()
	print no
	return True

def mkran4xi_allweights(file,NS,N,zmin=.43,zmax=.7,c='sci'):
	from random import random
	dirout = ''
	dir = ''
	if c == 'sci':
		dir = 'DR12/'
		dirout = 'DR12/weightfiles/'

	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	fstr = file.split('-')
	gf = open(dirout+'g'+fstr[0]+NS+fstr[1]+'_'+wz+'4xiweights.dat')
	zl = []
	for line in gf:
		ln = line.split()
		wi =  []
		for i in range(2,len(ln)):
			wi.append(float(ln[i]))
		zl.append(wi)
	fr = open(dir+'r'+fstr[0]+'_'+NS+fstr[1]+'_'+str(N)+'.dat')
	fo = open(dirout+'r'+fstr[0]+'_'+NS+fstr[1]+'_'+str(N)+wz+'4xiweights.dat','w')
	for line in fr:
		ln = line.split()
		indz = int(random()*len(zl))
		fo.write(ln[0]+' '+ln[1]+' ')
		for i in range(0,len(zl[indz])):
			fo.write(str(zl[indz][i])+ ' ')
		fo.write('\n')
		
	fo.close()
	return True


def mkgal4xi(sample,NS,version,cm='',zmin=.6,zmax=1.,c='sci',app='.fits',wm='',compmin=0,gmin=0,gmax=30,zpl=False,gri22='gri22',znudge=False,wmin=False,extc=False,depthc=False,depthextc=False,decmax=False):
	#note, zmin/zmax assume LRG sample these need to be change for QSO files
	from healpix import healpix, radec2thphi
	from optimize import fmin
	from random import random
	if wm == 'wstar' or wm == 'cpstar':
		wsys = np.loadtxt('allstars17.519.9Healpixall256.dat')
		b,m = findlinmb(sample,NS,version,'star',zmin,zmax,dir='')
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
	if wm == 'wdepthext':
		wsys = np.loadtxt('healdepthinm512.dat').transpose()	
		b,m = findlinmb(sample,NS,version,'depth',.8,2.2,wm='nosys'+gri22,dir='')
		be,me = findlinmb(sample,NS,version,'ext',.8,2.2,wm='wdepth'+gri22,dir='')
		wext = np.loadtxt('healSFD_r_'+str(256)+'_fullsky.dat')/2.751
	if wm == 'wgdepthextext':
		b,m = findlinmb(sample,NS,version,'IMAGE_DEPTH_EXT1',.8,2.2,dir='')
		be,me = findlinmb(sample,NS,version,'EB_MINUS_V-1',.8,2.2,wm='wgdepthext',dir='')
		print b,m,be,me
	if wm == 'wdepthgmag':
		#ml = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		ml,bl = fitallQSOvdepthNS(NS,version,sys='depth',wm='nosys')
		node = 22.45
		wsys = np.loadtxt('healdepthinm512.dat').transpose()
		h = healpix()	
	if wm == 'wdepthimag':
		#ml = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		ml,bl = fitallQSOvdepthNSi(NS,version,sys='depth',wm='nosys',gri22=gri22)
		node = 22.45
		wsys = np.loadtxt('healdepthinm512.dat').transpose()
		h = healpix()	

	if wm == 'wdepthimagext':
		#ml = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		ml,bl = fitallQSOvdepthNSi(NS,version,sys='depth',wm='nosys',gri22=gri22)
		node = 22.45
		wsys = np.loadtxt('healdepthinm512.dat').transpose()
		be,me = findlinmb(sample,NS,version,'ext',0.8,2.2,wm='wdepthimag'+gri22,dir='')
		wext = np.loadtxt('healSFD_r_'+str(256)+'_fullsky.dat')/2.751
	
		h = healpix()	

	if wm == 'wdepthgmagext':
		slpl = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		wsys = np.loadtxt('healdepthinm512.dat').transpose()	
		node = 22.45
		be,me = findlinmb(sample,NS,version,'ext',zmin,zmax,wm='wdepthgmag')
		wext = np.loadtxt('healSFD_r_'+str(256)+'_fullsky.dat')/2.751

	if wm == 'wext':
		wsys = np.loadtxt('healSFD_r_256_fullsky.dat')/2.751 #note, file is zipped in directory	
		#b,m = findlinmb(sample,NS,version,'depth',zmin,zmax)
		d = np.loadtxt('ngebossQSO_'+NS+version+'_mz0.9xz2.2fid256vext.dat').transpose()
		dt = d[1]
		e = d[2]	
		lf = linfit(d[0],dt,e)
		inl = np.array([1.,0])
		b,m = fmin(lf.chilin,inl)
		h = healpix()
		
	if c == 'sci': #AJR uses this define directory for machine he uses
		dir = dirsci
	f = fitsio.read(dir+cm+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat'+app) #read galaxy/quasar file
	no = 0
	gw = ''
	if gmax != 30:
		gw = 'gx'+str(gmax)
	dthw = ''
	if depthc:
		dthw += 'depthi'+str(depthc)
		dmap = np.loadtxt('healdepthinm512.dat').transpose()
		h = healpix()	
	if depthextc:
		dthw += 'depthexti'+str(depthextc)
		dmap = np.loadtxt('healdepthinm512.dat').transpose()
		dmape = np.loadtxt('healSFD_r_256_fullsky.dat')/2.751*1.698
		h = healpix()	
	if wmin:
		dthw += 'wmin'+str(wmin)
		dmap = np.loadtxt('healdepthinm512.dat').transpose()
		dmape = np.loadtxt('healSFD_r_256_fullsky.dat')/2.751
		b,m = findlinmb(sample,NS,version,'depth',.8,2.2,wm='nosys'+gri22,dir='')
		be,me = findlinmb(sample,NS,version,'ext',.8,2.2,wm='wdepth'+gri22,dir='')
		
		h = healpix()	

	if extc:
		dthw += 'ext'+str(extc)
		dmape = np.loadtxt('healSFD_r_256_fullsky.dat')/2.751
		h = healpix()
	if decmax:
		dthw += 'decx'+str(decmax)		
	if znudge:
		dthw += 'znudge'	
	if zpl:
		dthw += 'zpl'
	fo = open(dir+'geboss'+cm+sample+'_'+NS+version+'_mz'+str(zmin)+'xz'+str(zmax)+gw+gri22+dthw+wm+'4xi.dat','w')
	for i in range(0,len(f)):
		z = f[i]['Z']
		if zpl:
			if z == f[i]['Z_VI']:
				z = f[i]['Z_PL']
		if znudge:
			if len(str(z)) <= 5:
				z += .001*random()-.0005
		comp = f[i]['COMP']
		gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
		c = 1
		if gri22 == 'gri22':
			rm = f[i]['MODELMAG'][2]-f[i]['EXTINCTION'][2]
			im = f[i]['MODELMAG'][3]-f[i]['EXTINCTION'][3]
			if gm > 22 or rm > 22 or im > 22:
				c = 0
		if depthc:
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(512,th,phi)
			if dmap[pix2] < depthc:
				c = 0		
		if depthextc:
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(512,th,phi)
			pixe = h.ang2pix_nest(256,th,phi)
			if dmap[pix2]-dmape[pixe] < depthextc:
				c = 0		
		if wmin:
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(512,th,phi)
			pixe = h.ang2pix_nest(256,th,phi)
			ns = dmap[pix2]
			ws = (b+m*ns)
			ext = dmape[pixe]
			we = (be+me*ext)

			if ws*we < wmin:
				c = 0		

		if extc:
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(256,th,phi)
			if dmape[pix2] > extc:
				c = 0		
		if decmax:
			ra,dec = f[i]['RA'],f[i]['DEC']
			if dec > decmax:
				c = 0
		if z > zmin and z < zmax and comp > compmin and gm < gmax and c == 1:
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
			if wm == 'wstar' or wm == 'cpstar' or wm == 'wext':
				pix2 = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws
			if wm == 'wgdepthextext':
				sysw = f[i]['IMAGE_DEPTH'][1]
				sysw = luptm(sysw,1)-3.303*f[i]['EB_MINUS_V']
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*(1./(b+m*sysw))
				ext = f[i]['EB_MINUS_V']
				w = w*(1./(be+me*ext))

			if wm == 'wdepth' or wm == 'cpdepth':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				
			if wm == 'wdepthext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				pixe = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				ext = wext[pixe]
				we = 1./(be+me*ext)
				w = w*ws*we				
			if wm == 'wdepthgmag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if gm < 20.75:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 21.25*slp
					m = gm*slp+b
				if gm >= 21.75:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				
			if wm == 'wdepthimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if im < 20.25:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 19.75*slp
					m = im*slp+b
				if im >= 20.25 and im < 20.75:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.25*slp
					m = im*slp+b
				if im >= 20.75 and im < 21.25:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 20.75*slp
					m = im*slp+b
				if im >= 21.25:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.25*slp
					m = im*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				
			if wm == 'wdepthimagext':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if im < 20.25:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 19.75*slp
					m = im*slp+b
				if im >= 20.25 and im < 20.75:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.25*slp
					m = im*slp+b
				if im >= 20.75 and im < 21.25:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 20.75*slp
					m = im*slp+b
				if im >= 21.25:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.25*slp
					m = im*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				pix2 = h.ang2pix_nest(256,th,phi)
				ne = wext[pix2]
				we = 1./(me*ne+be)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws*we				
			if wm == 'wdepthgmagext':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if gm < 20.75:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 21.25*slp
					m = gm*slp+b
				if gm >= 21.75:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				pix2 = h.ang2pix_nest(256,th,phi)
				ne = wext[pix2]
				we = 1./(me*ne+be)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws*we				
			if wm == 'nw':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']	
			fo.write(str(f[i]['RA'])+' '+str(f[i]['DEC'])+' '+str(z)+' '+str(w)+'\n')
	fo.close()
	print no
	return True

def mkran4xi(sample,NS,version,cm='',N=0,wm='',zmin=.6,zmax=.1,comp = 'sci',gmax=30,zpl=False,gri22='',znudge=False,wmin=False,depthc=False,depthextc=False,extc=False,decmax=False):
	from random import random
	from healpix import healpix, radec2thphi
	if comp == 'sci':
		dir = dirsci 
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gw = ''
	print depthextc
	if gmax != 30:
		gw = 'gx'+str(gmax)
	if depthc:
		gri22 += 'depthi'+str(depthc)
		dmap = np.loadtxt('healdepthinm512.dat').transpose()	
		h = healpix()
	if depthextc:
		gri22 += 'depthexti'+str(depthextc)
		dmap = np.loadtxt('healdepthinm512.dat').transpose()
		dmape = np.loadtxt('healSFD_r_256_fullsky.dat')/2.751*1.698
		h = healpix()	
	if wmin:		
		dmap = np.loadtxt('healdepthinm512.dat').transpose()
		dmape = np.loadtxt('healSFD_r_256_fullsky.dat')/2.751
		b,m = findlinmb(sample,NS,version,'depth',.8,2.2,wm='nosys'+gri22,dir='')
		be,me = findlinmb(sample,NS,version,'ext',.8,2.2,wm='wdepth'+gri22,dir='')
		gri22 += 'wmin'+str(wmin)
		h = healpix()	

	if extc:
		gri22 += 'ext'+str(extc)
		dmape = np.loadtxt('healSFD_r_256_fullsky.dat')/2.751
		h = healpix()	
	if decmax:
		gri22 += 'decx'+str(decmax)
	if znudge:
		gri22 += 'znudge'
	if zpl:
		gri22 += 'zpl'		
	gf = np.loadtxt(dir+'geboss'+cm+sample+'_'+NS+version+'_'+wz+gw+gri22+wm+'4xi.dat').transpose()
	fr = np.loadtxt(dir+'reboss'+cm+sample+'_'+NS+version+'_'+str(N)+'.dat').transpose()
	fo = open(dir+'reboss'+cm+sample+'_'+NS+version+'_'+str(N)+wz+gw+gri22+wm+'4xi.dat','w')
	n = 0
	for i in range(0,len(fr[0])):
		c = 1
		if depthc:
			ra,dec = fr[0][i],fr[1][i]
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(512,th,phi)
			if dmap[pix2] < depthc:
				c = 0		
		if depthextc:
			ra,dec = fr[0][i],fr[1][i]
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(512,th,phi)
			pixe = h.ang2pix_nest(256,th,phi)
			if dmap[pix2]-dmape[pixe] < depthextc:
				c = 0		
		if wmin:
			ra,dec = fr[0][i],fr[1][i]
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(512,th,phi)
			pixe = h.ang2pix_nest(256,th,phi)
			ns = dmap[pix2]
			ws = (b+m*ns)
			ext = dmape[pixe]
			we = (be+me*ext)

			if ws*we < wmin:
				c = 0		

		if extc:
			ra,dec = fr[0][i],fr[1][i]
			th,phi = radec2thphi(ra,dec)
			pix2 = h.ang2pix_nest(256,th,phi)
			if dmape[pix2] > extc:
				c = 0		
		if decmax:
			ra,dec = fr[0][i],fr[1][i]
			if dec > 10:
				c = 0
		if c == 1:
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

def createalladfilesfb(sample,NS,version,cm='',nran=1,wm='',zmin=.6,zmax=1.,gmax=30,gri22='',zpl=False,znudge=False,wmin=False,depthc=False,depthextc=False,extc=False,decmax=False):
	#after defining jack-knifes, this makes all of the divided files and the job submission scripts
	#./suball.sh sends all of the jobs to the queue on the system I use
	mkgal4xi(sample,NS,version,cm=cm,wm=wm,zmin=zmin,zmax=zmax,gmax=gmax,gri22=gri22,zpl=zpl,znudge=znudge,wmin=wmin,depthc=depthc,depthextc=depthextc,extc=extc,decmax=decmax)
	gw = ''
	if gmax != 30:
		gw = 'gx'+str(gmax)

	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	sysw = gw
	sysw += gri22
	if depthc:
		sysw += 'depthi'+str(depthc)
	if depthextc:
		sysw += 'depthexti'+str(depthextc)
	if wmin:
		sysw += 'wmin'+str(wmin)	
	if extc:
		sysw += 'ext'+str(extc)
	if decmax:
		sysw += 'decx'+str(decmax)	
	if znudge:
		sysw += 'znudge'	
	if zpl:
		sysw += 'zpl'	
	sysw += wm
	gf = 'geboss'+cm+sample+'_'+NS+version+'_'+wz+sysw#gw+gri22+wm
	createSourcesrd_ad(gf)
	for i in range(0,20):
		createSourcesrd_adJack(gf,i,'eboss'+cm+sample+'_'+NS+version)

	rf = 'reboss'+cm+sample+'_'+NS+version+'_'
	for rann in range(0,nran):	
		print rann
		mkran4xi(sample,NS,version,cm=cm,N=rann,wm=wm,zmin=zmin,zmax=zmax,gmax=gmax,zpl=zpl,gri22=gri22,znudge=znudge,wmin=wmin,depthc=depthc,depthextc=depthextc,extc=extc,decmax=decmax)
		#mkran4xifit(sample,NS,version,N=rann,wm=wm,zmin=zmin,zmax=zmax)
		rfi = rf+str(rann)+wz+sysw#gw+gri22+wm
		createSourcesrd_ad(rfi)
		for i in range(0,20):
			createSourcesrd_adJack(rfi,i,'eboss'+cm+sample+'_'+NS+version)
	mksuball_nran_Dmufbfjack(rf,gf,nran,wr=wz+sysw)#gw+gri22+wm)
	return True

def ppxilcalc_LSD_bin(file,mom,NS='ngc',v='v1.6',bs=8,start=0,rmax=250,nranf=1,njack=20,fa='',md='EZmock_QSO',zmin=.8,zmax=2.2,pp=False):
	#mom get multiplied by two, so mom=1 is quadrupole
	mdr = md
	from numpy import zeros
	from time import time
	DDnl = []	
	DDnorml = 0
	DDnormt = 0
	DRnl = []
	DRnorml = 0
	DRnormt = 0
	RRnl = []
	RRnl0 = []
	nbin = rmax/bs
	nmubin = 100
	for i in range(0,rmax*nmubin):
		DDnl.append(0)
		DRnl.append(0)
		RRnl.append(0)
		RRnl0.append(0)
	RRnorml = 0
	RRnormt = 0
	pl = []
	print file,zmin,zmax,NS,v
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	rf0 = 'r'+mdr+'_'+NS+v+'_'
	gfr = 'g'+md+file+NS+v+wz
	for i in range(0,nmubin):
		mu = i/float(nmubin)+.5/float(nmubin)
		mub = int(mu*nmubin)
		#print mu
		pl.append((1.,P2(mu),P4(mu),P6(mu),P8(mu),mub))
	#print len(pl)
	fdp = open(gfr+gfr+'2ptdmu.dat').readlines()
	DDnormt += float(fdp[0])
	for k in range(1,len(fdp)):
		dp = float(fdp[k])
		DDnl[k-1] += dp

	for N in range(0,nranf):
#		if N != 2:
		#try:
		#	fdnp = open('r'+md+'_'+NS+'_'+str(N)+rb+file+'2ptdmu.dat').readlines()
		#except:
		fdnp = open(gfr+rf0+str(N)+wz+'2ptdmu.dat').readlines()
		DRnormt += float(fdnp[0])
		for k in range(1,len(fdp)):
			dr = float(fdnp[k])
			DRnl[k-1] += dr
		fr = open(rf0+str(N)+wz+'_all2ptdmu.dat').readlines()
		RRnormt += float(fr[0])
		for k in range(1,len(fdp)):
			rbin = int((k-1)/100.)
			rp = float(fr[k])
			RRnl[k-1] += rp

				
					
	xil = zeros((nbin),'f')
	#for i in range(0,nbin):
	#print DDnormt,DRnormt,RRnormt
	r = start+bs/2.
	for i in range(start,rmax,bs):
		xi = 0
		ddt = 0
		drt = 0
		rrt = 0
		for j in range(0,nmubin):
			rr0 = 0
			dd = 0
			dr = 0
			rr = 0

			for k in range(0,bs):
				bin = nmubin*(i+k)+j			
				if bin < len(RRnl):
					#if RRnl[bin] == 0:
					#	pass
						
					#else:
					dd += DDnl[bin]#*pl[bin][mom]
					rr += RRnl[bin]#*pl[bin][mom]
					dr += DRnl[bin]#*pl[bin][mom]#/RRnl[bin]
			ddt += dd
			rrt += rr
			drt += dr
			if rr != 0:		
				xi += pl[j][mom]/float(nmubin)*(dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr
		#print ddt/DDnormt,drt/DRnormt,rrt/RRnormt
		#if rr0 > 0:
		#	xi = (dd/DDnormt-2.*dr/DRnormt+rr/RRnormt)*RRnormt/rr0
		#else:
		#	xi = 0
		if i/bs < nbin:
			xil[i/bs] = xi
		if pp == True:
			#r = start+bs/2.+i*bs
			print r,sqrt(1./ddt)*r	
		r += bs	
	return xil


def ppxilcalc_LSDfjack_bs(sample,NS,version,jack,mom,zmin=.6,zmax=1.,wm='',bs=5,start=0,rmax=250,mumin=0,mumax=1.,nranf=1,njack=20,wf=False,wmu = ''):
	#finds xi, for no jack-knife, set jack = -1, otherwise can be used to calculate jack-knife xi for 0 <= jack < Njack
	from numpy import zeros
	from time import time
	#DDnl = zeros((nbin,njack),'f')
	rf = 'reboss'+sample+'_'+NS+version+'_'
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gf = 'geboss'+sample+'_'+NS+version+'_'+wz+wm
	fl = sample+'_'+NS+version+'_'+wz+wm
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
	if wf:
		fD = open(dirsci+'paircounts/Paircounts_'+fl+'.dat','w')
		#fDR = open('DRcounts'+file+'.dat','w')
		#fR = open('RRcounts'+file+'.dat','w')
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
	
	if wf:
		fD.write('#normalized paicounts; 100 dmu = 0.01 bins, 250 dr = 1mpc/h bins; every new line increases in r; every 100 lines increases in mu\n')
		fD.write('#r_center mu_center DD DR RR\n')
		#fD.write(str(DDnormt)+'\n')
		#fDR.write(str(DRnormt)+'\n')
		#fR.write(str(RRnormt)+'\n')
		for j in range(0,100):
			for i in range(0,rmax):
				fD.write(str(.5+i)+' '+str(.01*j+.005)+' '+str(DDnl[j+100*i]/DDnormt)+' '+str(DRnl[j+100*i]/DRnormt)+' '+str(RRnl[j+100*i]/RRnormt)+'\n')
				#fDR.write(str(DRnl[j+100*i])+' ')
				#fR.write(str(RRnl[j+100*i])+' ')
			#fD.write('\n')
			#fDR.write('\n')
			#fR.write('\n')
		fD.close()
		#fDR.close()
		#fR.close()
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


def ppxilfile_bs(sample,NS,version,mom,zmin=.6,zmax=1.,wm='',bs=5,start=0,rmax=250,nranf=1,njack=20,mumin=0,mumax=1.,wmu='',wf=False):
	#write out xi to a file, no jack-knife errors
	#for quadrupole, set mom = 1, for hexadecapole set mom = 2, ect. (odd moments not supported)
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gf = 'geboss'+sample+'_'+NS+version+'_'+wz+wm
	ans = ppxilcalc_LSDfjack_bs(sample,NS,version,-1,mom,zmin=zmin,zmax=zmax,wm=wm,mumin=mumin,mumax=mumax,bs=bs,start=start,rmax=rmax,nranf=nranf,wmu=wmu,wf=wf)
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

def doallQSOvdepth(v,wm='nosys',gri22='gri22'):
	#magl = [0,20.5,21.,21.5,22.,30]
	magl = [0,20.,20.5,21.,21.5,22.]
	for i in range(0,len(magl)-1):
		gmag = (magl[i],magl[i+1])
		ngvsys('QSO','N',v,'depth',22,23,512,.8,2.2,wm,gmag,gri22=gri22)
		ngvsys('QSO','S',v,'depth',22,23,512,.8,2.2,wm,gmag,gri22=gri22)
	return True

def doallQSOvdepthr(v,wm='nosys',gri22='gri22'):
	magl = [0,20.,20.5,21.,21.5,22.]
	for i in range(0,len(magl)-1):
		gmag = (magl[i],magl[i+1])
		ngvsys('QSO','N',v,'depth',22,23,512,.8,2.2,wm,rmag=gmag,gri22=gri22)
		ngvsys('QSO','S',v,'depth',22,23,512,.8,2.2,wm,rmag=gmag,gri22=gri22)
	return True

def doallQSOvdepthi(v,wm='nosys',gri22='gri22'):
	magl = [0,20.,20.5,21.,21.5,22.,30]
	for i in range(0,len(magl)-1):
		gmag = (magl[i],magl[i+1])
		ngvsys('QSO','N',v,'depth',22,23,512,.8,2.2,wm,imag=gmag,gri22=gri22)
		ngvsys('QSO','S',v,'depth',22,23,512,.8,2.2,wm,imag=gmag,gri22=gri22)
	return True


def fitallQSOvdepth(v,sys='depth',wm='nosys'):
	mnl = []
	bnl = []
	msl = []
	bsl = []
	magl = [0,20.5,21.,21.5,22.,30]
	for i in range(0,len(magl)-1):
		gmw = 'gm'+str(magl[i])+str(magl[i+1])
		bn,mn = findlinmb('QSO','N',v,sys,.8,2.2,wm+gmw)
		mnl.append(mn)
		bnl.append(bn)
		bs,ms = findlinmb('QSO','S',v,sys,.8,2.2,wm+gmw)
		msl.append(ms)
		bsl.append(bs)
	return mnl,bnl,msl,bsl	

def fitallQSOvdepthNS(NS,v,sys='depth',wm='nosys',dir=''):
	mnl = []
	bnl = []
	magl = [0,20.5,21.,21.5,22.,30]
	for i in range(0,len(magl)-1):
		gmw = 'gm'+str(magl[i])+str(magl[i+1])
		bn,mn = findlinmb('QSO',NS,v,sys,.8,2.2,wm+gmw,dir=dir)
		mnl.append(mn)
		bnl.append(bn)
	return mnl,bnl	

def fitallQSOvdepthNSi(NS,v,bnd='i',sys='depth',wm='nosys',gri22='',dir=''):
	mnl = []
	bnl = []
	magl = [0,20.,20.5,21.,21.5,22.]
	for i in range(0,len(magl)-1):
		gmw = bnd+'m'+str(magl[i])+str(magl[i+1])
		bn,mn = findlinmb('QSO',NS,v,sys,.8,2.2,wm+gmw+gri22,dir=dir)
		mnl.append(mn)
		bnl.append(bn)
	return mnl,bnl	

def fitallQSOvdepthNSu(NS,v,sys='depth',wm='nosys',gri22='',dir=''):
	mnl = []
	bnl = []
	magl = [0,20.5,21,21.5,30]
	for i in range(0,len(magl)-1):
		gmw = 'um'+str(magl[i])+str(magl[i+1])
		bn,mn = findlinmb('QSO',NS,v,sys,.8,2.2,wm+gmw+gri22,dir=dir)
		mnl.append(mn)
		bnl.append(bn)
	return mnl,bnl	


def fitallQSOvdepthNSgri(NS,v,sys='depth',wm='nosys',gri22='',dir=''):
	mnl = []
	bnl = []
	magl = [0,60,61.5,63,64.5,66]
	for i in range(0,len(magl)-1):
		gmw = 'gri'+str(magl[i])+str(magl[i+1])
		bn,mn = findlinmb('QSO',NS,v,sys,.8,2.2,wm+gmw+gri22,dir=dir)
		mnl.append(mn)
		bnl.append(bn)
	return mnl,bnl	

def mksysmapran(res):
	#make healpix systematic maps from eBOSSrandoms.ran.fits, in nested format
	from healpix import healpix, radec2thphi
	h = healpix()
	nc = 1+5+1+1+5+5 #one column for nran, 5 for skyflux, 1 for aimass, 1 for EB-V, 5 for depth, 5 for PSF
	nrow = 12*res*res #number of healpix pixels in full sky
	d = np.zeros((nrow,nc))
	f = fitsio.read(dirfits+'eBOSSrandoms.ran.fits')
	for i in range(0,len(f)):
		ra,dec = f[i]['RA'],f[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = h.ang2pix_nest(res,th,phi)
		d[p][0] += 1.
		for j in range(0,5):
			d[p][1+j] += f[i]['SKYFLUX'][j]
		d[p][6] += 	f[i]['AIRMASS']
		d[p][7] += 	f[i]['EB_MINUS_V']
		for j in range(0,5):
			d[p][8+j] += f[i]['IMAGE_DEPTH'][j]
		for j in range(0,5):
			d[p][13+j] += f[i]['PSF_FWHM'][j]
	print 'filled data table'
	fo = open(ebossdir+'healmap_eBOSS_nested_'+str(res)+'.dat','w')
	fo.write('# pix SKYFLUX_u SKYFLUX_g SKYFLUX_r SKYFLUX_i SKYFLUX_z AIRMASS EB_MINUS_V IMAGE_DEPTH_u IMAGE_DEPTH_g IMAGE_DEPTH_r IMAGE_DEPTH_i IMAGE_DEPTH_z PSF_FWHM_u PSF_FWHM_g PSF_FWHM_r PSF_FWHM_i PSF_FWHM_z\n')
	for i in range(0,nrow):
		if d[i][0] != 0:
			fo.write(str(i)+' ')
			for j in range(1,nc):
				fo.write(str(d[i][j]/d[i][0])+' ')
			fo.write('\n')
	fo.close()
	return True
		
		

def ngvsys(sampl,NS,ver,sys,sysmin,sysmax,res,zmin,zmax,wm='',umag=False,gmag=False,rmag=False,imag=False,umg=False,gri=False,gri22=''):
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
	if wm == 'wgdepthextext':
		b,m = findlinmb(sampl,NS,ver,'IMAGE_DEPTH_EXT1',zmin,zmax)
		be,me = findlinmb(sampl,NS,ver,'EB_MINUS_V-1',zmin,zmax,wm='wgdepthext')

	if wm == 'wstar':
		wsys = np.loadtxt(dirsys+'allstars17.519.9Healpixall256.dat')
		b,m = findlinmb(sampl,NS,ver,'star',zmin,zmax)
	if wm == 'wdepth':
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()	
		b,m = findlinmb(sampl,NS,ver,'depth',.8,2.2,wm='nosys'+gri22,dir=ebossdir)
	if wm == 'wdepthext':
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()	
		b,m = findlinmb(sampl,NS,ver,'depth',.8,2.2,wm='nosys'+gri22,dir=ebossdir)
		be,me = findlinmb(sampl,NS,ver,'ext',.8,2.2,wm='wdepth'+gri22,dir=ebossdir)
		wext = np.loadtxt(dirsys+'healSFD_r_'+str(256)+'_fullsky.dat')/2.751
	if wm == 'wdepthgmagext':
		slpl = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()	
		node = 22.45
		be,me = findlinmb(sampl,NS,ver,'ext',zmin,zmax,wm='wdepthgmag')
		wext = np.loadtxt(dirsys+'healSFD_r_'+str(256)+'_fullsky.dat')/2.751

	if wm == 'wdepthgmag':
		slpl = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		bl = [-1.1444474848040611, -3.2505181700338035, -5.8979701872596051, -15.205991475527554, -29.455524298695721]
		node = 22.45
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()	

	if wm == 'wdepthimag':
		#ml = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		ml,bl = fitallQSOvdepthNSi(NS,ver,sys='depth',wm='nosys',gri22=gri22,dir=ebossdir)
		node = 22.45
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()
		h = healpix()	
	if wm == 'wdepthumag':
		#ml = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		ml,bl = fitallQSOvdepthNSu(NS,ver,sys='depth',wm='nosys',gri22=gri22,dir=ebossdir)
		node = 22.45
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()
		h = healpix()	
	if wm == 'wdepthgrimag':
		#ml = [0.09558874587471361, 0.18952081376668517, 0.30767488133801107, 0.72263365380178768, 1.3575184365813411]
		ml,bl = fitallQSOvdepthNSgri(NS,ver,sys='depth',wm='nosys',gri22=gri22,dir=ebossdir)
		node = 22.45
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()
		h = healpix()	
	
	npo = 12*res**2
		
	if sys != 'star' and sys != 'ext' and sys != 'depth':
		#fsys = np.loadtxt(dirsys+'heal'+sys+'nm'+str(res)+'.dat')
		filesys = np.loadtxt(ebossdir+'healmap_eBOSS_nested_'+str(res)+'.dat').transpose()
		hd = open(ebossdir+'healmap_eBOSS_nested_'+str(res)+'.dat').readline().split()
		for i in range(0,len(hd)):
			if hd[i] == sys:
				col = i-1 #because of space between # and first text in file
				break
		print col
		fsys = np.zeros((npo))
		for i in range(0,len(filesys[0])):
			p = filesys[0][i]
			if sys.split('_')[0] == 'IMAGE':
				if sys.split('_')[-1] == 'i':
					fsys[p] = luptm(filesys[col][i],3)
				if sys.split('_')[-1] == 'g':
					fsys[p] = luptm(filesys[col][i],1)
				if sys.split('_')[-1] == 'r':
					fsys[p] = luptm(filesys[col][i],2)
				if sys.split('_')[-1] == 'z':
					fsys[p] = luptm(filesys[col][i],4)
				if sys.split('_')[-1] == 'u':
					fsys[p] = luptm(filesys[col][i],0)

			else:
				fsys[p] = filesys[col][i]
		
	if sys == 'ext':
		fsys = np.loadtxt(dirsys+'healSFD_r_'+str(res)+'_fullsky.dat')/2.751 #E(B-V)#*2.285/2.751 #r-band extinction
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

	#ml = []
	

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
		um = f[i]['MODELMAG'][0]-f[i]['EXTINCTION'][0]
		gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
		rm = f[i]['MODELMAG'][2]-f[i]['EXTINCTION'][2]
		im = f[i]['MODELMAG'][3]-f[i]['EXTINCTION'][3]
		c = 1
		if gri22 == 'gri22':
			if gm > 22 or rm > 22 or im > 22:
				c = 0
		if umg != False:
			um = f[i]['MODELMAG'][0]-f[i]['EXTINCTION'][0]
			if um-gm < umg[0] or um-gm > umg[1]:
				gc = False

		if gri != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if gm+rm+im < gri[0] or gm+rm+im > gri[1]:
				gc = False

		if umag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if um < umag[0] or um > umag[1]:
				gc = False

		if gmag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if gm < gmag[0] or gm > gmag[1]:
				gc = False
		if rmag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if rm < rmag[0] or rm > rmag[1]:
				gc = False
		if imag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if im < imag[0] or im > imag[1]:
				gc = False
		if z > zmin and z < zmax and c == 1 and gc:
			no += 1
			#w = 1.
			#if wm == '':
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT'] #weight in current catalog
			if wm == 'nosys':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP'] #compare to result with no systematic weights
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
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws				
			if wm == 'wdepthext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				pixe = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				ext = wext[pixe]
				we = 1./(be+me*ext)
				w = w*ws*we				
			if wm == 'wdepthgmag':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = w*ws
			if wm == 'wdepthimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if im < 20.25:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 19.75*slp
					m = im*slp+b
				if im >= 20.25 and im < 20.75:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.25*slp
					m = im*slp+b
				if im >= 20.75 and im < 21.25:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 20.75*slp
					m = im*slp+b
				if im >= 21.25:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.25*slp
					m = im*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				
			if wm == 'wdepthumag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if um < 20.5:
					m = ml[0]
					b = bl[0]
				if um >= 20.5 and um < 21:
					m = ml[1]
					b = bl[1]
				if um >= 21 and um < 21.5:
					m = ml[2]
					b = bl[2]
				if um >= 21.5:
					m = ml[3]
					b = bl[3]
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				
			if wm == 'wdepthgrimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				grim = gm+rm+im
# 				if grim < 60.75:
# 					slp = (ml[1]-ml[0])/.5
# 					b = ml[0] - 59.25*slp
# 					m = grim*slp+b
# 				if grim >= 60.75 and grim < 62.25:
# 					slp = (ml[2]-ml[1])/.5
# 					b = ml[1] - 60.75*slp
# 					m = grim*slp+b
# 				if grim >= 62.25 and grim < 63.75:
# 					slp = (ml[3]-ml[2])/.5
# 					b = ml[2] - 62.25*slp
# 					m = grim*slp+b
# 				if grim >= 63.75:
# 					slp = (ml[4]-ml[3])/.5
# 					b = ml[3] - 63.75*slp
# 					m = grim*slp+b
				if grim < 60.:
					m = ml[0]
					b = bl[0]
				if grim >= 60 and grim < 61.5:
					m = ml[1]
					b = bl[1]
				if grim >= 61.5 and grim < 63:
					m = ml[2]
					b = bl[2]
				if grim >= 63 and grim < 64.5:
					m = ml[3]
					b = bl[3]
				if grim >= 64.5:
					m = ml[4]
					b = bl[4]
					
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				

			if wm == 'wdepthgmagext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				pix2 = h.ang2pix_nest(256,th,phi)
				ne = wext[pix2]
				we = 1./(be+me*ne)
				w = w*ws*we
			if wm == 'wgdepthextext':
				sysw = f[i]['IMAGE_DEPTH'][1]
				sysw = luptm(sysw,1)-3.303*f[i]['EB_MINUS_V']
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*(1./(b+m*sysw))
				ext = f[i]['EB_MINUS_V']
				w = w*(1./(be+me*ext))

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
	if umg != False:
		wm += 'umg'+str(umg[0])+str(umg[1])
	if gri != False:
		wm += 'gri'+str(gri[0])+str(gri[1])
	if umag != False:
		wm += 'um'+str(umag[0])+str(umag[1])
	
	if gmag != False:
		wm += 'gm'+str(gmag[0])+str(gmag[1])
	if rmag != False:
		wm += 'rm'+str(rmag[0])+str(rmag[1])
	if imag != False:
		wm += 'im'+str(imag[0])+str(imag[1])
	fs = open(ebossdir+'n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+gri22+str(res)+'v'+sys+'.dat','w')
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
	fs.close()
	plotvssys(sampl,NS,ver,sys,sysmin,sysmax,res,zmin,zmax,wm+gri22)
	return xl,yl,el

def ngvsys_ran(sampl,NS,ver,sys,sysmin,sysmax,zmin,zmax,band=-1,wm='',umag=False,gmag=False,rmag=False,imag=False,umg=False,gri=False,gri22=''):
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
	if wm == 'wstar':
		wsys = np.loadtxt(dirsys+'allstars17.519.9Healpixall256.dat')
		b,m = findlinmb(sampl,NS,ver,'star',zmin,zmax)
	if wm == 'wgdepthext':
		b,m = findlinmb(sampl,NS,ver,'IMAGE_DEPTH_EXT1',zmin,zmax)
	if wm == 'wgdepthextext':
		b,m = findlinmb(sampl,NS,ver,'IMAGE_DEPTH_EXT1',.8,2.2)
		be,me = findlinmb(sampl,NS,ver,'EB_MINUS_V-1',.8,2.2,wm='wgdepthext')
		print b,m
		print be,me
	if wm == 'wgidepthext':
		b,m = findlinmb(sampl,NS,ver,'IMAGE_DEPTH_EXT1',zmin,zmax)
		bi,mi = findlinmb(sampl,NS,ver,'IMAGE_DEPTH_EXT3',zmin,zmax,wm='wgdepthext')
	binng = []
	binnr = []
	nsysbin = 10 #10 bins are being used
	for i in range(0,nsysbin):
		binng.append(0)
		binnr.append(0)

	
	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.ran.fits') #read galaxy/quasar file
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

	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.dat.fits') #read galaxy/quasar file
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
		c = 1
		if gri22 == 'gri22':
			if gm > 22 or rm > 22 or im > 22:
				c = 0
		if umg != False:
			um = f[i]['MODELMAG'][0]-f[i]['EXTINCTION'][0]
			if um-gm < umg[0] or um-gm > umg[1]:
				gc = False

		if gri != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if gm+rm+im < gri[0] or gm+rm+im > gri[1]:
				gc = False

		if umag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if um < umag[0] or um > umag[1]:
				gc = False

		if gmag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if gm < gmag[0] or gm > gmag[1]:
				gc = False
		if rmag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if rm < rmag[0] or rm > rmag[1]:
				gc = False
		if imag != False:
			#gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			if im < imag[0] or im > imag[1]:
				gc = False
		if z > zmin and z < zmax and c == 1 and gc:
			no += 1
			#w = 1.
			#if wm == '':
			
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT'] #weight in current catalog
			if wm == 'nosys':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP'] #compare to result with no systematic weights
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
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws				
			if wm == 'wdepthext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				pixe = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				ext = wext[pixe]
				we = 1./(be+me*ext)
				w = w*ws*we				
			if wm == 'wdepthgmag':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = w*ws
			if wm == 'wdepthimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if im < 20.25:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 19.75*slp
					m = im*slp+b
				if im >= 20.25 and im < 20.75:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.25*slp
					m = im*slp+b
				if im >= 20.75 and im < 21.25:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 20.75*slp
					m = im*slp+b
				if im >= 21.25:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.25*slp
					m = im*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				
			if wm == 'wdepthumag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if um < 20.5:
					m = ml[0]
					b = bl[0]
				if um >= 20.5 and um < 21:
					m = ml[1]
					b = bl[1]
				if um >= 21 and um < 21.5:
					m = ml[2]
					b = bl[2]
				if um >= 21.5:
					m = ml[3]
					b = bl[3]
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				
			if wm == 'wdepthgrimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				grim = gm+rm+im
# 				if grim < 60.75:
# 					slp = (ml[1]-ml[0])/.5
# 					b = ml[0] - 59.25*slp
# 					m = grim*slp+b
# 				if grim >= 60.75 and grim < 62.25:
# 					slp = (ml[2]-ml[1])/.5
# 					b = ml[1] - 60.75*slp
# 					m = grim*slp+b
# 				if grim >= 62.25 and grim < 63.75:
# 					slp = (ml[3]-ml[2])/.5
# 					b = ml[2] - 62.25*slp
# 					m = grim*slp+b
# 				if grim >= 63.75:
# 					slp = (ml[4]-ml[3])/.5
# 					b = ml[3] - 63.75*slp
# 					m = grim*slp+b
				if grim < 60.:
					m = ml[0]
					b = bl[0]
				if grim >= 60 and grim < 61.5:
					m = ml[1]
					b = bl[1]
				if grim >= 61.5 and grim < 63:
					m = ml[2]
					b = bl[2]
				if grim >= 63 and grim < 64.5:
					m = ml[3]
					b = bl[3]
				if grim >= 64.5:
					m = ml[4]
					b = bl[4]
					
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				

			if wm == 'wdepthgmagext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				pix2 = h.ang2pix_nest(256,th,phi)
				ne = wext[pix2]
				we = 1./(be+me*ne)
				w = w*ws*we
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
			if wm == 'wgidepthext':
				sysw = f[i]['IMAGE_DEPTH'][1]
				sysw = luptm(sysw,1)-extc[band]*f[1]['EB_MINUS_V']
				syswi = f[i]['IMAGE_DEPTH'][3]
				syswi = luptm(syswi,3)-extc[band]*f[3]['EB_MINUS_V']
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*(1./(b+m*sysw))
				w = w*(1./(bi+mi*syswi))
					
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
	if umg != False:
		wm += 'umg'+str(umg[0])+str(umg[1])
	if gri != False:
		wm += 'gri'+str(gri[0])+str(gri[1])
	if umag != False:
		wm += 'um'+str(umag[0])+str(umag[1])
	
	if gmag != False:
		wm += 'gm'+str(gmag[0])+str(gmag[1])
	if rmag != False:
		wm += 'rm'+str(rmag[0])+str(rmag[1])
	if imag != False:
		wm += 'im'+str(imag[0])+str(imag[1])
	fs = open(ebossdir+'n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+gri22+'v'+sys+str(band)+'.dat','w')
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
	print chin
	xl = np.array(xl)
	yl = np.array(yl)
	el = np.array(el)
	plotvssys_simp(xl,yl,el,sys)
	return xl,yl,el

def doallnsys4qsoplot():
	nsl = ['N','S']
	wml = ['nosys','wgdepthextext']
	sysl = ['SKYFLUX','AIRMASS','EB_MINUS_V','IMAGE_DEPTH_EXT','PSF_FWHM']
	bandl = [3,-1,-1,3,3]
	sysrl = [(4,20),(1,2),(0.01,0.15),(21.8,22.9),(.7,2)]
	#for i in range(0,len(sysl)):
	for i in range(3,4):	
		for ns in nsl:
			for wm in wml:
				print sysl[i],ns,wm
				ngvsys_ran('QSOsys',ns,'v1.6',sysl[i],sysrl[i][0],sysrl[i][1],.8,2.2,band=bandl[i],wm=wm)
	#for ns in nsl:
	#	for wm in wml:
	#		ngvsys('QSOsys',ns,'v1.6','star',30,300,256,.8,2.2,wm=wm)			
	return True


def findlinmb(sampl,NS,ver,sys,zmin,zmax,wm='nosys',res='512',dir=ebossdir):
	#finds linear fit parameters (depth or stellar density relationships hard-coded to expect given resolutions)
	from optimize import fmin
	res = ''
	if sys == 'star' or sys == 'ext':
		res = '256'
	if sys == 'depth':
		res = '512'
	d = np.loadtxt(dir+'n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+str(res)+'v'+sys+'.dat').transpose()	
	lf = linfit(d[0],d[1],d[2])
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	print lf.chilin((b0,m0))
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
	#sl = []
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
		C1 = 1.
		pk1 = C1**2.*pk1
		dk = ldk*k
		isigk = 0
		for i in range(0,len(nl)):
			isigk += vl[i]/(pk1**2.+2.*pk1/nl[i]) #sum inverse variance for each redshift shell
			#if r1==r2:
			#	isigk += nl[i]**2.*vl[i]*2.*pi*r1**2.*dr
		sigk = dk/(r1*r2)/isigk*(sin(k*r1)*sin(k*r2))#*k*k
		sumxi += sigk*8. #factor of 4 makes BOSS results ~correct
	sumxi = sumxi/(2.*pi**2.)#/vol
	if r1 == r2: #pure shot-noise piece has been separated out since this depends on the r bin size
		
	#	sumxi += 1./npp
	#	print r1,sqrt(1./npp)*r1
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
 			#npp += 1./3.*pi*((r1+dr/2.)**3.-(r1-dr/2.)**3.)*vt*nave**2.
			npp  += nave**2.*vt*2.*pi*r1**2.*dr
		sumxi += 1./npp
		print r1, r1*sqrt(sumxi), r1*sqrt(1./npp)
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

# def putallBAOqsomocks(sig=1,covmd='mock',bs=10):
# 	d = np.loadtxt('BAOmockQPM_QSOv1.0'+covmd+str(bs)+'.dat')
# 	ma = 0
# 	sa = 0
# 	siga = 0
# 	chia = 0
# 	n = 0
# 	print len(d)
# 	for i in range(0,len(d)):
# 		a = d[i]
# 		if sig == 1:
# 			s1b = float(a[1]),float(a[2])
# 		if sig == 2:
# 			s1b = float(a[3]),float(a[4])	
# 		if s1b[0] > .8 and s1b[1] < 1.2:
# 
# 			ma += a[0]
# 			sa += a[0]**2.
# 			siga += (float(a[2])-float(a[1]))/2.
# 			chia += a[-1]
# 			n += 1.
# 	ma = ma/n
# 	sa = sqrt(sa/n-ma**2.)
# 	siga = siga/n
# 	chia = chia/n
# 	return ma,sa,siga,chia,n

def putallBAOqsomocks(N=1000,sig=1,mock='qpm_qso',covmd='EZmock',bs=5,start=0,version='v1.6',mb='',Bp=0.4):
	ma = 0
	sa = 0
	siga = 0
	chia = 0
	n = 0
	bsst = str(bs)+'st'+str(start)
	for i in range(1,1+N):
		fl = ''
		if i < 1000:
			fl += '0'
		if i < 100:
			fl += '0'
		if i < 10:
			fl += '0'
		fl += str(i)

		a = sigreg_c12('BAOxichil'+mock+fl+version+covmd+mb+str(Bp)+bsst)
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
	

def xibao(sample,zmin,zmax,version='v1.6',wm='',bs=5,start=0,rmin=30,rmax=180.,md=1.,m=1.,mb='',Bp=.4,v='n',mockn='',covmd='EZmock',damp='3.04.08.0',Nmock=1000,template='Challenge_matterpower'):
	#does baofits, set mb='nobao' to do no BAO fit
	ebossdir = '/Users/ashleyross/eboss/'
	from baofit_pubtest import doxi_isolike
	from Cosmo import distance
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	bsst = str(bs)+'st'+str(start)
	if sample == 'lrg' or sample == 'QSO' or sample == 'QSOsys':
		dn = np.loadtxt(ebossdir+'xi0geboss'+sample+'_N'+version+'_'+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0geboss'+sample+'_S'+version+'_'+wz+wm+bsst+'.dat').transpose()
		#dn = (np.loadtxt(ebossdir+'xi0geboss'+sample+'_N'+version+'_'+wz+wm+bsst+'.dat').transpose()+np.loadtxt(ebossdir+'xi0geboss'+sample+'_N'+version+'_'+wz+wm+bsst+'_1.dat').transpose())/2.
		#ds = (np.loadtxt(ebossdir+'xi0geboss'+sample+'_S'+version+'_'+wz+wm+bsst+'.dat').transpose()+np.loadtxt(ebossdir+'xi0geboss'+sample+'_S'+version+'_'+wz+wm+bsst+'_1.dat').transpose())/2.
	if sample == 'QPM_QSO':	
		dn = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'NGC'+version+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'SGC'+version+wz+wm+bsst+'.dat').transpose()
	if sample == 'qpm_qso':	
		ebossdir=''
		dn = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'ngc'+version+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'sgc'+version+wz+wm+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
	if sample == 'aveQPM_QSO':	
		dn = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'NGC'+version+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'SGC'+version+wz+wm+bsst+'.dat').transpose()
	if sample == 'aveqpm_qso':	
		dn = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'ngc'+version+wz+wm+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0g'+sample+mockn+'sgc'+version+wz+wm+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
	if sample == 'aveQPMmockv1.6.1PZ':	
		dn = np.loadtxt(ebossdir+'xi0'+sample+'ngc'+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi0'+sample+'sgc'+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)

	if sample == 'QSOEZPZ':	
		dn = np.loadtxt(dirscio+'xi02EZmockNGC'+mockn+bsst+'.dat').transpose()
		ds = np.loadtxt(dirscio+'xi02EZmockSGC'+mockn+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
		ebossdir=''
	if sample == 'EZmock_QSO':	
		fl = 'g'+sample
		dn = np.loadtxt('xi0'+fl+mockn+'ngc'+version+wz+bsst+'.dat').transpose()
		ds = np.loadtxt('xi0'+fl+str(1000+int(mockn))+'sgc'+version+wz+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
		ebossdir=''

	if sample == 'QSOaveEZ':
		dn = np.loadtxt(ebossdir+'xiave0gEZmock_QSOngcv1.6mz0.8xz2.2'+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xiave0gEZmock_QSOsgcv1.6mz0.8xz2.2'+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
	if sample == 'QSOaveEZPZ':
		dn = np.loadtxt(ebossdir+'xiave0EZPZNGC'+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xiave0EZPZSGC'+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
	if sample == 'QSOaveQPMPZ':
		dn = np.loadtxt(ebossdir+'xiave0QPMmockv1.6PZngc'+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xiave0QPMmockv1.6PZsgc'+bsst+'.dat').transpose()
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
			
	if sample == 'lrg':
		wt = (ds[1]*1.8+1.3*dn[1])/3.1
	if sample == 'QSO' or sample == 'QPM_QSO' or sample == 'aveQPM_QSO' or sample == 'QSOsys':
		#if version == 'v1.2' or version == 'v1.5' or version == 'v1.3':
		wt = (ds[1]*.537+dn[1]*.716)/(.537+.716)
		#else:
		#	wt = (ds[1]+dn[1])/2.
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
		mod = np.loadtxt('BAOtemplates/xi0sm'+template+'0.4'+damp+'15.00.dat').transpose()[1]
	else:
		#mod = np.loadtxt('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()[1]
		if sample == 'QPM_QSO':
			mod = np.loadtxt('BAOtemplates/xi0'+template+'0.4'+damp+'15.00.dat').transpose()[1]
		else:
			mod = np.loadtxt('BAOtemplates/xi0'+template+'0.4'+damp+'15.00.dat').transpose()[1]


	csample = sample
	#facc = .716/(.537+.716)
	#print facc
	if sample == 'aveQPM_QSO' or sample == 'QPM_QSO' or sample == 'QSOaveEZ':
		if covmd == 'an':
			csample = 'QSO'
		if covmd == 'mock':
			csample = 'QPM_QSO'	
	if sample == 'QSO' and covmd == 'mock':
		csample = 'QPM_QSO'	
	if covmd == 'an':	
		cov = np.loadtxt(ebossdir+'covxiNS'+csample+'v1.6'+wz+str(float(bs))+'.dat')
	if covmd == 'mock':
		cov = np.loadtxt(ebossdir+'cov0'+csample+version+'NScomb'+str(bs)+'st'+str(start)+'.dat')
	cov2 = ''
	if covmd == 'EZmock':
		cov = np.loadtxt(ebossdir+'cov0EZmock_QSOv1.6ngc'+bsst+'.dat')
		cov2 = np.loadtxt(ebossdir+'cov0EZmock_QSOv1.6sgc'+bsst+'.dat')
	if covmd == 'QPMmock':
		cov = np.loadtxt(ebossdir+'cov0qpm_qsov1.6ngc'+str(bs)+'st0.dat')
		cov2 = np.loadtxt(ebossdir+'cov0qpm_qsov1.6sgc'+str(bs)+'st0.dat')
	if covmd == 'QPMmockPZ':
		cov = np.loadtxt(ebossdir+'covQPMmockv1.6.1PZ0ngc'+str(bs)+'st0.dat')
		cov2 = np.loadtxt(ebossdir+'covQPMmockv1.6.1PZ0sgc'+str(bs)+'st0.dat')

	cov = cov*m	
	if md != 1:
		for i in range(0,len(cov)):
			cov[i][i] = cov[i][i]*md
				
	
	chil = doxi_isolike(dd,cov,mod,rl,rmin=rmin,rmax=rmax,v=v,wo=sample+version+covmd+mb,Bp=Bp,cov2=cov2,Nmock=Nmock)
	fo = open(ebossdir+'BAOxichil'+sample+mockn+version+covmd+mb+str(Bp)+bsst+'.dat','w')
	print ebossdir+'BAOxichil'+sample+mockn+version+covmd+mb+str(Bp)+bsst+'.dat'
	for i in range(0,len(chil)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chil[i])+'\n')
	fo.close()
	a = sigreg_c12(ebossdir+'BAOxichil'+sample+mockn+version+covmd+mb+str(Bp)+bsst)
	#print a
	return a

def xibaoNS(sample,zmin,zmax,version='v1.5',wm='',bs=5,start=0,rmin=30,rmax=150.,md=1.,m=1.,mb='',Bp=.4,v='n',mockn='',covmd='mock'):
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
	if sample == 'QSOaveEZPZ':	
		dn = np.loadtxt(ebossdir+'xiave0EZPZNGC'+bsst+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xiave0EZPZSGC'+bsst+'.dat').transpose()

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
	cov2 = ''
	if covmd == 'EZmock':
		covN = np.loadtxt(ebossdir+'covEZPZ0NGC5st0.dat')
		covS = np.loadtxt(ebossdir+'covEZPZ0SGC5st0.dat')
				
	
	chiln = doxi_isolike(dn[1],covN,mod,rl,rmin=rmin,rmax=rmax,v=v,wo=sample+version+mb,Bp=Bp)
	fo = open(ebossdir+'BAOxichilNGC'+sample+mockn+version+mb+str(Bp)+'.dat','w')
	for i in range(0,len(chiln)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chiln[i])+'\n')
	fo.close()
	an = sigreg_c12(ebossdir+'BAOxichilNGC'+sample+mockn+version+mb+str(Bp))
	print an
	chils = doxi_isolike(ds[1],covS,mod,rl,rmin=rmin,rmax=rmax,v=v,wo=sample+version+mb,Bp=Bp)
	fo = open(ebossdir+'BAOxichilSGC'+sample+mockn+version+mb+str(Bp)+'.dat','w')
	for i in range(0,len(chils)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chils[i])+'\n')
	fo.close()
	als = sigreg_c12(ebossdir+'BAOxichilSGC'+sample+mockn+version+mb+str(Bp))
	print als
	chilt = np.array(chiln)+np.array(chils)
	fo = open(ebossdir+'BAOxichilNScomb'+sample+mockn+version+mb+str(Bp)+'.dat','w')
	for i in range(0,len(chilt)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chilt[i])+'\n')
	fo.close()
	a = sigreg_c12(ebossdir+'BAOxichilNScomb'+sample+mockn+version+mb+str(Bp))
	
	#print a
	return a

def mksubfile_mocks(ind):
	
	fo = open('sub'+str(ind)+'.sh','w')
	fl = ''
	if ind < 1000:
		fl += '0'
	if ind < 100:
		fl += '0'
	if ind < 10:
		fl += '0'
	fl += str(ind)
	#print fl
	fo.write('#!/bin/bash\n')
	fo.write('#$ -V -cwd\n')
	fo.write('. /etc/profile.d/modules.sh \n')
	fo.write('module add  apps/gcc/python/2.7.3 \n')
	#fo.write('/opt/gridware/apps/gcc/python/2.7.3/bin/python baofit.py '+file +' '+str(B)+' CMASS '+reg+' '+col+' '+tp +' '+str(bs)+' '+str(st)+'\n')
	fo.write('python xitools_eboss.py '+fl+'\n')
	fo.close()
	return True

	
def mksuball_mocks(start,N):
	fo = open('suball.sh','w')
	fo.write('#!/bin/bash\n')
	for i in range(start,start+N):
		mksubfile_mocks(i)
		#fo.write('qsub sub'+str(i)+col+str(bs)+'.sh -q sciama1.q\n')
		fo.write('qsub sub'+str(i)+'.sh \n')
		fo.write('sleep .2 \n')
	fo.close()
	return True


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
	for j in range(0,10000):
		z = .905
		l = []
		wl = []
		while z < 2.2:
			sigma = d.D(z)**p #testing how BAO scales with growth factor due to damping
			w = 1/sigma**2.
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
	#return vu,vw
	return sqrt(vu/vw) #factor by which sigma should change

def testmuweights():
	from Cosmo import distance
	from random import gauss
	#test difference in recovered signal to noise with and without weighting
	url = []
	wrl = []
	for j in range(0,10000):
		mu = .005
		l = []
		wl = []
		while mu < 1.:
			w = 1.+2.*mu**2.
			sigma = sqrt(1./w)
			for i in range(0,10):
				v = gauss(0,sigma)
				l.append(v)
				wl.append(w)
			mu += .01
		u = np.mean(l)
		wr = sum(np.array(l)*np.array(wl))/sum(wl)
		url.append(u)
		wrl.append(wr)
	print len(url),len(wrl),np.mean(url),np.mean(wrl)	
	vu = np.mean(np.array(url)**2.)-np.mean(url)**2.
	vw = np.mean(np.array(wrl)**2.)-np.mean(wrl)**2.
	#return vu,vw
	return sqrt(vu/vw) #factor by which sigma should change


def QSObaovsz(n=1.6e-05,dz=.1,zmin=.9,zmax=2.2,area=2000.):
	from baoerr import baoerr
	from Cosmo import distance
	d = distance(.3,.7)
	ol = []
	zl = []
	nl = []
	nb = int((zmax-zmin+dz*.01)/dz)
	for i in range(0,nb):
		z = zmin+dz/2.+dz*i
		b = 0.53+0.29*(1.+z)**2.
		errz = baoerr(z-dz/2.,z+dz/2.,dz,.001,area,n,b)
		print z,errz
		ol.append(1./errz**2.)
		zl.append(z)
		vol = (d.dc(z+dz/2.)**3.-d.dc(z-dz/2.)**3.) #volume of shell is proportional to this
		numbin = n*vol #number per zbin will be proportional to this
		nl.append(numbin)
	ivt = sum(ol)
	errt = sqrt(1./ivt)
	#print ol
	ol = np.array(ol)
	nl = np.array(nl)
	from matplotlib import pyplot as plt
	plt.plot(zl,ol**.5/nl,'k-')
	plt.show()
	return errt	

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

def compdecreg(sample='QSO',version='v1.6',NS='S',dir='/Users/ashleyross/fitsfiles/',zmin=.8,zmax=2.2,wm='',gri22=''):
	from healpix import healpix,radec2thphi
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat.fits')
	fr = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.ran.fits')
	ngh = 0
	nrh = 0
	ngl = 0
	nrl = 0
	ncph = 0
	ncpl = 0
	nzfh = 0
	nzfl = 0
	nwsh = 0
	nwsl = 0
	ngth = 0
	ngtl = 0
	nrth = 0
	nrtl = 0
	nzvih = 0
	nzvil = 0
	umgh = 0
	umgl = 0
	umgvh = 0
	umgvl = 0
	if wm == 'wdepthext':
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()	
		b,m = findlinmb(sample,NS,version,'depth',.8,2.2,wm='nosys'+gri22,dir=ebossdir)
		be,me = findlinmb(sample,NS,version,'ext',.8,2.2,wm='wdepth'+gri22,dir=ebossdir)
		wext = np.loadtxt(dirsys+'healSFD_r_'+str(256)+'_fullsky.dat')/2.751
		h = healpix()

	for i in range(0,len(f)):
		z = f[i]['Z']
		if z > zmin and z < zmax:
			um = f[i]['MODELMAG'][0]-f[i]['EXTINCTION'][0]
			gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT']
			wss = f[i]['WEIGHT_SYSTOT']
			if wm == 'wdepthext':
				th,phi = radec2thphi(f[i]['RA'],f[i]['DEC'])
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				pixe = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				ext = wext[pixe]
				we = 1./(be+me*ext)
				wss = ws*we
				w = w*wss			
			zv = 0
			if f[i]['Z_VI'] == z:
				zv = 1.
			if f[i]['DEC'] > 10:
				ngh += w
				ngth += 1.
				ncph += f[i]['WEIGHT_CP']
				nzfh += f[i]['WEIGHT_NOZ']
				nwsh += wss
				nzvih += zv
				umgh += um-gm
				umgvh += (um-gm)**2.
			else:
				ngl += w
				ngtl += 1.
				ncpl += f[i]['WEIGHT_CP']
				nzfl += f[i]['WEIGHT_NOZ']
				nwsl += wss
				nzvil += zv
				umgl += um-gm
				umgvl += (um-gm)**2.

	for i in range(0,len(fr)):
		z = fr[i]['Z']
		if z > zmin and z < zmax:
			if fr[i]['DEC'] > 10:
				nrh += fr[i]['WEIGHT_FKP']
				nrth += 1.
			else:
				nrl += fr[i]['WEIGHT_FKP']
				nrtl += 1.
	umgl = umgl/ngtl
	umgvl = sqrt(umgvl/ngtl-umgl*umgl)
	umgh = umgh/ngth
	umgvh = sqrt(umgvh/ngth-umgh*umgh)
	print umgl,umgh,umgvl,umgvh,ngl/nrtl,ngh/nrth,ngtl/nrl,ngth/nrh,ncph/ngth,ncpl/ngtl,nzfh/ngth,nzfl/ngtl,nwsh/ngth,nwsl/ngtl,nzvih/ngth,nzvil/ngtl
	return True

def compQSONS(sample='QSO',version='v1.6',dir='/Users/ashleyross/fitsfiles/',zmin=.8,zmax=2.2,wm='',gri22=''):
	from healpix import healpix,radec2thphi
	fn = fitsio.read(dir+'eboss_'+version+'-'+sample+'-N'+'-eboss_'+version+'.dat.fits')
	fs = fitsio.read(dir+'eboss_'+version+'-'+sample+'-S'+'-eboss_'+version+'.dat.fits')
	frn = fitsio.read(dir+'eboss_'+version+'-'+sample+'-N'+'-eboss_'+version+'.2.ran.fits')
	frs = fitsio.read(dir+'eboss_'+version+'-'+sample+'-S'+'-eboss_'+version+'.2.ran.fits')
	ngh = 0
	nrh = 0
	ngl = 0
	nrl = 0
	ncph = 0
	ncpl = 0
	nzfh = 0
	nzfl = 0
	nwsh = 0
	nwsl = 0
	ngth = 0
	ngtl = 0
	nrth = 0
	nrtl = 0
	nzvih = 0
	nzvil = 0
	umgh = 0
	umgl = 0
	umgvh = 0
	umgvl = 0
	if wm == 'wdepthext':
		wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()	
		b,m = findlinmb(sample,NS,version,'depth',.8,2.2,wm='nosys'+gri22,dir=ebossdir)
		be,me = findlinmb(sample,NS,version,'ext',.8,2.2,wm='wdepth'+gri22,dir=ebossdir)
		wext = np.loadtxt(dirsys+'healSFD_r_'+str(256)+'_fullsky.dat')/2.751
		h = healpix()

	for i in range(0,len(fn)):
		z = fn[i]['Z']
		if z > zmin and z < zmax:
			w = (fn[i]['WEIGHT_NOZ']+fn[i]['WEIGHT_CP']-1.)*fn[i]['WEIGHT_FKP']*fn[i]['WEIGHT_SYSTOT']
			wss = fn[i]['WEIGHT_SYSTOT']
			um = fn[i]['MODELMAG'][0]-fn[i]['EXTINCTION'][0]
			gm = fn[i]['MODELMAG'][1]-fn[i]['EXTINCTION'][1]

			if wm == 'wdepthext':
				th,phi = radec2thphi(fn[i]['RA'],fn[i]['DEC'])
				w = (fn[i]['WEIGHT_NOZ']+fn[i]['WEIGHT_CP']-1.)*fn[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				pixe = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				ext = wext[pixe]
				we = 1./(be+me*ext)
				wss = ws*we
				w = w*wss			
			zv = 0
			if fn[i]['Z_VI'] == z:
				zv = 1.
			ngh += w
			ngth += 1.
			ncph += fn[i]['WEIGHT_CP']
			nzfh += fn[i]['WEIGHT_NOZ']
			nwsh += wss
			nzvih += zv
			umgh += um-gm
			umgvh += (um-gm)**2.
	for i in range(0,len(fs)):
		z = fs[i]['Z']
		if z > zmin and z < zmax:
			w = (fs[i]['WEIGHT_NOZ']+fs[i]['WEIGHT_CP']-1.)*fs[i]['WEIGHT_FKP']*fs[i]['WEIGHT_SYSTOT']
			wss = fs[i]['WEIGHT_SYSTOT']
			um = fs[i]['MODELMAG'][0]-fs[i]['EXTINCTION'][0]
			gm = fs[i]['MODELMAG'][1]-fs[i]['EXTINCTION'][1]

			if wm == 'wdepthext':
				th,phi = radec2thphi(fs[i]['RA'],fs[i]['DEC'])
				w = (fs[i]['WEIGHT_NOZ']+fs[i]['WEIGHT_CP']-1.)*fs[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				pixe = h.ang2pix_nest(256,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				ext = wext[pixe]
				we = 1./(be+me*ext)
				wss = ws*we
				w = w*wss			
			zv = 0
			if fs[i]['Z_VI'] == z:
				zv = 1.
			ngl += w
			ngtl += 1.
			ncpl += fs[i]['WEIGHT_CP']
			nzfl += fs[i]['WEIGHT_NOZ']
			nwsl += wss
			nzvil += zv
			umgl += um-gm
			umgvl += (um-gm)**2.

	nrh = float(len(frn))
	nrl = float(len(frs))
	nrth = nrh
	nrtl = nrl
	#for i in range(0,len(fr)):
	#	z = fr[i]['Z']
	#	if z > zmin and z < zmax:
	#		if fr[i]['DEC'] > 10:
	#			nrh += fr[i]['WEIGHT_FKP']
	#			nrth += 1.
	#		else:
	#			nrl += fr[i]['WEIGHT_FKP']
	#			nrtl += 1.
	umgl = umgl/ngtl
	umgvl = sqrt(umgvl/ngtl-umgl*umgl)
	umgh = umgh/ngth
	umgvh = sqrt(umgvh/ngth-umgh*umgh)

	print ngl,nrtl,ngh,nrth
	print umgl,umgh,umgvl,umgvh,ngl/nrtl,ngh/nrth,ngtl/nrl,ngth/nrh,ncph/ngth,ncpl/ngtl,nzfh/ngth,nzfl/ngtl,nwsh/ngth,nwsl/ngtl,nzvih/ngth,nzvil/ngtl
	return True


def visinspectQSOcomp():
	f = fitsio.FITS('/Users/ashleyross/fitsfiles/DR14Q_v1_1.fits')
	nt = 0
	nvi = 0
	std = 0
	zi = []
	zp = []
	zdl = []
	dmin = -4.
	dmax = 4.
	sp = .01
	nzb = int((dmax-dmin+sp/10.)/sp)
	zdl = np.zeros((nzb))
	nb = 0
	for i in range(0,526359):
		tf = f[1]['EBOSS_TARGET1'][i][0]
		if tf & 2**10 > 0:
			nt += 1.
			zvi = f[1]['Z_VI'][i][0]
			if zvi > .01 and zvi < 10.:
				nvi += 1.
				std += (zvi-f[1]['Z_PIPE'][i][0])**2.
				zi.append(zvi)
				zp.append(f[1]['Z_PIPE'][i][0])
				zd = zvi - 	f[1]['Z_PIPE'][i][0]
				if zd > dmin and zd < dmax and abs(zd) > 0.01:
					ind = int((zd-dmin)/sp)
					zdl[ind] += 1.					
					nb += 1.
	print nb
	print nt
	print nvi,sqrt(std/nvi)
	zl = np.arange(dmin+sp/2.,dmax+sp/2.,.01)
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	import pylab
	pylab.ion()
	plt.clf()
	plt.minorticks_on()
	plt.plot(zl,zdl,'k-')
	#plt.save()
	return True
	#return zi,zp

### Below is for plots


def plotreg(sample='QSO',version='v1.6',NS='S',md='ran',zmin=.8,zmax=2.2,ramin=-180,ramax=180,decmin=-90,decmax=90,dir='/Users/ashleyross/fitsfiles/'):
	from matplotlib import pyplot as plt
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat.fits')
	fr = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.2.ran.fits')
	ral = []
	decl = []
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z > zmin and z < zmax:
			ra = f[i]['RA']
			if NS == 'S':
				if ra > 180:
					ra -= 360
			ral.append(ra)
			decl.append(f[i]['DEC'])
	rarl = []
	decrl = []
	for i in range(0,len(fr)):
		ra = fr[i]['RA']
		if NS == 'S':
			if ra > 180:
				ra -= 360
		rarl.append(ra)
		decrl.append(fr[i]['DEC'])
	plt.plot(rarl,decrl,'k,')
	plt.plot(ral,decl,'ro',markeredgecolor='r',markersize=1)
	plt.xlim(ramin,ramax)
	plt.ylim(decmin,decmax)
	plt.show()
	return True

def plotregetalam(sample='QSO',version='v1.6',NS='S',md='ran',ramin=-90,ramax=90,decmin=0,decmax=360,depthmin=0,dir='/Users/ashleyross/fitsfiles/',immin=0):
	from matplotlib import pyplot as plt
	from stripeDR8 import radec2le
	from healpix import radec2thphi, healpix
	h = healpix()
	f = fitsio.read(dir+'eboss_'+version+'-'+sample+'-'+NS+'-eboss_'+version+'.dat.fits')
	wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()
	ral = []
	decl = []
	for i in range(0,len(f)):
		ra = f[i]['RA']
		dec = f[i]['DEC']
		im = f[i]['IMATCH']
		th,phi = radec2thphi(ra,dec)
		p = h.ang2pix_nest(512,th,phi)		
		if wsys[p] > depthmin and im >= immin:
			eta,lam = radec2le(ra,dec)
			ral.append(eta)
			decl.append(lam)
	print max(ral),max(decl)	
	plt.plot(ral,decl,'k,')
	plt.xlim(ramin,ramax)
	plt.ylim(decmin,decmax)
	plt.show()
	return True


def plotdepthmap(ramin=-180,ramax=180,decmin=-20,decmax=90,size=1,smin=21):
	from healpix import radec2thphi, healpix,thphi2radec
	from matplotlib import pyplot as plt
	import matplotlib.cm as cm
	h = healpix()
	wsys = np.loadtxt(dirsys+'healdepthinm512.dat').transpose()
	ral = []
	decl = []
	wl = []
	for i in range(0,len(wsys)):
		if wsys[i] < 24:
			th,phi = h.pix2ang_nest(512,i)
			ra,dec = thphi2radec(th,phi)
			if ra > 180:
				ra = ra -360
			ral.append(ra)
			decl.append(dec)
			wl.append(wsys[i])
	map = plt.scatter(ral,decl,c=wl,s=size,cmap=cm.rainbow,lw=0,vmin=smin)
	cbar = plt.colorbar(map)
	plt.xlim(ramin,ramax)
	plt.ylim(decmin,decmax)
	plt.show()
	return True


def plotzhistplvi(zmin=.8,zmax=2.2,nb=100,mn=-.2,mx=.2,sig=.001,dir='/Users/ashleyross/fitsfiles/'):
	from matplotlib import pyplot as plt
	import matplotlib.mlab as mlab
	f = fitsio.read(dir+'sequels_v1.61-QSO-N-sequels_v1.61.dat.fits')
	hst = np.zeros((nb))
	nout = 0
	nt = 0
	std = 0
	stdno = 0
	ntno = 0
	for i in range(0,len(f)):
		if f[i]['Z_VI'] > 0 and f[i]['Z_PL'] > zmin and f[i]['Z_PL'] < zmax and f[i]['ZWARNING'] <= 0:
			df = f[i]['Z_VI']- f[i]['Z_PL']
			if df < mx and df > mn:
				stdno += (f[i]['Z_VI']- f[i]['Z_PL'])**2.
				ntno += 1.
			if df >= mx:
				df = mx-(mx-mn)/(10.*nb)
				nout += 1.
			if df <= mn:
				df = mn+(mx-mn)/(10.*nb)
				nout += 1.	
			nt += 1.
			std += (f[i]['Z_VI']- f[i]['Z_PL'])**2.	
			bn = int(nb*((df-mn)/(mx-mn)))
			hst[bn] += 1
	sp = (mx-mn)/float(nb)		
	xl = np.arange(mn+sp/2.,mx,sp)
	lzl = []
	stdc = 0
	norm = 0
	for i in range(0,len(xl)):
		lzi = pi*sqrt(stdno/ntno)*.17*(1.+(xl[i]/(sqrt(stdno/ntno)*.17))**2.)
		lzl.append(1./lzi)
		norm += 1./lzi
		stdc += xl[i]**2./lzi
	print nout,nt,ntno,sqrt(std/nt),sqrt(stdno/(ntno)),sqrt(stdc/norm)
	plt.plot(xl,hst/(sum(hst)*sp),'k-')
	#plt.plot(xl,mlab.normpdf(xl, 0, sqrt(stdno/(ntno))),'r-')
	plt.plot(xl,lzl,'r-')
	plt.plot(xl,mlab.normpdf(xl, 0, sig),'b-')
	plt.xlabel(r'$Z_{VI}-Z_{PL}$')
	plt.ylabel('N')
	plt.show()
	return True

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

def plotxiCLRGNScomb(bs='5st0',wm='',v='v1.5_IRt'):
	#Plots comparison between NGC and SGC clustering for and to theory LRGs, no stellar density correction
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xicLRGNScomb'+v+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.43.06.010.015.00.dat').transpose()
	#cov = np.loadtxt(ebossdir+'covxiNSlrgv0.8_IRcmz0.6xz1.010.dat')
	#et = []
	#for i in range(0,len(cov)):
	#	et.append(sqrt(cov[i][i]))
	#plt.fill_between(dt[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	ds = np.loadtxt(ebossdir+'xi0gebosscmass-lrg_S'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi0gebosscmass-lrg_N'+v+'_mz0.6xz1.0'+bs+'.dat').transpose()
	
	t = (ds[1]*.54+.63*dn[1])/(.63+.54)
	#plt.plot(ds[0],ds[0]**2.*ds[1],'b-s')
	#plt.plot(dn[0],dn[0]**2.*dn[1],'r-d')
	plt.plot(dn[0],dn[0]**2.*t,'k-o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*2.1,'k:')
	plt.xlim(20,200)
	plt.ylim(-49,120)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,180,'SGC',color='b')
	#plt.text(30,170,'NGC',color='r')
	plt.text(120,100,'N+S Combined',color='k')
	plt.title(r'Correlation function of v1.5 CMASS + eboss LRGs, 0.6 < z < 1.0')
	pp.savefig()
	pp.close()
	return True

def plotnzCLRG():
	ds = np.loadtxt(ebossdir+'nbar-cmass-eboss_v1.5_IRt-lrg-S-eboss_v1.5_IRt.dat').transpose()
	dn = np.loadtxt(ebossdir+'nbar-cmass-eboss_v1.5_IRt-lrg-N-eboss_v1.5_IRt.dat').transpose()
	nt = (ds[3]+dn[3])/2.
	dse = np.loadtxt(ebossdir+'nbar-eboss_v1.5_IRt-lrg-S-eboss_v1.5_IRt.dat').transpose()
	dne = np.loadtxt(ebossdir+'nbar-eboss_v1.5_IRt-lrg-N-eboss_v1.5_IRt.dat').transpose()
	nte = (dse[3]+dne[3])/2.
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	#pp = PdfPages(ebossdir+'nzcmassLRGv1.5.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.plot(ds[0],nt,'k-',ds[0],nte,'r-')
	plt.xlim(0.5,1.1)
	plt.ylim(0,.0003)
	plt.show()
	return True

def plotnzQSO():
	ds = np.loadtxt(ebossdir+'nbar-eboss_v1.5-QSO-S-eboss_v1.5.dat').transpose()
	dn = np.loadtxt(ebossdir+'nbar-eboss_v1.5-QSO-N-eboss_v1.5.dat').transpose()
	nt = (ds[3]+dn[3])/2.*1.e5
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'nzQSOv1.5.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.plot(ds[0],nt,'k-')
	plt.xlim(0.5,2.5)
	#plt.ylim(0,.0003)
	plt.ylabel(r'$10^{5}$ $n$ ($h^{3}$Mpc$^{-3}$)',size=16)
	plt.xlabel('redshift',size=16)

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

def plotxiQSONScompEZ(mom='0',bs='8st0',v='v1.62',mini=1,maxi=25,wm='',mumin=0,mumax=1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	muw = ''
	muwd = ''
	if mumin != 0:
		muw += 'mmu'+str(mumin)
		muwd += 'mum'+str(mumin)
	if mumax != 1:
		muw += 'xmu'+str(mumax)	
		muwd += 'mux'+str(mumax)
	norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSONScompEZ'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	aves = np.loadtxt(ebossdir+'xiave'+mom+'gEZmock_QSOsgcv1.6mz0.8xz2.2'+muw+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xiave'+mom+'gEZmock_QSOngcv1.6mz0.8xz2.2'+muw+bs+'.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov'+mom+'EZmock_QSOv1.6sgc'+bs+'.dat')[mini:maxi,mini:maxi]
	covn = np.loadtxt(ebossdir+'cov'+mom+'EZmock_QSOv1.6ngc'+bs+'.dat')[mini:maxi,mini:maxi]
	diffs = aves[1][mini:maxi]*norm**2.-ds[1][mini:maxi]
	facn = 1.
	facs = 1.
	if wm == 'gri22depthi22' or wm == 'gri22depthi22wdepthimag' or wm == 'gri22depthi22ext0.15wdepthimagext' or wm == 'gri22depthi22ext0.15wdepthext':
		facs = 47494/53693.
		facn = 68488/71576.0

	chis = np.dot(np.dot(diffs,np.linalg.pinv(covs)),diffs)*facs
	diffn = aven[1][mini:maxi]*norm**2.-dn[1][mini:maxi]
	chin = np.dot(np.dot(diffn,np.linalg.pinv(covn)),diffn)*facn
	diff = ds[1][mini:maxi]-dn[1][mini:maxi]
	cov = covn+covs
	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)
	print chis,chin,chi
	plt.plot(aven[0][mini:maxi],aven[0][mini:maxi]**2.*aven[1][mini:maxi]*norm**2.,'r--')
	plt.plot(aves[0][mini:maxi],aves[0][mini:maxi]**2.*aves[1][mini:maxi]*norm**2.,'b--')
	plt.errorbar(ds[0][mini:maxi]-.5,ds[0][mini:maxi]**2.*ds[1][mini:maxi],ds[0][mini:maxi]**2.*aves[2][mini:maxi]*norm**2.,fmt='bs')
	plt.errorbar(dn[0][mini:maxi]+.5,dn[0][mini:maxi]**2.*dn[1][mini:maxi],dn[0][mini:maxi]**2.*aven[2][mini:maxi]*norm**2.,fmt='rd')
	
	plt.xlim(ds[0][mini]-2.,ds[0][maxi]+2.)
	if mom == '0':
		plt.ylim(-80,100)
	else:
		plt.ylim(-200,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{'+mom+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == '0':
		plt.text(30,90,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color='b')
		plt.text(30,83,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color='r')
	else:
		plt.text(30,180,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color='b')
		plt.text(30,165,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color='r')
	
	#plt.text(30,160,'Combined',color='k')
	#plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True

def plotxiQSOcompEZ(mom='0',NS='ngc',bs='8st0',v='v1.62',mini=3,maxi=25,wm='',mumin=0,mumax=1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	muw = ''
	muwd = ''
	if mumin != 0:
		muw += 'mmu'+str(mumin)
		muwd += 'mum'+str(mumin)
	if mumax != 1:
		muw += 'xmu'+str(mumax)	
		muwd += 'mux'+str(mumax)
	norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSO'+NS+'compEZ'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	if NS == 'sgc':
		d = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	else:
		d = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	ave = np.loadtxt(ebossdir+'xiave'+mom+'gEZmock_QSO'+NS+'v1.6mz0.8xz2.2'+muw+bs+'.dat').transpose()
	cov = np.loadtxt(ebossdir+'cov'+mom+'EZmock_QSOv1.6'+NS+bs+'.dat')[mini:maxi,mini:maxi]
	diff = ave[1][mini:maxi]-d[1][mini:maxi]
	facn = 1.
	facs = 1.

	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)
	print chi
	plt.plot(ave[0][mini:maxi],ave[0][mini:maxi]**2.*ave[1][mini:maxi]*norm**2.,'k--')
	plt.errorbar(d[0][mini:maxi],d[0][mini:maxi]**2.*d[1][mini:maxi],d[0][mini:maxi]**2.*ave[2][mini:maxi]*norm**2.,fmt='ko')
	
	plt.xlim(d[0][mini]-2.,d[0][maxi]+2.)
	if mom == '0':
		plt.ylim(-80,100)
	else:
		plt.ylim(-200,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{'+mom+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == '0':
		plt.text(30,90,NS+r', $\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='b')
	else:
		plt.text(30,180,NS+r', $\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='b')
	
	#plt.text(30,160,'Combined',color='k')
	#plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True


def plotxiQSONScompQPM(mom='0',bs='8st0',v='v1.6',mini=1,maxi=25,wm='',mumin=0,mumax=1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	muw = ''
	muwd = ''
	if mumin != 0:
		muw += 'mmu'+str(mumin)
		muwd += 'mum'+str(mumin)
	if mumax != 1:
		muw += 'xmu'+str(mumax)	
		muwd += 'mux'+str(mumax)
	norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSONScompQPM'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	aves = np.loadtxt(ebossdir+'xi'+mom+'gaveqpm_qsosgcv1.6mz0.8xz2.2'+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xi'+mom+'gaveqpm_qsongcv1.6mz0.8xz2.2'+bs+'.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov'+mom+'qpm_qsov1.6sgc8st0.dat')[mini:maxi,mini:maxi]*norm**4.
	covn = np.loadtxt(ebossdir+'cov'+mom+'qpm_qsov1.6ngc8st0.dat')[mini:maxi,mini:maxi]*norm**4.
	#aves = np.loadtxt(ebossdir+'xi'+mom+'aveQPMmockv1.6.1PZsgc'+bs+'.dat').transpose()
	#aven = np.loadtxt(ebossdir+'xi'+mom+'aveQPMmockv1.6.1PZngc'+bs+'.dat').transpose()
	#covs = np.loadtxt(ebossdir+'covQPMmockv1.6.1PZ'+mom+'sgc8st0.dat')[mini:maxi,mini:maxi]*norm**4.
	#covn = np.loadtxt(ebossdir+'covQPMmockv1.6.1PZ'+mom+'ngc8st0.dat')[mini:maxi,mini:maxi]*norm**4.


	#covs = np.loadtxt(ebossdir+'covEZPZ'+mom+'SGC'+muw+bs+'.dat')[mini:maxi,mini:maxi]*norm**4.
	#covn = np.loadtxt(ebossdir+'covEZPZ'+mom+'NGC'+muw+bs+'.dat')[mini:maxi,mini:maxi]*norm**4.
	diffs = aves[1][mini:maxi]*norm**2.-ds[1][mini:maxi]
	facn = 1.
	facs = 1.
	if wm == 'gri22depthi22' or wm == 'gri22depthi22wdepthimag' or wm == 'gri22depthi22ext0.15wdepthimagext' or wm == 'gri22depthi22ext0.15wdepthext':
		facs = 47494/53693.
		facn = 68488/71576.0

	chis = np.dot(np.dot(diffs,np.linalg.pinv(covs)),diffs)*facs
	diffn = aven[1][mini:maxi]*norm**2.-dn[1][mini:maxi]
	chin = np.dot(np.dot(diffn,np.linalg.pinv(covn)),diffn)*facn
	print chis,chin
	plt.plot(aven[0][mini:maxi],aven[0][mini:maxi]**2.*aven[1][mini:maxi]*norm**2.,'r--')
	plt.plot(aves[0][mini:maxi],aves[0][mini:maxi]**2.*aves[1][mini:maxi]*norm**2.,'b--')
	plt.errorbar(ds[0][mini:maxi]-.5,ds[0][mini:maxi]**2.*ds[1][mini:maxi],ds[0][mini:maxi]**2.*aves[2][mini:maxi]*norm**2.,fmt='bs')
	plt.errorbar(dn[0][mini:maxi]+.5,dn[0][mini:maxi]**2.*dn[1][mini:maxi],dn[0][mini:maxi]**2.*aven[2][mini:maxi]*norm**2.,fmt='rd')
	
	plt.xlim(ds[0][mini]-2.,ds[0][maxi]+2.)
	if mom == '0':
		plt.ylim(-80,100)
	else:
		plt.ylim(-200,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{'+mom+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == '0':
		plt.text(30,90,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color='b')
		plt.text(30,83,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color='r')
	else:
		plt.text(30,180,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color='b')
		plt.text(30,165,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color='r')
	
	#plt.text(30,160,'Combined',color='k')
	#plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True


def plotxiQSONScompEZQPM(mom='0',bs='8st0',v='v1.62',mini=1,maxi=25,wm='',mumin=0,mumax=1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	muw = ''
	muwd = ''
	if mumin != 0:
		muw += 'mmu'+str(mumin)
		muwd += 'mum'+str(mumin)
	if mumax != 1:
		muw += 'xmu'+str(mumax)	
		muwd += 'mux'+str(mumax)
	norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSONScompEZQPM'+muw+v+wm+bs+'t.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	dt = (ds*.537+dn*.716)/(.537+.716)
	dsnw = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2nw'+muwd+bs+'.dat').transpose()
	dnnw = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2nw'+muwd+bs+'.dat').transpose()
	dtnw = (dsnw*.537+dnnw*.716)/(.537+.716)

	aves = np.loadtxt(ebossdir+'xiave'+mom+'gEZmock_QSOsgcv1.6mz0.8xz2.2'+muw+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xiave'+mom+'gEZmock_QSOngcv1.6mz0.8xz2.2'+muw+bs+'.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov'+mom+'EZmock_QSOv1.6sgc'+bs+'.dat')[mini:maxi,mini:maxi]
	covn = np.loadtxt(ebossdir+'cov'+mom+'EZmock_QSOv1.6ngc'+bs+'.dat')[mini:maxi,mini:maxi]
	
	avet = (aves*.537+aven*.716)/(.537+.716)
	covti = np.linalg.pinv(covs)+np.linalg.pinv(covn)
	#avesq = np.loadtxt(ebossdir+'xi'+mom+'gaveqpm_qsosgcv1.6mz0.8xz2.2'+bs+'.dat').transpose()
	#avenq = np.loadtxt(ebossdir+'xi'+mom+'gaveqpm_qsongcv1.6mz0.8xz2.2'+bs+'.dat').transpose()
	avesq = np.loadtxt(ebossdir+'xi'+mom+'aveQPMmockv1.6.1PZsgc'+bs+'.dat').transpose()
	avenq = np.loadtxt(ebossdir+'xi'+mom+'aveQPMmockv1.6.1PZngc'+bs+'.dat').transpose()
	
	aveqt = (avesq*.537+avenq*.716)/(.537+.716)
	covsq = np.loadtxt(ebossdir+'cov'+mom+'qpm_qsov1.6sgc8st0.dat')[mini:maxi,mini:maxi]*norm**4.
	covnq = np.loadtxt(ebossdir+'cov'+mom+'qpm_qsov1.6ngc8st0.dat')[mini:maxi,mini:maxi]*norm**4.
	#covsq = np.loadtxt(ebossdir+'covQPMmockv1.6.1PZ'+mom+'sgc8st0.dat')[mini:maxi,mini:maxi]*norm**4.
	#covnq = np.loadtxt(ebossdir+'covQPMmockv1.6.1PZ'+mom+'ngc8st0.dat')[mini:maxi,mini:maxi]*norm**4.

	covtqi = np.linalg.pinv(covsq)+np.linalg.pinv(covnq)
	diff = avet[1][mini:maxi]*norm**2.-dt[1][mini:maxi]
	diffq = aveqt[1][mini:maxi]*norm**2.-dt[1][mini:maxi]
	diffnw= avet[1][mini:maxi]*norm**2.-dtnw[1][mini:maxi]
	diffnwq = aveqt[1][mini:maxi]*norm**2.-dtnw[1][mini:maxi]
	diffd = dt[1][mini:maxi]-dtnw[1][mini:maxi]
	chi = np.dot(np.dot(diff,covti),diff)
	chiq = np.dot(np.dot(diffq,covtqi),diffq)
	chinw = np.dot(np.dot(diffnw,covti),diffnw)
	chinwq = np.dot(np.dot(diffnwq,covtqi),diffnwq)
	chidnw = np.dot(np.dot(diffd,covti),diffd)
	chidnwq = np.dot(np.dot(diffd,covtqi),diffd)
	print chi,chiq,chinw,chinwq,chidnw,chidnwq
	plt.plot(dt[0][mini:maxi],dt[0][mini:maxi]**2.*dt[1][mini:maxi],'k-',linewidth=2)
	plt.plot(dt[0][mini:maxi],dt[0][mini:maxi]**2.*dtnw[1][mini:maxi],'--',color='.5',linewidth=1)
	plt.errorbar(avet[0][mini:maxi]-.5,avet[0][mini:maxi]**2.*avet[1][mini:maxi],avet[0][mini:maxi]**2.*avet[2][mini:maxi]*norm**2.,fmt=':bs')
	plt.errorbar(avet[0][mini:maxi]+.5,aveqt[0][mini:maxi]**2.*aveqt[1][mini:maxi],aveqt[0][mini:maxi]**2.*aveqt[2][mini:maxi]*norm**2.,fmt=':rd')
	
	plt.xlim(dt[0][mini-1],dt[0][maxi])
	print dt[0][mini]
	if mom == '0':
		plt.ylim(-60,80)
	else:
		plt.ylim(-200,100)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{'+mom+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == '0':
		plt.text(20,-20,r'EZ, $\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini)+' ('+str(chinw)[:3]+')',color='b')
		plt.text(20,-27,r'QPM, $\chi^2$/dof ='+str(chiq)[:4]+'/'+str(maxi-mini)+' ('+str(chinwq)[:3]+')',color='r')
		#if wm == '':
		#	plt.text(30,76,'DR14 QSO sample, '+v+', 0.8<z<2.2, no extra cuts/weights')
		#else:
		plt.text(30,-34,'DR14 QSO sample')#, '+v+' no extra cuts/weights='+wm)
		plt.text(30,-41,r'DR14 QSO sample, no $w_{\rm sys}$',color='.5')
	else:
		plt.text(30,80,r'EZ, $\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='b')
		plt.text(30,70,r'QPM, $\chi^2$/dof ='+str(chiq)[:4]+'/'+str(maxi-mini),color='r')
	
	#plt.text(30,160,'Combined',color='k')
	#plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	print ebossdir
	return True


def plotxiQSONScompEZone(NS,mom='0',bs='8st0',v='v1.6',samp='QSOsys',zmin=.8,mini=3,maxi=25,wm='',mumin=0,mumax=1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	muw = ''
	muwd = ''
	if mumin != 0:
		muw += 'mmu'+str(mumin)
		muwd += 'mum'+str(mumin)
	if mumax != 1:
		muw += 'xmu'+str(mumax)	
		muwd += 'mux'+str(mumax)
	norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+samp+NS+'compEZ'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	d = np.loadtxt(ebossdir+'xi'+mom+'geboss'+samp+'_'+NS+v+'_mz'+str(zmin)+'xz2.2'+wm+muwd+bs+'.dat').transpose()
	ave = np.loadtxt(ebossdir+'xiave'+mom+'EZPZ'+NS+'GC'+muw+bs+'.dat').transpose()
	cov = np.loadtxt(ebossdir+'covEZPZ'+mom+NS+'GC'+muw+bs+'.dat')[mini:maxi,mini:maxi]*norm**4.
	diff = ave[1][mini:maxi]*norm**2.-d[1][mini:maxi]
	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)
	print chi
	plt.plot(ave[0][mini:maxi],ave[0][mini:maxi]**2.*ave[1][mini:maxi]*norm**2.,'k--')
	plt.errorbar(d[0][mini:maxi]-.5,d[0][mini:maxi]**2.*d[1][mini:maxi],d[0][mini:maxi]**2.*ave[2][mini:maxi]*norm**2.,fmt='ko')
	plt.xlim(d[0][mini]-2.,d[0][maxi]+2.)
	if mom == '0':
		plt.ylim(-80,100)
	else:
		plt.ylim(-200,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{'+mom+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == '0':
		plt.text(30,90,NS+'GC, '+r'$\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='k')
	else:
		plt.text(30,180,NS+'GC, '+r'$\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='k')
	
	#plt.text(30,160,'Combined',color='k')
	#plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True

def plotxiQSONScompEZonecuts(NS,mom='0',bs='8st0',v='v1.6',mini=3,maxi=25,wm='gri22depthi22ext0.15wdepthext',mumin=0,mumax=1,covm='qpm',zmin=.8):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	muw = ''
	muwd = ''
	if mumin != 0:
		muw += 'mmu'+str(mumin)
		muwd += 'mum'+str(mumin)
	if mumax != 1:
		muw += 'xmu'+str(mumax)	
		muwd += 'mux'+str(mumax)
	norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSO'+NS+'comp'+covm+'cut'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	fac = 1.
	if wm == 'gri22':
		if NS == 'S':
			fac = 51671/53693.
		if NS == 'N':
			fac = 68488/71576.0	
	if wm == 'gri22depthi22' or wm == 'gri22depthi22wdepthimag' or wm == 'gri22depthi22ext0.15wdepthimagext' or wm == 'gri22depthi22ext0.15wdepthext' or wm == 'gri22depthi22nw':
		if NS == 'S':
			fac = 47494/53693.
		if NS == 'N':
			fac = 68488/71576.0
	if wm == 'gri22depthi22decx10wdepthext':
		if NS == 'S':
			fac = 23170/53693.

	if wm == 'gri22depthi22wdepthext':
		if NS == 'S':
			fac = 47494/53693.
		if NS == 'N':
			fac = 68268/71576.0
	if wm == 'gri22depthi22ext0.14wdepthext':
		if NS == 'S':
			fac = 46798/53693.
		if NS == 'N':
			fac = 68264/71576.0
	if wm == 'gri22depthi22ext0.12wdepthext':
		if NS == 'S':
			fac = 46798/53693.
		if NS == 'N':
			fac = 68223/71576.0
	if wm == 'gri22depthi22ext0.1wdepthext':
		if NS == 'S':
			fac = 45076/53693.
		if NS == 'N':
			fac = 68168/71576.0
	if wm == 'gri22depthi22ext0.08wdepthext':
		if NS == 'S':
			fac = 41510/53693.
		if NS == 'N':
			fac = 67934/71576.0
				
	if wm == 'gri22depthi22.2ext0.15wdepthext':
		if NS == 'S':
			fac = 38511/53693.
		if NS == 'N':
			fac = 68488/71576.0
	if wm == 'gri22depthi22.1ext0.15wdepthext':
		if NS == 'S':
			fac = 43360/53693.
		if NS == 'N':
			fac = 68488/71576.0
	if wm == 'depthi22' or wm == 'depthi22wdepthext' or wm == 'depthi22znudge':
		if NS == 'S':
			fac = 49412/53693.
		if NS == 'N':
			fac = 71348/71576.0

	if wm == 'depthexti22' or wm == 'depthexti22wdepthext':
		if NS == 'S':
			fac = 45012/53693.
		if NS == 'N':
			fac = 70936/71576.0

	if wm == 'depthi22.1' or wm == 'depthi22.1wdepthext':
		if NS == 'S':
			fac = 45195/53693.
		if NS == 'N':
			fac = 71348/71576.0
	if wm == 'wmin0.9' or wm == 'wmin0.9wdepthext':
		if NS == 'S':
			fac = 43954/53693.
		if NS == 'N':
			fac = 69055/71576.0

	if wm == 'depthi21.9' or wm == 'depthi21.9wdepthext':
		if NS == 'S':
			fac = 52041/53693.
		if NS == 'N':
			fac = 70412/71576.0

					
	d = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_'+NS+v+'_mz'+str(zmin)+'xz2.2'+wm+muwd+bs+'.dat').transpose()
	dnc = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_'+NS+v+'_mz'+str(zmin)+'xz2.2'+muwd+bs+'.dat').transpose()
	if covm == 'EZ':
		ave = np.loadtxt(ebossdir+'xiave'+mom+'EZPZ'+NS+'GC'+muw+bs+'.dat').transpose()
		cov = np.loadtxt(ebossdir+'covEZPZ'+mom+NS+'GC'+muw+bs+'.dat')[mini:maxi,mini:maxi]*norm**4.
	if covm == 'qpm':
		NSM = 'sgc'
		if NS == 'N':
			NSM = 'ngc'
		ave = np.loadtxt(ebossdir+'xi'+mom+'gaveqpm_qso'+NSM+'v1.6mz0.8xz2.2'+bs+'.dat').transpose()
		cov = np.loadtxt(ebossdir+'cov'+mom+'qpm_qsov1.6'+NSM+bs+'.dat')[mini:maxi,mini:maxi]*norm**4.

	diff = ave[1][mini:maxi]*norm**2.-d[1][mini:maxi]
	diffnc = ave[1][mini:maxi]*norm**2.-dnc[1][mini:maxi]
	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)*fac
	chinc = np.dot(np.dot(diffnc,np.linalg.pinv(cov)),diffnc)
	print chi,chinc
	plt.plot(ave[0][mini:maxi],ave[0][mini:maxi]**2.*ave[1][mini:maxi]*norm**2.,'k--')
	plt.errorbar(d[0][mini:maxi],d[0][mini:maxi]**2.*d[1][mini:maxi],d[0][mini:maxi]**2.*ave[2][mini:maxi]*norm**2.,fmt='ko')
	plt.plot(d[0][mini:maxi],d[0][mini:maxi]**2.*dnc[1][mini:maxi],'rd')
	plt.xlim(d[0][mini]-2.,d[0][maxi]+2.)
	if mom == '0':
		plt.ylim(-80,100)
	else:
		plt.ylim(-200,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{'+mom+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == '0':
		plt.text(30,90,NS+'GC, '+r'after cut $\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='k')
	else:
		plt.text(30,180,NS+'GC, '+r'after cut $\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='k')
	if mom == '0':
		plt.text(30,82,NS+'GC, '+r'before cut $\chi^2$/dof ='+str(chinc)[:4]+'/'+str(maxi-mini),color='r')
	else:
		plt.text(30,164,NS+'GC, '+r'before cut $\chi^2$/dof ='+str(chinc)[:4]+'/'+str(maxi-mini),color='r')
	
	#plt.text(30,160,'Combined',color='k')
	#plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()

	return True

def plotxiQSONScomp(mom='0',bs='8st0',v='v1.6',mini=3,maxi=25,wm='gri22depthi22ext0.15wdepthext',mumin=0,mumax=1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	muw = ''
	muwd = ''
	if mumin != 0:
		muw += 'mmu'+str(mumin)
		muwd += 'mum'+str(mumin)
	if mumax != 1:
		muw += 'xmu'+str(mumax)	
		muwd += 'mux'+str(mumax)
	norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSONScomp'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	if wm == 'gri22depthi22' or wm == 'gri22depthi22wdepthimag' or wm == 'gri22depthi22ext0.15wdepthimagext' or wm == 'gri22depthi22ext0.15wdepthext':
		facs = 47494/53693.
		facn = 68488/71576.0
	if wm == 'gri22depthi22ext0.12wdepthext':
		facs = 46798/53693.
		facn = 68223/71576.0
	if wm == 'gri22depthi22ext0.1wdepthext':
		facs = 45076/53693.
		facn = 68168/71576.0
	if wm == 'gri22depthi22ext0.08wdepthext':
		facs = 41510/53693.
		facn = 67934/71576.0
				
	aves = np.loadtxt(ebossdir+'xiave'+mom+'EZPZ'+'SGC'+muw+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xiave'+mom+'EZPZ'+'NGC'+muw+bs+'.dat').transpose()				
	ds = (np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.9xz2.2'+wm+muwd+bs+'.dat').transpose()+np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.9xz2.2'+wm+muwd+bs+'_1.dat').transpose())/2.
	dn = (np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.9xz2.2'+wm+muwd+bs+'.dat').transpose()+np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.9xz2.2'+wm+muwd+bs+'_1.dat').transpose())/2.
	cov = np.loadtxt(ebossdir+'covEZPZ'+mom+'NGC'+muw+bs+'.dat')[mini:maxi,mini:maxi]/facn+np.loadtxt(ebossdir+'covEZPZ'+mom+'SGC'+muw+bs+'.dat')[mini:maxi,mini:maxi]/facs
	diff = ds[1][mini:maxi]-dn[1][mini:maxi]
	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)
	print chi
	#plt.plot(ave[0][mini:maxi],ave[0][mini:maxi]**2.*ave[1][mini:maxi]*norm**2.,'k--')
	plt.errorbar(ds[0][mini:maxi]-.5,ds[0][mini:maxi]**2.*ds[1][mini:maxi],ds[0][mini:maxi]**2.*aves[2][mini:maxi]*norm**2.,fmt='bo')
	plt.errorbar(dn[0][mini:maxi]+.5,dn[0][mini:maxi]**2.*dn[1][mini:maxi],dn[0][mini:maxi]**2.*aven[2][mini:maxi]*norm**2.,fmt='rd')
	plt.xlim(dn[0][mini]-2.,dn[0][maxi]+2.)
	if mom == '0':
		plt.ylim(-80,100)
	else:
		plt.ylim(-200,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{'+mom+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == '0':
		plt.text(30,90,'NGC/SGC '+r'$\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='k')
	else:
		plt.text(30,180,'NGC/SGC '+r'$\chi^2$/dof ='+str(chi)[:4]+'/'+str(maxi-mini),color='k')
	
	#plt.text(30,160,'Combined',color='k')
	#plt.title(r'Correlation function of v0.7 quasars, 0.9 < z < 2.2')
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
	plt.plot(dn1[0],dn1[0]**2.*t3,'y-o')
	plt.plot(dt[0],dt[0]**2.*dt[1]*1.3,'k:')
	plt.xlim(20,250)
	plt.ylim(-49,150)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,140,'w(depth,g)',color='b')
	plt.text(30,130,'w(depth)',color='r')
	plt.text(30,120,'w(depth,g); g<22',color='y')
	plt.title(r'Correlation function of quasars, v1.3, 0.9 < z < 2.2')
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

def plotxiQSONSbaofit(bs='8st0',v='v1.62',a='',wm='',maxi=25):
	#Plots comparison between QSO clustering and best-fit BAO theory
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONSbaofit'+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	if bs == '10st0':
		bsc = '10.0'
	if bs == '5st0':
		bsc = '5.0'
	dt = np.loadtxt(ebossdir+'ximodQSO'+v+a+'.dat').transpose()
	dtn = np.loadtxt(ebossdir+'ximodQSO'+v+a+'nobao.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov0EZmock_QSOv1.6sgc'+bs+'.dat')#[mini:maxi,mini:maxi]
	covn = np.loadtxt(ebossdir+'cov0EZmock_QSOv1.6ngc'+bs+'.dat')#[mini:maxi,mini:maxi]
	covi = np.linalg.pinv(covn)+np.linalg.pinv(covs)
	cov = np.linalg.pinv(covi)
	#cov = np.loadtxt('covxiNSQSO'+v+'mz0.9xz2.2'+bsc+'.dat')
	et = []
	for i in range(0,maxi):
		et.append(sqrt(cov[i][i]))
	#plt.fill_between(dt[0],dt[0]**2.*(dt[1]*2.-et),dt[0]**2.*(dt[1]*2.+et),color='0.75')
	#plt.plot(dt[0],dt[0]**2.*dt[1]*2.,'k:')
	#if v != 'v1.0':
	#	dsw = np.loadtxt('xi0gebossQSO_S'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	#	dnw = np.loadtxt('xi0gebossQSO_N'+v+'_mz0.9xz2.2wdepth'+bs+'.dat').transpose()
	#else:
	dsw = np.loadtxt(ebossdir+'xi0gebossQSO_S'+v+'_mz0.8xz2.2'+wm+bs+'.dat').transpose()
	dnw = np.loadtxt(ebossdir+'xi0gebossQSO_N'+v+'_mz0.8xz2.2'+wm+bs+'.dat').transpose()
	#fac = int(10/float(bsc))	
	wt = (dsw[1]*.537+dnw[1]*.716)/(.537+.716)
	plt.errorbar(dnw[0][:maxi],dnw[0][:maxi]**2.*wt[:maxi],dnw[0][:maxi]**2.*et,fmt='ko')
	plt.plot(dt[0],dt[0]**2.*dt[1],'k-')
	plt.plot(dt[0],dt[0]**2.*dtn[1],'k--')
	plt.xlim(20,190)
	plt.ylim(-35,80)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,71,r'$\alpha=1.044\pm0.042$',color='k',size=16)
	#plt.text(35,64,r'$\chi^2$/dof = 5.8/7',color='k',size=16)
	plt.text(30,71,r'$\alpha=0.986\pm0.040$',color='k',size=16)
	plt.text(35,64,r'$\chi^2$/dof = 11.4/13',color='k',size=16)
	#plt.title(r'BAO best-fit for v1.6 eboss QSOs, 0.9 < z < 2.2')
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

def plotvssys(sampl,NS,ver,sys,sysmin,sysmax,res,zmin,zmax,wm='',xlab=''):
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	from optimize import fmin
	pp = PdfPages(ebossdir+'n'+sampl+NS+ver+'v'+sys+wm+'.pdf')
	plt.clf()
	plt.minorticks_on()
	d = np.loadtxt(ebossdir+'n'+'geboss'+sampl+'_'+NS+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+str(res)+'v'+sys+'.dat').transpose()
	chin = sum((d[1]-1.)**2./d[2]**2.)
	print chin
	lf = linfit(d[0],d[1],d[2])
	inl = np.array([1.,0])
	b,m = fmin(lf.chilin,inl)
	chilin = sum((d[1]-(m*d[0]+b))**2./d[2]**2.)
	print chilin
	plt.errorbar(d[0],d[1],d[2],fmt='ko')
	ol = np.ones((len(d[0])))
	plt.plot(d[0],ol,'k:')
	plt.plot(d[0],m*d[0]+b,'k--')
	if xlab == '':
		plt.xlabel(sys,size=16)
	else:
		plt.xlabel(xlab,size=16)
		
	plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
	plt.ylim(.7,1.19)
	plt.text(min(d[0])+0.1*(max(d[0])-min(d[0])),1.1,r'$\chi^2$ null ='+str(chin)[:4],color='k')
	plt.text(min(d[0])+0.1*(max(d[0])-min(d[0])),1.08,r'$\chi^2$ lin ='+str(chilin)[:4],color='k')
	#plt.title(r'galaxy density vs. $i$-band depth for v0.7 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True

def plotvssys_simp(xl,yl,el,sys):
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	from optimize import fmin
	plt.clf()
	plt.minorticks_on()
	chin = sum((yl-1.)**2./el**2.)
	print chin
	lf = linfit(xl,yl,el)
	inl = np.array([1.,0])
	b,m = fmin(lf.chilin,inl)
	chilin = sum((yl-(m*xl+b))**2./el**2.)
	print chilin
	plt.errorbar(xl,yl,el,fmt='ko')
	ol = np.ones((len(el)))
	plt.plot(xl,ol,'k:')
	plt.plot(xl,m*xl+b,'k--')
	#if xlab == '':
	plt.xlabel(sys,size=16)
	#else:
	#	plt.xlabel(xlab,size=16)
		
	plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
	plt.ylim(.7,1.19)
	plt.text(min(xl)+0.1*(max(xl)-min(xl)),1.1,r'$\chi^2$ null ='+str(chin)[:4],color='k')
	plt.text(min(xl)+0.1*(max(xl)-min(xl)),1.08,r'$\chi^2$ lin ='+str(chilin)[:4],color='k')
	plt.show()
	#plt.title(r'galaxy density vs. $i$-band depth for v0.7 eboss QSOs, 0.9 < z < 2.2')
	#pp.savefig()
	#pp.close()
	return True


def sysplotsQSO6pan():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	from pylab import *
	from numpy import loadtxt as load
	import numpy as np
	pp = PdfPages(ebossdir+'nQSOdr14vsys.pdf')
	
	fig = plt.figure(figsize=(8.5,7))
	lskn = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2nosysvSKY_FLUX3.dat').transpose()
	lsks = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2nosysvSKY_FLUX3.dat').transpose()
	lskt = (lskn[-2]/lskn[-1]**2.+lsks[-2]/lsks[-1]**2.)/(1./lskn[-1]**2.+1./lsks[-1]**2.)
	lske = (1./(1./lskn[-1]**2.+1./lsks[-1]**2.))**.5
	chi = 0
	for i in range(0,len(lskt)):
		chi += (lskt[i]-1.)**2./lske[i]**2.
	print chi	
 	cskn = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2wgdepthextextvSKY_FLUX3.dat').transpose()
 	csks = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2wgdepthextextvSKY_FLUX3.dat').transpose()
 	cskt = (cskn[-2]/cskn[-1]**2.+csks[-2]/csks[-1]**2.)/(1./cskn[-1]**2.+1./csks[-1]**2.)
 	cske = (1./(1./cskn[-1]**2.+1./csks[-1]**2.))**.5
 	chi = 0
 	for i in range(0,len(cskt)):
 		chi += (cskt[i]-1.)**2./cske[i]**2.
 	print chi	

	leen = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2nosysvEB_MINUS_V-1.dat').transpose()
	les = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2nosysvEB_MINUS_V-1.dat').transpose()
	cen = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2wgdepthextextvEB_MINUS_V-1.dat').transpose()
	ces = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2wgdepthextextvEB_MINUS_V-1.dat').transpose()
	let = (leen[1]/leen[2]**2.+les[1]/les[2]**2.)/(1./leen[2]**2.+1./les[2]**2.)
	lee = (1./(1./leen[2]**2.+1./les[2]**2.))**.5
	chi = 0
	for i in range(0,len(let)):
		chi += (let[i]-1.)**2./lee[i]**2.
		#print chi

	print chi
		
	cet = (cen[1]/cen[2]**2.+ces[1]/ces[2]**2.)/(1./cen[2]**2.+1./ces[2]**2.)
	cee = (1./(1./cen[2]**2.+1./ces[2]**2.))**.5
	chi = 0
	#chim = 0
	for i in range(0,len(cet)):
		chi += (cet[i]-1.)**2./cee[i]**2.
		#print chi
	print chi

	lan = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2nosysvAIRMASS-1.dat').transpose()
	las = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2nosysvAIRMASS-1.dat').transpose()
	can = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2wgdepthextextvAIRMASS-1.dat').transpose()
	cas = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2wgdepthextextvAIRMASS-1.dat').transpose()
			
	lat = (lan[1]/lan[2]**2.+las[1]/las[2]**2.)/(1./lan[2]**2.+1./las[2]**2.)
	lae = (1./(1./lan[2]**2.+1./las[2]**2.))**.5
	chi = 0
	for i in range(0,len(lat)):
		chi += (lat[i]-1.)**2./lae[i]**2.
	print chi	
	cat = (can[1]/can[2]**2.+cas[1]/cas[2]**2.)/(1./can[2]**2.+1./cas[2]**2.)
	cae = (1./(1./can[2]**2.+1./cas[2]**2.))**.5
 	chi = 0
	for i in range(0,len(cat)):
		chi += (cat[i]-1.)**2./cae[i]**2.
	print chi	

	lsen = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2nosysvPSF_FWHM3.dat').transpose()
	lses = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2nosysvPSF_FWHM3.dat').transpose()
	csen = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2wgdepthextextvPSF_FWHM3.dat').transpose()
	cses = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2wgdepthextextvPSF_FWHM3.dat').transpose()
			
	lset = (lsen[1]/lsen[2]**2.+lses[1]/lses[2]**2.)/(1./lsen[2]**2.+1./lses[2]**2.)
	lsee = (1./(1./lsen[2]**2.+1./lses[2]**2.))**.5
	chi = 0
	for i in range(0,len(lset)):
		chi += (lset[i]-1.)**2./lsee[i]**2.
	print chi	
	cset = (csen[1]/csen[2]**2.+cses[1]/cses[2]**2.)/(1./csen[2]**2.+1./cses[2]**2.)
	csee = (1./(1./csen[2]**2.+1./cses[2]**2.))**.5
 	chi = 0
	for i in range(0,len(cat)):
		chi += (cset[i]-1.)**2./csee[i]**2.
	print chi	

	lstn = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2nosys256vstar.dat').transpose()
	lsts = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2nosys256vstar.dat').transpose()
	cstn = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2wgdepthextext256vstar.dat').transpose()
	csts = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2wgdepthextext256vstar.dat').transpose()
			
	lstt = (lstn[1]/lstn[2]**2.+lsts[1]/lsts[2]**2.)/(1./lstn[2]**2.+1./lsts[2]**2.)
	lste = (1./(1./lstn[2]**2.+1./lsts[2]**2.))**.5
	chi = 0
	for i in range(0,len(lstt)):
		chi += (lstt[i]-1.)**2./lste[i]**2.
	print chi	
	cstt = (cstn[1]/cstn[2]**2.+csts[1]/csts[2]**2.)/(1./cstn[2]**2.+1./csts[2]**2.)
	cste = (1./(1./cstn[2]**2.+1./csts[2]**2.))**.5
 	chi = 0
	for i in range(0,len(cstt)):
		chi += (cstt[i]-1.)**2./cste[i]**2.
	print chi	

	ldn = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2nosysvIMAGE_DEPTH_EXT3.dat').transpose()
	lds = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2nosysvIMAGE_DEPTH_EXT3.dat').transpose()
	cdn = load(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz2.2wgdepthextextvIMAGE_DEPTH_EXT3.dat').transpose()
	cds = load(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz2.2wgdepthextextvIMAGE_DEPTH_EXT3.dat').transpose()
			
	ldt = (ldn[1]/ldn[2]**2.+lds[1]/lds[2]**2.)/(1./ldn[2]**2.+1./lds[2]**2.)
	lde = (1./(1./ldn[2]**2.+1./lds[2]**2.))**.5
	chi = 0
	for i in range(0,len(ldt)):
		chi += (ldt[i]-1.)**2./lde[i]**2.
	print chi	
	cdt = (cdn[1]/cdn[2]**2.+cds[1]/cds[2]**2.)/(1./cdn[2]**2.+1./cds[2]**2.)
	cde = (1./(1./cdn[2]**2.+1./cds[2]**2.))**.5
 	chi = 0
	for i in range(0,len(cdt)):
		chi += (cdt[i]-1.)**2./cde[i]**2.
	print chi	

	
	ols = ones((len(lan[0])))
	#plt.ylabel('number density / average number density',size=13)
	ax = fig.add_subplot(2,3,1)
	ax.set_xlim(4,21)
	ax.set_ylim(.76,1.19)
	ax.minorticks_on()
# 	for tick in ax.xaxis.get_major_ticks():
# 		tick.label.set_fontsize(8)
# 	for tick in ax.yaxis.get_major_ticks():
# 		tick.label.set_fontsize(8)
	#ax.set_ylabel('number density / average number density',size=13)
	ax.set_xlabel('i-band sky background',size=13)
	ax.plot(lskn[0],ols,'k:')
	ax.plot(lskn[0]+.07,lskt,'--',color='steelblue') 	
	ax.errorbar(cskn[0]-.07,cskt,cske,fmt='s',color='firebrick',markersize=7,elinewidth=2)
	ax2 = fig.add_subplot(2,3,2,sharey=ax)
	ax2.minorticks_on()
# 	for tick in ax2.xaxis.get_major_ticks():
# 		tick.label.set_fontsize(8)
	ax2.plot(leen[0],ols,'k:')
 	ax2.plot(leen[0]+.001,let,'--',color='steelblue') 
 	ax2.errorbar(cen[0]-.001,cet,cee,fmt='s',color='firebrick',markersize=7,elinewidth=2)
  	ax2.set_xlabel('E(B-V)',size=14)
 	ax2.set_xlim(0.001,.152)
	ax2.set_ylim(.751,1.2)
	ax2.xaxis.set_ticks(np.arange(0.02,.152,.04))
	#setp( ax2.get_yticklabels(), visible=False)
	for ylabel_2 in ax2.axes.get_yticklabels():
		ylabel_2.set_visible(False)

	ax3 = fig.add_subplot(2,3,3,sharey=ax)
	ax3.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax3.set_xlabel('airmass',size=13)
	ax3.plot(lan[0],ols,'k:')
  	ax3.plot(lan[0]+.007,lat,'--',color='steelblue')
  	ax3.errorbar(can[0]-.007,cat,cae,fmt='s',color='firebrick',markersize=7,elinewidth=2)
	ax3.set_xlim(1.01,2)
	ax3.set_ylim(.751,1.2)
	for ylabel_i in ax3.axes.get_yticklabels():
		ylabel_i.set_visible(False)

	ax4 = fig.add_subplot(2,3,4)
	ax4.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax4.set_xlabel('i-band seeing',size=13)
	ax4.plot(lsen[0],ols,'k:')
  	ax4.plot(lsen[0]+.007,lset,'--',color='steelblue')
  	ax4.errorbar(csen[0]-.007,cset,csee,fmt='s',color='firebrick',markersize=7,elinewidth=2)
	ax4.set_xlim(0.7,1.99)
	ax4.set_ylim(.751,1.2)
	ax4.set_ylabel(r'                                          number density / average number density',size=13)
	ax5 = fig.add_subplot(2,3,5,sharey=ax4)
	ax5.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax5.set_xlabel(r'$N_{\rm star}/Nside=256$',size=13)
	ax5.plot(lstn[0],ols,'k:')
  	ax5.plot(lstn[0]+.007,lstt,'--',color='steelblue')
  	ax5.errorbar(cstn[0]-.007,cstt,cste,fmt='s',color='firebrick',markersize=7,elinewidth=2)
	ax5.set_xlim(31,299)
	ax5.set_ylim(.751,1.2)
	for ylabel_5 in ax5.axes.get_yticklabels():
		ylabel_5.set_visible(False)

	ax6 = fig.add_subplot(2,3,6,sharey=ax4)
	ax6.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax6.set_xlabel(r'$i$-band depth',size=13)
	ax6.plot(ldn[0],ols,'k:')
  	ax6.plot(ldn[0]+.007,ldt,'--',color='steelblue')
  	ax6.errorbar(cdn[0]-.007,cdt,cde,fmt='s',color='firebrick',markersize=7,elinewidth=2)
	ax6.set_xlim(21.81,22.9)
	ax6.set_ylim(.751,1.2)
	for ylabel_6 in ax6.axes.get_yticklabels():
		ylabel_6.set_visible(False)
	
	
	fig.subplots_adjust(wspace=0)
 	xl = [6]
 	yl = [1.14]
 	#ax.plot(xl,yl,'s',color='firebrick')
 	yl = [1.1]
 	#ax.plot(xl,yl,'d',color='steelblue')
 	#ax.text(7,1.13,'CMASS',fontsize=16,color='firebrick')
 	#ax.text(7,1.09,'LOWZ',fontsize=16,color='steelblue')
	cnr = [5,13.5,1.08,1.165]
	xl = [cnr[0],cnr[1]]
	yl = [cnr[2],cnr[2]]
	#ax.plot(xl,yl,'k-')
	yl = [cnr[3],cnr[3]]
	#ax.plot(xl,yl,'k-')
	xl = [cnr[0],cnr[0]]
	yl = [cnr[2],cnr[3]]
	#ax.plot(xl,yl,'k-')
	xl = [cnr[1],cnr[1]]
	#ax.plot(xl,yl,'k-')

	pp.savefig()
	pp.close()
	return True


def plotQSOgmagNSvsdepth(v='v1.5'):
	#plots N_QSO vs. i-band depth
	from optimize import fmin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'nQSONS'+v+'vdepth.pdf')
	plt.clf()
	plt.minorticks_on()
	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2nosysgm020.5512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2nosysgm020.5512vdepth.dat').transpose()
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
	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2nosysgm20.521.0512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2nosysgm20.521.0512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='rd')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'r--')

	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2nosysgm21.021.5512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2nosysgm21.021.5512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='bs')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'b--')

	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2nosysgm21.522.0512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2nosysgm21.522.0512vdepth.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0],dt,e,fmt='g^')
	lf = linfit(ds[0],dt,e)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(ds[0],b0+m0*ds[0],'g--')
	ds = np.loadtxt(ebossdir+'ngebossQSO_S'+v+'_mz0.9xz2.2nosysgm22.030512vdepth.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSO_N'+v+'_mz0.9xz2.2nosysgm22.030512vdepth.dat').transpose()
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
	plt.title(r'galaxy density vs. $i$-band depth for v1.5 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return bl,ml

def plotQSOgmagNSvsredshift(v='v1.6',wm='wgdepthextext'):
	#plots N_QSO vs. i-band depth
	from optimize import fmin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'nQSONS'+v+wm+'vredshift.pdf')
	plt.clf()
	plt.minorticks_on()
	
	
	ds = np.loadtxt(ebossdir+'ngebossQSOsys_Sv1.6_mz0.8xz1.15'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSOsys_Nv1.6_mz0.8xz1.15'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0]-.015,dt,e,fmt='ko')
	lf = linfit(ds[0],dt,e)
	chi0 = lf.chilin((1,0))
	
	ol = np.ones((len(ds[0])))
	plt.plot(ds[0],ol,'k--')
	
	ds = np.loadtxt(ebossdir+'ngebossQSOsys_Sv1.6_mz1.15xz1.5'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSOsys_Nv1.6_mz1.15xz1.5'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0]-.005,dt,e,fmt='rd')
	lf = linfit(ds[0],dt,e)
	chi1 = lf.chilin((1,0))

	ds = np.loadtxt(ebossdir+'ngebossQSOsys_Sv1.6_mz1.5xz1.85'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSOsys_Nv1.6_mz1.5xz1.85'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0]+.005,dt,e,fmt='bs')
	lf = linfit(ds[0],dt,e)
	chi2 = lf.chilin((1,0))
	
	ds = np.loadtxt(ebossdir+'ngebossQSOsys_Sv1.6_mz1.85xz2.2'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dn = np.loadtxt(ebossdir+'ngebossQSOsys_Nv1.6_mz1.85xz2.2'+wm+'vIMAGE_DEPTH_EXT3.dat').transpose()
	dt = (ds[1]/ds[2]**2.+dn[1]/dn[2]**2.)/(1./ds[2]**2.+1./dn[2]**2.)
	e = (1./(1./ds[2]**2.+1./dn[2]**2.))**.5	
	plt.errorbar(ds[0]+.015,dt,e,fmt='g^')
	lf = linfit(ds[0],dt,e)
	chi3 = lf.chilin((1,0))
	plt.text(22.1,.87,r'$0.8 < z < 1.15, \chi^2=$'+str(chi0)[:3],fontsize=18,color='k')
	plt.text(22.1,.85,r'$1.15 < z < 1.5, \chi^2=$'+str(chi1)[:3],fontsize=18,color='r')
	plt.text(22.1,.83,r'$1.5 < z < 1.85, \chi^2=$'+str(chi2)[:3],fontsize=18,color='b')
	plt.text(22.1,.81,r'$1.85 < z < 2.2, \chi^2=$'+str(chi3)[:3],fontsize=18,color='g')
	plt.xlabel(r'5$\sigma$ $i$-band depth (magnitudes)',size=16)
	plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
	plt.ylim(.8,1.2)
	#plt.title(r'galaxy density vs. $i$-band depth for v1.6 eboss QSOs')
	pp.savefig()
	pp.close()
	return True


def plotQSOgmagNSregvsdepth(NS,v='v1.5'):
	#plots N_QSO vs. i-band depth
	from optimize import fmin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'nQSO'+NS+v+'vdepth.pdf')
	plt.clf()
	plt.minorticks_on()
	dt = np.loadtxt(ebossdir+'ngebossQSO_'+NS+v+'_mz0.9xz2.2nosysgm020.5512vdepth.dat').transpose()
	plt.errorbar(dt[0],dt[1],dt[2],fmt='ko')
	lf = linfit(dt[0],dt[1],dt[2])
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl = []
	ml = []
	bl.append(b0)
	ml.append(m0)
	plt.plot(dt[0],b0+m0*dt[0],'k--')
	dt = np.loadtxt(ebossdir+'ngebossQSO_'+NS+v+'_mz0.9xz2.2nosysgm20.521.0512vdepth.dat').transpose()
	plt.errorbar(dt[0],dt[1],dt[2],fmt='rd')
	lf = linfit(dt[0],dt[1],dt[2])
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(dt[0],b0+m0*dt[0],'r--')

	dt = np.loadtxt(ebossdir+'ngebossQSO_'+NS+v+'_mz0.9xz2.2nosysgm21.021.5512vdepth.dat').transpose()
	plt.errorbar(dt[0],dt[1],dt[2],fmt='bs')
	lf = linfit(dt[0],dt[1],dt[2])
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(dt[0],b0+m0*dt[0],'b--')

	dt = np.loadtxt(ebossdir+'ngebossQSO_'+NS+v+'_mz0.9xz2.2nosysgm21.522.0512vdepth.dat').transpose()
	plt.errorbar(dt[0],dt[1],dt[2],fmt='g^')
	lf = linfit(dt[0],dt[1],dt[2])
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(dt[0],b0+m0*dt[0],'g--')
	dt = np.loadtxt(ebossdir+'ngebossQSO_'+NS+v+'_mz0.9xz2.2nosysgm22.030512vdepth.dat').transpose()
	plt.errorbar(dt[0],dt[1],dt[2],fmt='^',color='purple')
	lf = linfit(dt[0],dt[1],dt[2])
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	bl.append(b0)
	ml.append(m0)
	plt.plot(dt[0],b0+m0*dt[0],'--',color='purple')
	plt.text(22.6,.81,r'$g < 20.5$',fontsize=18,color='k')
	plt.text(22.6,.77,r'$20.5 < g < 21$',fontsize=18,color='r')
	plt.text(22.6,.73,r'$21 < g < 21.5$',fontsize=18,color='b')
	plt.text(22.6,.69,r'$21. 5 < g < 22$',fontsize=18,color='g')
	plt.text(22.6,.65,r'$g > 22$',fontsize=18,color='purple')
	plt.xlabel(r'5$\sigma$ $i$-band depth, averaged in Nside=512 pixels (magnitudes)',size=16)
	plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
	plt.ylim(.6,1.4)
	plt.title(r'galaxy density vs. $i$-band depth for v1.5 eboss ' +NS+'GC QSOs, 0.9 < z < 2.2')
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

def plotQSONSbaolike(v='v1.0',title='QSOs'):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONSbaolik'+v+'.pdf')
	plt.clf()
	plt.minorticks_on()
	db = np.loadtxt(ebossdir+'BAOxichilQSO'+v+'0.48st0.dat').transpose()
	dnb = np.loadtxt(ebossdir+'BAOxichilQSO'+v+'nobao0.48st0.dat').transpose()
	chim = min(db[1])
	ol = np.ones((len(db[0])))
	plt.plot(db[0],db[1]-chim,'k-',linewidth=4)
	plt.plot(db[0],dnb[1]-chim,'k--',linewidth=3)
	plt.plot(db[0],ol,'k:',linewidth=1)
	plt.text(0.825,1.1,r'$1\sigma$')
	plt.plot(db[0],ol*4,'k:',linewidth=1)
	plt.text(0.825,4.1,r'$2\sigma$')
	plt.plot(db[0],ol*9,'k:',linewidth=1)
	plt.text(0.825,9.1,r'$3\sigma$')
	#plt.xlim(20,165)
	plt.ylim(0,16)
	plt.xlabel(r'$\alpha_{\rm BAO}$',size=18)
	plt.ylabel(r'$\Delta\chi^2$',size=18)
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	#plt.title(r'BAO likelihood for '+title)
	pp.savefig()
	pp.close()
	return True
	
def plot3Dphotz(b=1.,mumax=.8,pimax=1000.):
	from matplotlib import pyplot as plt
	d = np.loadtxt('/Users/ashleyross/eBOSS/LRG_3d_ps_30t_1016.out').transpose()
	nb = int(len(d[0])/200)
	xil = []
	rl = []
	rpold = d[1][0]
	xi = 0
	n = 0 
	for i in range(0,len(d[0])):	
		pi = d[0][i]
		
		if pi < pimax:
			rp = d[1][i]
			mu = pi/sqrt(pi**2.+rp**2.)
			if mu < mumax:
				if rp == rpold:
					xi += d[2][i]
					n += 1.
				else:
					#print mu
					if n != 0:
						xil.append(xi/n)
						rl.append(rpold)
					xi = d[2][i]
					n = 0	
			rpold = rp
	if n != 0:
		xil.append(xi/n)
		rl.append(rp)
	print rl
	rl = np.array(rl)
	xil = np.array(xil)
	dt = np.loadtxt('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumax'+str(mumax)+'0.406.010.0combzsiglsp1.0.dat').transpose()
	plt.plot(rl,rl**2.*xil,'ko-',dt[0],dt[0]**2.*dt[1]*b,'r-')
	plt.xlim(0,160)
	plt.xlabel(r'$r_{\perp}$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$r^2_{\perp}\xi$ ($h^{-2}$Mpc$^2$)')
	plt.title(r'$\mu <$ '+str(mumax))
	plt.show()
	
	return True

def plotumgmap(sampl='QSO',NS='S',ver='v1.6',zmin=.8,zmax=2.2,ramin=-180,ramax=180,decmin=-90,decmax=90,res=64,wm='',gri22='',smin=-5,smax=5,size=1):
	from healpix import healpix,radec2thphi,thphi2radec
	h = healpix()

	npo = 12*res**2
	pixlg = []
	pixlug = []
	pixlr = []
	ral = []
	decl = []

	for i in range(0,npo):
		pixlg.append(0)
		pixlug.append(0)
		pixlr.append(0)
		th,phi = h.pix2ang_nest(res,i)
		ra,dec = thphi2radec(th,phi)
		if NS == 'S':
			if ra > 180:
				ra = ra -360
		ral.append(ra)
		decl.append(dec)


	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.ran.fits') #read galaxy/quasar file
	nr = 0
	for i in range (0,len(f)):
		ra,dec = f[i]['RA'],f[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(res,th,phi))
		pixlr[p] += 1.
		nr += 1.
	print nr
	sumr = 0
	sump = 0
	for i in range(0,len(pixlr)):
		if pixlr[i] > 0:
			sumr += pixlr[i]
			sump += 1.
	print sumr/sump
	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.dat.fits') #read galaxy/quasar file
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
		c = 1
		if gri22 == 'gri22':
			if gm > 22 or rm > 22 or im > 22:
				c = 0

		if z > zmin and z < zmax and c == 1 and gc:
			no += 1
			#w = 1.
			#if wm == '':
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT'] #weight in current catalog
			if wm == 'nosys':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP'] #compare to result with no systematic weights
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
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws				
			if wm == 'wdepthgmag':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = w*ws
			if wm == 'wdepthimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if im < 20.25:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 19.75*slp
					m = im*slp+b
				if im >= 20.25 and im < 20.75:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.25*slp
					m = im*slp+b
				if im >= 20.75 and im < 21.25:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 20.75*slp
					m = im*slp+b
				if gm >= 21.25:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.25*slp
					m = im*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				

			if wm == 'wdepthgmagext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				pix2 = h.ang2pix_nest(256,th,phi)
				ne = wext[pix2]
				we = 1./(be+me*ne)
				w = w*ws*we
			pixlg[p] += 1.*w
			pixlug[p] += 1.*w*(um-gm)
			zm += w*z
			nt += w
	pug = []
	ralug = []
	declug = []
	no = 0
	nor = 0
	for i in range(0,len(pixlg)):
		if pixlr[i] > .1*sumr/sump:#100.:
			if pixlg[i] == 0:
				print i,pixlr[i],ral[i],decl[i]
				pug.append(-.2)
				ralug.append(ral[i])
				declug.append(decl[i])
			else:
				pug.append(pixlug[i]/pixlg[i])
				ralug.append(ral[i])
				declug.append(decl[i])
				if pixlug[i]/pixlg[i] < 0.06 or pixlug[i]/pixlg[i] > .32:
					no += pixlg[i]
					nor += pixlr[i]		
		#else:
		#	pug.append(-99)
	print no,nor
	if nor > 0:
		print no/nor,sum(pixlg)/sum(pixlr)	
	from matplotlib import pyplot as plt
	import matplotlib.cm as cm
	map = plt.scatter(ralug,declug,c=pug,s=size,cmap=cm.rainbow,lw=0,edgecolors='none',vmin=smin,vmax=smax)
	cbar = plt.colorbar(map)
	#plt.xlim(ramin,ramax)
	#plt.ylim(decmin,decmax)
	plt.show()
	return True

def plotumgmapspix(sampl='QSO',NS='S',ver='v1.6',zmin=.8,zmax=2.2,ramin=-180,ramax=180,decmin=-90,decmax=90,res=16,wm='',gri22='',smin=-5,smax=5,size=1):
	from stripeDR8 import radec2le,ang2pix,pix2etalam,le2radec
	#h = healpix()

	npo = 36*13*res**2
	pixlg = []
	pixlug = []
	pixlr = []
	ral = []
	decl = []

	for i in range(0,npo):
		pixlg.append(0)
		pixlug.append(0)
		pixlr.append(0)
		eta,lam = pix2etalam(res,i)
		ra,dec = le2radec(lam,eta)
		if NS == 'S':
			if ra > 180:
				ra = ra -360
		ral.append(ra)
		decl.append(dec)


	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.ran.fits') #read galaxy/quasar file
	nr = 0
	for i in range (0,len(f)):
		ra,dec = f[i]['RA'],f[i]['DEC']
		lam,eta = radec2le(ra,dec)
		p = int(ang2pix(res,lam,eta))
		pixlr[p] += 1.
		nr += 1.
	print nr
	sumr = 0
	sump = 0
	for i in range(0,len(pixlr)):
		if pixlr[i] > 0:
			sumr += pixlr[i]
			sump += 1.
	print sumr/sump
	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.dat.fits') #read galaxy/quasar file
	no = 0
	zm = 0
	nt = 0
	umgl = []
	nugb = 0
	for i in range (0,len(f)):
		z = f[i]['Z']
		gc = True
		um = f[i]['MODELMAG'][0]-f[i]['EXTINCTION'][0]
		gm = f[i]['MODELMAG'][1]-f[i]['EXTINCTION'][1]
		rm = f[i]['MODELMAG'][2]-f[i]['EXTINCTION'][2]
		im = f[i]['MODELMAG'][3]-f[i]['EXTINCTION'][3]
		c = 1
		if gri22 == 'gri22':
			if gm > 22 or rm > 22 or im > 22:
				c = 0

		if z > zmin and z < zmax and c == 1 and gc:
			no += 1
			#w = 1.
			#if wm == '':
			ra,dec = f[i]['RA'],f[i]['DEC']
			lam,eta = radec2le(ra,dec)
			p = int(ang2pix(res,lam,eta))

			#th,phi = radec2thphi(ra,dec)
			#p = int(h.ang2pix_nest(res,th,phi))
			
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT'] #weight in current catalog
			if wm == 'nosys':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP'] #compare to result with no systematic weights
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
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws				
			if wm == 'wdepthgmag':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = w*ws
			if wm == 'wdepthimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if im < 20.25:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 19.75*slp
					m = im*slp+b
				if im >= 20.25 and im < 20.75:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.25*slp
					m = im*slp+b
				if im >= 20.75 and im < 21.25:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 20.75*slp
					m = im*slp+b
				if gm >= 21.25:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.25*slp
					m = im*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				

			if wm == 'wdepthgmagext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				pix2 = h.ang2pix_nest(256,th,phi)
				ne = wext[pix2]
				we = 1./(be+me*ne)
				w = w*ws*we
			if abs(um-gm) < 2.:
				pixlg[p] += 1.*w
				pixlug[p] += 1.*w*(um-gm)
				zm += w*z
				nt += w
				umgl.append(um-gm)
			else:
				nugb += 1.
				print z
	print nugb				
	pug = []
	ralug = []
	declug = []
	no = 0
	nor = 0
	for i in range(0,len(pixlg)):
		if pixlr[i] > .5*sumr/sump:#100.:
			if pixlg[i] == 0:
				print i,pixlr[i],ral[i],decl[i]
				pug.append(-.2)
				ralug.append(ral[i])
				declug.append(decl[i])
			else:
				pug.append(pixlug[i]/pixlg[i])
				ralug.append(ral[i])
				declug.append(decl[i])
				if pixlug[i]/pixlg[i] < 0.0 or pixlug[i]/pixlg[i] > .26:
					no += pixlg[i]
					nor += pixlr[i]		
		#else:
		#	pug.append(-99)
	print no,nor
	if nor > 0:
		print no/nor,sum(pixlg)/sum(pixlr)	
	from matplotlib import pyplot as plt
	import matplotlib.cm as cm
	plt.hist(umgl)
	plt.show()
	map = plt.scatter(ralug,declug,c=pug,s=size,cmap=cm.rainbow,lw=0,edgecolors='none',vmin=smin,vmax=smax)
	cbar = plt.colorbar(map)
	#plt.xlim(ramin,ramax)
	#plt.ylim(decmin,decmax)
	plt.show()
	plt.hist(pug)
	plt.show()
	return True


def plotmeanzmap(sampl='QSO',NS='S',ver='v1.6',zmin=.8,zmax=2.2,ramin=-180,ramax=180,decmin=-90,decmax=90,res=64,wm='',gri22='',smin=-5,smax=5,size=1):
	from healpix import healpix,radec2thphi,thphi2radec
	h = healpix()

	npo = 12*res**2
	pixlg = []
	pixlz = []
	pixlr = []
	ral = []
	decl = []

	for i in range(0,npo):
		pixlg.append(0)
		pixlz.append(0)
		pixlr.append(0)
		th,phi = h.pix2ang_nest(res,i)
		ra,dec = thphi2radec(th,phi)
		if NS == 'S':
			if ra > 180:
				ra = ra -360
		ral.append(ra)
		decl.append(dec)


	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.ran.fits') #read galaxy/quasar file
	nr = 0
	for i in range (0,len(f)):
		ra,dec = f[i]['RA'],f[i]['DEC']
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(res,th,phi))
		pixlr[p] += 1.
		nr += 1.
	print nr
	sumr = 0
	sump = 0
	for i in range(0,len(pixlr)):
		if pixlr[i] > 0:
			sumr += pixlr[i]
			sump += 1.
	print sumr/sump
	f = fitsio.read(dirfits+'eboss_'+ver+'-'+sampl+'-'+NS+'-eboss_'+ver+'.dat.fits') #read galaxy/quasar file
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
		c = 1
		if gri22 == 'gri22':
			if gm > 22 or rm > 22 or im > 22:
				c = 0

		if z > zmin and z < zmax and c == 1 and gc:
			no += 1
			#w = 1.
			#if wm == '':
			ra,dec = f[i]['RA'],f[i]['DEC']
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			
			w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*f[i]['WEIGHT_SYSTOT'] #weight in current catalog
			if wm == 'nosys':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP'] #compare to result with no systematic weights
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
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				ws = 1./(b+m*ns)
				w = w*ws				
			if wm == 'wdepthgmag':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = w*ws
			if wm == 'wdepthimag':
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
				if im < 20.25:
					slp = (ml[1]-ml[0])/.5
					b = ml[0] - 19.75*slp
					m = im*slp+b
				if im >= 20.25 and im < 20.75:
					slp = (ml[2]-ml[1])/.5
					b = ml[1] - 20.25*slp
					m = im*slp+b
				if im >= 20.75 and im < 21.25:
					slp = (ml[3]-ml[2])/.5
					b = ml[2] - 20.75*slp
					m = im*slp+b
				if gm >= 21.25:
					slp = (ml[4]-ml[3])/.5
					b = ml[3] - 21.25*slp
					m = im*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']*ws				

			if wm == 'wdepthgmagext':
				w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
				pix2 = h.ang2pix_nest(512,th,phi)
				ns = wsys[pix2]
					
				if gm < 20.75:
					slp = (slpl[1]-slpl[0])/.5
					b = slpl[0] - 20.25*slp
					m = gm*slp+b
				if gm >= 20.75 and gm < 21.25:
					slp = (slpl[2]-slpl[1])/.5
					b = slpl[1] - 20.75*slp
					m = gm*slp+b
				if gm >= 21.25 and gm < 21.75:
					slp = (slpl[3]-slpl[2])/.5
					b = slpl[2] - 21.25*slp
					m = gm*slp+b
					if m > slpl[3]:
						print m,hm,b,slp
				if gm >= 21.75:
					slp = (slpl[4]-slpl[3])/.5
					b = slpl[3] - 21.75*slp
					m = gm*slp+b
				bw = 1.-node*m
				ws = 1./(bw+m*ns)
				pix2 = h.ang2pix_nest(256,th,phi)
				ne = wext[pix2]
				we = 1./(be+me*ne)
				w = w*ws*we
			pixlg[p] += 1.*w
			pixlz[p] += 1.*w*z
			zm += w*z
			nt += w
	pz = []
	ralz = []
	declz = []
	no = 0
	nor = 0
	for i in range(0,len(pixlg)):
		if pixlr[i] > .1*sumr/sump:#100.:
			if pixlg[i] == 0:
				print i,pixlr[i],ral[i],decl[i]
				pz.append(.8)
				ralug.append(ral[i])
				declug.append(decl[i])
			else:
				pz.append(pixlz[i]/pixlg[i])
				ralz.append(ral[i])
				declz.append(decl[i])
				if pixlz[i]/pixlg[i] < 1.2 or pixlz[i]/pixlg[i] > 1.8:
					no += pixlg[i]
					nor += pixlr[i]		
		#else:
		#	pug.append(-99)
	print no,nor
	if nor > 0:
		print no/nor,sum(pixlg)/sum(pixlr)	
	from matplotlib import pyplot as plt
	import matplotlib.cm as cm
	map = plt.scatter(ralz,declz,c=pz,s=size,cmap=cm.rainbow,lw=0,edgecolors='none')#,vmin=smin,vmax=smax)
	cbar = plt.colorbar(map)
	#plt.xlim(ramin,ramax)
	#plt.ylim(decmin,decmax)
	plt.show()
	return True			
	
def BAOrelPlanck(wo='QSO4SH',xmax=2.,BOSS=False,BOSSDR12=True,MGS=True,wz=True,sdss=False,df6=True,QSODR14=True,LRGDR14=False,des=False,desy1=False,eboss=False,desi=False):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.axes as ax
	import numpy as np
	from Cosmo import *
	from numpy import loadtxt as load

	pp = PdfPages('BAOrelPlanck'+wo+'.pdf')
	plt.clf()
	fig=plt.figure()
	ax0=fig.add_subplot(1,1,1)
	ax.Axes.axis(ax0,'auto')
	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.ylim ( 0.911, 1.09 )
	pe = np.loadtxt('/Users/ashleyross/DR7VAC/DVordPlanck.txt').transpose()
	x = pe[0]
	#x=np.arange(0.02,2.,0.01)
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
 	print wigf
 	#ru = du.rs/df.rs
 	#rd = dd.rs/df.rs
 	#re = (ru-rd)/2.
 	#fo = open('planckdvrserr_om.dat','w')
 	for i in range(0,len(pe[1])):
 		dvrsf = df.dV(x[i])/df.rs
 		dd = delm*(dm.dV(x[i])/dm.rs/dvrsf-1.)/.01
		#fo.write(str(x[i])+' '+str(dd)+'\n')
		yu[i] = y[i]+pe[1][i]#+dd#sqrt(pe[1][i]**2.)#+re**2.)
		yd[i] = y[i]-pe[1][i]#dd#sqrt(pe[1][i]**2.)#+re**2.)
	plt.fill_between(x,yd,yu,color='0.75')
	#plt.fill_between(x,yd2,yu2,color='0.75')
	plt.plot(x,1.0*y,'k-')
	#plt.plot(x,1.0*yu,'k-')
	#plt.plot(x,yd,'k-')
	if df6:
		plt.text(0.11,.95,'6dFGS',fontsize=18,color='green')
		x6 = [0.1]
		y6 = [0.987]
		e6 = [0.045]
		plt.errorbar(x6, y6,e6,fmt='s',markersize=7,elinewidth=1.75,color='green')
	if BOSS:
		plt.text(0.3,.945,'BOSS DR11',fontsize=18)
		xl = [0.32,0.57]
		yl = [0.974,0.983]
		el = [0.020,0.01]
		plt.errorbar(xl, yl,el,fmt='o',markersize=7,elinewidth=1.75,color='k')
	if BOSSDR12:
		plt.text(0.44,.96,'BOSS DR12',fontsize=18,color='purple')
		xl = [0.38,0.61]
		yl = [0.999*.995,0.985*.996]
		el = [0.01,0.01]
		plt.errorbar(xl, yl,el,fmt='d',markersize=7,elinewidth=1.75,color='purple')

	#plt.text(0.26,.9,'LOWZ',fontsize=18)
	#plt.text(0.51,.92,'BOSS',fontsize=18)
	#plt.text(0.5,.9,'CMASS',fontsize=18)
	if MGS:
		plt.text(0.1,1.08,'SDSS MGS',fontsize=18,color='r')
		xl = [0.15]
		yl = [1.04*.9957]
		el = [0.037]
		plt.errorbar(xl,yl,el,fmt='D',markeredgecolor='r',markersize=7,elinewidth=1.75,color='r')
	if wz:
		plt.text(0.5,1.07,'WiggleZ',fontsize=18,color='.5')
		ywl = [1.061*wigf[0],1.065*wigf[1],1.039*wigf[2]]
		xwl = [.44,.6,.73] 
		ewl = [0.048,0.045,0.034]
		plt.errorbar(xwl, ywl,ewl,fmt='s',markersize=6,elinewidth=1.,color='.5')

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
		plt.text(1.2,.95,'DR14 quasars',fontsize=18,color='b')
		plt.text(1.2,.94,'(predicted)',fontsize=18,color='b')
		xl = [1.5]
		yl = [1.0]
		el = [0.039]
		plt.errorbar(xl,yl,el,fmt='o',markeredgecolor='b',markersize=7,elinewidth=1.75,color='b')

	if LRGDR14:
		plt.text(.76,1.025,'DR14 LRGs',fontsize=18,color='b')
		xl = [.75]
		yl = [1.005]
		el = [0.021]
		plt.errorbar(xl,yl,el,fmt='o',markeredgecolor='b',markersize=7,elinewidth=1.75,color='b')

	if des:
		plt.text(.95,1.015,'DES projected',fontsize=18,color='b')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.875,0.925]
		yld = [1.0,1.0]
		eld = [0.019,0.019]
		yu = [1.019,1.019]
		yd = [0.981,0.981]
		#plt.errorbar(xld,yld,eld,fmt='-s',markeredgecolor='b',markersize=7,elinewidth=1.75,color='b')
		plt.fill_between(xld,yd,yu,color='b')
	if desy1:
		plt.text(.85,1.015,'DES Y1',fontsize=18,color='b')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.8]
		yld = [1.02]
		eld = [0.04]
		yu = [1.019,1.019]
		yd = [0.981,0.981]
		plt.errorbar(xld,yld,eld,fmt='-s',markeredgecolor='b',markersize=7,elinewidth=1.75,color='b')
		#plt.fill_between(xld,yd,yu,color='b')
		plt.ylim(.9,1.1)
	if desi:
		plt.text(1.4,1.03,'DESI',fontsize=18,color='purple')		
		d = load('/Users/ashleyross/BAOpredictionscode/BAO_BigBOSS.dat').transpose()
		xl = d[0]
		yl = np.ones((len(d[0])))
		el = d[-1]
		plt.errorbar(xl,yl,el,fmt='-*',markeredgecolor='purple',markersize=7,elinewidth=1.75,color='purple')
	plt.xlim ( 0.0, xmax )
	
	plt.xlabel ('Redshift', fontsize=18)
	#plt.ylabel (r'$(D_{\rm V}/r_{\rm d})/(D_{\rm V}/r_{\rm d})_{\rm Planck}$', fontsize=18)
	plt.ylabel (r'Distance/Distance(Planck$\Lambda$CDM)', fontsize=18)
	pp.savefig()
	pp.close()
	return True


def test():
	return True

if __name__ == '__main__':
	#do bao fit on EZ mocks
	import sys
	nm = str(sys.argv[1])
	xibao('qpm_qso',0.8,2.2,version='v1.6',wm='',bs=8,start=0,rmin=30,rmax=180.,mb='',Bp=.4,v='n',mockn=nm,covmd='QPMmock',Nmock=400)
	
	  