from math import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import fitsio
from xitools_eboss import *

dirsci = '/mnt/lustre/ashleyr/eboss/' #where AJR puts eboss catalogs, change this to wherever you have put catalogs
dirsys = 'maps/' #change to local directory where ngalvsys from wiki was put, note star map and depth map included
dirfits = '/Users/ashleyross/fitsfiles/' #change to where your catalog files are
ebossdir = '/Users/ashleyross/eboss/' #where AJR puts correlation functions, writes out results
dirscio = '/mnt/lustre/ashleyr/eboss/mocks/'

def P2(mu):
	return .5*(3.*mu**2.-1.)
	
def P4(mu):
	return .125*(35.*mu**4.-30.*mu**2.+3.)


def mkjackf_elg(samp,v='v5_10_7',cm='',Njack=20):
	#defines jack-knifes
	ranHealp_elg(samp,v=v)
	mf = open('ranHeal_pix256eboss'+cm+samp+v+'_elg.dat').readlines()
	fo = open('jackhpixeboss'+cm+samp+v+'_elg'+str(Njack)+'.dat','w')
	for i in range(0,Njack-1):
		fo.write(str(mf[(len(mf)/Njack)*(i+1)].split()[0])+'\n')
	fo.close()
	return True

def ranHealp_elg(samp,v='v5_10_7',cm='',res=256,rad=''):
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
	d2 = False
	if samp == '21p22':
		f = fitsio.read(dirsci+'eboss21'+'.'+v+'.latest.rands.fits')	
		d2 = True
		f2 = fitsio.read(dirsci+'eboss22'+'.'+v+'.latest.rands.fits')	
	else:
		f = fitsio.read(dirsci+'eboss'+samp+'.'+v+'.latest.rands.fits')	
	for i in range(0,len(f)):
		ra,dec = f[i]['ra'],f[i]['dec']
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(res,th,phi))
		pixl[p] += 1.
	if d2:
		for i in range(0,len(f2)):
			ra,dec = f2[i]['ra'],f2[i]['dec']
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			pixl[p] += 1.

	fo = open('ranHeal_pix'+str(res)+'eboss'+cm+samp+v+'_elg.dat','w')
	for i in range(0,len(pixl)):
		if pixl[i] > 0:
			fo.write(str(i)+' '+str(pixl[i])+'\n')
	fo.close()


def mkgalELG4xi(zmin=.6,zmax=1.1,samp='21',v='v5_10_7',c='sci',app='.fits',compl=.8,compls=.7,fkp='fkp',wm='wstar'):
	from healpix import healpix, radec2thphi
	if c == 'sci': #AJR uses this define directory for machine he uses
		dir = dirsci
	if wm == 'wstar' or wm == 'cpstar':
		wsys = np.loadtxt('allstars17.519.9Healpixall256.dat')
		b,ms = findlinmb('ELG',samp,v,'star',zmin,zmax,dir='')
		h = healpix()
		print b,ms
	#print wm.split('_')[1]	
	iv = False
	if len(wm.split('_')) > 1:
		if wm.split('_')[1] == 'ivar':
			b,ms = findlinmb(samp,v,'',wm,zmin,zmax,dir='')
			iv = True
			print b,ms
			
	ffkp = np.loadtxt('nbarELG'+samp+v+'.dat').transpose()
	d2 = False
	if samp == '21p22':
		f = fitsio.read(dir+'eboss21'+'.'+v+'.latest'+app)
		d2 = True
		f2 = fitsio.read(dir+'eboss22'+'.'+v+'.latest'+app)
	else:
		f = fitsio.read(dir+'eboss'+samp+'.'+v+'.latest'+app) #read galaxy/quasar file
	fo = open(dir+'gebosselg_'+samp+v+'_mz'+str(zmin)+'xz'+str(zmax)+fkp+wm+'4xi.dat','w')
	n = 0
	mins = 100
	maxs = 0
	for i in range(0,len(f)):
		z = f[i]['Z']
		
		w =1.
		m = 0
		#if f[i]['dec'] < .5 and f[i]['ra'] > 350 and f[i]['ra'] < 355:
		#	m = 1
		if f[i]['sector_TSR'] < compl or f[i]['sector_SSR'] < compls:
			m = 1
		if z > zmin and z < zmax and m == 0 and f[i]['Z_reliable'] == True and f[i]['isdupl'] == False:# and f[i]['depth_ivar_r'] > 100 and f[i]['depth_ivar_g'] > 300 and f[i]['depth_ivar_z'] > 50:
			ra,dec = f[i]['ra'],f[i]['dec']
			th,phi = radec2thphi(ra,dec)
			if wm == 'wstar' or wm == 'cpstar' or wm == 'wext':
				pix2 = int(h.ang2pix_nest(256,th,phi))
				#print pix2
				ns = wsys[pix2]
				ws = 1./(b+ms*ns)
				if ws > maxs:
					maxs = ws
				if ws < mins:
					mins = ws
				w = w*ws
			else:	
				if iv:
					sysv = f[i][wm]
					if sysv > 0:
						sysv = -2.5*(log(5./sqrt(sysv),10.)-9)
						ws = 1./(b+sysv*ms)
						w = w*ws	
			zind = int(z/.01)
			if fkp == 'fkp':
				fkpw = ffkp[-1][zind]
				w = w*fkpw
			fo.write(str(f[i]['ra'])+' '+str(f[i]['dec'])+' '+str(z)+' '+str(w)+'\n')
			n += 1.
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
		
			w =1.
			m = 0
			#if f[i]['dec'] < .5 and f[i]['ra'] > 350 and f[i]['ra'] < 355:
			#	m = 1
			if f2[i]['sector_TSR'] < compl or f2[i]['sector_SSR'] < compls:
				m = 1
			if z > zmin and z < zmax and m == 0 and f2[i]['Z_reliable'] == True and f2[i]['isdupl'] == False:# and f[i]['depth_ivar_r'] > 100 and f[i]['depth_ivar_g'] > 300 and f[i]['depth_ivar_z'] > 50:
				ra,dec = f2[i]['ra'],f2[i]['dec']
				th,phi = radec2thphi(ra,dec)
				if wm == 'wstar' or wm == 'cpstar' or wm == 'wext':
					pix2 = int(h.ang2pix_nest(256,th,phi))
					#print pix2
					ns = wsys[pix2]
					ws = 1./(b+ms*ns)
					if ws > maxs:
						maxs = ws
					if ws < mins:
						mins = ws
					w = w*ws
				else:
					if iv:
						sysv = f2[i][wm]
						if sysv > 0:
							sysv = -2.5*(log(5./sqrt(sysv),10.)-9)
							ws = 1./(b+sysv*ms)
							w = w*ws	

				zind = int(z/.01)
				if fkp == 'fkp':
					fkpw = ffkp[-1][zind]
					w = w*fkpw
				fo.write(str(f2[i]['ra'])+' '+str(f2[i]['dec'])+' '+str(z)+' '+str(w)+'\n')
				n += 1.

	print n	
	print mins,maxs,wm	
	fo.close()
	return True		

def mkranELG4xi(samp='21',v='v5_10_7',zmin=.7,zmax=1.1,comp = 'sci',N=0,app='.fits',compl=.8,compls=.7,fkp='fkp',wm='wstar'):
	from random import random
	if comp == 'sci':
		dir = dirsci 
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	gf = np.loadtxt(dir+'gebosselg_'+samp+v+'_mz'+str(zmin)+'xz'+str(zmax)+fkp+wm+'4xi.dat').transpose()
	d2 = False
	if samp == '21p22':
		f = fitsio.read(dir+'eboss21'+'.'+v+'.latest.rands'+app)
		d2 = True
		f2 = fitsio.read(dir+'eboss22'+'.'+v+'.latest.rands'+app)
		ns = (len(f)+len(f2))/1000000
		print len(f)+len(f2),ns
		
	else:
		f = fitsio.read(dir+'eboss'+samp+'.'+v+'.latest.rands'+app)
		ns = len(f)/1000000
		print len(f),ns

	fo = open(dir+'rebosselg'+'_'+samp+v+'_'+str(N)+wz+fkp+wm+'4xi.dat','w')
	n = 0
	nw = 0
	#minc = N*10**6
	#maxc = (N+1)*10**6 #will become relevant once files are big enough
	
	#if len(f) < maxc:
	#	maxc = len(f)
	#for i in range(minc,maxc):
	for i in range(N,len(f),ns):
		indz = int(random()*len(gf[2]))
		z = gf[2][indz]
		wg = gf[3][indz]
		if wm == 'noSSR':
			w = wg*f[i]['sector_TSR']
		else:	
			w = wg*f[i]['sector_TSR']*f[i]['plate_SSR']
		m = 0
		#if f[i]['dec'] < .5 and f[i]['ra'] > 350 and f[i]['ra'] < 355:
		#	m = 1

		#if z > zmin and z < zmax:
		if f[i]['sector_TSR'] >= compl and f[i]['sector_SSR'] >= compls:# and f[i]['depth_ivar_r'] > 100 and f[i]['depth_ivar_g'] > 300 and f[i]['depth_ivar_z'] > 50:
			fo.write(str(f[i]['ra'])+' '+str(f[i]['dec'])+' '+str(z)+' '+str(w)+'\n')
			n += 1.
			nw += w
	if d2:
		for i in range(N,len(f2),ns):
			indz = int(random()*len(gf[2]))
			z = gf[2][indz]
			wg = gf[3][indz]
			if wm == 'noSSR':
				w = wg*f2[i]['sector_TSR']
			else:	
				w = wg*f2[i]['sector_TSR']*f2[i]['plate_SSR']
			m = 0
			#if f[i]['dec'] < .5 and f[i]['ra'] > 350 and f[i]['ra'] < 355:
			#	m = 1

			#if z > zmin and z < zmax:
			if f2[i]['sector_TSR'] >= compl and f2[i]['sector_SSR'] >= compls:# and f[i]['depth_ivar_r'] > 100 and f[i]['depth_ivar_g'] > 300 and f[i]['depth_ivar_z'] > 50:
				fo.write(str(f2[i]['ra'])+' '+str(f2[i]['dec'])+' '+str(z)+' '+str(w)+'\n')
				n += 1.
				nw += w

	print n,nw #just helps to know things worked properly
	print n/10000.*ns #area in sq degrees
	fo.close()
	return True


def mkodens(reg,zmin,zmax,v='4',res=256,tol=.2,app='.fits'):
	#dir = 'output_v4/'
	from healpix import radec2thphi
	import healpy as hp
	np = 12*res*res
	gl = []
	rl = []
	for i in range(0,np):
		gl.append(0)
		rl.append(0)
	f = fitsio.read(dirsci+'eBOSS_ELG'+'_clustering_'+reg+'_v'+v+'.dat.fits')
	n = 0
	nw = 0
	nnan = 0
	mins = 100
	maxs = 0
	for i in range(0,len(f)):
		z = f[i]['Z']
		
		w =1.
		m = 0
		#if f[i]['dec'] < .5 and f[i]['ra'] > 350 and f[i]['ra'] < 355:
		#	m = 1
		
		if z > zmin and z < zmax: 
			ra,dec = f[i]['RA'],f[i]['DEC']
			#if rec == '_rec':
			#	w = f[i]['WEIGHT_SYSTOT']
			#else:	
			w = f[i]['WEIGHT_SYSTOT']*f[i]['WEIGHT_CP']*f[i]['WEIGHT_NOZ']

	for i in range(0,len(f)):
		z = f[i]['Z']
		if z > zmin and z < zmax:		
			ra,dec = f[i]['RA'],f[i]['DEC']
			w = f[i]['WEIGHT_SYSTOT']*f[i]['WEIGHT_CP']*f[i]['WEIGHT_NOZ']
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi)
			gl[p] += w
	f = fitsio.read(dirsci+'eBOSS_ELG'+'_clustering_'+reg+'_v'+v+'.ran.fits')
	for i in range(0,len(f)):
		z = f[i]['Z']	
		if z > zmin and z < zmax:		
			ra,dec = f[i]['RA'],f[i]['DEC']
			w = f[i]['WEIGHT_SYSTOT']*f[i]['WEIGHT_CP']*f[i]['WEIGHT_NOZ']
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi)
			rl[p] += w

	avet = sum(gl)/sum(rl)	

	ng = 0
	np = 0 
	for i in range(0,len(gl)):
		if rl[i] > tol:
			ng += gl[i]
			np += rl[i]
	print ng,np,sum(gl),sum(rl)
	ave = ng/np
	print ave,avet
	fo = open('galelg'+reg+v+str(zmin)+str(zmax)+str(res)+'odenspczw.dat','w')
	ft = open('galelg'+reg+v+str(zmin)+str(zmax)+str(res)+'rdodens.dat','w')
	no = 0		
	for i in range(0,len(rl)):
		if rl[i] > tol:
			th,phi = hp.pix2ang(res,i)
			sra = sin(phi)
			cra = cos(phi)
			sdec = sin(-1.*(th-pi/2.))
			cdec = cos(-1.*(th-pi/2.))
			od = gl[i]/(ave*rl[i]) -1.
			#od = gl[i]
			fo.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(od)+' '+str(rl[i])+'\n')
			ft.write(str(th)+' '+str(phi)+' '+str(od)+' '+str(rl[i])+'\n')
	print no
	#ft.close()
	fo.close()
	return True

def mkodens_BOSS(samp,reg,zmin,zmax,res=256,tol=.2,app='.fits'):
	#dir = 'output_v4/'
	dirb = '/mnt/lustre/ashleyr/BOSS/'
	from healpix import radec2thphi
	import healpy as hp
	np = 12*res*res
	gl = []
	rl = []
	for i in range(0,np):
		gl.append(0)
		rl.append(0)
	f = fitsio.read(dirb+'galaxy_DR12v5_'+samp+'_'+reg+'.fits.gz')
	n = 0
	nw = 0
	nnan = 0
	mins = 100
	maxs = 0
	for i in range(0,len(f)):
		z = f[i]['Z']
		
		w =1.
		m = 0
		#if f[i]['dec'] < .5 and f[i]['ra'] > 350 and f[i]['ra'] < 355:
		#	m = 1
		
		if z > zmin and z < zmax: 
			ra,dec = f[i]['RA'],f[i]['DEC']
			#if rec == '_rec':
			#	w = f[i]['WEIGHT_SYSTOT']
			#else:	
			w = (f[i]['WEIGHT_SYSTOT']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_NOZ']
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi)
			gl[p] += w

	f = fitsio.read(dirb+'random0_DR12v5_'+samp+'_'+reg+'.fits.gz')
	for i in range(0,len(f)):
		z = f[i]['Z']	
		if z > zmin and z < zmax:		
			ra,dec = f[i]['RA'],f[i]['DEC']
			w = 1.
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi)
			rl[p] += w

	avet = sum(gl)/sum(rl)	

	ng = 0
	np = 0 
	for i in range(0,len(gl)):
		if rl[i] > tol:
			ng += gl[i]
			np += rl[i]
	print ng,np,sum(gl),sum(rl)
	ave = ng/np
	print ave,avet
	fo = open('gal'+samp+reg+str(zmin)+str(zmax)+str(res)+'odenspczw.dat','w')
	ft = open('gal'+samp+reg+str(zmin)+str(zmax)+str(res)+'rdodens.dat','w')
	no = 0		
	for i in range(0,len(rl)):
		if rl[i] > tol:
			th,phi = hp.pix2ang(res,i)
			sra = sin(phi)
			cra = cos(phi)
			sdec = sin(-1.*(th-pi/2.))
			cdec = cos(-1.*(th-pi/2.))
			od = gl[i]/(ave*rl[i]) -1.
			#od = gl[i]
			fo.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(od)+' '+str(rl[i])+'\n')
			ft.write(str(th)+' '+str(phi)+' '+str(od)+' '+str(rl[i])+'\n')
	print no
	#ft.close()
	fo.close()
	return True


def plotcrosscor(reg):
	import numpy as np
	from matplotlib import pyplot as plt
	d1 = np.loadtxt('galelg'+reg+'40.60.7256galelg'+reg+'41.01.12562ptPixc.dat').transpose()
	d2 = np.loadtxt('galelg'+reg+'40.60.7256galelg'+reg+'40.91.02562ptPixc.dat').transpose()
	d3 = np.loadtxt('galelg'+reg+'40.60.7256galelg'+reg+'40.80.92562ptPixc.dat').transpose()
	d4 = np.loadtxt('galelg'+reg+'40.80.9256galelg'+reg+'41.01.12562ptPixc.dat').transpose()
	d5 = np.loadtxt('galelg'+reg+'40.70.8256galelg'+reg+'41.01.12562ptPixc.dat').transpose()
	plt.plot(d1[0],d1[1],d2[0],d2[1],d3[0],d3[1],d4[0],d4[1],d5[0],d5[1])
	plt.legend(('0.60.7x1.01.1','0.60.7x0.91.0','0.60.7x0.80.9','0.80.9x1.01.1','0.70.8x1.01.1'))
	plt.ylim(-0.02,0.02)
	plt.show()
	return True


def plotodens(chunk,res=256,zmin=.6,zmax=1.1,v='v5_10_7',nranmin=40,cmap=cm.coolwarm,vmax=10):
	from healpix import thphi2radec
	try:
		d = np.loadtxt(ebossdir+'galelg'+chunk+v+str(zmin)+str(zmax)+str(res)+'rdodens.dat').transpose()
	except:	
		mkodens(chunk,zmin,zmax,v,res=res)
		d = np.loadtxt(ebossdir+'galelg'+chunk+v+str(zmin)+str(zmax)+str(res)+'rdodens.dat').transpose()
	rl = []
	dl = []
	cl = []
	for i in range(0,len(d[0])):
		if d[-1][i] >= nranmin:
			th,phi = d[0][i],d[1][i]
			ra,dec = thphi2radec(th,phi)
			rl.append(ra)
			dl.append(dec)
			cl.append(d[2][i])
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	if vmax > max(cl):
		vmax = max(cl)
	map = plt.scatter(rl,dl,c=cl,s=10,cmap=cmap,lw=0,vmax=vmax)#,vmin=-1,vmax=5,lw=0)
	#cbar = plt.colorbar(map)
	#plt.xlabel('right ascension J2000 (degrees)',fontsize=24)
	#plt.ylabel('declination J2000 (degrees)',fontsize=24)
	plt.colorbar(map,orientation='horizontal')#,ticks=[.5, .6, .7, .8, .9, 1.])
	#plt.drawmeridians(arange(0,360,20),linewidth=.2,fontsize=20)
	#plt.savefig(ebossdir+'elg'+chunk+comp+'.png')
	plt.show()
	#pp.savefig()
	#pp.close()
	return True

def plotELGcomp_full(ver='test',reg='SGC',cmap=cm.coolwarm,comp='COMP_BOSS'):
	#remember!!! export PATH=$PATH:/Users/ashleyross/mangle2.2/bin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	import fitsio
	import matplotlib.cm as cm
	
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	f = fitsio.read(dirfits +'ELG'+reg+ver+'full.ran.fits')
	for i in range(0,f.size):
		if f[i]['RA'] > 180:
			f[i]['RA'] -= 360
	#w = (f['RA'] > 180)
	#f[w]['RA'] -= 360
	map = plt.scatter(f['RA'],f['DEC'],c=f[comp],s=.1,cmap=cmap,vmin=.5,lw=0)
	#cbar = plt.colorbar(map)
	plt.xlabel('right ascension J2000 (degrees)',fontsize=24)
	plt.ylabel('declination J2000 (degrees)',fontsize=24)
	plt.colorbar(map,orientation='horizontal',ticks=[.5, .6, .7, .8, .9, 1.])
	#plt.drawmeridians(arange(0,360,20),linewidth=.2,fontsize=20)
	plt.savefig(ebossdir+'elg'+reg+ver+comp+'.png')
	#plt.show()
	#pp.savefig()
	#pp.close()
	return True

def plotELGcomp_AR(chunk,comp='sector_TSR',cmap=cm.coolwarm,v='v5_10_7'):
	#remember!!! export PATH=$PATH:/Users/ashleyross/mangle2.2/bin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	import fitsio
	import matplotlib.cm as cm
	
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	f = fitsio.read(dirfits +'ELG'+'.'+v+'.latest.rands.fits')
	w = (f['ra'] <1e6) # just something stupid
	if chunk != 'ALL':
		w = (f['chunk'] == 'eboss'+str(chunk))
	for i in range(0,f.size):
		if f[i]['ra'] > 180:
			f[i]['ra'] -= 360
	
	map = plt.scatter(f[w]['ra'],f[w]['dec'],c=f[w][comp],s=.1,cmap=cmap,vmin=.5,lw=0)
	#cbar = plt.colorbar(map)
	plt.xlabel('right ascension J2000 (degrees)',fontsize=24)
	plt.ylabel('declination J2000 (degrees)',fontsize=24)
	plt.colorbar(map,orientation='horizontal',ticks=[.5, .6, .7, .8, .9, 1.])
	#plt.drawmeridians(arange(0,360,20),linewidth=.2,fontsize=20)
	plt.savefig(ebossdir+'elg'+chunk+comp+'.png')
	#plt.show()
	#pp.savefig()
	#pp.close()
	return True


def plotELGcomp_simp(ver='test',reg='SGC',sys='cpgdepth',cmap=cm.coolwarm,comp='COMP'):
	#remember!!! export PATH=$PATH:/Users/ashleyross/mangle2.2/bin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	import fitsio
	import matplotlib.cm as cm
	
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	f = fitsio.read(dirfits +'ELG'+reg+ver+sys+'.ran.fits')
	map = plt.scatter(f['RA'],f['DEC'],c=f[comp],s=.1,cmap=cmap,vmin=.5,lw=0)
	#cbar = plt.colorbar(map)
	plt.xlabel('right ascension J2000 (degrees)',fontsize=24)
	plt.ylabel('declination J2000 (degrees)',fontsize=24)
	plt.colorbar(map,orientation='horizontal',ticks=[.5, .6, .7, .8, .9, 1.])
	#plt.drawmeridians(arange(0,360,20),linewidth=.2,fontsize=20)
	plt.savefig(ebossdir+'elg'+reg+ver+comp+'.png')
	#plt.show()
	#pp.savefig()
	#pp.close()
	return True


def plotELGgalchunk(chunk,v='v5_10_7',zmin=.6,zmax=1.1,compl=.5,compls=0):
	#remember!!! export PATH=$PATH:/Users/ashleyross/mangle2.2/bin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	import fitsio
	import matplotlib.cm as cm
	
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	#f = fitsio.read(dirfits +'eboss'+chunk+'.'+v+'.latest.fits')
	f = fitsio.read(dirfits +'ELG.'+v+'.latest.fits')
	ral = []
	decl = []
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			if z > zmin and z < zmax and f[i]['chunk'] == chunk and f[i]['Z_reliable']==True and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['isdupl'] == False:
				ral.append(f[i]['ra'])
				decl.append(f[i]['dec'])

	#f = fitsio.read(dirfits +'eboss'+chunk+'.'+v+'.latest.rands.fits')
	f = fitsio.read(dirfits +'ELG.'+v+'.latest.rands.fits')
	rarl = []
	decrl = []
	for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['chunk'] == chunk:
				rarl.append(f[i]['ra'])
				decrl.append(f[i]['dec'])

	plt.plot(rarl,decrl,'k,')
	plt.plot(ral,decl,'ro',markersize=.4)			

	plt.xlabel('right ascension J2000 (degrees)',fontsize=24)
	plt.ylabel('declination J2000 (degrees)',fontsize=24)
	plt.show()
	#pp.savefig()
	#pp.close()
	return True

def plotELGgalchunk_zsplit(chunk,v='v5_10_7',zmin=.6,zmax=1.1,compl=.8,compls=.7,zspl=.8):
	#remember!!! export PATH=$PATH:/Users/ashleyross/mangle2.2/bin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	import fitsio
	import matplotlib.cm as cm
	
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	f = fitsio.read(dirfits +'eboss'+chunk+'.'+v+'.latest.fits')
	ral = []
	decl = []
	ral2 = []
	decl2 = []

	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			if z > zmin and z < zmax and f[i]['Z_reliable']==True and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['isdupl'] == False:
				if z < zspl:
					ral.append(f[i]['ra'])
					decl.append(f[i]['dec'])
				else:
					ral2.append(f[i]['ra'])
					decl2.append(f[i]['dec'])
						
	plt.plot(ral2,decl2,'ro',markersize=.4)
	plt.plot(ral,decl,'ko',markersize=.4)			
				

	plt.xlabel('right ascension J2000 (degrees)',fontsize=24)
	plt.ylabel('declination J2000 (degrees)',fontsize=24)
	plt.show()
	#pp.savefig()
	#pp.close()
	return True



def plotNbarELG(ver='v5_10_7',nocut='nocut'):
	from matplotlib import pyplot as plt
	d1 = np.loadtxt(ebossdir+'nbarELG21'+ver+nocut+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'nbarELG22'+ver+nocut+'.dat').transpose()
	d3 = np.loadtxt(ebossdir+'nbarELG23'+ver+nocut+'.dat').transpose()
	plt.plot(d1[0],d1[1]*10000,'k-')
	plt.plot(d1[0],d2[1]*10000,'r-')
	plt.plot(d1[0],d3[1]*10000,'b-')
	plt.xlim(0.5,1.2)
	plt.xlabel('redshift')
	plt.ylabel(r'number density (10$^4$ [$h$/Mpc]$^3$)')
	plt.text(1.1,4.6,'chunk 21',color='k')
	plt.text(1.1,4.4,'chunk 22',color='r')
	plt.text(1.1,4.2,'chunk 23',color='b')
	plt.show()
	return True


def mkNbarELG(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,snrm=-1,snrx=1000):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0
	if chunk == '21p22':
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['plate_rSN2'] > snrm and f[i]['plate_rSN2'] < snrx:
				sumr += f[i]['sector_TSR']
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f2)):
			if f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls and f2[i]['plate_rSN2'] > snrm and f2[i]['plate_rSN2'] < snrx:
				sumr += f2[i]['sector_TSR']

		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:
		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['plate_rSN2'] > snrm and f[i]['plate_rSN2'] < snrx:
				sumr += f[i]['sector_TSR']

		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	area = (sumr)/10000.
	print 'effective area is '+str(area)
	no = 0
	fo = open(ebossdir+'nbarELG'+chunk+ver+'.dat','w')
	nb = int(zmax/sp)
	zl = []
	for i in range(0,nb):
		zl.append(0)
	sum = 0
	sumw = 0
	sumt = 0
	ndup = 0	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
			
			wfczss = 1./(f[i]['plate_SSR'])
			
			sum += wfczss
			if z > zmin and z < zmax and f[i]['Z_reliable']==True and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['isdupl'] == False and f[i]['plate_rSN2'] > snrm and f[i]['plate_rSN2'] < snrx:
				sumt += 1.
				sumw += wfczss
				zl[zind] += wfczss
			if f[i]['isdupl']:
				ndup += 1.	
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				wfczss = 1./(f2[i]['plate_SSR'])
			
				sum += wfczss
				if z > zmin and z < zmax and f2[i]['Z_reliable']==True and f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls and f2[i]['isdupl'] == False and f2[i]['plate_rSN2'] > snrm and f2[i]['plate_rSN2'] < snrx:
					sumt += 1.
					sumw += wfczss
					zl[zind] += wfczss
				if f2[i]['isdupl']:
					ndup += 1.	
			
	print ndup
	print sumw/sumt
	#areaf = sumt/sumw
	#area = area*areaf
	#print 'effective area is '+str(area)
	vl = []
	veffl = []
	nl = []
	for i in range(0,len(zl)):
		zlo = i*sp
		zh = (i+1)*sp
		v = area/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		vl.append(v)
		nbarz =  zl[i]/v	
		nl.append(nbarz)
		veff = v*(nbarz*P0/(1.+nbarz*P0))**2.
		veffl.append(veff)

	
	f = 1.
	zpl = []
	nbl = []
	nbwl = []
	wtot = 0
	zm = 0
	nw = 0
	nnw = 0
	for i in range(0,nb):
		z = sp/2.+sp*i
		zpl.append(z)
		fo.write(str(z)+' '+str(nl[i])+' '+str(zl[i])+' '+str(vl[i])+' '+str(veffl[i])+' '+str(1./(1.+zl[i]/vl[i]/f*P0))+'\n')
		if z > .6 and z < 1.1:
			zm += z*1./(1.+zl[i]/vl[i]/f*P0)
			wtot += 1./(1.+zl[i]/vl[i]/f*P0)
			nw += zl[i]/(1.+zl[i]/vl[i]/f*P0)
			nnw += zl[i]
		nbl.append(zl[i]/vl[i]/f)
		nbwl.append(zl[i]/vl[i]/f*1./(1.+zl[i]/vl[i]/f*P0))
	fo.close()
	print sum,zm/wtot,nw,nnw
	plt.plot(zpl,nbl)
	plt.show()
	plt.plot(zpl,nbwl)
	plt.show()
	return True

def plotnzplate(chunk,plate,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,snrm=-1,snrx=1000):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0

	f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	nb = int(zmax/sp)
	zl = []
	for i in range(0,nb):
		zl.append(0)
	sum = 0
	sumw = 0
	sumt = 0
	ndup = 0
	nlow = 0	
	for i in range(0,len(f)):
		if f[i]['PLATE'] == plate:
			z = f[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				
				sum += 1.
				if z > zmin and z < zmax and f[i]['Z_reliable']==True and f[i]['isdupl'] == False:
					sumt += 1.
					zl[zind] += 1.
					if z > 0.6 and z < 0.75:
						nlow += 1.

				if f[i]['isdupl']:
					ndup += 1.	
					
	print ndup,sum,sumt,nlow,nlow/sumt

def platevsfraclow(chunk='23',platemin=9546,platemax=9631,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,snrm=-1,snrx=1000):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0

	f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	nb = int(zmax/sp)
	zl = []
	for i in range(0,nb):
		zl.append(0)
	sum = 0
	sumw = 0
	sumt = 0
	ndup = 0
	nlow = 0
	nll = np.zeros((platemax-platemin+1))
	nlt = np.zeros((platemax-platemin+1))	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
		
			
			sum += 1.
			if z > zmin and z < zmax and f[i]['Z_reliable']==True and f[i]['isdupl'] == False:
				sumt += 1.
				plind = f[i]['PLATE']-platemin
				nlt[plind] += 1.
				if z > 0.6 and z < 0.75:
					nll[plind] += 1.
					nlow += 1.

			if f[i]['isdupl']:
				ndup += 1.	
					
	print ndup,sum,sumt,nlow,nlow/sumt

	pl = []
	for i in range(platemin,platemax+1):
		pl.append(i)
	plt.plot(pl,nll/nlt)
	plt.xlabel('Plate')
	plt.ylabel(r'$N(0.6 < z < 0.75)/N(0.6 < z < 1.1)$')
	plt.title('Fraction of low redshift ELGs vs. plate, Chunk 23')
	plt.show()
	return True

def platevsdepth(chunk='23',platemin=9546,platemax=9631,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,snrm=-1,snrx=1000):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0

	f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	nb = int(zmax/sp)
	zl = []
	for i in range(0,nb):
		zl.append(0)
	sum = 0
	sumw = 0
	sumt = 0
	ndup = 0
	nlow = 0
	nll = np.zeros((platemax-platemin+1))
	nlt = np.zeros((platemax-platemin+1))	
	for i in range(0,len(f)):
		sumt += 1.
		plind = f[i]['PLATE']-platemin
		nlt[plind] += 1.
		nll[plind] += f[i]['depth_ivar_g']
					
#	print ndup,sum,sumt,nlow,nlow/sumt

	pl = []
	for i in range(platemin,platemax+1):
		pl.append(i)
	plt.plot(pl,nll/nlt)
	plt.xlabel('Plate')
	plt.ylabel('mean ivar g depth')
	#plt.title('Fraction of low redshift ELGs vs. plate, Chunk 23')
	plt.show()
	return True

def platevsdepthr(chunk='23',platemin=9546,platemax=9631,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,snrm=-1,snrx=1000):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0

	f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	nb = int(zmax/sp)
	zl = []
	for i in range(0,nb):
		zl.append(0)
	sum = 0
	sumw = 0
	sumt = 0
	ndup = 0
	nlow = 0
	nll = np.zeros((platemax-platemin+1))
	nlt = np.zeros((platemax-platemin+1))	
	for i in range(0,len(f)):
		sumt += 1.
		plind = f[i]['PLATE']-platemin
		nlt[plind] += 1.
		nll[plind] += f[i]['depth_ivar_r']
					
#	print ndup,sum,sumt,nlow,nlow/sumt

	pl = []
	for i in range(platemin,platemax+1):
		pl.append(i)
	plt.plot(pl,nll/nlt)
	plt.xlabel('Plate')
	plt.ylabel('mean ivar r depth')
	#plt.title('Fraction of low redshift ELGs vs. plate, Chunk 23')
	plt.show()
	return True

def platevsdepthz(chunk='23',platemin=9546,platemax=9631,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,snrm=-1,snrx=1000):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0

	f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	nb = int(zmax/sp)
	zl = []
	for i in range(0,nb):
		zl.append(0)
	sum = 0
	sumw = 0
	sumt = 0
	ndup = 0
	nlow = 0
	nll = np.zeros((platemax-platemin+1))
	nlt = np.zeros((platemax-platemin+1))	
	for i in range(0,len(f)):
		sumt += 1.
		plind = f[i]['PLATE']-platemin
		nlt[plind] += 1.
		nll[plind] += f[i]['depth_ivar_z']
					
#	print ndup,sum,sumt,nlow,nlow/sumt

	pl = []
	for i in range(platemin,platemax+1):
		pl.append(i)
	plt.plot(pl,nll/nlt)
	plt.xlabel('Plate')
	plt.ylabel('mean ivar z depth')
	#plt.title('Fraction of low redshift ELGs vs. plate, Chunk 23')
	plt.show()
	return True


def mkNbarELG_nocut(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0
	if chunk == '21p22':
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				sumr += f[i]['sector_TSR']
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f2)):
			if f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls:
				sumr += f2[i]['sector_TSR']

		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:
		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				sumr += f[i]['sector_TSR']

		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	area = (sumr)/10000.
	print 'effective area is '+str(area)
	no = 0
	fo = open(ebossdir+'nbarELG'+chunk+ver+'nocut.dat','w')
	nb = int(zmax/sp)
	zl = []
	for i in range(0,nb):
		zl.append(0)
	sum = 0
	sumw = 0
	sumt = 0
	ndup = 0	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
			
			wfczss = 1.#/(f[i]['plate_SSR'])
			
			sum += wfczss
			if z > zmin and z < zmax and f[i]['CLASS'] == 'GALAXY' and f[i]['isdupl'] == False:#f[i]['Z_reliable']==True and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['isdupl'] == False:
				sumt += 1.
				sumw += wfczss
				zl[zind] += wfczss
			if f[i]['isdupl']:
				ndup += 1.	
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				wfczss = 1.#/(f2[i]['plate_SSR'])
			
				sum += wfczss
				if z > zmin and z < zmax and f2[i]['CLASS'] == 'GALAXY' and f2[i]['isdupl'] == False:#f2[i]['Z_reliable']==True and f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls and f2[i]['isdupl'] == False:
					sumt += 1.
					sumw += wfczss
					zl[zind] += wfczss
				if f2[i]['isdupl']:
					ndup += 1.	
			
	print ndup
	print sumw/sumt
	#areaf = sumt/sumw
	#area = area*areaf
	#print 'effective area is '+str(area)
	vl = []
	veffl = []
	nl = []
	for i in range(0,len(zl)):
		zlo = i*sp
		zh = (i+1)*sp
		v = area/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		vl.append(v)
		nbarz =  zl[i]/v	
		nl.append(nbarz)
		veff = v*(nbarz*P0/(1.+nbarz*P0))**2.
		veffl.append(veff)

	
	f = 1.
	zpl = []
	nbl = []
	nbwl = []
	wtot = 0
	zm = 0
	nw = 0
	nnw = 0
	for i in range(0,nb):
		z = sp/2.+sp*i
		zpl.append(z)
		fo.write(str(z)+' '+str(nl[i])+' '+str(zl[i])+' '+str(vl[i])+' '+str(veffl[i])+' '+str(1./(1.+zl[i]/vl[i]/f*P0))+'\n')
		if z > .6 and z < 1.1:
			zm += z*1./(1.+zl[i]/vl[i]/f*P0)
			wtot += 1./(1.+zl[i]/vl[i]/f*P0)
			nw += zl[i]/(1.+zl[i]/vl[i]/f*P0)
			nnw += zl[i]
		nbl.append(zl[i]/vl[i]/f)
		nbwl.append(zl[i]/vl[i]/f*1./(1.+zl[i]/vl[i]/f*P0))
	fo.close()
	print sum,zm/wtot,nw,nnw
	plt.plot(zpl,nbl)
	plt.show()
	plt.plot(zpl,nbwl)
	plt.show()
	return True

def ELG_SSR(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0
	if chunk == '21p22':
		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:

		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	ngal = 0
	nq = 0
	nqz = 0
	nstar = 0
	ngoodgal = 0
	nSSR = 0
	ntot = 0
	ngoodtot = 0
	ngz = 0
	nqzz = 0
	nqf = 0
	nqfz = 0
	nbg = 0
	nbgz = 0
	for i in range(0,len(f)):
		if f[i]['isdupl'] == False:
			ntot += 1.
			nSSR += f[i]['sector_SSR']
			if f[i]['Z_reliable']:
				ngoodtot += 1.
				if 0.6 < f[i]['Z'] < 1.1:
					ngz += 1.
				
			if f[i]['CLASS'] == 'GALAXY':
				ngal += 1.
				if f[i]['Z_reliable']:
					ngoodgal += 1.
				else:
					nbg += 1.	
					if 0.6 < f[i]['Z'] < 1.1:
						nbgz += 1.
			if f[i]['CLASS'] == 'QSO':
				nq += 1.
				if f[i]['Z_reliable']:
					nqz += 1.
					if 0.6 < f[i]['Z'] < 1.1:
						nqzz += 1.
				else:
					nqf += 1.
					if 0.6 < f[i]['Z'] < 1.1:	
						nqfz += 1.	

			if f[i]['CLASS'] == 'STAR':
				nstar += 1.
			
	print ngal,nq,nstar,ngoodgal,nqz,ntot,ngal+nq+nstar
	print ngoodgal/ngal,nSSR/ntot,ngoodtot/ntot,ngoodgal/ntot,ngoodgal/(ngal+nq),(ngoodgal+nqz)/(ntot),nqz/nq,nq/ntot,ngz/ngoodtot,nqzz/nqz,nqfz/nqf,nbgz/nbg
	return True

def ELG_SSRvsSNR(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr = 0
	if chunk == '21p22':
		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:

		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	
	pm = min(f['PLATE'])
	px = max(f['PLATE'])
	pl = []
	pltot = []
	plsnr = []
	for i in range(pm,px+1):
		pl.append(0)
		pltot.append(0)
		plsnr.append(0)
	ngal = 0
	nq = 0
	nqz = 0
	nstar = 0
	ngoodgal = 0
	nSSR = 0
	ntot = 0
	for i in range(0,len(f)):
		if f[i]['isdupl'] == False:
			ntot += 1.
			plate = f[i]['PLATE']
			if f[i]['CLASS'] == 'GALAXY':
				ngal += 1.
				pltot[plate-pm] += 1.
				plsnr[plate-pm] = f[i]['plate_rSN2']
				if f[i]['Z_reliable']:
					pl[plate-pm] += 1.
					ngoodgal += 1.
					nSSR += f[i]['sector_SSR']
			if f[i]['CLASS'] == 'QSO':
				nq += 1.
				if f[i]['Z_reliable']:
					nqz += 1.
			if f[i]['CLASS'] == 'STAR':
				nstar += 1.
			
	print ngal,nq,nstar,ngoodgal,nqz,ntot,ngal+nq+nstar
	print ngoodgal/ngal,nSSR/ngoodgal,ngoodgal/ntot,ngoodgal/(ngal+nq),(ngoodgal+nqz)/(ngal+nq-nqz)
	pl = np.array(pl)
	pltot = np.array(pltot)
	plt.plot(plsnr,pl/pltot,'ko')
	plt.show()
	return True


def mkNbarELG_splitSSR(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,split=.8,zm=''):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr1 = 0
	sumr2 = 0
	if chunk == '21p22':
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				if f[i]['plate_SSR'] > split:
					sumr1 += f[i]['sector_TSR']
				else:
					sumr2 += f[i]['sector_TSR']	
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f2)):
			if f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls:
				if f2[i]['plate_SSR'] > split:
					sumr1 += f2[i]['sector_TSR']
				else:
					sumr2 += f2[i]['sector_TSR']	

		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:
		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				if f[i]['plate_SSR'] > split:
					sumr1 += f[i]['sector_TSR']
				else:
					sumr2 += f[i]['sector_TSR']	

		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	area1 = (sumr1)/10000.
	area2 = (sumr2)/10000.
	print 'effective areas are '+str(area1) +' and '+str(area2)
	no = 0
	#fo = open(ebossdir+'nbarELG'+chunk+ver+'.dat','w')
	nb = int(zmax/sp)
	zl1 = []
	zl2 = []
	for i in range(0,nb):
		zl1.append(0)
		zl2.append(0)
	sum = 0
	sumw = 0
	sumt = 0	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
			if zm == 'allz':
				wfczss = 1.
			else:
				wfczss = 1./(f[i]['plate_SSR'])
			
			sum += wfczss
			kz = 0
			if f[i]['Z_reliable']==True or zm == 'allz':
				kz = 1
			if z > zmin and z < zmax and kz == 1 and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				sumt += 1.
				sumw += wfczss
				if f[i]['plate_SSR'] > split:
					zl1[zind] += wfczss
				else:
					zl2[zind] += wfczss	
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				if zm == 'allz':
					wfczss = 1.
				else:
					wfczss = 1./(f2[i]['plate_SSR'])
			
				sum += wfczss
				kz = 0
				if f2[i]['Z_reliable']==True or zm == 'allz':
					kz = 1
			
				sum += wfczss
				if z > zmin and z < zmax and kz == 1 and f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls:
					sumt += 1.
					sumw += wfczss
					if f2[i]['plate_SSR'] > split:
						zl1[zind] += wfczss
					else:
						zl2[zind] += wfczss	
			
	print sumw/sumt
	#areaf = sumt/sumw
	#area = area*areaf
	#print 'effective area is '+str(area)
	vl = []
	veffl = []
	nl = []
	nl2 = []
	zpl = []
	for i in range(0,len(zl1)):
		z = sp/2.+sp*i
		zpl.append(z)

		zlo = i*sp
		zh = (i+1)*sp
		v = area1/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		v2 = area2/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		vl.append(v)
		nbarz =  zl1[i]/v	
		nl.append(nbarz)
		nbarz2 =  zl2[i]/v2	
		nl2.append(nbarz2)
		veff = v*(nbarz*P0/(1.+nbarz*P0))**2.
		veffl.append(veff)

	
	f = 1.
	
	nbl = []
	nbwl = []
	wtot = 0
	zm = 0
	nw = 0
	nnw = 0
# 	for i in range(0,nb):
# 		fo.write(str(z)+' '+str(nl[i])+' '+str(zl[i])+' '+str(vl[i])+' '+str(veffl[i])+' '+str(1./(1.+zl[i]/vl[i]/f*P0))+'\n')
# 		if z > .6 and z < 1.1:
# 			zm += z*1./(1.+zl[i]/vl[i]/f*P0)
# 			wtot += 1./(1.+zl[i]/vl[i]/f*P0)
# 			nw += zl[i]/(1.+zl[i]/vl[i]/f*P0)
# 			nnw += zl[i]
# 		nbl.append(zl[i]/vl[i]/f)
# 		nbwl.append(zl[i]/vl[i]/f*1./(1.+zl[i]/vl[i]/f*P0))
# 	fo.close()
# 	print sum,zm/wtot,nw,nnw
	plt.plot(zpl,nl,zpl,nl2)
	plt.show()
	return True

def mkNbarELG_splitdepthr(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,split=220):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr1 = 0
	sumr2 = 0
	if chunk == '21p22':
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				if f[i]['decam_depth_r'] > split:
					sumr1 += f[i]['sector_TSR']
				else:
					sumr2 += f[i]['sector_TSR']	
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f2)):
			if f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls:
				if f2[i]['decam_depth_r'] > split:
					sumr1 += f2[i]['sector_TSR']
				else:
					sumr2 += f2[i]['sector_TSR']	

		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:
		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.rands.fits')
		for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				if f[i]['decam_depth_r'] > split:
					sumr1 += f[i]['sector_TSR']
				else:
					sumr2 += f[i]['sector_TSR']	

		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	area1 = (sumr1)/10000.
	area2 = (sumr2)/10000.
	print 'effective areas are '+str(area1) +' and '+str(area2)
	no = 0
	#fo = open(ebossdir+'nbarELG'+chunk+ver+'.dat','w')
	nb = int(zmax/sp)
	zl1 = []
	zl2 = []
	for i in range(0,nb):
		zl1.append(0)
		zl2.append(0)
	sum = 0
	sumw = 0
	sumt = 0	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
			
			wfczss = 1./(f[i]['plate_SSR'])
			
			sum += wfczss
			if z > zmin and z < zmax and f[i]['Z_reliable']==True and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				sumt += 1.
				sumw += wfczss
				if f[i]['decam_depth'][2] > split:
					zl1[zind] += wfczss
				else:
					zl2[zind] += wfczss	
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				wfczss = 1./(f2[i]['plate_SSR'])
			
				sum += wfczss
				if z > zmin and z < zmax and f2[i]['Z_reliable']==True and f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls:
					sumt += 1.
					sumw += wfczss
					if f2[i]['decam_depth'][2] > split:
						zl1[zind] += wfczss
					else:
						zl2[zind] += wfczss	
			
	print sumw/sumt
	#areaf = sumt/sumw
	#area = area*areaf
	#print 'effective area is '+str(area)
	vl = []
	veffl = []
	nl = []
	nl2 = []
	zpl = []
	for i in range(0,len(zl1)):
		z = sp/2.+sp*i
		zpl.append(z)

		zlo = i*sp
		zh = (i+1)*sp
		v = area1/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		v2 = area2/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		vl.append(v)
		nbarz =  zl1[i]/v	
		nl.append(nbarz)
		nbarz2 =  zl2[i]/v2	
		nl2.append(nbarz2)
		veff = v*(nbarz*P0/(1.+nbarz*P0))**2.
		veffl.append(veff)

	
	f = 1.
	
	nbl = []
	nbwl = []
	wtot = 0
	zm = 0
	nw = 0
	nnw = 0
# 	for i in range(0,nb):
# 		fo.write(str(z)+' '+str(nl[i])+' '+str(zl[i])+' '+str(vl[i])+' '+str(veffl[i])+' '+str(1./(1.+zl[i]/vl[i]/f*P0))+'\n')
# 		if z > .6 and z < 1.1:
# 			zm += z*1./(1.+zl[i]/vl[i]/f*P0)
# 			wtot += 1./(1.+zl[i]/vl[i]/f*P0)
# 			nw += zl[i]/(1.+zl[i]/vl[i]/f*P0)
# 			nnw += zl[i]
# 		nbl.append(zl[i]/vl[i]/f)
# 		nbwl.append(zl[i]/vl[i]/f*1./(1.+zl[i]/vl[i]/f*P0))
# 	fo.close()
# 	print sum,zm/wtot,nw,nnw
	nl = np.array(nl)
	nl2 = np.array(nl2)
	plt.plot(zpl,nl*1e4,zpl,nl2*1e4)
	plt.xlabel('redshift')
	plt.ylabel(r'$10^4 n(z) (h/Mpc)^3$')
	plt.legend(labels=['decam_depth_r>'+str(split),'decam_depth_r<='+str(split)])
	plt.title('Chunk '+chunk)
	plt.show()
	return True


def mkNbarELG_splitplate(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=.8,compls=.7,split=9590):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr1 = 0
	sumr2 = 0
	if chunk == '21p22':

		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:
		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	#fo = open(ebossdir+'nbarELG'+chunk+ver+'.dat','w')
	nb = int(zmax/sp)
	zl1 = []
	zl2 = []
	for i in range(0,nb):
		zl1.append(0)
		zl2.append(0)
	#sum = 0
	sumw = 0
	sumt = 0	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
			
			wfczss = 1./(f[i]['plate_SSR'])
			
			#sum += wfczss
			if z > zmin and z < zmax and f[i]['Z_reliable']==True and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				sumt += 1.
				sumw += wfczss
				if f[i]['PLATE'] > split:
					zl1[zind] += wfczss
				else:
					zl2[zind] += wfczss	
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				wfczss = 1./(f2[i]['plate_SSR'])
			
				#sum += wfczss
				if z > zmin and z < zmax and f2[i]['Z_reliable']==True and f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls:
					sumt += 1.
					sumw += wfczss
					if f2[i]['PLATE'] > split:
						zl1[zind] += wfczss
					else:
						zl2[zind] += wfczss	
			
	print sumw/sumt
	#areaf = sumt/sumw
	#area = area*areaf
	#print 'effective area is '+str(area)
	vl = []
	veffl = []
	nl = []
	nl2 = []
	zpl = []
	for i in range(0,len(zl1)):
		z = sp/2.+sp*i
		zpl.append(z)


	
	f = 1.
	
	nbl = []
	nbwl = []
	wtot = 0
	zm = 0
	nw = 0
	nnw = 0
# 	for i in range(0,nb):
# 		fo.write(str(z)+' '+str(nl[i])+' '+str(zl[i])+' '+str(vl[i])+' '+str(veffl[i])+' '+str(1./(1.+zl[i]/vl[i]/f*P0))+'\n')
# 		if z > .6 and z < 1.1:
# 			zm += z*1./(1.+zl[i]/vl[i]/f*P0)
# 			wtot += 1./(1.+zl[i]/vl[i]/f*P0)
# 			nw += zl[i]/(1.+zl[i]/vl[i]/f*P0)
# 			nnw += zl[i]
# 		nbl.append(zl[i]/vl[i]/f)
# 		nbwl.append(zl[i]/vl[i]/f*1./(1.+zl[i]/vl[i]/f*P0))
# 	fo.close()
# 	print sum,zm/wtot,nw,nnw
	plt.plot(zpl,zl1/sum(zl1),zpl,zl2/sum(zl2))
	plt.xlabel('redshift')
	plt.ylabel('Normalized redshift histogram')
	plt.legend(labels=['plate>9590','plate<=9590'])
	plt.show()
	return True

def nzELG_splitfailures(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=0,compls=0):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr1 = 0
	sumr2 = 0
	if chunk == '21p22':

		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:
		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	#fo = open(ebossdir+'nbarELG'+chunk+ver+'.dat','w')
	nb = int(zmax/sp)
	zl1 = []
	zl2 = []
	for i in range(0,nb):
		zl1.append(0)
		zl2.append(0)
	#sum = 0
	sumw = 0
	sumt = 0	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
			
			wfczss = 1.#/(f[i]['plate_SSR'])
			
			#sum += wfczss
			if z > zmin and z < zmax and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
				sumt += 1.
				sumw += wfczss
				if f[i]['Z_reliable']==True:
					zl1[zind] += wfczss
				else:
					zl2[zind] += wfczss	
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				wfczss = 1.#/(f2[i]['plate_SSR'])
			
				#sum += wfczss
				if z > zmin and z < zmax and f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls:
					sumt += 1.
					sumw += wfczss
					if f2[i]['Z_reliable']==True:
						zl1[zind] += wfczss
					else:
						zl2[zind] += wfczss	
			
	print sumw/sumt
	#areaf = sumt/sumw
	#area = area*areaf
	#print 'effective area is '+str(area)
	vl = []
	veffl = []
	nl = []
	nl2 = []
	zpl = []
	for i in range(0,len(zl1)):
		z = sp/2.+sp*i
		zpl.append(z)


	
	f = 1.
	
	nbl = []
	nbwl = []
	wtot = 0
	zm = 0
	nw = 0
	nnw = 0
# 	for i in range(0,nb):
# 		fo.write(str(z)+' '+str(nl[i])+' '+str(zl[i])+' '+str(vl[i])+' '+str(veffl[i])+' '+str(1./(1.+zl[i]/vl[i]/f*P0))+'\n')
# 		if z > .6 and z < 1.1:
# 			zm += z*1./(1.+zl[i]/vl[i]/f*P0)
# 			wtot += 1./(1.+zl[i]/vl[i]/f*P0)
# 			nw += zl[i]/(1.+zl[i]/vl[i]/f*P0)
# 			nnw += zl[i]
# 		nbl.append(zl[i]/vl[i]/f)
# 		nbwl.append(zl[i]/vl[i]/f*1./(1.+zl[i]/vl[i]/f*P0))
# 	fo.close()
# 	print sum,zm/wtot,nw,nnw
	zl1 = np.array(zl1)
	zl2 = np.array(zl2)
	plt.plot(zpl,zl1/sum(zl1),zpl,zl2/sum(zl2))
	plt.xlabel('redshift')
	plt.ylabel('Normalized redshift histogram')
	plt.legend(labels=['reliable','not reliable'])
	plt.show()
	return True

def nzELG_splitSNR(chunk,ver='v5_10_7',sp=0.01,zmin=0.1,zmax=1.5,P0=5000.,compl=0,compls=0,splitv=20,plus=2):
	from Cosmo import distance
	d = distance(.31,.69)
	from matplotlib import pyplot as plt
	d2 = False
	sumr1 = 0
	sumr2 = 0
	if chunk == '21p22':

		
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		d2 = True

	else:
		f = fitsio.read(dirfits+'eboss'+chunk+'.'+ver+'.latest.fits')
	no = 0
	#fo = open(ebossdir+'nbarELG'+chunk+ver+'.dat','w')
	nb = int(zmax/sp)
	zl1 = []
	zl2 = []
	for i in range(0,nb):
		zl1.append(0)
		zl2.append(0)
	#sum = 0
	sumw = 0
	sumt = 0
	ssr1 = 0
	ssr2 = 0
	n1 = 0
	n2 = 0	
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			zind = int(z/sp)
			
			wfczss = 1.#/(f[i]['plate_SSR'])
			
			#sum += wfczss
			if z > zmin and z < zmax and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['Z_reliable']==True:
				sumt += 1.
				sumw += wfczss
				if f[i]['plate_rSN2'] > splitv+plus:
					zl1[zind] += wfczss
					ssr1 += f[i]['sector_SSR']
					n1 += 1.
				if f[i]['plate_rSN2'] < splitv-plus:
					zl2[zind] += wfczss	
					ssr2 += f[i]['sector_SSR']
					n2 += 1.
	if d2:
		for i in range(0,len(f2)):
			z = f2[i]['Z']
			if z < zmax:
				zind = int(z/sp)
			
				wfczss = 1.#/(f2[i]['plate_SSR'])
			
				#sum += wfczss
				if z > zmin and z < zmax and f2[i]['sector_TSR'] > compl and f2[i]['sector_SSR'] > compls and f2[i]['Z_reliable']==True:
					sumt += 1.
					sumw += wfczss
					if f2[i]['plate_rSN2'] > splitv+plus:
						zl1[zind] += wfczss
						ssr1 += f2[i]['sector_SSR']
						n1 += 1.
					if f2[i]['plate_rSN2'] < splitv-plus:
						zl2[zind] += wfczss	
						ssr2 += f2[i]['sector_SSR']
						n2 += 1.
			
	print sumw/sumt
	fac = (ssr1/n1)/(ssr2/n2)
	#areaf = sumt/sumw
	#area = area*areaf
	#print 'effective area is '+str(area)
	vl = []
	veffl = []
	nl = []
	nl2 = []
	zpl = []
	for i in range(0,len(zl1)):
		z = sp/2.+sp*i
		zpl.append(z)


	
	f = 1.
	
	nbl = []
	nbwl = []
	wtot = 0
	zm = 0
	nw = 0
	nnw = 0
# 	for i in range(0,nb):
# 		fo.write(str(z)+' '+str(nl[i])+' '+str(zl[i])+' '+str(vl[i])+' '+str(veffl[i])+' '+str(1./(1.+zl[i]/vl[i]/f*P0))+'\n')
# 		if z > .6 and z < 1.1:
# 			zm += z*1./(1.+zl[i]/vl[i]/f*P0)
# 			wtot += 1./(1.+zl[i]/vl[i]/f*P0)
# 			nw += zl[i]/(1.+zl[i]/vl[i]/f*P0)
# 			nnw += zl[i]
# 		nbl.append(zl[i]/vl[i]/f)
# 		nbwl.append(zl[i]/vl[i]/f*1./(1.+zl[i]/vl[i]/f*P0))
# 	fo.close()
# 	print sum,zm/wtot,nw,nnw
	zl1 = np.array(zl1)
	zl2 = np.array(zl2)
	diffl = zl1/sum(zl1)*fac-zl2/sum(zl2)
	print fac
	plt.plot(zpl,zl1/sum(zl1),zpl,zl2/sum(zl2),zpl,diffl/sum(diffl))
	plt.xlabel('redshift')
	plt.ylabel('Normalized redshift histogram')
	plt.legend(labels=['SNR > '+str(splitv+plus),'SNR < '+str(splitv-plus)])
	plt.show()
	return True


def mkcov_mockELG_EZ(reg,bs=8,mom=0,N=1000,rec='_recon',v='v4'):
	if bs == 5:
		#dir = ('/Users/ashleyross/eBOSS/ELG_EZmock_clustering/2PCF_ELG'+rec+'/')
		dir = (dirsci+'EZmockELG'+v+'/')
		nbin=40
	if bs == 8:
		dir = ('/Users/ashleyross/eBOSS/ELGv4_2PCF_bin8/2PCF/')
		#dir = ('/Users/ashleyross/eBOSS/ELG_EZmock_clustering/2PCF_ELG'+rec+'_bin8/')
		nbin=25
	xiave = np.zeros((nbin))
	cov = np.zeros((nbin,nbin))

	Ntot = 0
	fac = 1.
	#if reg == 'SGC':
	#	fac = 1.4
	for i in range(1,1+N):
		zer = ''
		if i < 1000:
			zer += '0'
		if i < 100:
			zer += '0'
		if i < 10:
			zer += '0'
		#try:		
		xiave += np.loadtxt(dir+'2PCF_EZmock'+rec+'_eBOSS_ELG_'+reg+'_'+v+'_z0.6z1.1_'+zer+str(i)+'.dat').transpose()[1+mom/2]*fac
		Ntot += 1.
		#except:
		#	print i
	print Ntot		
	xiave = xiave/float(Ntot)
	for i in range(1,1+N):
		zer = ''
		if i < 1000:
			zer += '0'
		if i < 100:
			zer += '0'
		if i < 10:
			zer += '0'
		#try:		
		xii = np.loadtxt(dir+'2PCF_EZmock'+rec+'_eBOSS_ELG_'+reg+'_'+v+'_z0.6z1.1_'+zer+str(i)+'.dat').transpose()[1+mom/2]*fac
		for j in range(0,nbin):
			xij = xii[j]
			for k in range(0,nbin):
				xik = xii[k]
				cov[j][k] += (xij-xiave[j])*(xik-xiave[k])
#		except:
#			print i
	cov = cov/float(Ntot)					
	fo = open('xiave'+rec+str(mom)+reg+'ELG_EZ'+v+str(bs)+'st0.dat','w')
	errl = []
	for i in range(0,nbin):
		fo.write(str(bs/2.+bs*i)+ ' '+str(xiave[i])+ ' '+str(sqrt(cov[i][i]))+'\n')
		errl.append(sqrt(cov[i][i]))
	fo.close()	
	fo = open('cov'+rec+str(mom)+reg+'ELG_EZ'+v+str(bs)+'st0.dat','w')
	
	for i in range(0,nbin):
		for j in range(0,nbin):
			fo.write(str(cov[i][j])+' ')
		fo.write('\n')
	fo.close()
# 	from matplotlib import pyplot as plt
# 	
# 	if reg != 'comb':
# 		if rec == '_recon':
# 			rec = '_rec'
# 
# 		d = np.loadtxt('/Users/ashleyross/eBOSS/xi'+str(mom)+'gebossELG_'+reg+'3'+rec+'_mz0.6xz1.1fkp'+str(bs)+'st0.dat').transpose()
# 		print(len(d[0]),len(errl),len(d[1]))
# 		plt.errorbar(d[0][:nbin],d[0][:nbin]**2.*d[1][:nbin],d[0][:nbin]**2.*errl,fmt='ko')
# 		plt.plot(d[0][:nbin],d[0][:nbin]**2.*xiave,'r-')
# 		plt.show()
# 	else:
# 		if rec == '_recon':
# 			rec = '_rec'
# 		dn = np.loadtxt('/Users/ashleyross/eBOSS/xi'+str(mom)+'gebossELG_NGC3'+rec+'_mz0.6xz1.1fkp'+str(bs)+'st0.dat').transpose()
# 		ds = np.loadtxt('/Users/ashleyross/eBOSS/xi'+str(mom)+'gebossELG_SGC3'+rec+'_mz0.6xz1.1fkp'+str(bs)+'st0.dat').transpose()
# 		d = (dn*.8+1.*ds)/1.8
# 		plt.errorbar(d[0][:nbin],d[0][:nbin]**2.*d[1][:nbin],d[0][:nbin]**2.*errl,fmt='ko')
# 		plt.plot(d[0][:nbin],d[0][:nbin]**2.*xiave,'r-')
# 		plt.show()
		
	return True


def ngvsys_ELG(sampl,ver,sys,sysmin,sysmax,zmin,zmax,band=-1,compl=.5,wm='',umag=False,gmag=False,rmag=False,imag=False,umg=False,gri=False,gri22=''):
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
	d2 = False
	stl = []
	wstl = []
	errl = []
	binng = []
	binnr = []
	nsysbin = 10 #10 bins are being used
	for i in range(0,nsysbin):
		binng.append(0)
		binnr.append(0)

	if sampl == '21p22':
		f = fitsio.read(dirfits+'eboss21.'+ver+'.latest.rands.fits')
		f2 = fitsio.read(dirfits+'eboss22.'+ver+'.latest.rands.fits')
		d2 = True
	else:
		f = fitsio.read(dirfits+'eboss'+sampl+'.'+ver+'.latest.rands.fits') #read galaxy/quasar file
	nr = 0
	sysm = float(nsysbin)/(sysmax-sysmin)
	bs = 0
	bsr = 0
	#extc = [4.239,3.303,2.285,1.698,1.263]
	for i in range (0,len(f)):
		if f[i]['sector_TSR'] > compl:
			w = f[i]['sector_TSR']*f[i]['plate_SSR']
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
			if sys.split('_')[1] == 'ivar':
				if sysv > 0:
					sysv = -2.5*(log(5./sqrt(sysv),10.)-9)
			bins = int((sysv-sysmin)*sysm)
			if bins >= 0 and bins < nsysbin:
				binnr[bins] += w
			else:
				bsr += w

			nr += w
	if d2:
		for i in range (0,len(f2)):
			if f2[i]['sector_TSR'] > compl:
				w = f2[i]['sector_TSR']*f2[i]['plate_SSR']
				if band == -1:
					sysv = f2[i][sys]
				else:
					if sys == 'IMAGE_DEPTH_EXT':
						sysv = f2[i]['IMAGE_DEPTH'][band]
						sysv = luptm(sysv,band)-extc[band]*f[i]['EB_MINUS_V']
					else:
						sysv = f2[i][sys][band]	 
						if sys == 'IMAGE_DEPTH':
							sysv = luptm(sysv,band)
				if sys.split('_')[1] == 'ivar':
					if sysv > 0:
						sysv = -2.5*(log(5./sqrt(sysv),10.)-9)

				bins = int((sysv-sysmin)*sysm)
				if bins >= 0 and bins < nsysbin:
					binnr[bins] += w
				else:
					bsr += w

				nr += w
	print min(f[sys]),max(f[sys])
	print nr
	ffkp = np.loadtxt(ebossdir+'nbarELG22v5_10_7.dat').transpose()
	if sampl == '21p22':
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		
	else:
		f = fitsio.read(dirfits+'eboss'+sampl+'.'+ver+'.latest.fits') #read galaxy/quasar file
	
	no = 0
	zm = 0
	nt = 0
	avebs = 0
	sysstr = sys.split('_')
	
	for i in range (0,len(f)):
		z = f[i]['Z']
		c = 1
		gc = True
		if z > zmin and z < zmax and f[i]['Z_reliable'] == True and c == 1 and gc and f[i]['sector_TSR'] > compl and f[i]['isdupl'] == False:
			no += 1
			w = 1.
			#if wm == '':
			zind = int(z/.01)
			wfkp =ffkp[-1][zind]
			w = w*wfkp
			
			if sysstr[0] == 'decam' and sysstr[1] == 'depth':
				sys = 'decam_depth'
				if sysstr[2] == 'g':
					band = 1
				if sysstr[2] == 'r':
					band = 2
				if sysstr[2] == 'z':
					band = 4
			if band == -1:
				sysv = f[i][sys]
			else:
				sysv = f[i][sys][band]
				if sys == 'IMAGE_DEPTH_EXT':
					sysv = f[i]['IMAGE_DEPTH'][band]
					sysv = luptm(sysv,band)-extc[band]*f[i]['EB_MINUS_V']
				else:
					sysv = f[i][sys][band]	 
					if sys == 'IMAGE_DEPTH':
						sysv = luptm(sysv,band)
			if sys.split('_')[1] == 'ivar':
				if sysv > 0:
				
					sysv = -2.5*(log(5./sqrt(sysv),10.)-9)

			bins = int((sysv-sysmin)*sysm)
					
			if bins >= 0 and bins < nsysbin:
				binng[bins] += 1.*w
			else:
				bs += w #count numbers outside of sysmin/sysmax
				avebs += w*sysv
			zm += w*z
			nt += w
	if d2:
		for i in range (0,len(f2)):
			z = f2[i]['Z']
			c = 1
			gc = True
			if z > zmin and z < zmax and f2[i]['Z_reliable'] == True and c == 1 and gc and f2[i]['sector_TSR'] > compl and f2[i]['isdupl'] == False:
				no += 1
				w = 1.
				#if wm == '':
				zind = int(z/.01)
				wfkp =ffkp[-1][zind]
				w = w*wfkp
				if band == -1:
					sysv = f2[i][sys]
				else:
					if sys == 'IMAGE_DEPTH_EXT':
						sysv = f2[i]['IMAGE_DEPTH'][band]
						sysv = luptm(sysv,band)-extc[band]*f[i]['EB_MINUS_V']
					else:
						sysv = f2[i][sys][band]	 
						if sys == 'IMAGE_DEPTH':
							sysv = luptm(sysv,band)

				if sys.split('_')[1] == 'ivar':
					if sysv > 0:
						sysv = -2.5*(log(5./sqrt(sysv),10.)-9)

				bins = int((sysv-sysmin)*sysm)
					
				if bins >= 0 and bins < nsysbin:
					binng[bins] += 1.*w
				else:
					bs += w #count numbers outside of sysmin/sysmax
					avebs += w*sysv
				zm += w*z
				nt += w
			
	if bs > 0:
		avebs = avebs/bs
	print avebs,sysmin		
	print 'total number, weighted number'
	print no,nt
	print 'mean redshift'
	print zm/nt

	print 'total number of randoms/objects '+str(nr)+'/'+str(nt)
	print 'number of randoms/objects outside tested range '+str(bsr)+'/'+str(bs)			
	ave = nt/nr
	print 'average number of objects per random is '+ str(ave)
	fs = open(ebossdir+'n'+'geboss'+sampl+'_'+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'v'+sys+'.dat','w')
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

def ngvsys_ELGhp(sampl,ver,sys,sysmin,sysmax,zmin,zmax,nside=256,band=-1,compl=.5,wm='',umag=False,gmag=False,rmag=False,imag=False,umg=False,gri=False,gri22=''):
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
	from healpix import radec2thphi,ang2pix_ring
	mapf = fitsio.read(ebossdir+'ELG_systhp.hp'+str(nside)+'.fits')
	mapl = np.zeros((12*nside*nside))
	for i in range(0,len(mapf)):
		pix = mapf[i]['hpind']
		mapl[pix] = mapf[i][sys]
	sysm = min(mapf[sys])
	sysx = max(mapf[sys])
	if sysm > sysmin:
		sysmin = sysm
	if sysx < sysmax:
		sysmax = sysx
	print sysmin,sysmax	
	d2 = False
	stl = []
	wstl = []
	errl = []
	binng = []
	binnr = []
	nsysbin = 10 #10 bins are being used
	for i in range(0,nsysbin):
		binng.append(0)
		binnr.append(0)

	if sampl == '21p22':
		f = fitsio.read(dirfits+'eboss21.'+ver+'.latest.rands.fits')
		f2 = fitsio.read(dirfits+'eboss22.'+ver+'.latest.rands.fits')
		d2 = True
	else:
		f = fitsio.read(dirfits+'eboss'+sampl+'.'+ver+'.latest.rands.fits') #read galaxy/quasar file
	nr = 0
	sysm = float(nsysbin)/(sysmax-sysmin)
	bs = 0
	bsr = 0
	#extc = [4.239,3.303,2.285,1.698,1.263]
	for i in range (0,len(f)):
		if f[i]['sector_TSR'] > compl:
			w = f[i]['sector_TSR']*f[i]['plate_SSR']
			ra,dec = f[i]['ra'],f[i]['dec']
			th,phi = radec2thphi(ra,dec)
			pix = ang2pix_ring(nside,th,phi)
			sysv = mapl[pix]			
			bins = int((sysv-sysmin)*sysm)
			if bins >= 0 and bins < nsysbin:
				binnr[bins] += w
			else:
				bsr += w

			nr += w
	if d2:
		for i in range (0,len(f2)):
			if f2[i]['sector_TSR'] > compl:
				w = f2[i]['sector_TSR']*f2[i]['plate_SSR']
				ra,dec = f2[i]['ra'],f2[i]['dec']
				th,phi = radec2thphi(ra,dec)
				pix = ang2pix_ring(nside,th,phi)
				sysv = mapl[pix]			

				bins = int((sysv-sysmin)*sysm)
				if bins >= 0 and bins < nsysbin:
					binnr[bins] += w
				else:
					bsr += w

				nr += w
	print nr
	ffkp = np.loadtxt(ebossdir+'nbarELG22v5_10_7.dat').transpose()
	if sampl == '21p22':
		f = fitsio.read(dirfits+'eboss21'+'.'+ver+'.latest.fits')
		f2 = fitsio.read(dirfits+'eboss22'+'.'+ver+'.latest.fits')
		
	else:
		f = fitsio.read(dirfits+'eboss'+sampl+'.'+ver+'.latest.fits') #read galaxy/quasar file
	
	no = 0
	zm = 0
	nt = 0
	avebs = 0
	sysstr = sys.split('_')
	
	for i in range (0,len(f)):
		z = f[i]['Z']
		c = 1
		gc = True
		if z > zmin and z < zmax and f[i]['Z_reliable'] == True and c == 1 and gc and f[i]['sector_TSR'] > compl and f[i]['isdupl'] == False:
			no += 1
			w = 1.
			#if wm == '':
			zind = int(z/.01)
			wfkp =ffkp[-1][zind]
			w = w*wfkp
			
			ra,dec = f[i]['ra'],f[i]['dec']
			th,phi = radec2thphi(ra,dec)
			pix = ang2pix_ring(nside,th,phi)
			sysv = mapl[pix]			

			bins = int((sysv-sysmin)*sysm)
					
			if bins >= 0 and bins < nsysbin:
				binng[bins] += 1.*w
			else:
				bs += w #count numbers outside of sysmin/sysmax
				avebs += w*sysv
			zm += w*z
			nt += w
	if d2:
		for i in range (0,len(f2)):
			z = f2[i]['Z']
			c = 1
			gc = True
			if z > zmin and z < zmax and f2[i]['Z_reliable'] == True and c == 1 and gc and f2[i]['sector_TSR'] > compl and f2[i]['isdupl'] == False:
				no += 1
				w = 1.
				#if wm == '':
				zind = int(z/.01)
				wfkp =ffkp[-1][zind]
				w = w*wfkp
				ra,dec = f2[i]['ra'],f2[i]['dec']
				th,phi = radec2thphi(ra,dec)
				pix = ang2pix_ring(nside,th,phi)
				sysv = mapl[pix]			

				bins = int((sysv-sysmin)*sysm)
					
				if bins >= 0 and bins < nsysbin:
					binng[bins] += 1.*w
				else:
					bs += w #count numbers outside of sysmin/sysmax
					avebs += w*sysv
				zm += w*z
				nt += w
			
	if bs > 0:
		avebs = avebs/bs
	print avebs,sysmin		
	print 'total number, weighted number'
	print no,nt
	print 'mean redshift'
	print zm/nt

	print 'total number of randoms/objects '+str(nr)+'/'+str(nt)
	print 'number of randoms/objects outside tested range '+str(bsr)+'/'+str(bs)			
	ave = nt/nr
	print 'average number of objects per random is '+ str(ave)
	fs = open(ebossdir+'n'+'geboss'+sampl+'_'+ver+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'v'+sys+'.dat','w')
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



def nzra(reg='SGC',zb=0.01):
	from matplotlib import pyplot as plt
	import fitsio
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	fr = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.ran.fits')
	offl = [0,0,0.0,0.0]
	cl = ['r','b','purple','green']
	hl = []
	w = ( (f['Z'] > 0.6) & (f['Z'] < 1.1))
	#hist = np.histogram((0.112*f[w]['rz']-f[w]['gr']),bins=20,normed=False)
	#hist = np.histogram((-f[w]['rz']+.218*f[w]['gr']),bins=20,normed=False)	
	hist = np.histogram(f['Z'],bins=20,normed=False)	
	hbins = hist[1]
	print(hbins)
	histm = hist[0]/float(len(fr))
	binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
	xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
	for i in range(0,2):
		ramin = i*25
		ramax = (i+1)*25
		w = ((f['RA'] > ramin) & (f['RA'] < ramax) & (f['Z'] > 0.6) & (f['Z'] < 1.1))
		wr = ((fr['RA'] > ramin) & (fr['RA'] < ramax))
		#print(len(f[w]))
		#plt.hist(f[w]['Z'],bins=20,normed=False,histtype='step')
		hist = np.histogram(f[w]['Z'],bins=hbins,normed=False)
		#hist = np.histogram((0.112*f[w]['rz']-f[w]['gr']),bins=hbins,normed=False)
		#hist = np.histogram((-f[w]['rz']+.218*f[w]['gr']),bins=hbins,normed=False)
		#plt.hist((-f[w]['rz']+.218*f[w]['gr'])/float(len(fr[wr])),bins=20,normed=False,histtype='step')
		
		
		xl = xlm+offl[i]
		#print(hist[1])
		#print(xl)
		#print(len(xl),len(hist[0]))
		plt.plot(xl,hist[0]/float(len(fr[wr])),color=cl[i])
		print(ramin,ramax,len(f[w]),len(f[w])/float(len(fr[wr])))
		#print(sum(hist[0])/float(len(fr[wr])))
		#print( hist[0]/float(len(fr[wr])))
		hl.append(hist[0]/float(len(fr[wr])))
		#wc = (w & (-f['rz']+.218*f['gr'] < -0.581))
		#print(len(f[wc])/float(len(fr[wr])))
	plt.show()
	#plt.plot(xlm,hl[3]-hl[1])
	#plt.plot(xlm,hl[2]-hl[1])
	#plt.show()
	#zl = np.arange(0.605,1.1,.01)
	#nzl = np.zeros(50)
	
def nzdr(reg='SGC',zb=0.01):
	from matplotlib import pyplot as plt
	import fitsio
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	fr = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.ran.fits')
	drl = np.unique(f['decals_dr'])
	offl = [0,0,0.0,0.0]
	cl = ['r','b','purple','green']
	hl = []
	w = ( (f['Z'] > 0.6) & (f['Z'] < 1.1))
	#hist = np.histogram((0.112*f[w]['rz']-f[w]['gr']),bins=20,normed=False)
	#hist = np.histogram((-f[w]['rz']+.218*f[w]['gr']),bins=20,normed=False)	
	hist = np.histogram(f['Z'],bins=20,normed=False)	
	hbins = hist[1]
	print(hbins)
	histm = hist[0]/float(len(fr))
	binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
	xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
	for i in range(0,len(drl)):
		w = (f['decals_dr'] == drl[i])
		wr = (fr['decals_dr'] == drl[i])
		hist = np.histogram(f[w]['Z'],bins=hbins,normed=False)
		
		
		xl = xlm
		plt.plot(xl,hist[0]/float(len(fr[wr])),color=cl[i])
		print(len(f[w]),len(f[w])/float(len(fr[wr])))
		hl.append(hist[0]/float(len(fr[wr])))
	plt.show()
	#plt.plot(xlm,hl[3]-hl[1])
	#plt.plot(xlm,hl[2]-hl[1])
	#plt.show()
	#zl = np.arange(0.605,1.1,.01)
	#nzl = np.zeros(50)

def mknznorm(reg='NGC',dz=0.001):
	from matplotlib import pyplot as plt
	import fitsio
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	zl = np.arange(.6,1.11,dz)
	zhist = np.histogram(f['Z'],bins=zl,normed=True)
	print(sum(zhist[0])*dz)
	fo = open('nznormELG'+reg+'.dat','w')
	for i in range(0,len(zl)-1):
		zb = zl[i] +dz/2.
		fo.write(str(zb)+' '+str(zhist[0][i])+'\n')
	fo.close()
	return True

def mkwrp(reg,dr=1.,rmin=1.,rmax=300.,zeff=.85):		
	from Cosmo import distance
	from photoztools import calcwxi
	d = distance(.31,.69)
	rz = d.dc(zeff) #in radians theta = rp/rz
	rp = rmin
	fo = open('/Users/ashleyross/eBOSS/xirpELG'+reg+'mod.dat','w')
	while rp < rmax:
		th = rp/rz*180./pi
		w = calcwxi('/Users/ashleyross/eBOSS/nznormELG'+reg,th,cosm='BOSS')
		fo.write(str(rp)+' '+str(w)+'\n')
		rp += dr
		print rp
	fo.close()
	return True
 
def calc_w024(r,reg='SGC',dmu=0.01):
	d = np.loadtxt('/Users/ashleyross/eBOSS/xirpELG'+reg+'mod.dat').transpose()
	rpl = d[0]
	mu = dmu/2.
	xi0 = 0
	xi2 = 0
	xi4 = 0
	while mu < 1.:
		rp = sqrt(1.-mu)*r
		if rp < 1.:
			xirp = d[1][0]
		else:
			bl = int(rp/1.)
			bh = bl + 1
			fac = rp/1.-bl
			xirp = d[1][bh]*fac+(1.-fac)*d[1][bl]
		xi0 += dmu*xirp
		xi2 += dmu*5/2.*P2(mu)*xirp
		xi4 += dmu*9/2.*P4(mu)*xirp
		mu += dmu
	return xi0,xi2,xi4	

def mkw024(reg='SGC',dmu=0.01,rmin=10,rmax=300,dr=1.):
	d = np.loadtxt('/Users/ashleyross/eBOSS/xirpELG'+reg+'mod.dat').transpose()
	fo = open('wrp024ELG'+reg+'.dat','w')
	rpl = d[0]
	r = rmin
	
	while r < rmax:
		mu = dmu/2.
		xi0 = 0
		xi2 = 0
		xi4 = 0
	
		while mu < 1.:
			rp = sqrt(1.-mu)*r
			if rp < 1.:
				xirp = d[1][0]
			else:
				bl = int(rp/1.)
				if bl >= len(d[1])-2:
					bh = bl
					print(r,rp,mu)
				else:
					bh = bl + 1
				fac = rp/1.-bl
				xirp = d[1][bh]*fac+(1.-fac)*d[1][bl]
			xi0 += dmu*xirp
			xi2 += dmu*5/2.*P2(mu)*xirp
			xi4 += dmu*9/2.*P4(mu)*xirp
			mu += dmu
		fo.write(str(r)+' '+str(xi0)+' '+str(xi2)+' '+str(xi4)+'\n')
		#print r
		r += dr
	fo.close()	
	return True

def mkw024_data(reg='SGC',dmu=0.01,bs=8,rmax=250):
	d = np.loadtxt('/Users/ashleyross/eBOSS/xirprpgebossELG_'+reg+'4_mz0.6xz1.1fkp1st0.dat').transpose()
	fo = open('/Users/ashleyross/eBOSS/wrp024ELG'+reg+'_data'+str(bs)+'st0.dat','w')
	rpl = d[0]
	nb = int(rmax/bs)
	for i in range(0,nb):
		mu = dmu/2.
		rbmin = float(i*bs)
		rbmax = rbmin+bs
		xi0 = 0
		xi2 = 0
		xi4 = 0
	
		while mu < 1.:
			rpmin = sqrt(1.-mu)*rbmin-.5
			if rpmin < 4.:
				print(rbmin,mu)
			rpmax = sqrt(1.-mu)*rbmax+.5
			w = (rpl > rpmin) & (rpl <= rpmax)
			#print(mu,rpmin,rpmax,rbmin,rbmax,sum(w))
			xirp = np.mean(d[1][w])
			xi0 += dmu*xirp
			xi2 += dmu*5/2.*P2(mu)*xirp
			xi4 += dmu*9/2.*P4(mu)*xirp
			mu += dmu
		r = i*bs+bs/2.
		fo.write(str(r)+' '+str(xi0)+' '+str(xi2)+' '+str(xi4)+'\n')
		#print r
		
	fo.close()	
	return True



def nzdepth(reg='NGC',gmax=30,rzmin=0,rzmax=2,grmin=0,grmax=2,sldmin=0,rmax=30,zmax=30):
	from matplotlib import pyplot as plt
	import fitsio
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	#fr = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.ran.fits')
	#drl = np.unique(f['decals_dr'])
	offl = [0,0,0.0,0.0]
	cl = ['r','b','purple','green']
	hl = []
	w = ( (f['Z'] > 0.6) & (f['Z'] < 1.1))
	#hist = np.histogram((0.112*f[w]['rz']-f[w]['gr']),bins=20,normed=False)
	#hist = np.histogram((-f[w]['rz']+.218*f[w]['gr']),bins=20,normed=False)	
	hist = np.histogram(f['Z'],bins=20,normed=False)	
	hbins = hist[1]
	print(hbins)
	
	binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
	xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
	chunkl = np.unique(f['chunk'])
	for chunk in chunkl:
		print(chunk)
		if reg == 'NGC':
			sld = f['rz']-0.637*f['gr']
		if reg == 'SGC':
			sld = f['rz']-0.218*f['gr']	
		rmag = -1.*(f['gr'] - f['g'])
		zmag = -1.*(f['rz']-rmag)	
		plt.hist(sld,bins='auto')
		plt.show()
		w = ((f['galdepth_g']+f['galdepth_r']*2+f['galdepth_z']*10 < 700) & (f['g'] < gmax) & (sld > sldmin) & (f['gr'] > grmin) & (f['gr'] < grmax) & (f['rz'] > rzmin) & (f['rz'] < rzmax) & (f['chunk'] == chunk) & (rmag < rmax) & (zmag < zmax))
		hist = np.histogram(f[w]['Z'],bins=hbins,normed=True)
		xl = xlm
		plt.plot(xl,hist[0],color=cl[0])
		print(len(f[w]))
		hl.append(hist[0])
		w = ((f['galdepth_g']+f['galdepth_r']*2+f['galdepth_z']*10 > 700) & (f['g'] < gmax) & (sld > sldmin) & (f['gr'] > grmin) & (f['gr'] < grmax) & (f['rz'] > rzmin) & (f['rz'] < rzmax) & (f['chunk'] == chunk) & (rmag < rmax) & (zmag < zmax))
		hist = np.histogram(f[w]['Z'],bins=hbins,normed=True)
		xl = xlm
		plt.plot(xl,hist[0],color=cl[1])
		print(len(f[w]))
		hl.append(hist[0])
	
		plt.show()
	#plt.plot(xlm,hl[3]-hl[1])
	#plt.plot(xlm,hl[2]-hl[1])
	#plt.show()
	#zl = np.arange(0.605,1.1,.01)
	#nzl = np.zeros(50)

def plotcol(reg='NGC',gmax=30):
	from matplotlib import pyplot as plt
	import fitsio
	
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	chunkl = np.unique(f['chunk'])
	for chunk in chunkl:
		w = ((f['Z'] > 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(f[w]['rz'],f[w]['gr'],'ko')

		w = ((f['Z'] < 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(f[w]['rz'],f[w]['gr'],'ro')
		plt.show()
		w = ((f['Z'] > 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(f[w]['g'],f[w]['gr'],'ko')

		w = ((f['Z'] < 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(f[w]['g'],f[w]['gr'],'ro')
		plt.show()
		w = ((f['Z'] > 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(f[w]['g'],f[w]['rz'],'ko')

		w = ((f['Z'] < 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(f[w]['g'],f[w]['rz'],'ro')
		plt.show()

		r = -1.*(f['gr'] - f['g'])
		plt.hist(r,bins='auto')
		plt.show()
		
		w = ((f['Z'] > 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(r[w],f[w]['gr'],'ko')

		w = ((f['Z'] < 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(r[w],f[w]['gr'],'ro')
		plt.show()

		zmag = -1.*(f['rz']-r)
		plt.hist(zmag,bins='auto')
		plt.xlabel('z band magnitude')
		plt.show()
		
		w = ((f['Z'] > 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(zmag[w],f[w]['rz'],'ko')

		w = ((f['Z'] < 0.7) & (f['chunk'] == chunk) & (f['g'] < gmax))
		plt.plot(zmag[w],f[w]['rz'],'ro')
		plt.show()

	
def gmaghist(reg='SGC'):
	from matplotlib import pyplot as plt
	import fitsio
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	fr = fitsio.read(dir+'eBOSS_ELG_clustering_'+reg+'_v4.ran.fits')
	drl = np.unique(f['decals_dr'])
	print(drl)
	offl = [0,0,0.0,0.0]
	cl = ['r','b','purple','green']
	hl = []
	w = ( (f['Z'] > 0.6) & (f['Z'] < 1.1))
	#hist = np.histogram((0.112*f[w]['rz']-f[w]['gr']),bins=20,normed=False)
	#hist = np.histogram((-f[w]['rz']+.218*f[w]['gr']),bins=20,normed=False)	
	hist = np.histogram(f['g'],bins=20,normed=False)	
	hbins = hist[1]
	print(hbins)
	histm = hist[0]/float(len(fr))
	binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
	xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
	plt.plot(xlm,hist[0])
	plt.show()
	for i in range(0,len(drl)):
		w = (f['decals_dr'] == drl[i])
		wr = (fr['decals_dr'] == drl[i])
		hist = np.histogram(f[w]['g'],bins=hbins,normed=False)
		
		
		xl = xlm
		plt.plot(xl,hist[0]/float(len(fr[wr])),color=cl[i])
		print(len(f[w]),len(f[w])/float(len(fr[wr])))
		hl.append(hist[0]/float(len(fr[wr])))
	plt.show()
	#plt.plot(xlm,hl[3]-hl[1])
	#plt.plot(xlm,hl[2]-hl[1])
	#plt.show()
	#zl = np.arange(0.605,1.1,.01)
	#nzl = np.zeros(50)
	


def bricks():
	from matplotlib import pyplot as plt
	import fitsio
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_full_ALL_v4.dat.fits')
	fr = fitsio.read(dir+'eBOSS_ELG_full_ALL_v4.ran.fits')
	bricks = np.unique(f['brickname'])
	fo = open('brickstats.dat','w')
	bd = []
	bzf = []
	mfrac = 1
	xfrac = 0
	mbz = 1
	xbz = 0
	n = 0
	for brick in bricks:
		w = (f['brickname'] == brick)
		wr = (fr['brickname'] == brick)
		g = f[w]
		r = fr[wr]
		if len(r) > 0:
			frac = len(g)/float(len(r))
			bd.append(frac)
			wfail = (g['IMATCH'] == 7)
			gzf = g[wfail]
			fbz = len(gzf)/float(len(g))
			bzf.append(fbz)
			if frac > xfrac:
				xfrac = frac
				print(brick,xfrac,n)
			if frac < mfrac:
				mfrac = frac
				print(brick,mfrac,n)
			if fbz > xbz:
				xbz = fbz
				print(brick,xbz,n)
			if fbz < mbz:
				mbz = fbz
				print(brick,mbz,n)
		else:
			print(brick,len(g),len(r))
		fo.write(brick+' '+str(len(g))+' '+str(len(r))+' '+str(len(gzf))+'\n')			
		n += 1	
	fo.close()
	#plt.hist(bd,bins=30)
	#plt.show()
	#plt.hist(bzf,bins=30)
	#plt.show()
	return bd,bzf

def brick_zmag():
	from matplotlib import pyplot as plt
	import fitsio
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/4/'
	f = fitsio.read(dir+'eBOSS_ELG_full_ALL_v4.dat.fits')
	fr = fitsio.read(dir+'eBOSS_ELG_full_ALL_v4.ran.fits')
	bricks = np.unique(f['brickname'])
	fo = open('brickstats_zmag.dat','w')
	bd = []
	bzf = []
	mfrac = 1
	xfrac = 0
	mbz = 1
	xbz = 0
	n = 0
	for brick in bricks:
		w = (f['brickname'] == brick)
		gals = f[w]
		rmag = -1.*(gals['gr'] - gals['g'])
		zmag = -1.*(gals['rz']-rmag)	
		zm = np.mean(zmag)
		bd.append(zm)
		fo.write(brick+' '+str(len(gals))+' '+str(zm)+'\n')			
		n += 1	
	fo.close()
	plt.hist(bd,bins='auto')
	plt.show()
	#plt.hist(bzf,bins=30)
	#plt.show()
	return True
	#return bd,bzf

def brick_zmag_ff():
	from matplotlib import pyplot as plt
	d = np.loadtxt('brickstats_zmag.dat',dtype={'names':('brick','ngal','meanz'),'formats':('S8','f4','f4')})
	plt.hist(d['meanz'],bins='auto')
	plt.show()
	return True

def exttest(reg='SGC'):
	import healpy as hp
	from healpix import radec2thphi
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'
	f = fitsio.read(dir+'4/eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	thphi = radec2thphi(f['RA'],f['DEC'])	
	r = hp.Rotator(coord=['C','G'])
	thphiG = r(thphi[0],thphi[1])
	pix = hp.ang2pix(1024,thphiG[0],thphiG[1])
	emap = fitsio.read('ebv_lhd.hpx.fits')['EBV']
	ebvl = []
	for i in range(0,len(pix)):
		ebvl.append(emap[pix[i]])
	ebvl = np.array(ebvl)
	de = ebvl-f['ebv']
	ral = f['RA']
	w = (ral > 180)
	ral[w] -= 360
	plt.scatter(ral,f['DEC'],c=de*(3.214-2.165),edgecolors='face')
	plt.colorbar()
	plt.title(r'$\Delta g-r$ based on Lenz et al. E(B-V)')
	plt.show()
	plt.savefig(reg+'gr_lenz.png')	
	plt.scatter(ral,f['DEC'],c=de*(2.165-1.211),edgecolors='face')
	plt.colorbar()
	plt.title(r'$\Delta r-z$ based on Lenz et al. E(B-V)')
	plt.show()	
	plt.savefig(reg+'rz_lenz.png')
	return True

def debv(reg='SGC',debv=0.01):
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'
	f = fitsio.read(dir+'4/eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
	dgr = debv*(3.214-2.165)
	drz = debv*(2.165-1.211)
	dg = debv*3.214
	n = 0
	ng = 0
	nc1 = 0
	nc2 = 0
	nc3 = 0
	nc4 = 0
	for i in range(0,len(f)):
		gn = f[i]['g']-dg
		grn = f[i]['gr']-dgr
		rzn = f[i]['rz']-drz
		k = 1
		if reg == 'SGC':
			if gn > 22.825:
				k = 0
				ng += 1
			if grn < -0.068*rzn +0.457:
				k = 0
				nc1 += 1
			if grn > 0.112*rzn+0.773:
				k = 0
				nc2 += 1
			if rzn < 0.218*grn+0.571:
				k = 0
				nc3 += 1
			if rzn > -0.555*grn+1.901:
				k = 0
				nc4 += 1
		if reg == 'NGC':
			if gn > 22.9:
				k = 0
				ng += 1
			if grn < -0.068*rzn +0.457:
				k = 0
				nc1 += 1
			if grn > 0.112*rzn+0.773:
				k = 0
				nc2 += 1
			if rzn < 0.637*grn+0.399:
				k = 0
				nc3 += 1
			if rzn > -0.555*grn+1.901:
				k = 0
				nc4 += 1
						
		n += k
	print(n/float(len(f)),n,ng,nc1,nc2,nc3,nc4)
	return True

def ngalvdebv():
	import healpy as hp
	from healpix import radec2thphi
	from optimize import fmin
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'
	regl = ['SGC','NGC']
	cl  =['b','r']
	for k in range(0,2):
		reg = regl[k]
		f = fitsio.read(dir+'4/eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
		thphi = radec2thphi(f['RA'],f['DEC'])	
		r = hp.Rotator(coord=['C','G'],deg=False)
		thphiG = r(thphi[0],thphi[1])
		pix = hp.ang2pix(1024,thphiG[0],thphiG[1])
		emap = fitsio.read('ebv_lhd.hpx.fits')['EBV']
		ebvl = []
		for i in range(0,len(pix)):
			ebvl.append(emap[pix[i]])
		ebvl = np.array(ebvl)
		de = ebvl-f['ebv']
		#denan = de[np.isnan(de)]
		wts=f[np.isfinite(de)]['WEIGHT_SYSTOT']
		w = (np.isfinite(de) & (de > -0.01) & (de < 0.015))
		denan = de[~w]
		de = de[w]
		
		be = [-0.06,-0.02,-0.01,-0.007,-0.005,-0.003,-0.001,0.001,0.005,0.01,0.02]
		hist = np.histogram(de,bins=10,normed=False)#,weights=wts)	
		hbins = hist[1]
		print(hbins)
		print(hist[0])
		binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
		#xlm = []
		#for i in range(0,len(be)-1):
		#	xlm.append((be[i]+be[i+1])/2.) 
		xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
	
		fr = fitsio.read(dir+'4/eBOSS_ELG_clustering_'+reg+'_v4.ran.fits')
		nr = len(fr)/float(len(f))
		print(nr)
		thphir = radec2thphi(fr['RA'],fr['DEC'])	
		thphiG = r(thphir[0],thphir[1])
		pixr = hp.ang2pix(1024,thphiG[0],thphiG[1])
		ebvlr = []
		for i in range(0,len(pixr)):
			ebvlr.append(emap[pixr[i]])
		ebvlr = np.array(ebvlr)
		der = ebvlr-fr['ebv']
		#dernan = der[np.isnan(der)]
		wr = (np.isfinite(der) & (der > -0.01) & (der < 0.015))
		dernan = der[~wr]
		der = der[wr]
		
		histr = np.histogram(der,bins=hbins,normed=False)
		nr = len(der)/float(len(de))
		print(nr)
		
		print(histr[0])
		print(hist[0]/histr[0]*nr)
		print('nan ratio is',len(denan)/float(len(dernan))*nr,len(denan),len(dernan)/float(len(fr)))
		fo = open('nELG'+reg+'vsdEBV.dat','w')
		for i in range(0,len(xlm)):
			fo.write(str(xlm[i])+' '+str(hist[0][i]/histr[0][i].astype('float')*nr)+' '+str(sqrt(hist[0][i])/histr[0][i].astype('float')*nr)+'\n')
		fo.close()
		lf = linfit(xlm,hist[0]/histr[0].astype('float')*nr,np.sqrt(hist[0])/histr[0].astype('float')*nr)
		inl = np.array([1.,-10.])
		b0,m0 = fmin(lf.chilin,inl)
		print(b0,m0)
		print( lf.chilin((b0,m0)))
		plt.plot(xlm,b0+m0*xlm,'--',color=cl[k])
		plt.errorbar(xlm,hist[0]/histr[0].astype('float')*nr,np.sqrt(hist[0])/histr[0].astype('float')*nr,color=cl[k])

	#plt.xlim(-0.05,0.02)
	plt.xlabel(r'$\Delta$E(B-V) (Lenz et al - SFD)')
	plt.ylabel(r'$n_{\rm gal}/\langle n_{\rm gal} \rangle$')
	#plt.text(-0.04,1.1,'NGC',color='r')
	#plt.text(-0.04,1.08,'SGC',color='b')
	plt.show()
	#plt.savefig('nELGvsdeltaEBV.png')
	return True

def nzsplitdebv():
	import healpy as hp
	from healpix import radec2thphi
	from optimize import fmin
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'
	regl = ['SGC','NGC']
	cl  =['b','r']
	for k in range(0,2):
		reg = regl[k]
		f = fitsio.read(dir+'4/eBOSS_ELG_clustering_'+reg+'_v4.dat.fits')
		thphi = radec2thphi(f['RA'],f['DEC'])	
		r = hp.Rotator(coord=['C','G'],deg=False)
		thphiG = r(thphi[0],thphi[1])
		pix = hp.ang2pix(1024,thphiG[0],thphiG[1])
		emap = fitsio.read('ebv_lhd.hpx.fits')['EBV']
		ebvl = []
		for i in range(0,len(pix)):
			ebvl.append(emap[pix[i]])
		ebvl = np.array(ebvl)
		de = ebvl-f['ebv']
		#denan = de[np.isnan(de)]
		wts=f[np.isfinite(de)]['WEIGHT_SYSTOT']
		if reg == 'NGC':
			w = (np.isfinite(de) & (de > -0.01) & (de < 0.00) & (f['chunk']=='eboss23'))
		else:
			w = (np.isfinite(de) & (de > -0.01) & (de < 0.00))
		
		hist = np.histogram(f[w]['Z'],bins=20,normed=True)#,weights=wts)	
		hbins = hist[1]
		print(hbins)
		print(hist[0])
		binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
		xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
		print(len(f[w]))
		plt.plot(xlm,hist[0],color=cl[0])
		if reg == 'NGC':
			w = (np.isfinite(de) & (de > 0.00) & (f['chunk']=='eboss23'))
		else:
			w = (np.isfinite(de) & (de > 0.00))	
		hist = np.histogram(f[w]['Z'],bins=hbins,normed=True)
		plt.plot(xlm,hist[0],color=cl[1])
		print(len(f[w]))
		plt.show()
	return True


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


def brickanalysis(ebossdir=ebossdir):
	d = np.loadtxt(ebossdir+'brickstats-dr3.dat',dtype={'names':('brick','ntar','nran','nzfail','ra','dec','nexp_g','nexp_r','nexp_z','nobjs'),'formats':('S8','f4','f4','f4','f4','f4','f4','f4','f4','f4')})
	print(min(d['ntar']),min(d['nran']),min(d['nzfail']))
	w = (d['nran'] > 0)
	d = d[w]
	nexp_t = d['nexp_g']+d['nexp_r']+d['nexp_z']
	print(max(nexp_t))
	print(max(d['ntar']),max(d['nran']),max(d['nzfail']),max(d['nexp_g']),max(d['nexp_r']),max(d['nexp_z']),max(d['nobjs']),min(d['nobjs']))	
	print(min(d['ntar']/d['nran']),max(d['ntar']/d['nran']))
	print(min(d['nzfail']/d['ntar']),max(d['nzfail']/d['ntar']))
	from matplotlib import pyplot as plt
	plt.plot(d['ntar']/d['nran'],d['nran'],'ko')
	plt.show()
	wr = (d['nran'] > 100)
	print(len(d[wr])/float(len(d)))
	drancut = d[wr]
	plt.hist(np.log(drancut['ntar']/drancut['nran']),bins=30)
	plt.show()
	plt.hist(drancut['nzfail']/drancut['ntar'],bins=30)
	plt.show()
	wzf = (drancut['nzfail']/drancut['ntar'] < 0.3)
	print(len(drancut[wzf])/float(len(drancut)))
	ral = []
	decl =[]
	for ng in range(1,22):
		w = (d['nexp_g'] ==ng)
		if len(d[w])>0:
			print(ng,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
	for nr in range(1,20):
		w = (d['nexp_r'] ==nr)
		if len(d[w])>0:
			print(nr,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
	for nz in range(1,17):
		w = (d['nexp_z'] ==nz)
		if len(d[w])>0:
			print(nz,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
	ntl = []
	densl = []
	for nt in range(3,47):
		w = (nexp_t ==nt)
		if len(d[w])>0:
			print(nt,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
			#ntl.append(nt)
			#densl.append(np.mean(d[w]['ntar']/d[w]['nran']))
	#plt.plot(ntl,densl)
	#plt.show()
	dgz = d['nexp_g']-d['nexp_z']
	print('difference between number of g and z exposures')
	print(min(dgz),max(dgz))
	for ndgz in range(-5,18):
		w = (dgz ==ndgz)
		if len(d[w])>0:
			print(ndgz,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
			ntl.append(ndgz)
			densl.append(np.mean(d[w]['ntar']/d[w]['nran']))
	plt.plot(ntl,densl)
	plt.show()

	ntl = []
	densl = []		
	drz = d['nexp_r']-d['nexp_z']
	print('difference between number of r and z exposures')
	print(min(drz),max(drz))
	for ndrz in range(-10,18):
		w = (drz ==ndrz)
		if len(d[w])>0:
			print(ndrz,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
			ntl.append(ndrz)
			densl.append(np.mean(d[w]['ntar']/d[w]['nran']))
	plt.plot(ntl,densl)
	plt.show()

	ntl = []
	densl = []
	dgr = d['nexp_g']-d['nexp_r']
	print('difference between number of g and r exposures')
	print(min(dgr),max(dgr))
	for ndgr in range(-12,18):
		w = (dgr ==ndgr)
		if len(d[w])>0:
			print(ndgr,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
			ntl.append(ndgr)
			densl.append(np.mean(d[w]['ntar']/d[w]['nran']))
	plt.plot(ntl,densl)
	plt.show()




	nx = max(d['nobjs'])
	for i in range(0,20):
		mx = (i+1)/20.*nx
		mm = i/20.*nx
		w = ((d['nobjs'] > mm) & (d['nobjs'] < mx))	
		if sum(d[w]['nran']) > 0:
			print(mx,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran']))
	w = ((dgz < -2) | (dgz > 5))
	#fo = open('brickzexpg7.dat','w')
	dg = d[w]
	#for i in range(0,len(dg)):
	#	fo.write(dg[i]['brick']+'\n')
	#fo.close()
	print(sum(d[w]['nran']),sum(d[w]['ntar']))
	plt.plot(d[w]['ra'],d[w]['dec'],'bo')
	plt.xlim(0,50)
	#plt.show()
	
	w = (dgr > 4)
	#fo = open('brickzexpg7.dat','w')
	dg = d[w]
	#for i in range(0,len(dg)):
	#	fo.write(dg[i]['brick']+'\n')
	#fo.close()
	print(sum(d[w]['nran']),sum(d[w]['ntar']))
	plt.plot(d[w]['ra'],d[w]['dec'],'ro')
	plt.xlim(0,50)
	plt.show()
	
	w = ((dgz < -2) | (dgz > 5) | (dgr > 4))
	print(sum(d[w]['nran']),sum(d[w]['ntar']))
	fo = open('brickexpdiffext.dat','w')
	for brick in d[w]['brick']:
		fo.write(brick+'\n')
	fo.close()
	
# 	f = fitsio.read('survey-bricks-dr3.fits.gz')
# 	for brick in drancut[~wzf]['brick']:
# 		if np.isin(brick,f['brickname']):
# 			wbrickd = (f['brickname'] == brick)
# 			brickd = f[wbrickd]
# 			ral.append(brickd['ra'])
# 			decl.append(brickd['dec'])
# 	print(len(ral),len(drancut[~wzf]))
# 	plt.plot(ral,decl,'ko')
# 	plt.xlim(0,50)
# 	plt.ylim(-5,5)
# 	plt.show()
# 
# 	we = (drancut['ntar']/drancut['nran'] > 0.065)
# 	print(len(drancut[we])/float(len(drancut)))
# 	ral = []
# 	decl =[]
# 	f = fitsio.read('survey-bricks-dr3.fits.gz')
# 	for brick in drancut[we]['brick']:
# 		if np.isin(brick,f['brickname']):
# 			wbrickd = (f['brickname'] == brick)
# 			brickd = f[wbrickd]
# 			ral.append(brickd['ra'])
# 			decl.append(brickd['dec'])
# 	print(len(ral),len(drancut[we]))
# 	plt.plot(ral,decl,'ko')
# 	plt.xlim(0,50)
# 	plt.ylim(-5,5)
# 
# 	plt.show()

def matchbrick(ebossdir=ebossdir):
	d = np.loadtxt(ebossdir+'brickstats.dat',dtype={'names':('brick','ntar','nran','nzfail'),'formats':('S8','f4','f4','f4')})
	f = fitsio.read('survey-bricks-dr3.fits.gz')
	fo = open('brickstats-dr3.dat','w')
	n = 0
	for i in range(0,len(d)):
		brick = d[i]['brick']
		if np.isin(brick,f['brickname']):
			wbrickd = (f['brickname'] == brick)
			brickd = f[wbrickd]
			ra = brickd['ra'][0]
			dec = brickd['dec'][0]
			nexp_g = brickd['nexp_g'][0]
			nexp_r = brickd['nexp_r'][0]
			nexp_z = brickd['nexp_z'][0]
			nobjs = brickd['nobjs'][0]
		else:
			ra = 0
			dec = 0
			nexp_g = 0
			nexp_r = 0
			nexp_z = 0
			nobs = 0
		fo.write(brick+' '+str(d[i]['ntar'])+' '+str(d[i]['nran'])+' '+str(d[i]['nzfail'])+' '+str(ra)+' '+str(dec)+' '+str(nexp_g)+' '+str(nexp_r)+' '+str(nexp_z)+' '+str(nobjs)+'\n')
		print(i)
	fo.close()
	return True
		
def plotvssys_simp(xl,yl,el,sys):
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	from optimize import fmin
	from xitools_eboss import linfit
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
	#plt.ylim(.7,1.19)
	plt.text(min(xl)+0.1*(max(xl)-min(xl)),1.1,r'$\chi^2$ null ='+str(chin)[:4],color='k')
	plt.text(min(xl)+0.1*(max(xl)-min(xl)),1.08,r'$\chi^2$ lin ='+str(chilin)[:4],color='k')
	plt.show()
	#plt.title(r'galaxy density vs. $i$-band depth for v0.7 eboss QSOs, 0.9 < z < 2.2')
	#pp.savefig()
	#pp.close()
	return True


def plotxi2ELG(reg,bs='8st0',v='test',wm='fkpcpgdepth',wm1='fkpgdepth',zmin=.6,zmax=1.1,l1='',l2=''):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi2ELG'+reg+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	d1 = np.loadtxt(ebossdir+'xi2gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm1+bs+'.dat').transpose()

	d2 = np.loadtxt(ebossdir+'xi2gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi2Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	#plt.plot(dt[0],dt[0]*dt[1]*.9,'k:')
	plt.plot(d1[0],d1[0]*(d1[1]))
	plt.plot(d2[0],d2[0]*(d2[1]))
	plt.xlim(10,200)
	#plt.ylim(-30,60)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s \xi_2(s)$ ($h^{-1}$Mpc)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of ELGs in, '+reg+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		plt.legend(labels=['model','eBOSS'])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True


def plotxiELGcomb(reg = 'ALL',mom=0,bs='8st0',v='test',rec='',zmin=.7,zmax=1.1,l1='',l2='',modplot=True):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcomb'+reg+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = 1.2
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss21'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = 3.1
	d2 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss22'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a3 = 2.5
	d3 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss23'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a5 = 1.4
	d5 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss25'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	if reg == 'ALL':
		dc = (d1*a1+d2*a2+d3*a3+d5*a5)/(a1+a2+a3+a5)
	if reg == 'SGC':
		dc = (d1*a1+d2*a2)/(a1+a2)
	if reg == 'NGC':
		dc = (d3*a3+d5*a5)/(a3+a5)		
	if rec == '_rec':
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.593.02.04.015.01.0.dat').transpose()
	else:
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*.9,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*.9,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		plt.plot(dc[0],dc[0]**2.*dc[1])
	if mom ==2:	
		plt.plot(dc[0],dc[0]*dc[1])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-30,60)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	if rec == '_rec':
		plt.title(r'post-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	else:
		plt.title(r'pre-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		if modplot:
			if reg == 'ALL':
				plt.legend(labels=['model','eboss21+eboss22+eboss23+eboss25'])
			if reg == 'NGC':
				plt.legend(labels=['model','eboss23+eboss25'])
			if reg == 'SGC':
				plt.legend(labels=['model','eboss21+eboss22'])
		else:
			plt.legend(labels=[wm1,wm])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxiELGcombNS(mom=0,bs='8st0',v='test',rec='',zmin=.6,zmax=1.1,l1='',l2='',comb='',modplot=True):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcombNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .558
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .479
	dc = (d1*a1+d2*a2)/(a1+a2)
	if rec == '_rec':
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.593.02.04.015.01.0.dat').transpose()
	else:
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*.75,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*.75,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		plt.plot(dc[0],dc[0]**2.*dc[1])
	if mom ==2:	
		plt.plot(dc[0],dc[0]*dc[1])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-15,45)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-1.5,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	if rec == '_rec':
		plt.title(r'post-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	else:
		plt.title(r'pre-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		if modplot:
			#if reg == 'ALL':
			plt.legend(labels=['model','NGC+SGC'])
			#if reg == 'NGC':
			#	plt.legend(labels=['model','eboss23+eboss25'])
			#if reg == 'SGC':
			#	plt.legend(labels=['model','eboss21+eboss22'])
		else:
			plt.legend(labels=[wm1,wm])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxiELGcompNS(mom=0,bs='8st0',v='test',rec='',zmin=.6,zmax=1.1,l1='',l2='',modplot=True,comb=''):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcompNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = 1.2+3.1
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = 2.5+1.4
	dc = (d1*a1+d2*a2)/(a1+a2)
	if rec == '_rec':
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.593.02.04.015.01.0.dat').transpose()
	else:
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*.75,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*.75,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		plt.plot(d2[0],d2[0]**2.*d2[1])
		plt.plot(d1[0],d1[0]**2.*d1[1])
	if mom ==2:	
		plt.plot(d2[0],d2[0]*d2[1])
		plt.plot(d1[0],d1[0]*d1[1])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-15,45)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-1.5,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	if rec == '_rec':
		plt.title(r'post-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	else:
		plt.title(r'pre-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		if modplot:
			#if reg == 'ALL':
			plt.legend(labels=['model','NGCcomb','SGCcomb'])
			#if reg == 'NGC':
			#	plt.legend(labels=['model','eboss23+eboss25'])
			#if reg == 'SGC':
			#	plt.legend(labels=['model','eboss21+eboss22'])
		else:
			plt.legend(labels=['NGCcomb','SGCcomb'])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxiELGcompth(mom=0,reg='SGC',bs='8st0',v='4',rec='',zmin=.6,zmax=1.1,thfac=.75,l1='',l2='',wm='',modplot=True,comb='',angfac=2.,md='subang'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcompNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+bs+'.dat').transpose()
	ave = np.loadtxt(ebossdir+'xiave'+str(mom)+reg+'ELG_EZv'+v+bs+'.dat').transpose()
	rl = d1[0][:len(ave[0])]
	xil = d1[1][:len(ave[0])]
	if md == 'subang':
		print('subtracting angular clustering with factor '+str(angfac))
		wp = np.loadtxt(ebossdir+'wrp024ELG'+reg+'_data'+bs+'.dat').transpose()
		xil = xil-angfac*wp[mom/2+1][:len(ave[0])]
	if rec == '_rec':
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.593.02.04.015.01.0.dat').transpose()
	else:
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if wm == 'nosysshuff' or md == 'subang':
		dt = dt - angfac*np.loadtxt('/Users/ashleyross/eBOSS/wrp024ELG'+reg+'.dat').transpose()[mom/2+1]
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*thfac,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*thfac,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		
		plt.errorbar(d1[0][:len(ave[0])],d1[0][:len(ave[0])]**2.*xil,d1[0][:len(ave[0])]**2.*ave[2][:len(ave[0])],fmt='ko')
	if mom ==2:	
		
		plt.errorbar(d1[0][:len(ave[0])],d1[0][:len(ave[0])]*xil,d1[0][:len(ave[0])]*ave[2][:len(ave[0])],fmt='ko')

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-15*thfac/.75,45*thfac/.75)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-1.5,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	if rec == '_rec':
		plt.title(r'post-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	else:
		plt.title(r'pre-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	plt.show()
	return True


def plotxiELGcomp(mom=0,bs='8st0',v='test',rec='',zmin=.7,zmax=1.1,l1='',l2='',modplot=False):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcomp'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = 1.2
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss21'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = 3.1
	d2 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss22'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a3 = 2.5
	d3 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss23'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a5 = 1.4
	d5 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_eboss25'+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	dc = (d1*a1+d2*a2+d3*a3+d5*a5)/(a1+a2+a3+a5)
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if modplot:
		plt.plot(dt[0],dt[0]**2.*dt[1]*.9,'k:')
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom == 0:
		plt.plot(d1[0],d1[0]**2.*d1[1],d2[0],d2[0]**2.*d2[1],d3[0],d3[0]**2.*d3[1],d5[0],d5[0]**2.*d5[1])
	if mom == 2:
		plt.plot(d1[0],d1[0]*d1[1],d2[0],d2[0]*d2[1],d3[0],d3[0]*d3[1],d5[0],d5[0]*d5[1])
		
	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-30,60)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	if rec == '_rec':
		plt.title(r'post-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	else:
		plt.title(r'pre-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	#if l1 == '':
	#	if modplot:
	plt.legend(labels=['eboss21','eboss22','eboss23','eboss25'])
	#	else:
	#		plt.legend(labels=[wm1,wm])
	#else:
	#	plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True


def plotxiELG(reg,mom=0,bs='8st0',v='1',wm='fkp',wm1=False,zmin=.7,zmax=1.1,l1='',l2='',modplot=True):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELG'+reg+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	if wm1:
		d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm1+bs+'.dat').transpose()
	
	d2 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if modplot:
		plt.plot(dt[0],dt[0]**2.*dt[1]*.9,'k:')
	if wm1:
		plt.plot(d1[0],d1[0]**2.*(d1[1]))
	plt.plot(d2[0],d2[0]**2.*(d2[1]))
	plt.xlim(10,200)
	if mom == 0:
		plt.ylim(-30,60)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of ELGs in, '+reg+' '+str(zmin)+' < z < '+str(zmax))
	
	#if l1 == '':
	#	if modplot:
	#		plt.legend(labels=['model',wm1,wm])
	#	else:
	#		plt.legend(labels=[wm1,wm])
	#else:
	#	plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True
	pp = PdfPages(ebossdir+'xi0ELG'+chunk+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'xr.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.plot(dt[0],dt[0]*dt[1]/1.44/1.035,'k:')
	plt.plot(d1[0],d1[0]*(d1[1]))
	plt.xlim(0,200)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s\xi(s)$ ($h^{-1}$Mpc)',size=16)
	plt.title(r'Correlation function of ELGs in chunk, '+chunk+' '+str(zmin)+' < z < '+str(zmax))
	plt.legend(labels=['model','eBOSS'])

	pp.savefig()
	pp.close()

	return True

def plotxiELG_wedge(reg,mom=0,bs='8st0',v='4',wm='fkp',wm1=False,zmin=.6,zmax=1.1,l1='',l2='',modplot=True):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	#pp = PdfPages(ebossdir+'xi'+str(mom)+'ELG'+reg+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'_5mu.pdf')
	plt.clf()
	plt.minorticks_on()
	#gebossELG_NGC4_mz0.6xz1.1fkpmux0.2
	d0 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'mux0.2'+bs+'.dat').transpose()
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'mum0.2mux0.4'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'mum0.4mux0.6'+bs+'.dat').transpose()
	d3 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'mum0.6mux0.8'+bs+'.dat').transpose()
	d4 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'mum0.8'+bs+'.dat').transpose()
	plt.plot(d0[0],d0[0]**2.*(d0[1]))
	plt.plot(d0[0],d0[0]**2.*(d1[1]))
	plt.plot(d0[0],d0[0]**2.*(d2[1]))
	plt.plot(d0[0],d0[0]**2.*(d3[1]))
	plt.plot(d0[0],d0[0]**2.*(d4[1]))
	plt.xlim(10,200)
	plt.legend((r'$\mu<0.2$',r'$0.2 < \mu < 0.4$',r'$0.4 < \mu < 0.6$',r'$0.4 < \mu < 0.8$',r'$ \mu > 0.8$'))
	plt.ylim(-100,100)
	#if mom == 0:
	#	plt.ylim(-30,60)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of ELGs in $\mu$ bins, '+reg+' '+str(zmin)+' < z < '+str(zmax))
	
	#if l1 == '':
	#	if modplot:
	#		plt.legend(labels=['model',wm1,wm])
	#	else:
	#		plt.legend(labels=[wm1,wm])
	#else:
	#	plt.legend(labels=[l1,l2])
	plt.savefig(ebossdir+'xi'+str(mom)+'ELG'+reg+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'_5mu.png')
	#pp.close()
	return True


def plotxiELGr(reg,bs='8st0',v='test',wm='fkpcpgdepth',wm1='fkpgdepth',zmin=.6,zmax=1.1,l1='',l2='',modplot=True):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi0ELG'+reg+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'r.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	d1 = np.loadtxt(ebossdir+'xi0gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm1+bs+'.dat').transpose()

	d2 = np.loadtxt(ebossdir+'xi0gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if modplot:
		plt.plot(dt[0],dt[0]*dt[1]*.9,'k:')
	plt.plot(d1[0],d1[0]*(d1[1]))
	plt.plot(d2[0],d2[0]*(d2[1]))
	plt.xlim(10,200)
	#plt.ylim(-30,60)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s\xi(s)$ ($h^{-1}$Mpc)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of ELGs in, '+reg+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		if modplot:
			plt.legend(labels=['model',wm1,wm])
		else:
			plt.legend(labels=[wm1,wm])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True



def plotxiELGrec(reg='eboss21',bs='8st0',v='test',wm='fkp',zmin=.6,zmax=1.1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi0ELGrec'+reg+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	d1 = np.loadtxt(ebossdir+'xi0gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi0gebossELG_'+reg+v+'_rec_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.593.04.07.015.00.dat').transpose()
	dtrec = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.593.02.04.015.01.0.dat').transpose()
	plt.plot(dt[0],dt[0]**2.*dt[1]*.9,'k:')
	plt.plot(dt[0],dt[0]**2.*dtrec[1]*.9,'k--')
	plt.plot(d1[0],d1[0]**2.*(d1[1]))
	plt.plot(d2[0],d2[0]**2.*(d2[1]))
	plt.xlim(10,200)
	plt.ylim(-30,60)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of ELGs in '+reg+' '+str(zmin)+' < z < '+str(zmax))
	plt.legend(labels=['model','eBOSS'])

	pp.savefig()
	pp.close()


def plotxiELGv21p2223(bs='8st0',v='v5_10_7',wm='fkp',zmin=.75,zmax=1.1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi0ELG21p2223'+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	d1 = np.loadtxt(ebossdir+'xi0gebosselg_21p22'+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi0gebosselg_23'+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = 466.
	a2 = 165.
	wave = (d1[1]*a1+d2[1]*a2)/(a1+a2)
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	plt.plot(dt[0],dt[0]**2.*dt[1]*.9,'k:')
	plt.plot(d1[0],d1[0]**2.*(d1[1]))
	plt.plot(d1[0],d1[0]**2.*(d2[1]))
	plt.plot(d1[0],d1[0]**2.*wave)
	plt.xlim(10,200)
	plt.ylim(-30,60)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of ELGs , '+str(zmin)+' < z < '+str(zmax))
	plt.legend(labels=['model','chunks 21+22','chunk 23','weighted average'])

	pp.savefig()
	pp.close()

	return True

def plotxiELGNScompEZ(mom='0',bs='5st0',v='4',mini=4,maxi=40,wm='fkp',mumin=0,mumax=1,bfs=1.,bfn=1.,covv='v4',zr='0.6xz1.1',mur=''):
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
	norm = 1.
	#norm = 1./float(mumax-mumin)
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGNScompEZ'+zr+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	#dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	ds = np.loadtxt(ebossdir+'xi'+mom+'gebossELG_SGC'+v+'_mz'+zr+wm+muwd+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi'+mom+'gebossELG_NGC'+v+'_mz'+zr+wm+muwd+bs+'.dat').transpose()

	aves = np.loadtxt(ebossdir+'xiave'+mom+'SGCELG_EZ'+covv+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xiave'+mom+'NGCELG_EZ'+covv+bs+'.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov'+mom+'SGCELG_EZ'+covv+bs+'.dat')[mini:maxi,mini:maxi]#*.857/1.4
	covn = np.loadtxt(ebossdir+'cov'+mom+'NGCELG_EZ'+covv+bs+'.dat')[mini:maxi,mini:maxi]#/2.
	diffs = bfs*aves[1][mini:maxi]*norm**2.-ds[1][mini:maxi]
	facn = 1.
	facs = 1.

	chis = np.dot(np.dot(diffs,np.linalg.pinv(covs)),diffs)*facs
	diffn = bfn*aven[1][mini:maxi]*norm**2.-dn[1][mini:maxi]
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

def xiELGNSdiffEZ(mom='0',bs='8st0',v='4',mini=3,maxi=25,wm='fkp',mumin=0,mumax=1,bfs=1.,bfn=1.):
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
	#pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGNScompEZ'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	#dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	ds = np.loadtxt(ebossdir+'xi'+mom+'gebossELG_SGC'+v+'_mz0.6xz1.1'+wm+muwd+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi'+mom+'gebossELG_NGC'+v+'_mz0.6xz1.1'+wm+muwd+bs+'.dat').transpose()

	aves = np.loadtxt(ebossdir+'xiave'+mom+'SGCELG_EZ'+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xiave'+mom+'NGCELG_EZ'+bs+'.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov'+mom+'SGCELG_EZ'+bs+'.dat')[mini:maxi,mini:maxi]#*.857/1.4
	covn = np.loadtxt(ebossdir+'cov'+mom+'NGCELG_EZ'+bs+'.dat')[mini:maxi,mini:maxi]#/2.
	diffns = bfn*dn[1][mini:maxi]-ds[1][mini:maxi]
	facn = 1.
	facs = 1.

	chins = np.dot(np.dot(diffns,np.linalg.pinv(covs+covn)),diffns)*facs
	print chins,maxi-mini,ds[0][mini],ds[0][maxi]
	return True


def plotxiELGNScompQPM(mom='0',bs='8st0',v='3',mini=3,maxi=25,wm='fkp',mumin=0,mumax=1,bfs=1.,bfn=1.,rec='_rec'):
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
	pp = PdfPages(ebossdir+'xi'+rec+str(mom)+'ELGNScompQSO'+muw+v+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	#dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	ds = np.loadtxt(ebossdir+'xi'+mom+'gebossELG_SGC'+v+rec+'_mz0.6xz1.1'+wm+muwd+bs+'.dat').transpose()
	dn = np.loadtxt(ebossdir+'xi'+mom+'gebossELG_NGC'+v+rec+'_mz0.6xz1.1'+wm+muwd+bs+'.dat').transpose()

	aves = np.loadtxt(ebossdir+'xiave'+rec+mom+'SGCELG_MV'+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xiave'+rec+mom+'NGCELG_MV'+bs+'.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov'+rec+mom+'SGCELG_MV'+bs+'.dat')[mini:maxi,mini:maxi]#*.857/1.4
	covn = np.loadtxt(ebossdir+'cov'+rec+mom+'NGCELG_MV'+bs+'.dat')[mini:maxi,mini:maxi]#/2.
	diffs = bfs*aves[1][mini:maxi]*norm**2.-ds[1][mini:maxi]
	facn = 1.
	facs = 1.

	chis = np.dot(np.dot(diffs,np.linalg.pinv(covs)),diffs)*facs
	diffn = bfn*aven[1][mini:maxi]*norm**2.-dn[1][mini:maxi]
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

def plotxiELGNSbaofit(bs='5st0',v='4',a='',rec='',wm='fkp',mini=10,maxi=30,mom='0',covv='v4'):
	#Plots comparison between QSO clustering and best-fit BAO theory
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiELGNSbaofit'+v+rec+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dts = np.loadtxt(ebossdir+'ximodELGSGC'+v+rec+a+'.dat').transpose()
	dts_nb = np.loadtxt(ebossdir+'ximodELGSGC'+v+rec+a+'nobao.dat').transpose()
	dts_iso = dts[1]-dts[2]
	dtn = np.loadtxt(ebossdir+'ximodELGNGC'+v+rec+a+'.dat').transpose()
	dtn_nb = np.loadtxt(ebossdir+'ximodELGNGC'+v+rec+a+'nobao.dat').transpose()
	dtn_iso = dtn[1]-dtn[2]
	if rec == '_rec':
		covs = np.loadtxt(ebossdir+'cov_recon'+mom+'SGCELG'+bs+'.dat')
		covn = np.loadtxt(ebossdir+'cov_recon'+mom+'NGCELG'+bs+'.dat')
	if rec == '':	
		covs = np.loadtxt(ebossdir+'cov'+mom+'SGCELG_EZ'+bs+'.dat')
		covn = np.loadtxt(ebossdir+'cov'+mom+'NGCELG_EZ'+bs+'.dat')
	
	covi = np.linalg.pinv(covn)+np.linalg.pinv(covs)
	cov = np.linalg.pinv(covi)
	#cov = np.loadtxt('covxiNSQSO'+v+'mz0.9xz2.2'+bsc+'.dat')
	et = []
	ets = []
	etn = []
	for i in range(0,maxi):
		et.append(sqrt(cov[i][i]))
		etn.append(sqrt(covn[i][i]))
		ets.append(sqrt(covs[i][i]))
	etn = np.array(etn)[mini:maxi]
	ets = np.array(ets)[mini:maxi]	
	et = np.array(et)
	dtt_iso = (dtn_iso/etn**2.+dts_iso/ets**2.)/(1./etn**2.+1./ets**2.)	
	dsw = np.loadtxt(ebossdir+'xi0gebossELG_SGC'+v+rec+'_mz0.6xz1.1'+wm+bs+'.dat').transpose()
	dsw_iso = dsw[1][mini:maxi]-dts[2]
	dnw = np.loadtxt(ebossdir+'xi0gebossELG_NGC'+v+rec+'_mz0.6xz1.1'+wm+bs+'.dat').transpose()
	dnw_iso = dnw[1][mini:maxi]-dtn[2]
	ddt_iso = (dnw_iso/etn**2.+dsw_iso/ets**2.)/(1./etn**2.+1./ets**2.)
	#print(ddt_iso)
	plt.errorbar(dnw[0][mini:maxi],ddt_iso*1.e3,et[mini:maxi]*1.e3,fmt='ko')
	plt.plot(dts[0],dtt_iso*1.e3,'k-')
	#plt.plot(dt[0],dt[0]**2.*dtn[1],'k--')
	plt.xlim(20,190)
	#plt.ylim(-35,80)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$10^3$($\xi(s)-\xi_{no BAO}$)',size=16)
	#plt.text(30,71,r'$\alpha=1.044\pm0.042$',color='k',size=16)
	#plt.text(35,64,r'$\chi^2$/dof = 5.8/7',color='k',size=16)
	#plt.text(30,71,r'$\alpha=0.986\pm0.040$',color='k',size=16)
	covi = np.linalg.pinv(cov[mini:maxi,mini:maxi])
	diff = (ddt_iso-dtt_iso)
	chi2 = np.dot(diff,(np.dot(diff,covi)))
	print chi2
	plt.text(120,1.00,r'$\chi^2$/dof = '+str(round(chi2,1))+'/'+str(maxi-mini-5),color='k',size=16)
	#plt.title(r'BAO best-fit for v1.6 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True


def plotELGNSbaolike(v='3',p='3',Bp='0.4',rec='',bs=''):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiELGNSbaolik'+v+rec+'.pdf')
	plt.clf()
	plt.minorticks_on()
	
	db = np.loadtxt(ebossdir+'BAOxichilNScombELG'+v+rec+Bp+bs+'.dat').transpose()
	a = db[0]
	ax = sigreg_c12l(db[1])
	print ax
	alph = (ax[2]+ax[1])/2.
	err = (ax[2]-ax[1])/2.
	dnb = np.loadtxt(ebossdir+'BAOxichilNScombELG'+v+rec+'nobao'+Bp+bs+'.dat').transpose()[1]
	chim = min(db[1])
	ol = np.ones((len(db[0])))
	plt.plot(db[0],db[1]-chim,'-',color='k',linewidth=2)
	plt.plot(db[0],dnb-chim,'--',color='k',linewidth=1)
	plt.plot(db[0],ol,'k:',linewidth=1)
	plt.text(0.825,1.1,r'$1\sigma$')
	plt.plot(db[0],ol*4,'k:',linewidth=1)
	plt.text(0.825,4.1,r'$2\sigma$')
	plt.plot(db[0],ol*9,'k:',linewidth=1)
	plt.text(0.825,9.1,r'$3\sigma$')
	plt.ylim(0,16)
	plt.xlabel(r'$\alpha_{\rm BAO}$',size=18)
	plt.ylabel(r'$\Delta\chi^2$',size=18)
	plt.text(.9,10,r'$\alpha=$'+str(round(alph,3))+r'$\pm$'+str(round(err,3)))
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	#plt.title(r'BAO likelihood for '+title)
	if rec == '_rec':
		plt.title(r'ELGs, $\xi$, post-recon')
	if rec == '':
		plt.title(r'ELGs, $\xi$, pre-recon')
	pp.savefig()
	pp.close()
	return True

if __name__ == "__main__":
	#brickanalysis(ebossdir=ebossdir)
	#matchbrick(ebossdir='')
	nzra()