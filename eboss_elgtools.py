from math import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import fitsio

dirsci = '/mnt/lustre/ashleyr/eboss/' #where AJR puts eboss catalogs, change this to wherever you have put catalogs
dirsys = 'maps/' #change to local directory where ngalvsys from wiki was put, note star map and depth map included
dirfits = '/Users/ashleyross/fitsfiles/' #change to where your catalog files are
ebossdir = '/Users/ashleyross/eboss/' #where AJR puts correlation functions, writes out results
dirscio = '/mnt/lustre/ashleyr/eboss/mocks/'


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


def mkodens(chunk,zmin,zmax,v='v5_10_7',res=256,tol=.2,app='.fits',compl=.8,compls=.7):
	#dir = 'output_v4/'
	from healpix import radec2thphi
	import healpy as hp
	np = 12*res*res
	gl = []
	rl = []
	for i in range(0,np):
		gl.append(0)
		rl.append(0)
	f = fitsio.read(dirfits +'eboss'+chunk+'.'+v+'.latest.fits')
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z > zmin and z < zmax and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['isdupl'] == False: #and f[i]['Z_reliable']==True:		
			ra,dec=f[i]['ra'],f[i]['dec']
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi)
			gl[p] += 1.
	f = fitsio.read(dirfits +'eboss'+chunk+'.'+v+'.latest.rands.fits')
	for i in range(0,len(f)):
		if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:		
			ra,dec=f[i]['ra'],f[i]['dec']
			th,phi = radec2thphi(ra,dec)
			p = hp.ang2pix(res,th,phi)
			rl[p] += f[i]['sector_TSR']#*f[i]['sector_SSR']

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
	fo = open(ebossdir+'galelg'+chunk+v+str(zmin)+str(zmax)+str(res)+'odenspczw.dat','w')
	ft = open(ebossdir+'galelg'+chunk+v+str(zmin)+str(zmax)+str(res)+'rdodens.dat','w')
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

def plotELGcomp(chunk,comp='plate_SSR',cmap=cm.coolwarm,v='v5_10_7'):
	#remember!!! export PATH=$PATH:/Users/ashleyross/mangle2.2/bin
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	import fitsio
	import matplotlib.cm as cm
	
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	f = fitsio.read(dirfits +'eboss'+chunk+'.'+v+'.latest.rands.fits')
	map = plt.scatter(f['ra'],f['dec'],c=f[comp],s=.1,cmap=cmap,vmin=.5,lw=0)
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

def plotELGgalchunk(chunk,v='v5_10_7',zmin=.6,zmax=1.1,compl=.8,compls=.7):
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
	for i in range(0,len(f)):
		z = f[i]['Z']
		if z < zmax:
			if z > zmin and z < zmax and f[i]['Z_reliable']==True and f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls and f[i]['isdupl'] == False:
				ral.append(f[i]['ra'])
				decl.append(f[i]['dec'])

	f = fitsio.read(dirfits +'eboss'+chunk+'.'+v+'.latest.rands.fits')
	rarl = []
	decrl = []
	for i in range(0,len(f)):
			if f[i]['sector_TSR'] > compl and f[i]['sector_SSR'] > compls:
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

def plotxiELG(chunk='21p22',bs='8st0',v='v5_10_7',wm='fkp',zmin=.7,zmax=1.1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi0ELG'+chunk+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi0gebossELG_SGCtest'+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	plt.plot(dt[0],dt[0]**2.*dt[1]*.9,'k:')
	plt.plot(d1[0],d1[0]**2.*(d1[1]))
	plt.plot(d2[0],d2[0]**2.*(d2[1]))
	plt.xlim(10,200)
	plt.ylim(-30,60)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	plt.ylabel(r'$s^2\xi(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'Correlation function of ELGs in chunk, '+chunk+' '+str(zmin)+' < z < '+str(zmax))
	plt.legend(labels=['model','eBOSS'])

	pp.savefig()
	pp.close()
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

def plotxiELGrec(reg='SGC',bs='8st0',v='test',wm='fkp',zmin=.6,zmax=1.1):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi0ELGrec'+reg+v+wm+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	d1 = np.loadtxt(ebossdir+'xi0gebossELG_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi0gebossELG_'+reg+v+'rec_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	dt = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.593.04.07.015.00.dat').transpose()
	dtrec = np.loadtxt('BAOtemplates/xi0Challenge_matterpower0.593.02.04.015.01.0.dat').transpose()
	plt.plot(dt[0],dt[0]**2.*dt[1]*.9+.3*dt[0]-.001*dt[0]**2.,'k:')
	plt.plot(dt[0],dt[0]**2.*dtrec[1]*.9+.3*dt[0]-.001*dt[0]**2.,'k--')
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

