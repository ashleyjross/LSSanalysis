from math import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import rc
from matplotlib import rcParams

import fitsio
#from xitools_eboss import *
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)


dirsci = '/mnt/lustre/ashleyr/eboss/' #where AJR puts eboss catalogs, change this to wherever you have put catalogs
dirsys = 'maps/' #change to local directory where ngalvsys from wiki was put, note star map and depth map included
dirfits = '/Users/ashleyross/fitsfiles/' #change to where your catalog files are
ebossdir = '/Users/ashleyross/Dropbox/eboss/' #where AJR puts correlation functions, writes out results
dirscio = '/mnt/lustre/ashleyr/eboss/mocks/'

def txt2html(filen):
	f = open(filen+'.txt').readlines()
	fo = open(filen+'htmlformat.txt','w')
	for i in range(0,len(f)):
		ln = f[i].split(',')
		fo.write('<tr><td>'+ln[0]+'</td><td>'+ln[1]+'</td><td>'+ln[2]+'</td></tr>\n')
	fo.close()
	return True

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
		print(( b,ms))
	#print wm.split('_')[1]	
	iv = False
	if len(wm.split('_')) > 1:
		if wm.split('_')[1] == 'ivar':
			b,ms = findlinmb(samp,v,'',wm,zmin,zmax,dir='')
			iv = True
			print(( b,ms))
			
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

	print( n)	
	print(( mins,maxs,wm))	
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
		print(( len(f)+len(f2),ns))
		
	else:
		f = fitsio.read(dir+'eboss'+samp+'.'+v+'.latest.rands'+app)
		ns = len(f)/1000000
		print(( len(f),ns))

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

	print(( n,nw)) #just helps to know things worked properly)
	print(( n/10000.*ns)) #area in sq degrees)
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
	print(( ng,np,sum(gl),sum(rl)))
	ave = ng/np
	print(( ave,avet))
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
	print( no)
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
	print(( ng,np,sum(gl),sum(rl)))
	ave = ng/np
	print(( ave,avet))
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
	print( no)
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
	print(( 'effective area is '+str(area)))
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
			
	print (ndup)
	print((sumw/sumt))
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
	print((sum,zm/wtot,nw,nnw))
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
					
	print((ndup,sum,sumt,nlow,nlow/sumt))

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
					
	print((ndup,sum,sumt,nlow,nlow/sumt))

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
	print(('effective area is '+str(area)))
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
			
	print( ndup)
	print((sumw/sumt))
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
	print(( sum,zm/wtot,nw,nnw))
	plt.plot(zpl,nbl)
	plt.show()
	plt.plot(zpl,nbwl)
	plt.show()
	return True

def LRGzfail(reg,v='test'):
	d = np.loadtxt(ebossdir+'eBOSS_LRG_fibeff_'+reg+'_v'+v+'.dat').transpose()
	gl = np.zeros(500)
	tl = np.zeros(500)
	for i in range(0,len(d[0])):
		fib = d[0][i]
		ind = int(fib/2)-1
		tl[ind] += d[2][i]
		gl[ind] += d[1][i]
	effl = gl/tl
	fib = np.arange(1.5,1000,2)	
	plt.plot(fib,effl,'ko')
	plt.show()	

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
			
	print(( ngal,nq,nstar,ngoodgal,nqz,ntot,ngal+nq+nstar))
	print(( ngoodgal/ngal,nSSR/ntot,ngoodtot/ntot,ngoodgal/ntot,ngoodgal/(ngal+nq),(ngoodgal+nqz)/(ntot),nqz/nq,nq/ntot,ngz/ngoodtot,nqzz/nqz,nqfz/nqf,nbgz/nbg))
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
			
	print(( ngal,nq,nstar,ngoodgal,nqz,ntot,ngal+nq+nstar))
	print(( ngoodgal/ngal,nSSR/ngoodgal,ngoodgal/ntot,ngoodgal/(ngal+nq),(ngoodgal+nqz)/(ngal+nq-nqz)))
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
	print(( 'effective areas are '+str(area1) +' and '+str(area2)))
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
			
	print(( sumw/sumt))
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
	print(( 'effective areas are '+str(area1) +' and '+str(area2)))
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
			
	print(( sumw/sumt))
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
			
	print(( sumw/sumt))
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
			
	print(( sumw/sumt))
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
			
	print(( sumw/sumt))
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
	print( fac)
	plt.plot(zpl,zl1/sum(zl1),zpl,zl2/sum(zl2),zpl,diffl/sum(diffl))
	plt.xlabel('redshift')
	plt.ylabel('Normalized redshift histogram')
	plt.legend(labels=['SNR > '+str(splitv+plus),'SNR < '+str(splitv-plus)])
	plt.show()
	return True



def calcxi_mockEZ(num,reg='SGC',samp='ELG',ver='7',pcmass=False,bs=8,mom=0,mumin=0,mumax=1,start=0,rec='',mupow=0,zmin=0.6,sys=''):
	#dir = '/mnt/lustre/ashleyr/eboss/EZmock'+samp+'v'+str(ver)+'/' #sciama
	dir = '/global/cscratch1/sd/ajross/ebossxi/'
	dirout = dir+ samp+'_v'+str(ver)+'/'
	muw = ''
	zw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	if mupow != 0:
		muw += 'mup'+str(mupow)		
	if samp == 'ELG':
	
		if zmin != 0.6:
			zw += 'zmin'+str(zmin)
		zmax = 1.1
	if samp == 'LRG' or samp == 'LRGpCMASS':
		zmin = 0.6
		zmax = 1.0
	if samp == 'QSO':
		zmin = 0.8
		zmax = 2.2
	predir = '/global/homes/z/zhaoc/cscratch/EZmock/2PCF/'+samp+'v'+ver+sys+'_rmu/z'+str(zmin)+'z'+str(zmax)+'/2PCF/'
	#predir = dir + samp+'_v'+str(ver)+'/2PCF/'
	recdir = '/global/homes/z/zhaoc/cscratch/EZmock/2PCF/'+samp+'v'+ver+rec+sys+'_rmu/z'+str(zmin)+'z'+str(zmax)+'/2PCF/'
	#shuffdir = dir + 'prerec_shuf/2PCF/'

	#if pcmass:
	#	samp += 'pCMASS'	
	#ddnorm = af['ddnorm'][0][num]/2.
	#drnorm = af['drnorm'][0][num]/2.
	#ddnorm = 1.
	#drnorm = .1
	#rrnorm = .01/4.
	#print(ddnorm,drnorm,ddnorm/drnorm)
	print((samp,zmin,zmax))
	zer = ''
	if num < 1000:
		zer += '0'
	if num < 100:
		zer += 	'0'
	if num < 10:
		zer += '0'
	fn = '2PCF_EZmock_eBOSS_'+samp+'_'+reg+'_v'+str(ver)+'_z'+str(zmin)+'z'+str(zmax)+'_'+zer+str(num)
	#if rec == '':
	#	af = fitsio.read(dir+'2PCF_ELGv4_'+reg+'_merge.fits')
	#	dd = af['dd'].transpose()[num]*ddnorm
	#	dr = af['dr'].transpose()[num]*drnorm
	#*rrnorm
	#normr = (np.loadtxt(predir+'2PCF_EZmock_eBOSS_'+samp+'_'+reg+'_v'+str(ver)+'_z'+str(zmin)+'z'+str(zmax)+'.rr').transpose()[-2]/rr)[0]
	if rec == '':
		
		dd = np.loadtxt(predir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(predir+fn+'.dr').transpose()[-1]#*drnorm
		rr = np.loadtxt(predir+fn+'.rr').transpose()[-1]
		#normd = (np.loadtxt(predir+fn+'.dd').transpose()[-2]/dd)[1000]
		#normdr = (np.loadtxt(predir+fn+'.dr').transpose()[-2]/dr)[0]
		#print(normdr/normd,normr/normdr)

	if rec == '_rec' or rec == '_rec':
		
		#fn += '_rec'
		dd = np.loadtxt(recdir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(recdir+fn+'.ds').transpose()[-1]#*drnorm
		ss = np.loadtxt(recdir+fn+'.ss').transpose()[-1]#*rrnorm	
		rr = np.loadtxt(predir+fn+'.rr').transpose()[-1]	

	if rec == 'shuff':		
		fn += '_rec'
		dd = np.loadtxt(shuffdir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(shuffdir+fn+'.ds').transpose()[-1]#*drnorm
		ss = np.loadtxt(shuffdir+fn+'.ss').transpose()[-1]#*rrnorm
		#normd = (np.loadtxt(shuffdir+fn+'.dd').transpose()[-2]/dd)[1000]
		#normdr = (np.loadtxt(shuffdir+fn+'.ds').transpose()[-2]/dr)[0]
		#norms = (np.loadtxt(shuffdir+fn+'.ss').transpose()[-2]/ss)[0]
		#print(normdr/normd,norms/normdr,norms/normr)
	
	#if subt:
	#if rec == 'rec':
	#	wrp = np.loadtxt(dir+'xi/wrpELG'+reg+'EZmock'+str(num-1)+muw+'1st0.dat').transpose()
	#else:
	#if rec != 'shuff' and samp != 'LRGpCMASS':
	#	wrp = np.loadtxt(dir+'xi/wrpELG'+reg+'EZmock'+rec+str(num)+muw+'1st0.dat').transpose()
	
	
	nb = (200-start)//bs
	xil = np.zeros(nb)
	xil2 = np.zeros(nb)
	xil4 = np.zeros(nb)
	wl = np.zeros(nb)
	wl2 = np.zeros(nb)
	wl4 = np.zeros(nb)

	nmub = 120
	dmu = 1./float(nmub)
	mubm = 0
	if mumin != 0:
		mubm = int(mumin*nmub)
	mubx = nmub
	if mumax != 1:
		mubx = int(mumax*nmub)
		print((mubx,mubx/nmub))
	for i in range(start,nb*bs+start,bs):
		xib = 0
		xib2 = 0
		xib4 = 0
		ddt = 0
		drt = 0
		rrt = 0
		sst = 0
		w = 0
		w2 = 0
		w4 = 0
		mut = 0
		rmin = i
		rmax = rmin+bs
		for m in range(mubm,mubx):
			ddb = 0
			drb = 0
			rrb = 0
			ssb = 0
			mu = m/float(nmub) + 0.5/float(nmub)
			for b in range(0,bs):
				bin = nmub*(i+b)+m
				if bin < 24000:
					ddb += dd[bin]
					drb += dr[bin]
					rrb += rr[bin]
					if rec == '_rec' or rec == 'shuff':
						ssb += ss[bin]
						sst += ss[bin]
				ddt += dd[bin]
				drt += dr[bin]
				rrt += rr[bin]
			if rec == '_rec' or rec == 'shuff':
				xi = (ddb-2.*drb+ssb)/rrb
			else:		
				xi = (ddb-2.*drb+rrb)/rrb
			#if subt:
# 			if rec != 'shuff' and samp != 'LRGpCMASS':
# 				mum = mu #m/float(nmub)
# 				mux = mu #mum+dmu
# 				rpmin = rmin*sqrt(1.-mux**2.)
# 				rpmax = rmax*sqrt(1.-mum**2.)
# 				xip = 0
# 				wr = ((wrp[0] >rpmin) & (wrp[0] < rpmax))
# 				if len(wrp[1][wr] > 0):
# 					xip = np.mean(wrp[1][wr])
# 					mut += dmu
# 				else:
# 					#print(rpmin,rpmax)	
# 					xip = 0
# 					
# 				w += xip*dmu
# 				w2 += xip*dmu*P2(mu)*5.
# 				w4 += xip*dmu*P4(mu)*9.

			xib += xi*dmu*(mu**mupow)
			xib2 += xi*dmu*P2(mu)*5.
			xib4 += xi*dmu*P4(mu)*9.		
		xil[i//bs] = xib
		xil2[i//bs] = xib2
		xil4[i//bs] = xib4
		if rec != 'shuff':
			wl[i//bs] = w
			wl2[i//bs] = w2
			wl4[i//bs] = w4
		#if rec != '':	
		#	print(ddt/sst,drt/sst,sst/rrt,rrt)
		#else:
		#	print(ddt/rrt,drt/rrt,rrt)
		#print(i,w,w2,w4,mut)	
	#print xil
	#if subt:
	#	muw += 'subt'	
	fo = open(dirout+'xi024'+samp+'_v'+ver+reg+zw+'EZmock'+rec+sys+str(num)+muw+str(bs)+'st'+str(start)+'.dat','w')
	for i in range(0,len(xil)):
		r = bs/2.+i*bs+start
		fo.write(str(r)+' '+str(xil[i])+' '+str(xil2[i])+' '+str(xil4[i])+' '+str(wl[i])+' '+str(wl2[i])+' '+str(wl4[i])+'\n')
	fo.close()
	return True

def putallOR(reg,rec,bs=5,samp='ELG',start=0):
	dir = '/global/cscratch1/sd/ajross/ebossxi/ELG_v7/'
	losl = ['x','y','z']
	indl = []
	for i in range(11,65):
		for los in losl:
			try:
				calcxi_mockOR(i,los=los,bs=bs,rec=rec,reg=reg)
				indl.append(i)
				print((i,los))
			except:
				pass
	print((len(indl)))
	avel = np.loadtxt(dir+'xi024ELG_v7'+reg+'ORmock'+rec+losl[0]+str(indl[0])+str(bs)+'st'+str(start)+'.dat')
	avel += np.loadtxt(dir+'xi024ELG_v7'+reg+'ORmock'+rec+losl[1]+str(indl[0])+str(bs)+'st'+str(start)+'.dat')
	avel += np.loadtxt(dir+'xi024ELG_v7'+reg+'ORmock'+rec+losl[2]+str(indl[0])+str(bs)+'st'+str(start)+'.dat')
	for i in range(1,len(indl)):
		for los in losl:
			avel += np.loadtxt(dir+'xi024ELG_v7'+reg+'ORmock'+rec+los+str(indl[i])+str(bs)+'st'+str(start)+'.dat')
	avel = avel/(3.*len(indl))
	np.savetxt('xiave024ELGORmock'+reg+rec+str(bs)+'st'+str(start)+'.dat',avel)
	#fo = open('xiave024ELGORmock'+reg+rec+str(bs+'st'+str(start)+'.dat','w')	
	return True	

def calcxi_mockOR(num,los='x',reg='SGC',samp='ELG',ver='7',pcmass=False,bs=8,mom=0,mumin=0,mumax=1,start=0,rec='',mupow=0):
	#dir = '/mnt/lustre/ashleyr/eboss/EZmock'+samp+'v'+str(ver)+'/' #sciama
	dir = '/global/cscratch1/sd/ajross/ebossxi/'
	dirout = dir+ samp+'_v'+str(ver)+'/'
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	if mupow != 0:
		muw += 'mup'+str(mupow)		
	if samp == 'ELG':
		zmin = 0.6
		zmax = 1.1
	if samp == 'LRG' or samp == 'LRGpCMASS':
		zmin = 0.6
		zmax = 1.0
	if samp == 'QSO':
		zmin = 0.8
		zmax = 2.2
	predir = '/global/homes/z/zhaoc/cscratch/EZmock/2PCF/OuterRim_cutsky/2PCF/OuterRimCutsky_mocks_vAvila_4/glos-'+los+'/'
	recdir = '/global/homes/z/zhaoc/cscratch/EZmock/2PCF/OuterRim_cutsky/2PCF/OuterRimCutsky_mocks_vAvila_4_rec/glos-'+los+'/'
	print((samp,zmin,zmax))
	fn = 'OuterRim_ELG_clustering_'+reg+'_vAvila_4'+rec+'_b132_n2410_f0.14_vb0.00.'+str(num)
	fnp = 'OuterRim_ELG_clustering_'+reg+'_vAvila_4'+'_b132_n2410_f0.14_vb0.00.'+str(num)
	if rec == '':
		
		dd = np.loadtxt(predir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(predir+fn+'.dr').transpose()[-1]#*drnorm
		rr = np.loadtxt(predir+fn+'.rr').transpose()[-1]

	if rec == '_rec':
		
		dd = np.loadtxt(recdir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(recdir+fn+'.ds').transpose()[-1]#*drnorm
		ss = np.loadtxt(recdir+fn+'.ss').transpose()[-1]#*rrnorm	
		rr = np.loadtxt(predir+fnp+'.rr').transpose()[-1]	

	
	
	nb = (200-start)//bs
	xil = np.zeros(nb)
	xil2 = np.zeros(nb)
	xil4 = np.zeros(nb)

	nmub = 120
	dmu = 1./float(nmub)
	mubm = 0
	if mumin != 0:
		mubm = int(mumin*nmub)
	mubx = nmub
	if mumax != 1:
		mubx = int(mumax*nmub)
		print((mubx,mubx/nmub))
	for i in range(start,nb*bs+start,bs):
		xib = 0
		xib2 = 0
		xib4 = 0
		ddt = 0
		drt = 0
		rrt = 0
		sst = 0
		mut = 0
		rmin = i
		rmax = rmin+bs
		for m in range(mubm,mubx):
			ddb = 0
			drb = 0
			rrb = 0
			ssb = 0
			mu = m/float(nmub) + 0.5/float(nmub)
			for b in range(0,bs):
				bin = nmub*(i+b)+m
				if bin < 24000:
					ddb += dd[bin]
					drb += dr[bin]
					rrb += rr[bin]
					if rec == '_rec' or rec == 'shuff':
						ssb += ss[bin]
						sst += ss[bin]
				ddt += dd[bin]
				drt += dr[bin]
				rrt += rr[bin]
			if rec == '_rec' or rec == 'shuff':
				xi = (ddb-2.*drb+ssb)/rrb
			else:		
				xi = (ddb-2.*drb+rrb)/rrb

			xib += xi*dmu*(mu**mupow)
			xib2 += xi*dmu*P2(mu)*5.
			xib4 += xi*dmu*P4(mu)*9.		
		xil[i//bs] = xib
		xil2[i//bs] = xib2
		xil4[i//bs] = xib4

	fo = open(dirout+'xi024'+samp+'_v'+ver+reg+'ORmock'+rec+los+str(num)+muw+str(bs)+'st'+str(start)+'.dat','w')
	for i in range(0,len(xil)):
		r = bs/2.+i*bs+start
		fo.write(str(r)+' '+str(xil[i])+' '+str(xil2[i])+' '+str(xil4[i])+'\n')
	fo.close()
	return True


def calcxi_dataCZ(reg='SGC',samp='ELG',ver='7',pcmass=False,bs=8,mom=0,mupow=0,mumin=0,mumax=1,start=0,rec='',zmin=0.6):
	#dir = '/mnt/lustre/ashleyr/eboss/EZmock'+samp+'v'+str(ver)+'/' #sciama
	dir = '/global/cscratch1/sd/ajross/ebossxi/'
	dirout = ''#dir+ samp+'_v'+str(ver)+'/'
	muw = ''
	zw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	if mupow != 0:
		muw += 'mup'+str(mupow)		
	if samp == 'ELG':
		#zmin = 0.6
		if zmin != 0.6:
			zw += 'zmin'+str(zmin)
		zmax = 1.1
	if samp == 'LRG' or samp == 'LRGpCMASS':
		zmin = 0.6
		zmax = 1.0
	if samp == 'QSO':
		zmin = 0.8
		zmax = 2.2
	if samp == 'ELG':	
		indir = '/global/homes/z/zhaoc/cscratch/EZmock/ELG/data/clustering/z'+str(zmin)+'z1.1/'	

	#if pcmass:
	#	samp += 'pCMASS'	
	#ddnorm = af['ddnorm'][0][num]/2.
	#drnorm = af['drnorm'][0][num]/2.
	#ddnorm = 1.
	#drnorm = .1
	#rrnorm = .01/4.
	#print(ddnorm,drnorm,ddnorm/drnorm)
	print((samp,zmin,zmax))
	fn = '2PCF_eBOSS_'+samp+'_clustering_'+reg+'_v'+ver+rec+'_z'+str(zmin)+'z'+str(zmax)
	fnnorec = '2PCF_eBOSS_'+samp+'_clustering_'+reg+'_v'+ver+'_z'+str(zmin)+'z'+str(zmax)
	#if rec == '':
	#	af = fitsio.read(dir+'2PCF_ELGv4_'+reg+'_merge.fits')
	#	dd = af['dd'].transpose()[num]*ddnorm
	#	dr = af['dr'].transpose()[num]*drnorm
	#*rrnorm
	#normr = (np.loadtxt(predir+'2PCF_EZmock_eBOSS_'+samp+'_'+reg+'_v'+str(ver)+'_z'+str(zmin)+'z'+str(zmax)+'.rr').transpose()[-2]/rr)[0]
	if rec == '':
		
		dd = np.loadtxt(indir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(indir+fn+'.dr').transpose()[-1]#*drnorm
		rr = np.loadtxt(indir+fn+'.rr').transpose()[-1]
		#normd = (np.loadtxt(predir+fn+'.dd').transpose()[-2]/dd)[1000]
		#normdr = (np.loadtxt(predir+fn+'.dr').transpose()[-2]/dr)[0]
		#print(normdr/normd,normr/normdr)

	if rec == '_rec':
		
		#fn += '_rec'
		dd = np.loadtxt(indir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(indir+fn+'.ds').transpose()[-1]#*drnorm
		ss = np.loadtxt(indir+fn+'.ss').transpose()[-1]#*rrnorm	
		rr = np.loadtxt(indir+fnnorec+'.rr').transpose()[-1]	

	if rec == 'shuff':		
		fn += '_rec'
		dd = np.loadtxt(shuffdir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(shuffdir+fn+'.ds').transpose()[-1]#*drnorm
		ss = np.loadtxt(shuffdir+fn+'.ss').transpose()[-1]#*rrnorm
		#normd = (np.loadtxt(shuffdir+fn+'.dd').transpose()[-2]/dd)[1000]
		#normdr = (np.loadtxt(shuffdir+fn+'.ds').transpose()[-2]/dr)[0]
		#norms = (np.loadtxt(shuffdir+fn+'.ss').transpose()[-2]/ss)[0]
		#print(normdr/normd,norms/normdr,norms/normr)
	
	#if subt:
	#if rec == 'rec':
	#	wrp = np.loadtxt(dir+'xi/wrpELG'+reg+'EZmock'+str(num-1)+muw+'1st0.dat').transpose()
	#else:
	#if rec != 'shuff' and samp != 'LRGpCMASS':
	#	wrp = np.loadtxt(dir+'xi/wrpELG'+reg+'EZmock'+rec+str(num)+muw+'1st0.dat').transpose()
	
	
	nb = (200-start)//bs
	xil = np.zeros(nb)
	xil2 = np.zeros(nb)
	xil4 = np.zeros(nb)
	wl = np.zeros(nb)
	wl2 = np.zeros(nb)
	wl4 = np.zeros(nb)

	nmub = 120
	dmu = 1./float(nmub)
	mubm = 0
	if mumin != 0:
		mubm = int(mumin*nmub)
	mubx = nmub
	if mumax != 1:
		mubx = int(mumax*nmub)
	for i in range(start,nb*bs+start,bs):
		xib = 0
		xib2 = 0
		xib4 = 0
		ddt = 0
		drt = 0
		rrt = 0
		sst = 0
		w = 0
		w2 = 0
		w4 = 0
		mut = 0
		rmin = i
		rmax = rmin+bs
		for m in range(mubm,mubx):
			ddb = 0
			drb = 0
			rrb = 0
			ssb = 0
			mu = m/float(nmub) + 0.5/float(nmub)
			for b in range(0,bs):
				bin = nmub*(i+b)+m
				if bin < 24000:
					ddb += dd[bin]
					drb += dr[bin]
					rrb += rr[bin]
					if rec == '_rec' or rec == 'shuff':
						ssb += ss[bin]
						sst += ss[bin]
				ddt += dd[bin]
				drt += dr[bin]
				rrt += rr[bin]
			if rec == '_rec' or rec == 'shuff':
				xi = (ddb-2.*drb+ssb)/rrb
			else:		
				xi = (ddb-2.*drb+rrb)/rrb
			#if subt:
# 			if rec != 'shuff' and samp != 'LRGpCMASS':
# 				mum = mu #m/float(nmub)
# 				mux = mu #mum+dmu
# 				rpmin = rmin*sqrt(1.-mux**2.)
# 				rpmax = rmax*sqrt(1.-mum**2.)
# 				xip = 0
# 				wr = ((wrp[0] >rpmin) & (wrp[0] < rpmax))
# 				if len(wrp[1][wr] > 0):
# 					xip = np.mean(wrp[1][wr])
# 					mut += dmu
# 				else:
# 					#print(rpmin,rpmax)	
# 					xip = 0
# 					
# 				w += xip*dmu
# 				w2 += xip*dmu*P2(mu)*5.
# 				w4 += xip*dmu*P4(mu)*9.

			xib += xi*dmu*(mu**mupow)
			xib2 += xi*dmu*P2(mu)*5.
			xib4 += xi*dmu*P4(mu)*9.		
		xil[i//bs] = xib
		xil2[i//bs] = xib2
		xil4[i//bs] = xib4
		if rec != 'shuff':
			wl[i//bs] = w
			wl2[i//bs] = w2
			wl4[i//bs] = w4
		#if rec != '':	
		#	print(ddt/sst,drt/sst,sst/rrt,rrt)
		#else:
		#	print(ddt/rrt,drt/rrt,rrt)
		#print(i,w,w2,w4,mut)	
	#print xil
	#if subt:
	#	muw += 'subt'	
	fo = open(dirout+'xi024'+samp+'_v'+ver+reg+zw+'CZdata'+rec+muw+str(bs)+'st'+str(start)+'.dat','w')
	for i in range(0,len(xil)):
		r = bs/2.+i*bs+start
		fo.write(str(r)+' '+str(xil[i])+' '+str(xil2[i])+' '+str(xil4[i])+' '+str(wl[i])+' '+str(wl2[i])+' '+str(wl4[i])+'\n')
	fo.close()
	return True


def calcxi_mockELGEZcompshuff(num,reg='SGC',bs=8,mom=0,mumin=0,mumax=1,start=0):
	dir = '/mnt/lustre/ashleyr/eboss/EZmockELGv4/'
	predir = dir + 'prerec/2PCF/'
	recdir = dir + 'recon/2PCF/'
	shuffdir = dir + 'prerec_shuf/2PCF/'
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	
	zer = ''
	if num < 1000:
		zer += '0'
	if num < 100:
		zer += 	'0'
	if num < 10:
		zer += '0'
	fn = '2PCF_EZmock_eBOSS_ELG_'+reg+'_v4_z0.6z1.1_'+zer+str(num)
	#if rec == '':
	#	af = fitsio.read(dir+'2PCF_ELGv4_'+reg+'_merge.fits')
	#	dd = af['dd'].transpose()[num]*ddnorm
	#	dr = af['dr'].transpose()[num]*drnorm
	rr = np.loadtxt(predir+'2PCF_EZmock_eBOSS_ELG_'+reg+'_v4_z0.6z1.1.rr').transpose()[-1]#*rrnorm
	normr = (np.loadtxt(predir+'2PCF_EZmock_eBOSS_ELG_'+reg+'_v4_z0.6z1.1.rr').transpose()[-2]/rr)[0]
		
	dd = np.loadtxt(predir+fn+'.dd').transpose()[-1]#*ddnorm
	dr = np.loadtxt(predir+fn+'.dr').transpose()[-1]#*drnorm
	#normd = (np.loadtxt(predir+fn+'.dd').transpose()[-2]/dd)[1000]
	#normdr = (np.loadtxt(predir+fn+'.dr').transpose()[-2]/dr)[0]
	#print(normdr/normd,normr/normdr)


	fn += '_rec'
	dds = np.loadtxt(shuffdir+fn+'.dd').transpose()[-1]#*ddnorm
	ds = np.loadtxt(shuffdir+fn+'.ds').transpose()[-1]#*drnorm
	ss = np.loadtxt(shuffdir+fn+'.ss').transpose()[-1]#*rrnorm
	#normd = (np.loadtxt(shuffdir+fn+'.dd').transpose()[-2]/dd)[1000]
	#normdr = (np.loadtxt(shuffdir+fn+'.ds').transpose()[-2]/dr)[0]
	#norms = (np.loadtxt(shuffdir+fn+'.ss').transpose()[-2]/ss)[0]
	#	print(normdr/normd,norms/normdr,norms/normr)
	
	#if subt:
	#if rec == 'rec':
	#	wrp = np.loadtxt(dir+'xi/wrpELG'+reg+'EZmock'+str(num-1)+muw+'1st0.dat').transpose()
	#else:
	
	
	nb = (200-start)/bs

	nmub = 120
	dmu = 1./float(nmub)
	mubm = 0
	if mumin != 0:
		mubm = mumin*nmub
	mubx = nmub
	if mumax != 1:
		mubx = mumax*nmub
	for i in range(start,nb*bs+start,bs):
		xib = 0
		xib2 = 0
		xib4 = 0
		ddt = 0
		drt = 0
		rrt = 0
		ddst = 0
		dst = 0
		sst = 0
		for m in range(mubm,mubx):
			ddb = 0
			ddsb = 0
			drb = 0
			dsb = 0
			rrb = 0
			ssb = 0
			mu = m/float(nmub) + 0.5/float(nmub)
			for b in range(0,bs):
				bin = nmub*(i+b)+m
				if bin < 24000:
					sst += ss[bin]
					ddt += dd[bin]
					ddst += dds[bin]
					drt += dr[bin]
					dst += ds[bin]
					rrt += rr[bin]
		print((ddt/rrt,ddst/rrt,drt/rrt,dst/rrt,sst/rrt,rrt))
	return True


def calcwrp_mockELGEZ(num,reg='SGC',bs=1,mom=0,mumin=0,mumax=1,start=0,rec=''):
	dir = '/mnt/lustre/ashleyr/eboss/EZmockELGv4/'
	predir = dir + 'prerec/2PCF/'
	recdir = dir + 'recon/2PCF/'
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	
	#ddnorm = af['ddnorm'][0][num]/2.
	#drnorm = af['drnorm'][0][num]/2.
	#ddnorm = 1.
	#drnorm = .1
	#rrnorm = .01/4.
	#print(ddnorm,drnorm,ddnorm/drnorm)
	zer = ''
	if num < 1000:
		zer += '0'
	if num < 100:
		zer += 	'0'
	if num < 10:
		zer += '0'
	fn = '2PCF_EZmock_eBOSS_ELG_'+reg+'_v4_z0.6z1.1_'+zer+str(num)
	if rec == '':
		
		dd = np.loadtxt(predir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(predir+fn+'.dr').transpose()[-1]#*drnorm

	if rec == 'rec':
		
		fn += '_rec'
		dd = np.loadtxt(recdir+fn+'.dd').transpose()[-1]#*ddnorm
		dr = np.loadtxt(recdir+fn+'.ds').transpose()[-1]#*drnorm
		ss = np.loadtxt(recdir+fn+'.ss').transpose()[-1]#*rrnorm		
	
	rr = np.loadtxt(dir+'2PCF_EZmock_eBOSS_ELG_'+reg+'_v4_z0.6z1.1.rr').transpose()[-1]#*rrnorm
	nb = 200/bs
	ddl = np.zeros(nb)
	drl = np.zeros(nb)
	rrl = np.zeros(nb)
	ssl = np.zeros(nb)
	nmub = 120
	dmu = 1./float(nmub)
	mubm = 0
	if mumin != 0:
		mubm = mumin*nmub
	mubx = nmub
	if mumax != 1:
		mubx = mumax*nmub
	for i in range(0,200):
		xib = 0
		ddt = 0
		drt = 0
		rrt = 0
		for m in range(mubm,mubx):
			bin = nmub*i+m
			ddb = dd[bin]
			drb = dr[bin]
			rrb = rr[bin]
				#ddt += dd[bin]
				#drt += dr[bin]
				#rrt += rr[bin]
			rp = (i+.5)*sqrt(1.-(m/float(nmub)+.5/float(nmub))**2.)
			brp = int(rp)
			ddl[brp] += ddb
			drl[brp] += drb
			rrl[brp] += rrb 
			if rec == 'rec':
				ssl[brp] += ss[bin]		
	#print xil
	#print(ddl,drl,rrl)
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	fo = open(dir+'xi/wrpELG'+reg+'EZmock'+rec+str(num)+muw+str(bs)+'st'+str(start)+'.dat','w')
	for i in range(0,nb):
		r = bs/2.+i*bs
		if rec == '':
			xi = (ddl[i]-2.*drl[i]+rrl[i])/rrl[i]
		if rec == 'rec':
			xi = (ddl[i]-2.*drl[i]+rrl[i])/ssl[i]	
		fo.write(str(r)+' '+str(xi)+'\n')
	fo.close()
	return True


def mkcov_mockELG_EZ(reg,samp='ELG',pcmass=False,bs=8,mom=0,N=1000,start=0,mumin=0,mumax=1,angfac=0,rec='',zmin=0.6):
	dir = '/mnt/lustre/ashleyr/eboss/EZmock'+samp+'v4/xi/'
	if pcmass:
		samp += 'pCMASS'
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	zw = ''
	if samp == 'ELG' and zmin != 0.6:
		zw += 'zmin'+str(zmin)	
	bsst = str(bs)+'st'+str(start)
	nbin = (200-start)/bs
	xiave = np.zeros((nbin))
	cov = np.zeros((nbin,nbin))

	Ntot = 0
	fac = 1.
	for i in range(1,N+1):
		nr = str(i)
		#print( i)
		#try:
		if mom != 'rp':
			xii = np.loadtxt(dir+'xi024'+samp+reg+zw+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
			xiave += xii[mom/2+1]-angfac*xii[mom/2+4]
		else:
			xii = np.loadtxt(dir+'wrpELG'+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
			xiave += xii[1]	
# 		except:
# 			if mom != 'rp':
# 				calcxi_mockELGEZ(i+1,reg=reg,bs=bs,mom=mom,mumin=mumin,mumax=mumax,start=start,rec=rec)
# 				xii = np.loadtxt(dir+'xi024'+samp+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
# 				xiave += xii[mom/2+1]-angfac*xii[mom/2+4]
# 			print(i,mom)
		Ntot += 1.
		#except:
		#	print i
		#print( i)
	print( Ntot)		
	xiave = xiave/float(Ntot)
	for i in range(1,N+1):
		nr = str(i)
		if mom != 'rp':
			xii = np.loadtxt(dir+'xi024'+samp+reg+zw+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()[mom/2+1]
			xiit = np.loadtxt(dir+'xi024'+samp+reg+zw+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()[mom/2+4]
		else:
			xii = np.loadtxt(dir+'wrpELG'+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()[1]
			xiit = xii
		for j in range(0,nbin):
			xij = xii[j]-angfac*xiit[j]
			for k in range(0,nbin):
				xik = xii[k]-angfac*xiit[k]
				cov[j][k] += (xij-xiave[j])*(xik-xiave[k])
#		except:
#			print i
	cov = cov/float(Ntot)					
	fo = open('xiave'+str(mom)+reg+zw+samp+'_EZ'+rec+muw+'angfac'+str(angfac)+bsst+'.dat','w')
	errl = []
	for i in range(0,nbin):
		fo.write(str(bs/2.+bs*i+start)+ ' '+str(xiave[i])+ ' '+str(sqrt(cov[i][i]))+'\n')
		errl.append(sqrt(cov[i][i]))
	fo.close()	
	fo = open('cov'+str(mom)+reg+zw+samp+'_EZ'+rec+muw+'angfac'+str(angfac)+bsst+'.dat','w')
	
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

def mkcov_mockELG_EZ_mu(reg,mu,samp='ELG',pcmass=False,bs=8,N=1000,start=0,mumin=0,mumax=1,rec=''):
	'''
	makes the covariance for a particular mu, reconstructed from xi0,xi2,xi4
	'''
	dir = '/mnt/lustre/ashleyr/eboss/EZmock'+samp+'v4/xi/'
	if pcmass:
		samp += 'pCMASS'
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	bsst = str(bs)+'st'+str(start)
	nbin = (200-start)/bs
	xiave = np.zeros((nbin))
	cov = np.zeros((nbin,nbin))

	Ntot = 0
	fac = 1.
	for i in range(1,N+1):
		nr = str(i)
		#print( i)
		#try:
		xii = np.loadtxt(dir+'xi024'+samp+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
		xiave += xii[1]+P2(mu)*xii[2]+P4(mu)*xii[3]
		Ntot += 1.
		#except:
		#	print i
		#print( i)
	print( Ntot)		
	xiave = xiave/float(Ntot)
	for i in range(1,N+1):
		nr = str(i)
		xiif = np.loadtxt(dir+'xi024'+samp+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
		xii = xiif[1]+P2(mu)*xiif[2]+P4(mu)*xiif[3]
		for j in range(0,nbin):
			xij = xii[j]
			for k in range(0,nbin):
				xik = xii[k]
				cov[j][k] += (xij-xiave[j])*(xik-xiave[k])
#		except:
#			print i
	cov = cov/float(Ntot)					
	fo = open('xiavemu'+str(mu)+reg+samp+'_EZ'+rec+muw+bsst+'.dat','w')
	errl = []
	for i in range(0,nbin):
		fo.write(str(bs/2.+bs*i+start)+ ' '+str(xiave[i])+ ' '+str(sqrt(cov[i][i]))+'\n')
		errl.append(sqrt(cov[i][i]))
	fo.close()	
	fo = open('covmu'+str(mu)+reg+samp+'_EZ'+rec+muw+bsst+'.dat','w')
	
	for i in range(0,nbin):
		for j in range(0,nbin):
			fo.write(str(cov[i][j])+' ')
		fo.write('\n')
	fo.close()


def mkcov_mock02_EZ(reg,ver=5,samp='ELG',pcmass=False,bs=8,mom=0,N=1000,rmax=200,start=0,mumin=0,mumax=1,angfac=0,md='me',rec=''):
	#dir = '/mnt/lustre/ashleyr/eboss/EZmock'+samp+'v'+str(ver)+'/'
	#if md == 'me':
	#	dir += '/xi/'
	#if md == 'cz' and rec != 'rec':
	#	dir += 'prerec/2PCF/'	
	#if pcmass:
	#	samp += 'pCMASS'
	
	dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ver)+'/'
	
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	bsst = str(bs)+'st'+str(start)
	nbin = (rmax-start)//bs
	xiave = np.zeros((nbin*2))
	cov = np.zeros((nbin*2,nbin*2))

	Ntot = 0
	fac = 1.
	if samp == 'QSO':
		zmin = 0.8
		zmax = 2.2
	for i in range(1,N+1):
		nr = str(i)
		#print( i)
		#try:
		if md == 'me':
			xii = np.loadtxt(dir+'xi024'+samp+'_v'+str(ver)+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
		if md == 'cz':
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dir+'2PCF_EZmock_eBOSS_'+samp+'_'+reg+'_v'+str(ver)+'_z'+str(zmin)+'z'+str(zmax)+'_'+mockn+'.dat').transpose()
		if md == 'jh':
			dirx = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys/'
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dirx+'ez_qsov7_0.31_'+str.lower(reg)+'_xi4m_wsys_ds5_'+mockn+'.dat').transpose()
		
		xic = np.concatenate((xii[1],xii[2]),axis=None)#-angfac*xii[mom/2+4] need to add this if desired
		if len(xic) != nbin*2:
			print('ERROR, wrong array size')
			return('ERROR!')
		xiave += xic
# 		except:
# 			if mom != 'rp':
# 				calcxi_mockELGEZ(i+1,reg=reg,bs=bs,mom=mom,mumin=mumin,mumax=mumax,start=start,rec=rec)
# 				xii = np.loadtxt(dir+'xi024'+samp+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
# 				xiave += xii[mom/2+1]-angfac*xii[mom/2+4]
# 			print(i,mom)
		Ntot += 1.
		#except:
		#	print i
		#print( i)
	print( Ntot)		
	xiave = xiave/float(Ntot)
	for i in range(1,N+1):
		nr = str(i)
		if md == 'me':
			xii = np.loadtxt(dir+'xi024'+samp+'_v'+str(ver)+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
		if md == 'cz':
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)		
			xii = np.loadtxt(dir+'2PCF_EZmock_eBOSS_'+samp+'_'+reg+'_v'+str(ver)+'_z'+str(zmin)+'z'+str(zmax)+'_'+mockn+'.dat').transpose()
		if md == 'jh':
			dirx = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys/'
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dirx+'ez_qsov7_0.31_'+str.lower(reg)+'_xi4m_wsys_ds5_'+mockn+'.dat').transpose()

		xic = np.concatenate((xii[1],xii[2]),axis=None)#-angfac*xii[mom/2+4] need to add this if desired
		for j in range(0,nbin*2):
			xij = xic[j]#-angfac*xiit[j]
			for k in range(0,nbin*2):
				xik = xic[k]#-angfac*xiit[k]
				cov[j][k] += (xij-xiave[j])*(xik-xiave[k])
#		except:
#			print i
	cov = cov/float(Ntot)					
	fo = open('xiave02'+reg+samp+'_v'+str(ver)+'_EZ'+rec+muw+bsst+'.dat','w')
	errl = []
	for i in range(0,nbin*2):
		if i < nbin:
			r = bs/2.+bs*i+start
		else:
			r = bs/2.+bs*(i-nbin)+start
		fo.write(str(r)+ ' '+str(xiave[i])+ ' '+str(sqrt(cov[i][i]))+'\n')
		errl.append(sqrt(cov[i][i]))
	fo.close()	
	fo = open('cov02'+reg+samp+'_v'+str(ver)+'_EZ'+rec+muw+bsst+'.dat','w')
	
	for i in range(0,nbin*2):
		for j in range(0,nbin*2):
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

def mkcov_mock_EZ(reg,mom=0,ver=5,samp='QSO',pcmass=False,bs=5,N=1000,maxr=200,start=0,mupow=0,mumin=0,mumax=1,angfac=0,md='me',rec='',zmin=0.6,sys=''):
	#dir = '/mnt/lustre/ashleyr/eboss/EZmock'+samp+'v'+str(ver)+'/'
	#if md == 'me':
	#	dir += '/xi/'
	dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ver)+'/'
	if md == 'cz' and rec != 'rec':
		dir += 'prerec/2PCF/'	
	if mom == 0:
		momi = 1
	if mom == 2:
		momi = 2	
	if pcmass:
		samp += 'pCMASS'
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(round(mumin,2))
	if mumax != 1:
		muw += 'mumax'+str(round(mumax,2))
	if mupow != 0:
		muw += 'mup'+str(mupow)		
	zw = ''
	if samp == 'ELG' and zmin != 0.6:
		zw += 'zmin'+str(zmin)	
	bsst = str(bs)+'st'+str(start)
	nbin = (maxr-start)//bs
	xiave = np.zeros((nbin))
	cov = np.zeros((nbin,nbin))

	Ntot = 0
	fac = 1.
	if samp == 'QSO':
		zmin = 0.8
		zmax = 2.2
	for i in range(1,N+1):
		nr = str(i)
		#print( i)
		#try:
		if md == 'me':
			xii = np.loadtxt(dir+'xi024'+samp+'_v'+str(ver)+reg+zw+'EZmock'+rec+sys+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
		if md == 'cz':
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dir+'2PCF_EZmock_eBOSS_'+samp+'_'+reg+'_v'+str(ver)+'_z'+str(zmin)+'z'+str(zmax)+'_'+mockn+'.dat').transpose()
		if md == 'jh':
			dirx = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys/'
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dirx+'ez_qsov7_0.31_'+str.lower(reg)+'_xi4m_wsys_ds5_'+mockn+'.dat').transpose()
		if md == 'jb':
			dirx = '/project/projectdirs/eboss/bautista/mocks/lrgs/EZmocks/LRG_v7_syst/multipoles/'
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dirx+'EZmock_eBOSS_LRGpCMASS_'+reg+'_v7'+rec+'_'+mockn+'_5mpc.mul').transpose()
		
		#xic = np.concatenate((xii[1],xii[2]),axis=None)#-angfac*xii[mom/2+4] need to add this if desired
		xic = xii[momi]
		xiave += xic
# 		except:
# 			if mom != 'rp':
# 				calcxi_mockELGEZ(i+1,reg=reg,bs=bs,mom=mom,mumin=mumin,mumax=mumax,start=start,rec=rec)
# 				xii = np.loadtxt(dir+'xi024'+samp+reg+'EZmock'+rec+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
# 				xiave += xii[mom/2+1]-angfac*xii[mom/2+4]
# 			print(i,mom)
		Ntot += 1.
		#except:
		#	print i
		#print( i)
	print( Ntot)		
	xiave = xiave/float(Ntot)
	for i in range(1,N+1):
		nr = str(i)
		if md == 'me':
			xii = np.loadtxt(dir+'xi024'+samp+'_v'+str(ver)+reg+zw+'EZmock'+rec+sys+nr+muw+str(bs)+'st'+str(start)+'.dat').transpose()
		if md == 'cz':
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)		
			xii = np.loadtxt(dir+'2PCF_EZmock_eBOSS_'+samp+'_'+reg+'_v'+str(ver)+'_z'+str(zmin)+'z'+str(zmax)+'_'+mockn+'.dat').transpose()
		if md == 'jh':
			dirx = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys/'
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dirx+'ez_qsov7_0.31_'+str.lower(reg)+'_xi4m_wsys_ds5_'+mockn+'.dat').transpose()
		if md == 'jb':
			dirx = '/project/projectdirs/eboss/bautista/mocks/lrgs/EZmocks/LRG_v7_syst/multipoles/'
			zer = ''
			if i < 10:
				zer += '0'
			if i < 100:
				zer += '0'
			if i < 1000:
				zer += 	'0'
			mockn = zer+str(i)	
			xii = np.loadtxt(dirx+'EZmock_eBOSS_LRGpCMASS_'+reg+'_v7'+rec+'_'+mockn+'_5mpc.mul').transpose()


		xic = xii[momi]
		#xic = np.concatenate((xii[1],xii[2]),axis=None)#-angfac*xii[mom/2+4] need to add this if desired
		for j in range(0,nbin):
			xij = xic[j]#-angfac*xiit[j]
			for k in range(0,nbin):
				xik = xic[k]#-angfac*xiit[k]
				cov[j][k] += (xij-xiave[j])*(xik-xiave[k])
#		except:
#			print i
	cov = cov/float(Ntot)					
	fo = open('xiave'+str(mom)+reg+zw+samp+'_v'+str(ver)+'_EZ'+rec+sys+muw+bsst+'.dat','w')
	errl = []
	for i in range(0,nbin):
		if i < nbin:
			r = bs/2.+bs*i+start
		else:
			r = bs/2.+bs*(i-nbin)+start
		fo.write(str(r)+ ' '+str(xiave[i])+ ' '+str(sqrt(cov[i][i]))+'\n')
		errl.append(sqrt(cov[i][i]))
	fo.close()	
	fo = open('cov'+str(mom)+reg+zw+samp+'_v'+str(ver)+'_EZ'+rec+sys+muw+bsst+'.dat','w')
	
	for i in range(0,nbin):
		for j in range(0,nbin):
			fo.write(str(cov[i][j])+' ')
		fo.write('\n')
	fo.close()
		
	return True



def compLRGBAOstats():
	xl = np.arange(-5,5,.1)
	norml = 1./sqrt(2.*np.pi)*np.exp(-1.*xl**2/2.)
	ap_pre, dap_pre, ap_pos, dap_pos = np.loadtxt('/Users/ashleyross/Dropbox/eboss/ap-pre-post.txt', unpack=1)
	at_pre, dat_pre, at_pos, dat_pos = np.loadtxt('/Users/ashleyross/Dropbox/eboss/at-pre-post.txt', unpack=1)
	#fap = np.loadtxt('/Users/ashleyross/Dropbox/eboss/ap-pre-post.txt').transpose()
	#fat = np.loadtxt('/Users/ashleyross/Dropbox/eboss/at-pre-post.txt').transpose()
	map = np.mean(ap_pos)
	msap = np.mean(dap_pos)
	print((map,msap,np.std(ap_pos)))
	pnorm = (ap_pos-map)/dap_pos
	plt.hist(pnorm,bins=30,normed=True)
	plt.plot(xl,norml,'k--')
	plt.show()
	mat = np.mean(at_pos)
	msat = np.mean(dat_pos)
	print((mat,msat,np.std(at_pos)))
	tnorm = (at_pos-mat)/dat_pos
	plt.hist(tnorm,bins=30,normed=True)
	plt.plot(xl,norml,'k--')
	plt.show()
	w = dap_pos < 0.02#0.8*dat_pos
	print((len(dap_pos[w])))
	print((np.mean(ap_pos[w]),np.mean(at_pos[w]),np.mean(dap_pos[w]),np.mean(dat_pos[w]),np.std(ap_pos[w]),np.std(at_pos[w])))
	plt.hist(pnorm[w],bins=10,normed=True)
	plt.plot(xl,norml,'k--')
	plt.show()
	
	return True

def putallBAOmocks(N=1000,sig=1,sigtest=.04,reg='NScombf',samp='ELGEZ',mock1=1,nmock=1000,bs=8,start=0,version='4',mb='',Bp='0.4',mumin=0,mumax=1,rec='',damp='0.5933.058.5',chitest=20):
	ma = 0
	sa = 0
	siga = 0
	chia = 0
	n = 0
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)

	bsst = str(bs)+'st'+str(start)
	ng = 0
	sg = 0
	ag = 0
	errg = 0
	wf = version+rec+damp+mb+muw+bsst
	wfnb = version+rec+damp+mb+muw
	nchi = 0
	fo = open('BAOfits'+samp+reg+wf+'.dat','w')
	for i in range(mock1,mock1+nmock):
		fl = ''
# 		if i < 1000:
# 			fl += '0'
# 		if i < 100:
# 			fl += '0'
# 		if i < 10:
# 			fl += '0'
		fl += str(i)
		if start != 'comb':
			a = sigreg_c12(dirsci+'EZmockELGv4/BAOfits/BAOxichil'+reg+samp+fl+wf)
		else:
			c1 = np.loadtxt(dirsci+'EZmockELGv4/BAOfits/BAOxichil'+reg+samp+fl+wfnb+'8st0.dat').transpose()
			c2 = np.loadtxt(dirsci+'EZmockELGv4/BAOfits/BAOxichil'+reg+samp+fl+wfnb+'8st2.dat').transpose()		
			c3 = np.loadtxt(dirsci+'EZmockELGv4/BAOfits/BAOxichil'+reg+samp+fl+wfnb+'8st4.dat').transpose()
			c4 = np.loadtxt(dirsci+'EZmockELGv4/BAOfits/BAOxichil'+reg+samp+fl+wfnb+'8st6.dat').transpose()
			ct = (c1+c2+c3+c4)/4.
			a = sigreg_c12(ct,md='l')		
		fo.write(str(a[0])+' '+str((a[2]-a[1])/2.)+'\n')
		if sig == 1:
			s1b = float(a[1]),float(a[2])
		if sig == 2:
			s1b = float(a[3]),float(a[4])	
		if s1b[0] > .8 and s1b[1] < 1.2:

			ma += a[0]
			sa += a[0]**2.
			sigone = (float(a[2])-float(a[1]))/2.
			if sigone < sigtest and abs((a[0]-1.)/sigone)<3.:
				ng += 1.
				sg += a[0]**2.
				ag += a[0]
				errg += sigone
			siga += (float(a[2])-float(a[1]))/2.
			chia += a[-1]
			n += 1.
		if a[-1] > chitest:
			nchi += 1	
	ma = ma/n
	sa = sqrt(sa/n-ma**2.)
	siga = siga/n
	chia = chia/n
	ag = ag/ng
	sg = sqrt(sg/ng-ag**2.)
	fo.close()
	return ma,sa,siga,chia,n,ng,ag,sg,errg/ng,nchi

	
def mkcov_mockELG_EZ_old(reg,bs=8,mom=0,N=1000,rec='_recon',v='v4'):
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
	print( Ntot)		
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
	print(( min(f[sys]),max(f[sys])))
	print( nr)
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
	print(( avebs,sysmin	))	
	print( 'total number, weighted number')
	print(( no,nt))
	print ('mean redshift')
	print((zm/nt))

	print(('total number of randoms/objects '+str(nr)+'/'+str(nt)))
	print(('number of randoms/objects outside tested range '+str(bsr)+'/'+str(bs)))			
	ave = nt/nr
	print(('average number of objects per random is '+ str(ave)))
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
	print(chin)
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
	print((sysmin,sysmax))	
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
	print (nr)
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
	#print avebs,sysmin		
	#print 'total number, weighted number'
# 	print no,nt
# 	print 'mean redshift'
# 	print zm/nt
# 
# 	print 'total number of randoms/objects '+str(nr)+'/'+str(nt)
# 	print 'number of randoms/objects outside tested range '+str(bsr)+'/'+str(bs)			
	ave = nt/nr
#	print 'average number of objects per random is '+ str(ave)
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
#	print chin
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
		print((ramin,ramax,len(f[w]),len(f[w])/float(len(fr[wr]))))
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
		print((len(f[w]),len(f[w])/float(len(fr[wr]))))
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
	print((sum(zhist[0])*dz))
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
		print (rp)
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
		rp = sqrt(1.-mu**2.)*r
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
	fo = open('/Users/ashleyross/eBOSS/wrp024ELG'+reg+'.dat','w')
	rpl = d[0]
	r = rmin
	
	while r < rmax:
		mu = dmu/2.
		xi0 = 0
		xi2 = 0
		xi4 = 0
	
		while mu < 1.:
			rp = sqrt(1.-mu**2.)*r
			if rp < 1.:
				xirp = d[1][0]
			else:
				bl = int(rp/1.)
				if bl >= len(d[1])-2:
					bh = bl
					print((r,rp,mu))
				else:
					bh = bl + 1
				fac = rp/1.-bl
				xirp = d[1][bh]*fac+(1.-fac)*d[1][bl]
			xi0 += dmu*xirp
			xi2 += dmu*5.*P2(mu)*xirp
			xi4 += dmu*9.*P4(mu)*xirp
			mu += dmu
		fo.write(str(r)+' '+str(xi0)+' '+str(xi2)+' '+str(xi4)+'\n')
		#print r
		r += dr
	fo.close()	
	return True

def mkw024_data(reg='SGC',md='mockave',mockn='',dmu=0.01,bs=8,rmax=200,zeff=0.85,rec=''):	
	from Cosmo import distance
	dz = distance(.31,.69)
	dm = dz.dc(zeff)*pi/180.
	if md == 'data':
		d = np.loadtxt('/Users/ashleyross/Dropbox/eBOSS/xirp0gebossELG_'+reg+'4_mz0.6xz1.1fkp1st0.dat').transpose()
	#d = np.loadtxt('/Users/ashleyross/eBOSS/wthgebossELG_'+reg+'4_mz0.6xz1.1fkp.dat').transpose()
		fo = open('/Users/ashleyross/Dropbox/eBOSS/wrp024ELG'+reg+'_data'+str(bs)+'st0.dat','w')
	if md == 'mockave':
		d = np.loadtxt('/Users/ashleyross/Dropbox/eBOSS/xiaverp'+reg+'ELG_EZangfac01st0.dat').transpose()
		fo = open('/Users/ashleyross/Dropbox/eBOSS/wrp024ELG'+reg+'_EZmockave'+str(bs)+'st0.dat','w')
	if md == 'mock':
		d = np.loadtxt('/Users/ashleyross/Dropbox/eBOSS/wrpELG'+reg+'EZmock'+rec+str(mockn)+'1st0.dat').transpose()
		fo = open('/Users/ashleyross/Dropbox/eBOSS/wrp024ELG'+reg+'_EZmock'+str(mockn)+str(bs)+'st0.dat','w')
	rpl = d[0]#*dm
	nb = int(rmax/bs)
	for i in range(0,nb):
		mu = dmu/2.
		rbmin = float(i*bs)
		rbmax = rbmin+bs
		xi0 = 0
		xi2 = 0
		xi4 = 0
		mutot = 0
		while mu < 1.:
			mum = mu#-dmu/2.*.999
			mux = mu#+dmu/2.*.999
			#print(mux)
			rpmin = sqrt(1.-mux**2.)*rbmin
			#if rpmin < 4.:
			#	print(rbmin,mu)
			rpmax = sqrt(1.-mum**2.)*rbmax
			w = (rpl > rpmin) & (rpl <= rpmax)
			if sum(w) == 0:
				#print(rpmin,rpmax,rbmin,mu)
				#xirp = d[1][int(rpmin/bs)]
				xirp = 0
			else:
			#print(mu,rpmin,rpmax,rbmin,rbmax,sum(w))
				xirp = np.mean(d[1][w])
				mutot += dmu
			xi0 += dmu*xirp
			xi2 += dmu*5.*P2(mu)*xirp
			xi4 += dmu*9.*P4(mu)*xirp
			mu += dmu
		r = i*bs+bs/2.
		fo.write(str(r)+' '+str(xi0)+' '+str(xi2)+' '+str(xi4)+'\n')
		print((r,mutot))
		
	fo.close()	
	return True

def mkwmubin_data(reg='SGC',dmu=0.01,bs=8,rmax=250,mumin=0,mumax=1):
	d = np.loadtxt('/Users/ashleyross/Dropbox/eBOSS/xirp0gebossELG_'+reg+'4_mz0.6xz1.1fkp1st0.dat').transpose()
	gf = ''
	if mumin != 0:
		gf += 'mum'+str(mumin)
	if mumax != 1.:
		gf += 'mux'+str(mumax)

	fo = open('/Users/ashleyross/Dropbox/eBOSS/wrp0'+gf+'ELG'+reg+'_data'+str(bs)+'st0.dat','w')
	rpl = d[0]
	nb = int(rmax/bs)
	for i in range(0,nb):
		mu = mumin + dmu/2.
		rbmin = float(i*bs)
		rbmax = rbmin+bs
		xi0 = 0
		mutot = 0
		while mu < mumax:
			rpmin = sqrt(1.-mu**2.)*rbmin
			#if rpmin < 4.:
			#	print(rbmin,mu)
			rpmax = sqrt(1.-mu**2.)*rbmax
			w = (rpl > rpmin) & (rpl <= rpmax)
			if sum(w) == 0:
				#print(rpmin,rpmax,rbmin,mu)
				#xirp = d[1][int(rpmin/bs)]
				xirp = 0
			else:
			#print(mu,rpmin,rpmax,rbmin,rbmax,sum(w))
				xirp = np.mean(d[1][w])
				mutot += dmu
			xi0 += dmu*xirp
			mu += dmu
		r = i*bs+bs/2.
		fo.write(str(r)+' '+str(xi0/mutot)+'\n')
		print((r,mutot))
		
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
		print((len(f[w])))
		hl.append(hist[0])
		w = ((f['galdepth_g']+f['galdepth_r']*2+f['galdepth_z']*10 > 700) & (f['g'] < gmax) & (sld > sldmin) & (f['gr'] > grmin) & (f['gr'] < grmax) & (f['rz'] > rzmin) & (f['rz'] < rzmax) & (f['chunk'] == chunk) & (rmag < rmax) & (zmag < zmax))
		hist = np.histogram(f[w]['Z'],bins=hbins,normed=True)
		xl = xlm
		plt.plot(xl,hist[0],color=cl[1])
		print((len(f[w])))
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
		print((len(f[w]),len(f[w])/float(len(fr[wr]))))
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
				print((brick,xfrac,n))
			if frac < mfrac:
				mfrac = frac
				print((brick,mfrac,n))
			if fbz > xbz:
				xbz = fbz
				print((brick,xbz,n))
			if fbz < mbz:
				mbz = fbz
				print((brick,mbz,n))
		else:
			print((brick,len(g),len(r)))
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
	print((n/float(len(f)),n,ng,nc1,nc2,nc3,nc4))
	return True

def ngalvdcol(res=128,mins=-.1,maxs=.1):
	import healpy as hp
	from healpix import radec2thphi
	from optimize import fmin
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'
	map = np.loadtxt('dmagDES_LSpix'+str(res)+'.dat').transpose()
	pixlgr = np.ones(res*res*12)*99
	pixlrz = np.ones(res*res*12)*99
	pixlg = np.ones(res*res*12)*99
	pixlr = np.ones(res*res*12)*99
	pixlz = np.ones(res*res*12)*99
	for i in range(0,len(map[0])):
		pix = int(map[0][i])
		pixlgr[pix] = map[1][i]-map[2][i]
		pixlrz[pix] = map[2][i]-map[3][i]
		pixlg[pix] = map[1][i]
		pixlr[pix] = map[2][i]
		pixlz[pix] = map[3][i]
	pixlgrz = pixlgr+pixlrz
	w = (pixlgr < 1)
	print((sum(w)))
	print((min(pixlgr[w]),max(pixlgr[w])))	
	f = fitsio.read(dir+'4/eBOSS_ELG_clustering_SGC_v4.dat.fits')
	thphi = radec2thphi(f['RA'],f['DEC'])	
	pix = hp.ang2pix(res,thphi[0],thphi[1])
	ng = np.zeros(res*res*12)
	for i in range(0,len(pix)):
		pixel = pix[i]
		ng[pixel] += f[i]['WEIGHT_CP']*f[i]['WEIGHT_NOZ']*f[i]['WEIGHT_FKP']
	fr = fitsio.read(dir+'4/eBOSS_ELG_clustering_SGC_v4.ran.fits')
	rann = len(fr)/sum(f['WEIGHT_FKP']*f['WEIGHT_CP']*f['WEIGHT_NOZ'])
	print((rann,sum(f['WEIGHT_FKP']*f['WEIGHT_CP']*f['WEIGHT_NOZ'])))
	thphir = radec2thphi(fr['RA'],fr['DEC'])	
	pixr = hp.ang2pix(res,thphir[0],thphir[1])
	nr = np.zeros(res*res*12)	
	for i in range(0,len(pixr)):
		pixel = pixr[i]
		nr[pixel] += fr[i]['WEIGHT_SYSTOT']
	nbin = 10
	hgr = np.zeros(nbin)
	hgr_ran = np.zeros(10)
	hrz = np.zeros(10)
	hrz_ran = np.zeros(10)
	binsize = (maxs-mins)/float(nbin)
	for i in range(0,res*res*12):
		if abs(pixlgr[i]) < maxs:
			bin = int((pixlgr[i]-mins)/binsize)	 	
			hgr[bin] += ng[bin]
			hgr_ran[bin] += nr[bin]
	print(hgr)
	print(hgr_ran)
	#xlm = np.arange(mins+binsize/2.,maxs,binsize)
	#plt.plot(xlm,hgr/hgr_ran*rann,'ko-')
	#plt.show()		
	histgrz = np.histogram(pixlgrz,bins=10,normed=False,weights=ng,range=(-0.02,.045))
	hbins = histgrz[1]
	histgrz_ran = np.histogram(pixlgrz,bins=hbins,normed=False,weights=nr,range=(-0.02,.045))	
	
	print(hbins)
	print((histgrz[0],sum(histgrz[0])))
	binsize = (histgrz[1][-1]-histgrz[1][0])/float(len(histgrz[1])-1)
	xlm = np.arange(histgrz[1][0]+binsize/2.,histgrz[1][-1]*1.001,binsize)
	plt.errorbar(xlm,histgrz[0]/histgrz_ran[0]*rann,np.sqrt(histgrz[0])/histgrz_ran[0]*rann,fmt='ko')
	plt.xlabel(r'$\Delta$ (g-r)+$\Delta$ (r-z)')
	print((sum((histgrz[0]/histgrz_ran[0]*rann-1.)**2./(np.sqrt(histgrz[0])/histgrz_ran[0]*rann)**2.)))
	plt.show()


	histgr = np.histogram(pixlgr,bins=10,normed=False,range=(-.02,.02),weights=ng)
	hbins = histgr[1]
	histgr_ran = np.histogram(pixlgr,bins=hbins,normed=False,range=(-.02,.02),weights=nr)	
	
	print(hbins)
	print((histgr[0],sum(histgr[0])))
	binsize = (histgr[1][-1]-histgr[1][0])/float(len(histgr[1])-1)
	xlm = np.arange(histgr[1][0]+binsize/2.,histgr[1][-1]*1.001,binsize)
	plt.errorbar(xlm,histgr[0]/histgr_ran[0]*rann,np.sqrt(histgr[0])/histgr_ran[0]*rann,fmt='ko')
	plt.xlabel(r'$\Delta$ (g-r)')
	print((sum((histgr[0]/histgr_ran[0]*rann-1.)**2./(np.sqrt(histgr[0])/histgr_ran[0]*rann)**2.)))
	plt.show()

	histg = np.histogram(pixlg,bins=10,normed=False,range=(0.01,.04),weights=ng)
	hbins = histg[1]
	histg_ran = np.histogram(pixlg,bins=hbins,normed=False,range=(0.01,.04),weights=nr)	
	
	print(hbins)
	print((histg[0],sum(histg[0])))
	binsize = (histg[1][-1]-histg[1][0])/float(len(histg[1])-1)
	xlm = np.arange(histg[1][0]+binsize/2.,histg[1][-1]*1.001,binsize)
	plt.errorbar(xlm,histg[0]/histg_ran[0]*rann,np.sqrt(histg[0])/histg_ran[0]*rann,fmt='ko')
	plt.xlabel(r'$\Delta$ (g)')
	print((sum((histg[0]/histg_ran[0]*rann-1.)**2./(np.sqrt(histg[0])/histg_ran[0]*rann)**2.)))
	plt.show()


	histrz = np.histogram(pixlrz,bins=10,normed=False,range=(-.015,.035),weights=ng)
	hbins = histrz[1]
	histrz_ran = np.histogram(pixlrz,bins=hbins,normed=False,range=(-.015,.035),weights=nr)	
	
	print(hbins)
	print((histrz[0],sum(histrz[0])))
	binsize = (histrz[1][-1]-histrz[1][0])/float(len(histrz[1])-1)
	xlm = np.arange(histrz[1][0]+binsize/2.,histrz[1][-1]*1.001,binsize)
	plt.errorbar(xlm,histrz[0]/histrz_ran[0]*rann,np.sqrt(histrz[0])/histrz_ran[0]*rann,fmt='ko')
	plt.xlabel(r'$\Delta$ (r-z)')
	print((sum((histrz[0]/histrz_ran[0]*rann-1.)**2./(np.sqrt(histrz[0])/histrz_ran[0]*rann)**2.)))
	plt.show()

	histr = np.histogram(pixlr,bins=10,normed=False,range=(0.005,.04),weights=ng)
	hbins = histr[1]
	histr_ran = np.histogram(pixlr,bins=hbins,normed=False,range=(0.005,.04),weights=nr)	
	
	print(hbins)
	print((histr[0],sum(histr[0])))
	binsize = (histr[1][-1]-histr[1][0])/float(len(histr[1])-1)
	xlm = np.arange(histr[1][0]+binsize/2.,histr[1][-1]*1.001,binsize)
	plt.errorbar(xlm,histr[0]/histr_ran[0]*rann,np.sqrt(histr[0])/histr_ran[0]*rann,fmt='ko')
	plt.xlabel(r'$\Delta$ (r)')
	print((sum((histr[0]/histr_ran[0]*rann-1.)**2./(np.sqrt(histr[0])/histr_ran[0]*rann)**2.)))
	plt.show()

	histz = np.histogram(pixlz,bins=10,normed=False,range=(-0.01,.035),weights=ng)
	hbins = histz[1]
	histz_ran = np.histogram(pixlz,bins=hbins,normed=False,range=(-0.01,.035),weights=nr)	
	
	print(hbins)
	print((histz[0],sum(histz[0])))
	binsize = (histz[1][-1]-histz[1][0])/float(len(histz[1])-1)
	xlm = np.arange(histz[1][0]+binsize/2.,histz[1][-1]*1.001,binsize)
	plt.errorbar(xlm,histz[0]/histz_ran[0]*rann,np.sqrt(histz[0])/histz_ran[0]*rann,fmt='ko')
	plt.xlabel(r'$\Delta$ (z)')
	print((sum((histz[0]/histz_ran[0]*rann-1.)**2./(np.sqrt(histz[0])/histz_ran[0]*rann)**2.)))
	plt.show()


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
		print((hist[0]))
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
		
		print((histr[0]))
		print((hist[0]/histr[0]*nr))
		print(('nan ratio is',len(denan)/float(len(dernan))*nr,len(denan),len(dernan)/float(len(fr))))
		fo = open('nELG'+reg+'vsdEBV.dat','w')
		for i in range(0,len(xlm)):
			fo.write(str(xlm[i])+' '+str(hist[0][i]/histr[0][i].astype('float')*nr)+' '+str(sqrt(hist[0][i])/histr[0][i].astype('float')*nr)+'\n')
		fo.close()
		lf = linfit(xlm,hist[0]/histr[0].astype('float')*nr,np.sqrt(hist[0])/histr[0].astype('float')*nr)
		inl = np.array([1.,-10.])
		b0,m0 = fmin(lf.chilin,inl)
		print((b0,m0))
		print(( lf.chilin((b0,m0))))
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

def ngalvH2():
	import healpy as hp
	from healpix import radec2thphi
	from optimize import fmin
	from astropy.io import fits
	H2map = fits.open('NHI_HPX.fits.gz')[1].data['NHI']
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
		#emap = fitsio.read('ebv_lhd.hpx.fits')['EBV']
		ebvl = []
		for i in range(0,len(pix)):
			ebvl.append(log(H2map[pix[i]]))
		ebvl = np.array(ebvl)
		de = ebvl#-f['ebv']
		#denan = de[np.isnan(de)]
		#wts=f[np.isfinite(de)]['WEIGHT_SYSTOT']
		w = (np.isfinite(de))# & (de > -0.01) & (de < 0.015))
		denan = de[~w]
		de = de[w]
		wts = f[w]['WEIGHT_SYSTOT']
		#be = [-0.06,-0.02,-0.01,-0.007,-0.005,-0.003,-0.001,0.001,0.005,0.01,0.02]
		hist = np.histogram(de,bins=10,normed=False,weights=wts)	
		hbins = hist[1]
		print(hbins)
		print((hist[0]))
		binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
		#xlm = []
		#for i in range(0,len(be)-1):
		#	xlm.append((be[i]+be[i+1])/2.) 
		xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
	
		fr = fitsio.read(dir+'4/eBOSS_ELG_clustering_'+reg+'_v4.ran.fits')
		#nr = len(fr)/float(len(f))
		nr = sum(fr['WEIGHT_SYSTOT'])/sum(f['WEIGHT_SYSTOT'])
		print(nr)
		thphir = radec2thphi(fr['RA'],fr['DEC'])	
		thphiG = r(thphir[0],thphir[1])
		pixr = hp.ang2pix(1024,thphiG[0],thphiG[1])
		ebvlr = []
		for i in range(0,len(pixr)):
			ebvlr.append(log(H2map[pixr[i]]))
		ebvlr = np.array(ebvlr)
		der = ebvlr#-fr['ebv']
		#dernan = der[np.isnan(der)]
		wr = (np.isfinite(der))# & (der > -0.01) & (der < 0.015))
		dernan = der[~wr]
		der = der[wr]
		wts = fr[wr]['WEIGHT_SYSTOT']
		histr = np.histogram(der,bins=hbins,normed=False,weights=wts)
		#nr = len(der)/float(len(de))
		#print(nr)
		
		print((histr[0]))
		print((hist[0]/histr[0]*nr))
		#print('nan ratio is',len(denan)/float(len(dernan))*nr,len(denan),len(dernan)/float(len(fr)))
		fo = open('nELG'+reg+'vsHII.dat','w')
		chinull = 0
		for i in range(0,len(xlm)):
			fo.write(str(xlm[i])+' '+str(hist[0][i]/histr[0][i].astype('float')*nr)+' '+str(sqrt(hist[0][i])/histr[0][i].astype('float')*nr)+'\n')
			diff = hist[0][i]/histr[0][i].astype('float')*nr-1.
			err = sqrt(hist[0][i])/histr[0][i].astype('float')*nr
			print((diff,err))
			chinull += (diff)**2./(err**2.)
		fo.close()
		print(chinull)
		lf = linfit(xlm,hist[0]/histr[0].astype('float')*nr,np.sqrt(hist[0])/histr[0].astype('float')*nr)
		inl = np.array([1.,-10.])
		b0,m0 = fmin(lf.chilin,inl)
		print((b0,m0))
		print(( lf.chilin((b0,m0))))
		plt.plot(xlm,b0+m0*xlm,'--',color=cl[k])
		plt.errorbar(xlm,hist[0]/histr[0].astype('float')*nr,np.sqrt(hist[0])/histr[0].astype('float')*nr,color=cl[k])

	#plt.xlim(-0.05,0.02)
	plt.xlabel(r'HII')
	plt.ylabel(r'$n_{\rm gal}/\langle n_{\rm gal} \rangle$')
	plt.text(hbins[1],1.1,'NGC',color='r')
	plt.text(hbins[1],1.08,'SGC',color='b')
	#plt.show()
	plt.savefig('nELGvsHII.png')
	return True


def nQSOvdebv(zmin=1.5):
	import healpy as hp
	from healpix import radec2thphi
	from optimize import fmin
	dir = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'
	regl = ['SGC','NGC']
	cl  =['b','r']
	emap = fitsio.read('ebv_lhd.hpx.fits')['EBV']
	SFD = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/maps/SDSSimageprop_Nside512.fits')['EBV']
	for k in range(0,2):
		reg = regl[k]
		f = fitsio.read(dir+'4/eBOSS_QSO_clustering_'+reg+'_v4.dat.fits')
		wz = (f['Z'] > zmin)
		thphi = radec2thphi(f[wz]['RA'],f[wz]['DEC'])	
		pixsfd = hp.ang2pix(512,thphi[0],thphi[1])
		sfdl = []
		for i in range(0,len(pixsfd)):
			sfdl.append(SFD[pixsfd[i]])
		sfdl = np.array(sfdl)		
		r = hp.Rotator(coord=['C','G'],deg=False)
		thphiG = r(thphi[0],thphi[1])
		pix = hp.ang2pix(1024,thphiG[0],thphiG[1])		
		ebvl = []
		for i in range(0,len(pix)):
			ebvl.append(emap[pix[i]])
		ebvl = np.array(ebvl)
		de = ebvl-sfdl
		#denan = de[np.isnan(de)]
		wts=f[wz][np.isfinite(de)]['WEIGHT_SYSTOT']
		w = (np.isfinite(de))# & (de > -0.02) & (de < 0.015))
		denan = de[~w]
		de = de[w]
		#plt.plot(sfdl[w],de,'ko')
		#plt.show()
		#plt.plot(ebvl[w],de,'ko')
		#plt.show()
		
		be = [-0.06,-0.02,-0.01,-0.007,-0.005,-0.003,-0.001,0.001,0.005,0.01,0.02]
		hist = np.histogram(de,bins=10,normed=False)#,weights=wts)	
		hbins = hist[1]
		print(hbins)
		print((hist[0]))
		binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
		#xlm = []
		#for i in range(0,len(be)-1):
		#	xlm.append((be[i]+be[i+1])/2.) 
		xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
	
		fr = fitsio.read(dir+'4/eBOSS_QSO_clustering_'+reg+'_v4.ran.fits')
		nr = len(fr)/float(len(f))
		print(nr)
		thphir = radec2thphi(fr['RA'],fr['DEC'])	
		pixsfdr = hp.ang2pix(512,thphir[0],thphir[1])
		sfdlr = []
		for i in range(0,len(pixsfdr)):
			sfdlr.append(SFD[pixsfdr[i]])
		sfdlr = np.array(sfdlr)		


		thphiG = r(thphir[0],thphir[1])
		pixr = hp.ang2pix(1024,thphiG[0],thphiG[1])
		ebvlr = []
		for i in range(0,len(pixr)):
			ebvlr.append(emap[pixr[i]])
		ebvlr = np.array(ebvlr)
		der = ebvlr-sfdlr
		#dernan = der[np.isnan(der)]
		wr = (np.isfinite(der))# & (der > -0.02) & (der < 0.015))
		dernan = der[~wr]
		der = der[wr]
		rap = fr[wr]['RA']
		if reg == 'SGC':
			wra = (rap >180)
			rap[wra] -= 360
		plt.scatter(rap,fr[wr]['DEC'],c=der,edgecolors='face',marker='.')
		plt.colorbar()
		plt.show()
		histr = np.histogram(der,bins=hbins,normed=False)
		nr = len(der)/float(len(de))
		print(nr)
		
		print((histr[0]))
		print((hist[0]/histr[0]*nr))
		print(('nan ratio is',len(denan)/float(len(dernan))*nr,len(denan),len(dernan)/float(len(fr))))
		fo = open('nQSO'+reg+'vsdEBV.dat','w')
		for i in range(0,len(xlm)):
			fo.write(str(xlm[i])+' '+str(hist[0][i]/histr[0][i].astype('float')*nr)+' '+str(sqrt(hist[0][i])/histr[0][i].astype('float')*nr)+'\n')
		fo.close()
		lf = linfit(xlm,hist[0]/histr[0].astype('float')*nr,np.sqrt(hist[0])/histr[0].astype('float')*nr)
		inl = np.array([1.,-10.])
		b0,m0 = fmin(lf.chilin,inl)
		print((b0,m0))
		print(( lf.chilin((b0,m0))))
		plt.plot(xlm,b0+m0*xlm,'--',color=cl[k])
		plt.errorbar(xlm,hist[0]/histr[0].astype('float')*nr,np.sqrt(hist[0])/histr[0].astype('float')*nr,color=cl[k])
		plt.show()
	#plt.xlim(-0.05,0.02)
	plt.xlabel(r'$\Delta$E(B-V) (Lenz et al - SFD)')
	plt.ylabel(r'$n_{\rm QSO}/\langle n_{\rm gal} \rangle$')
	#plt.text(-0.04,1.1,'NGC',color='r')
	#plt.text(-0.04,1.08,'SGC',color='b')
	plt.show()
	#plt.savefig('nELGvsdeltaEBV.png')
	return True

def readNOAOls(file):
	f = open(file)
	f.readline()
	a = f.readlines()
	ral = []
	decl = []
	for i in range(0,len(a)):
		ln = a[i].split(',')
		try:
			ra = float(ln[0])
			dec = float(ln[1])
			tp = ln[-1].strip('\n')
			if tp == 'PSF':
				ral.append(ra)
				decl.append(dec)
		except:
			pass	
	from matplotlib import pyplot as plt
	plt.plot(ral,decl,',')
	plt.show()

def readNOAOdes(file):
	f = open(file)
	f.readline()
	a = f.readlines()
	ral = []
	decl = []
	for i in range(0,len(a)):
		ln = a[i].split(',')
		try:
			ra = float(ln[0])
			dec = float(ln[1])
			ral.append(ra)
			decl.append(dec)
		except:
			pass	
	print((len(ral)))
	from matplotlib import pyplot as plt
	plt.plot(ral,decl,',')
	plt.show()


def put4dr7():
	from astropy.table import Table
	#'ra,dec,dered_mag_g,dered_mag_r,dered_mag_z,type\n'
	ral = []
	decl = []
	gl = []
	rl = []
	zl = []
	for i in range(0,6):
		file = '/Users/ashleyross/eBOSS/dr7g1821_'+str(i)+'.txt'
		f = open(file)
		f.readline()
		a = f.readlines()
		for i in range(0,len(a)):
			ln = a[i].split(',')
			try:
				ra = float(ln[0])
				dec = float(ln[1])
				g = float(ln[2])
				r = float(ln[3])
				z = float(ln[4])
				tp = ln[-1].strip('\n')
				if tp == 'PSF':
					ral.append(ra)
					decl.append(dec)
					gl.append(g)
					rl.append(r)
					zl.append(z)
			except:
				pass
	print((len(ral)))
	tab = Table([ral,decl,gl,rl,zl],names=('RA','DEC','dered_mag_g','dered_mag_r','dered_mag_z'))	
	tab.write('decalsDR7_g1821_ELGSGC_PSF.fits', format='fits', overwrite=True) 			
	from matplotlib import pyplot as plt
	plt.plot(ral,decl,',')
	plt.show()
	
def put4des():	
	#'ra,dec,wavg_mag_psf_g_dered,wavg_mag_psf_r_dered,wavg_mag_psf_z_dered\n'\n'
	ral = []
	decl = []
	gl = []
	rl = []
	zl = []
	for i in range(0,6):
		file = '/Users/ashleyross/eBOSS/DESdr1g1821_'+str(i)+'.txt'
		f = open(file)
		f.readline()
		a = f.readlines()
		for i in range(0,len(a)):
			ln = a[i].split(',')
			try:
				ra = float(ln[0])
				dec = float(ln[1])
				g = float(ln[2])
				r = float(ln[3])
				z = float(ln[4])
				ral.append(ra)
				decl.append(dec)
				gl.append(g)
				rl.append(r)
				zl.append(z)
			except:
				pass
	print((len(ral)))
	tab = Table([ral,decl,gl,rl,zl],names=('RA','DEC','wavg_mag_psf_g_dered','wavg_mag_psf_r_dered','wavg_mag_psf_z_dered'))	
	tab.write('desDR1_g1821_SGCELG_wavgpsf.fits', format='fits', overwrite=True) 			
	from matplotlib import pyplot as plt
	plt.plot(ral,decl,',')
	plt.show()

def matchDESdecals(angle=1/3600.,ramin=0,ramax=50): 
	des = fitsio.read('/Users/ashleyross/eBOSS/desDR1_g1821_SGCELG_wavgpsf.fits')
	ls = fitsio.read('/Users/ashleyross/eBOSS/decalsDR7_g1821_SGCELG_PSF.fits')
	#w = (ls['dered_mag_g'] < 20.9)
	w = ((ls['RA'] > ramin) & (ls['RA'] < ramax))
	ls = ls[w]
	c1 = SkyCoord(ra=des['RA']*u.degree, dec=des['DEC']*u.degree)     
	c2 = SkyCoord(ra=ls['RA']*u.degree, dec=ls['DEC']*u.degree)
	idx, d2d, d3d = c1.match_to_catalog_sky(c2)  
	w = d2d.value <= angle
	idx[~w] = -1

	idx1 = np.where(w)[0]
	idx2 = idx[idx>-1] 
	distance = d2d.value[w]
	print((len(des[idx1]),len(ls[idx2])))
	
	mdes = des[idx1]
	mls = ls[idx2]
	
	tab = Table([mdes['RA'],mdes['DEC'],mdes['wavg_mag_psf_g_dered'],mdes['wavg_mag_psf_r_dered'],mdes['wavg_mag_psf_z_dered'],mls['dered_mag_g'],mls['dered_mag_r'],mls['dered_mag_z']],names=('RA','DEC','wavg_mag_psf_g_dered','wavg_mag_psf_r_dered','wavg_mag_psf_z_dered','dered_mag_g','dered_mag_r','dered_mag_z'))	
	tab.write('matched_lsdesDR1_g1821_SGCELG_PSF.fits', format='fits', overwrite=True) 			
	

	from matplotlib import pyplot as plot
	plt.plot(des[idx1]['wavg_mag_psf_g_dered'],ls[idx2]['dered_mag_g'],'k,')
	#xl = [18,21]
	#yl = [18,21]
	#plt.plot(xl,yl,'r--')
	plt.hist(des[idx1]['wavg_mag_psf_g_dered']-ls[idx2]['dered_mag_g'],bins='auto',range=(-.1,.1))
	plt.xlim(-.12,.12)
	print((np.median(des[idx1]['wavg_mag_psf_g_dered']-ls[idx2]['dered_mag_g'])))
	plt.show()
	
	#plt.plot(des[idx1]['wavg_mag_psf_r_dered'],ls[idx2]['dered_mag_r'],'k,')
	#xl = [15,23]
	#yl = [15,23]
	#plt.plot(xl,yl,'r--')
	#plt.xlim(15,23)
	#plt.ylim(15,23)
	plt.hist(des[idx1]['wavg_mag_psf_r_dered']-ls[idx2]['dered_mag_r'],bins='auto',range=(-.1,.1))
	plt.xlim(-.12,.12)
	plt.xlabel('r diff')
	diff = des[idx1]['wavg_mag_psf_r_dered']-ls[idx2]['dered_mag_r']
	diff = diff[~np.isnan(diff)]
	print((np.median(diff)))
	plt.show()
	
	#plt.plot(des[idx1]['wavg_mag_psf_z_dered'],ls[idx2]['dered_mag_z'],'k,')
	#xl = [15,23]
	#yl = [15,23]
	#plt.plot(xl,yl,'r--')
	#plt.xlim(15,23)
	#plt.ylim(15,23)
	plt.hist(des[idx1]['wavg_mag_psf_z_dered']-ls[idx2]['dered_mag_z'],bins='auto',range=(-.1,.1))
	plt.xlim(-.12,.12)
	plt.xlabel('z diff')
	diff = des[idx1]['wavg_mag_psf_z_dered']-ls[idx2]['dered_mag_z']
	diff = diff[~np.isnan(diff)]
	print((np.median(diff)))
	plt.show()
	plt.clf()
	return True

def mkdiffplots_desls_hp(res=128,psize=50.):
	import healpy as hp
	from healpix import radec2thphi

	f = fitsio.read('/Users/ashleyross/eBOSS/matched_lsdesDR1_g1821_SGCELG_PSF.fits')

	thphi = radec2thphi(f['RA'],f['DEC'])	
	pix = hp.ang2pix(res,thphi[0],thphi[1])
	ral = np.zeros(max(pix)+1)
	decl = np.zeros(max(pix)+1)
	gl = np.zeros(max(pix)+1)
	rl = np.zeros(max(pix)+1)
	zl = np.zeros(max(pix)+1)
	nl = np.zeros(max(pix)+1)
	#pixl = np.zeros(max(pix)+1)
	for i in range(0,len(pix)):
		dg = f[i]['wavg_mag_psf_g_dered']-f[i]['dered_mag_g']
		dr = f[i]['wavg_mag_psf_r_dered']-f[i]['dered_mag_r']
		dz = f[i]['wavg_mag_psf_z_dered']-f[i]['dered_mag_z']
		
		if 0*dg == 0 and 0*dr == 0 and 0*dz == 0:
			if abs(dg) < 0.2 and abs(dr) < 0.2 and abs(dz) < 0.2:
				pixel = pix[i]
				nl[pixel] += 1.
				ra = f[i]['RA']
				if ra > 180:
					ra -= 360
				ral[pixel] += ra
				decl[pixel] += f[i]['DEC']
				gl[pixel] += dg
				rl[pixel] += dr
				zl[pixel] += dz
	print((sum(nl)))
	
	ral = ral/nl
	decl = decl/nl
	gl = gl/nl
	rl = rl/nl
	zl = zl/nl
	grl = gl-rl
	rzl = rl-zl			
	fo = open('dmagDES_LSpix'+str(res)+'.dat','w')
	for i in range(0,len(ral)):
		if nl[i] > 0:
			fo.write(str(i)+' '+str(gl[i])+' '+str(rl[i])+' '+str(zl[i])+'\n')
	fo.close()	
	sz = np.ones(len(ral))*psize
	plt.scatter(ral,decl,c=gl,vmin=-.1,vmax=.1,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$g$ DES - $g$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$; Nside ='+str(res))
	plt.savefig('/Users/ashleyross/eBOSS/deltag_chunk22_DESDECaLS_hp'+str(res)+'.png',dpi='figure')
	plt.clf()

	plt.scatter(ral,decl,c=gl,vmin=-.1,vmax=.1,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$r$ DES - $r$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$; Nside ='+str(res))
	plt.savefig('/Users/ashleyross/eBOSS/deltar_chunk22_DESDECaLS_hp'+str(res)+'.png',dpi='figure')
	plt.clf()

	plt.scatter(ral,decl,c=zl,vmin=-.1,vmax=.1,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$z$ DES - $z$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$; Nside ='+str(res))
	plt.savefig('/Users/ashleyross/eBOSS/deltaz_chunk22_DESDECaLS_hp'+str(res)+'.png',dpi='figure')
	plt.clf()

	plt.scatter(ral,decl,c=grl,vmin=-.05,vmax=.05,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$g-r$ DES - $g-r$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$; Nside ='+str(res))
	plt.savefig('/Users/ashleyross/eBOSS/deltagr_chunk22_DESDECaLS_hp'+str(res)+'.png',dpi='figure')
	plt.clf()

	plt.scatter(ral,decl,c=rzl,vmin=-.05,vmax=.05,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$r-z$ DES - $r-z$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$; Nside ='+str(res))
	plt.savefig('/Users/ashleyross/eBOSS/deltarz_chunk22_DESDECaLS_hp'+str(res)+'.png',dpi='figure')
	plt.clf()

	return True

def mkdiffplots_desls():
	f = fitsio.read('/Users/ashleyross/eBOSS/matched_lsdesDR1_g1821_ra050_absdec5_PSF.fits')
	sz = np.ones(len(f['RA']))*.1
	dg = f['wavg_mag_psf_g_dered']-f['dered_mag_g']
	plt.scatter(f['RA'],f['DEC'],c=dg,vmin=-.1,vmax=.1,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$g$ DES - $g$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$)')
	plt.savefig('/Users/ashleyross/eBOSS/deltag_chunk22_DESDECaLS.png',dpi='figure')
	plt.clf()
	
	dr = f['wavg_mag_psf_r_dered']-f['dered_mag_r']
	plt.scatter(f['RA'],f['DEC'],c=dr,vmin=-.1,vmax=.1,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$r$ DES - $r$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$)')
	plt.savefig('/Users/ashleyross/eBOSS/deltar_chunk22_DESDECaLS.png',dpi='figure')
	plt.clf()
	
	dz = f['wavg_mag_psf_z_dered']-f['dered_mag_z']
	plt.scatter(f['RA'],f['DEC'],c=dz,vmin=-.1,vmax=.1,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$z$ DES - $z$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$)')
	plt.savefig('/Users/ashleyross/eBOSS/deltaz_chunk22_DESDECaLS.png',dpi='figure')
	plt.clf()
	
	dgr = dg-dr
	plt.scatter(f['RA'],f['DEC'],c=dgr,vmin=-.05,vmax=.05,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$g-r$ DES - $g-r$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$)')
	plt.savefig('/Users/ashleyross/eBOSS/deltagr_chunk22_DESDECaLS.png',dpi='figure')
	plt.clf()
	
	drz = dr-dz
	plt.scatter(f['RA'],f['DEC'],c=drz,vmin=-.05,vmax=.05,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$r-z$ DES - $r-z$ DECaLS (PSF magnitudes DECaLS $18 < g < 21$)')
	plt.savefig('/Users/ashleyross/eBOSS/deltarz_chunk22_DESDECaLS.png',dpi='figure')
	plt.clf()
	return True

	
def DEBVchunk22():
	import healpy as hp
	from healpix import radec2thphi

	f = fitsio.read('/Users/ashleyross/eBOSS/matched_lsdesDR1_g1821_ra050_absdec5_PSF.fits')
	sz = np.ones(len(f['RA']))*.1
	emap = fitsio.read('maps/ebv_lhd.hpx.fits')['EBV']
	SFD = fitsio.read('maps/SDSSimageprop_Nside512.fits')['EBV']
	
	thphi = radec2thphi(f['RA'],f['DEC'])	
	pixsfd = hp.ang2pix(512,thphi[0],thphi[1])
	sfdl = []
	for i in range(0,len(pixsfd)):
		sfdl.append(SFD[pixsfd[i]])
	sfdl = np.array(sfdl)		
	r = hp.Rotator(coord=['C','G'],deg=False)
	thphiG = r(thphi[0],thphi[1])
	pix = hp.ang2pix(1024,thphiG[0],thphiG[1])		
	ebvl = []
	for i in range(0,len(pix)):
		ebvl.append(emap[pix[i]])
	ebvl = np.array(ebvl)
	de = ebvl-sfdl
	grd = (3.303-2.285)*de
	rzd = (2.285-1.263)*de
	plt.scatter(f['RA'],f['DEC'],c=grd,vmin=-.05,vmax=.05,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$\Delta$ $g-r$ from SFD/Lenz $\Delta$ E(B-V)')
	plt.savefig('/Users/ashleyross/eBOSS/deltagr_chunk22_EBV.png',dpi='figure')
	plt.clf()
	plt.scatter(f['RA'],f['DEC'],c=rzd,vmin=-.05,vmax=.05,marker=',',s=sz)
	plt.colorbar()
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(r'$\Delta$ $r-z$ from SFD/Lenz $\Delta$ E(B-V)')
	plt.savefig('/Users/ashleyross/eBOSS/deltarz_chunk22_EBV.png',dpi='figure')
	plt.clf()
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
		print((hist[0]))
		binsize = (hist[1][-1]-hist[1][0])/float(len(hist[1])-1)
		xlm = np.arange(hist[1][0]+binsize/2.,hist[1][-1]*1.001,binsize)
		print((len(f[w])))
		plt.plot(xlm,hist[0],color=cl[0])
		if reg == 'NGC':
			w = (np.isfinite(de) & (de > 0.00) & (f['chunk']=='eboss23'))
		else:
			w = (np.isfinite(de) & (de > 0.00))	
		hist = np.histogram(f[w]['Z'],bins=hbins,normed=True)
		plt.plot(xlm,hist[0],color=cl[1])
		print((len(f[w])))
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
	print((min(d['ntar']),min(d['nran']),min(d['nzfail'])))
	w = (d['nran'] > 0)
	d = d[w]
	nexp_t = d['nexp_g']+d['nexp_r']+d['nexp_z']
	print((max(nexp_t)))
	print((max(d['ntar']),max(d['nran']),max(d['nzfail']),max(d['nexp_g']),max(d['nexp_r']),max(d['nexp_z']),max(d['nobjs']),min(d['nobjs'])))	
	print((min(d['ntar']/d['nran']),max(d['ntar']/d['nran'])))
	print((min(d['nzfail']/d['ntar']),max(d['nzfail']/d['ntar'])))
	from matplotlib import pyplot as plt
	plt.plot(d['ntar']/d['nran'],d['nran'],'ko')
	plt.show()
	wr = (d['nran'] > 100)
	print((len(d[wr])/float(len(d))))
	drancut = d[wr]
	plt.hist(np.log(drancut['ntar']/drancut['nran']),bins=30)
	plt.show()
	plt.hist(drancut['nzfail']/drancut['ntar'],bins=30)
	plt.show()
	wzf = (drancut['nzfail']/drancut['ntar'] < 0.3)
	print((len(drancut[wzf])/float(len(drancut))))
	ral = []
	decl =[]
	for ng in range(1,22):
		w = (d['nexp_g'] ==ng)
		if len(d[w])>0:
			print((ng,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
	for nr in range(1,20):
		w = (d['nexp_r'] ==nr)
		if len(d[w])>0:
			print((nr,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
	for nz in range(1,17):
		w = (d['nexp_z'] ==nz)
		if len(d[w])>0:
			print((nz,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
	ntl = []
	densl = []
	for nt in range(3,47):
		w = (nexp_t ==nt)
		if len(d[w])>0:
			print((nt,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
			#ntl.append(nt)
			#densl.append(np.mean(d[w]['ntar']/d[w]['nran']))
	#plt.plot(ntl,densl)
	#plt.show()
	dgz = d['nexp_g']-d['nexp_z']
	print('difference between number of g and z exposures')
	print((min(dgz),max(dgz)))
	for ndgz in range(-5,18):
		w = (dgz ==ndgz)
		if len(d[w])>0:
			print((ndgz,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
			ntl.append(ndgz)
			densl.append(np.mean(d[w]['ntar']/d[w]['nran']))
	plt.plot(ntl,densl)
	plt.show()

	ntl = []
	densl = []		
	drz = d['nexp_r']-d['nexp_z']
	print('difference between number of r and z exposures')
	print((min(drz),max(drz)))
	for ndrz in range(-10,18):
		w = (drz ==ndrz)
		if len(d[w])>0:
			print((ndrz,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
			ntl.append(ndrz)
			densl.append(np.mean(d[w]['ntar']/d[w]['nran']))
	plt.plot(ntl,densl)
	plt.show()

	ntl = []
	densl = []
	dgr = d['nexp_g']-d['nexp_r']
	print('difference between number of g and r exposures')
	print((min(dgr),max(dgr)))
	for ndgr in range(-12,18):
		w = (dgr ==ndgr)
		if len(d[w])>0:
			print((ndgr,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
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
			print((mx,np.mean(d[w]['ntar']/d[w]['nran']),sum(d[w]['nran'])))
	w = ((dgz < -2) | (dgz > 5))
	#fo = open('brickzexpg7.dat','w')
	dg = d[w]
	#for i in range(0,len(dg)):
	#	fo.write(dg[i]['brick']+'\n')
	#fo.close()
	print((sum(d[w]['nran']),sum(d[w]['ntar'])))
	plt.plot(d[w]['ra'],d[w]['dec'],'bo')
	plt.xlim(0,50)
	#plt.show()
	
	w = (dgr > 4)
	#fo = open('brickzexpg7.dat','w')
	dg = d[w]
	#for i in range(0,len(dg)):
	#	fo.write(dg[i]['brick']+'\n')
	#fo.close()
	print((sum(d[w]['nran']),sum(d[w]['ntar'])))
	plt.plot(d[w]['ra'],d[w]['dec'],'ro')
	plt.xlim(0,50)
	plt.show()
	
	w = ((dgz < -2) | (dgz > 5) | (dgr > 4))
	print((sum(d[w]['nran']),sum(d[w]['ntar'])))
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




def sigreg_2dme(file,spar=.006,spat=.003,min=.8,max=1.2):
	#find the confidence region from the chi2 grid found in the module below
	dir = ''
	f = open(file+'.dat').readlines()
	sumt = 0
	nb1 = int((max-min)/spar)	
	nb2 = int((max-min)/spat)
	pl1 = []
	pl2 = []
	flar = []
	flat = []
	for i in range(0,nb1):
		pl1.append(0)
	for i in range(0,nb2):	
		pl2.append(0)
	for i in range(0,nb1):
		flar.append(min+i*spar+spar/2.)
	for i in range(0,nb2):
		flat.append(min+i*spat+spat/2.)
	pmax = 0	
	chimin = 1000
	corr = 0
	for i in range(0,len(f)):
		ln = f[i].split()
		if len(ln) == 3:
			chi = float(ln[2])
			if chi < chimin:
				chimin = chi
			p = exp(-.5*chi)
			a1 = float(ln[0])
			a2 = float(ln[1])
			if p > pmax:
				pmax = p
				armax = a1
				atmax = a2
			corr += a1*a2*p
			ind1 = int((a1-min)/spar)
			pl1[ind1] += p		
			ind2 = int((a2-min)/spat)
			pl2[ind2] += p
			sumt += p	
	corr = corr/sumt
	sum = 0
	sumc = 0
	zero68 = 0
	zero95 = 0
	sumca = 0
	zero68a = 0
	zero95a = 0
	thmax = 0
	thmin = 100
	ml = 0
	sig1 = .682
	sig2 = .95
	min2 = (1.-sig2)/2.
	min1 = (1.-sig1)/2.
	max1 = min1+sig1
	max2 = min2+sig2
	pmax = 0
	s = 0
	#ofn = float(f[0].split()[0])
	ofn = 0
	sumf = 0
	sumar = 0
	sumat = 0
	for i in range(0,len(pl1)):
		sumar += pl1[i]
	for i in range(0,len(pl2)):
		sumat += pl2[i]
	for i in range(0,len(flar)):
		fn = flar[i]
		od = sum#/sumt
		pb = pl1[i]/sumar
		sum += pb
		sumf += fn*pb
		if pb > pmax:
			pmax = pb
			imax = i
			fmax = fn
		d = sum#/sumt	
		if d > min2 and s == 0:
			abs1 = abs(od-min2)
			abs2 = abs(d-min2)
			fn2d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =1
		if sum > min1 and s == 1:
			abs1 = abs(od-min1)
			abs2 = abs(d-min1)
			fn1d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =2
		if sum > max1 and s == 2:
			abs1 = abs(od-max1)
			abs2 = abs(d-max1)
			fn1u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =3
		if sum > max2 and s == 3:
			abs1 = abs(od-max2)
			abs2 = abs(d-max2)
			fn2u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =4
		ofn = fn
	if imax != 0 and imax != len(flar)-1:
		pu = pl1[imax+1]/sumar
		pd = pl1[imax-1]/sumar
		dpu = pmax-pu
		dpd = pmax-pd
	
		fd = (float(flar[imax-1])+fmax)/2.
		fu = (float(flar[imax+1])+fmax)/2.	
		am = (fu/dpu+fd/dpd)/(1./dpu+1./dpd)
	else:
		am = float(flar[imax])
	a1b = sumf/sum
	err1 = (fn1u-fn1d)/2.		

	sum = 0
	sumc = 0
	zero68 = 0
	zero95 = 0
	sumca = 0
	zero68a = 0
	zero95a = 0
	thmax = 0
	thmin = 100
	ml = 0
	sig1 = .682
	sig2 = .95
	min2 = (1.-sig2)/2.
	min1 = (1.-sig1)/2.
	max1 = min1+sig1
	max2 = min2+sig2
	pmax = 0
	s = 0
	#ofn = float(f[0].split()[0])
	sumf = 0
	for i in range(0,len(flat)):
		fn = flat[i]
		od = sum#/sumt
		pb = pl2[i]/sumat
		sum += pb
		sumf += fn*pb
		if pb > pmax:
			pmax = pb
			imax = i
			fmax = fn
		d = sum#/sumt	
		if d > min2 and s == 0:
			abs1 = abs(od-min2)
			abs2 = abs(d-min2)
			fn2d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =1
		if sum > min1 and s == 1:
			abs1 = abs(od-min1)
			abs2 = abs(d-min1)
			fn1d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =2
		if sum > max1 and s == 2:
			abs1 = abs(od-max1)
			abs2 = abs(d-max1)
			fn1u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =3
		if sum > max2 and s == 3:
			abs1 = abs(od-max2)
			abs2 = abs(d-max2)
			fn2u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =4
		ofn = fn
	if imax != 0 and imax != len(flat)-1:
		pu = pl2[imax+1]/sumat
		pd = pl2[imax-1]/sumat
		dpu = pmax-pu
		dpd = pmax-pd
	
		fd = (float(flat[imax-1])+fmax)/2.
		fu = (float(flat[imax+1])+fmax)/2.	
		am = (fu/dpu+fd/dpd)/(1./dpu+1./dpd)
	else:
		am = float(flat[imax])
	a2b = sumf/sum
	err2 = (fn1u-fn1d)/2.		


	#print sumf/sum,am
	corr = corr-a1b*a2b
	return a1b,err1,a2b,err2,chimin,corr,corr/(err1*err2)

def sigreg_2dEZ(file):
	import numpy as np
	d = np.loadtxt(file).transpose()
	chi2 = d[-1]
	prob = np.exp(-0.5*chi2)
	pnorm = np.sum(prob)
	mar = np.sum(prob*d[0])/pnorm
	map = np.sum(prob*d[1])/pnorm
	sr = np.sum(prob*d[0]**2.)/pnorm
	sp = np.sum(prob*d[1]**2.)/pnorm
	sigr = sqrt(sr-mar**2.)
	sigp = sqrt(sp-map**2.)
	crp = np.sum(prob*d[0]*d[1])/pnorm-mar*map
	w = (chi2-np.min(chi2)) < 1
	indmin = np.argmin(chi2)
	print(d[0][indmin],d[1][indmin])
	print(np.max(abs(d[0][w]-d[0][indmin])),np.max(d[0][indmin]-abs(d[0][w])))
	print(np.max(abs(d[1][w]-d[1][indmin])),np.max(d[1][indmin]-abs(d[1][w])))
	return mar,sigr,map,sigp,np.min(chi2),crp,crp/(sigr*sigp)


def sigreg_2dmea(arr,spar=.006,spat=.003,min=.8,max=1.2):
	#find the confidence region from the chi2 grid found in the module below
	sumt = 0
	nb1 = int((max-min)/spar)	
	nb2 = int((max-min)/spat)
	pl1 = []
	pl2 = []
	flar = []
	flat = []
	for i in range(0,nb1):
		pl1.append(0)
	for i in range(0,nb2):	
		pl2.append(0)
	for i in range(0,nb1):
		flar.append(min+i*spar+spar/2.)
	for i in range(0,nb2):
		flat.append(min+i*spat+spat/2.)
	pmax = 0	
	chimin = 1000
	corr = 0
	arrt = arr.transpose()
	chil = arrt[2]
	chimin = np.min(chil)
	minind = np.argmin(chil)
	armax = arrt[0][minind]
	atmax = arrt[1][minind]
	pl = np.exp(-0.5*chil)
	corr = np.sum(arrt[0]*arrt[1]*pl)
	sumt = np.sum(pl)
	corr = corr/sumt
	for i in range(0,len(arr)):
		a1 = arrt[0][i]
		a2 = arrt[1][i]
		ind1 = int((a1-min)/spar)
		pl1[ind1] += pl[i]		
		ind2 = int((a2-min)/spat)
		pl2[ind2] += pl[i]
			
	sumar = np.sum(pl1)
	sumat = np.sum(pl2)
	pl1 = np.array(pl1)
	pl2 = np.array(pl2)
	flar = np.array(flar)
	flat = np.array(flat)
	a1m = np.sum(flar*pl1)/np.sum(pl1)
	a2m = np.sum(flat*pl2)/np.sum(pl2)
	s1 = sqrt(np.sum(flar**2.*pl1)/np.sum(pl1)-a1m**2.)
	s2 = sqrt(np.sum(flat**2.*pl2)/np.sum(pl2)-a2m**2.)
	return a1m,a2m,s1,s2,chimin,corr-a1m*a2m,(corr-a1m*a2m)/s1/s2
	sum = 0
	sumc = 0
	zero68 = 0
	zero95 = 0
	sumca = 0
	zero68a = 0
	zero95a = 0
	thmax = 0
	thmin = 100
	ml = 0
	sig1 = .682
	sig2 = .95
	min2 = (1.-sig2)/2.
	min1 = (1.-sig1)/2.
	max1 = min1+sig1
	max2 = min2+sig2
	pmax = 0
	s = 0
	#ofn = float(f[0].split()[0])
	ofn = 0
	sumf = 0
	sumar = 0
	sumat = 0
	for i in range(0,len(pl1)):
		sumar += pl1[i]
	for i in range(0,len(pl2)):
		sumat += pl2[i]
	for i in range(0,len(flar)):
		fn = flar[i]
		od = sum#/sumt
		pb = pl1[i]/sumar
		sum += pb
		sumf += fn*pb
		if pb > pmax:
			pmax = pb
			imax = i
			fmax = fn
		d = sum#/sumt	
		if d > min2 and s == 0:
			abs1 = abs(od-min2)
			abs2 = abs(d-min2)
			fn2d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =1
		if sum > min1 and s == 1:
			abs1 = abs(od-min1)
			abs2 = abs(d-min1)
			fn1d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =2
		if sum > max1 and s == 2:
			abs1 = abs(od-max1)
			abs2 = abs(d-max1)
			fn1u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =3
		if sum > max2 and s == 3:
			abs1 = abs(od-max2)
			abs2 = abs(d-max2)
			fn2u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =4
		ofn = fn
	if imax != 0 and imax != len(flar)-1:
		pu = pl1[imax+1]/sumar
		pd = pl1[imax-1]/sumar
		dpu = pmax-pu
		dpd = pmax-pd
	
		fd = (float(flar[imax-1])+fmax)/2.
		fu = (float(flar[imax+1])+fmax)/2.	
		am = (fu/dpu+fd/dpd)/(1./dpu+1./dpd)
	else:
		am = float(flar[imax])
	a1b = sumf/sum
	err1 = (fn1u-fn1d)/2.		

	sum = 0
	sumc = 0
	zero68 = 0
	zero95 = 0
	sumca = 0
	zero68a = 0
	zero95a = 0
	thmax = 0
	thmin = 100
	ml = 0
	sig1 = .682
	sig2 = .95
	min2 = (1.-sig2)/2.
	min1 = (1.-sig1)/2.
	max1 = min1+sig1
	max2 = min2+sig2
	pmax = 0
	s = 0
	#ofn = float(f[0].split()[0])
	sumf = 0
	for i in range(0,len(flat)):
		fn = flat[i]
		od = sum#/sumt
		pb = pl2[i]/sumat
		sum += pb
		sumf += fn*pb
		if pb > pmax:
			pmax = pb
			imax = i
			fmax = fn
		d = sum#/sumt	
		if d > min2 and s == 0:
			abs1 = abs(od-min2)
			abs2 = abs(d-min2)
			fn2d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =1
		if sum > min1 and s == 1:
			abs1 = abs(od-min1)
			abs2 = abs(d-min1)
			fn1d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =2
		if sum > max1 and s == 2:
			abs1 = abs(od-max1)
			abs2 = abs(d-max1)
			fn1u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =3
		if sum > max2 and s == 3:
			abs1 = abs(od-max2)
			abs2 = abs(d-max2)
			fn2u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =4
		ofn = fn
	if imax != 0 and imax != len(flat)-1:
		pu = pl2[imax+1]/sumat
		pd = pl2[imax-1]/sumat
		dpu = pmax-pu
		dpd = pmax-pd
	
		fd = (float(flat[imax-1])+fmax)/2.
		fu = (float(flat[imax+1])+fmax)/2.	
		am = (fu/dpu+fd/dpd)/(1./dpu+1./dpd)
	else:
		am = float(flat[imax])
	a2b = sumf/sum
	err2 = (fn1u-fn1d)/2.		


	#print sumf/sum,am
	corr = corr-a1b*a2b
	return a1b,err1,a2b,err2,chimin,corr,corr/(err1*err2)


def sigreg_2dJB(tmin,tmax,rmin,rmax):
	#find the confidence region from the chi2 grid found in the module below
	dir = ''
	file='/Users/ashleyross/Dropbox/eBOSS/eBOSS_LRGpCMASS_clustering_vDR16_COMB_rec_5mpc_shift0-rmin50.0-rmax150.0-bb-quad.at.ap.scan2d'
	f = open(file).readlines()
	sumt = 0
	flar = np.linspace(rmin,rmax)
	spar = flar[1]-flar[0]
	flat = np.linspace(tmin,tmax)
	spat = flat[1]-flat[0]
	nb1 = len(flar)	
	nb2 = len(flat)
	pl1 = []
	pl2 = []
	
	#flat = []
	for i in range(0,nb1):
		pl1.append(0)
	for i in range(0,nb2):	
		pl2.append(0)
	#for i in range(0,nb1):
	#	flar.append(min+i*spar+spar/2.)
	#for i in range(0,nb2):
	#	flat.append(min+i*spat+spat/2.)
	pmax = 0	
	chimin = 1000
	corr = 0
	for i in range(0,len(f)):
		ln = f[i].split()
		if len(ln) == 3:
			chi = float(ln[2])
			if chi < chimin:
				chimin = chi
			p = exp(-.5*chi)
			a1 = float(ln[0])
			a2 = float(ln[1])
			if p > pmax:
				pmax = p
				atmax = a1
				armax = a2
			corr += a1*a2*p
			ind1 = int((a1-tmin)/spat)
			#print(a1,ind1,tmin)
			pl1[ind1] += p		
			ind2 = int((a2-rmin)/spar)
			pl2[ind2] += p
			sumt += p	
	corr = corr/sumt
	sum = 0
	sumc = 0
	zero68 = 0
	zero95 = 0
	sumca = 0
	zero68a = 0
	zero95a = 0
	thmax = 0
	thmin = 100
	ml = 0
	sig1 = .682
	sig2 = .95
	min2 = (1.-sig2)/2.
	min1 = (1.-sig1)/2.
	max1 = min1+sig1
	max2 = min2+sig2
	pmax = 0
	s = 0
	#ofn = float(f[0].split()[0])
	ofn = 0
	sumf = 0
	sumar = 0
	sumat = 0
	for i in range(0,len(pl1)):
		sumat += pl1[i]
	for i in range(0,len(pl2)):
		sumar += pl2[i]
	for i in range(0,len(flar)):
		fn = flar[i]
		od = sum#/sumt
		pb = pl2[i]/sumar
		sum += pb
		sumf += fn*pb
		if pb > pmax:
			pmax = pb
			imax = i
			fmax = fn
		d = sum#/sumt	
		if d > min2 and s == 0:
			abs1 = abs(od-min2)
			abs2 = abs(d-min2)
			fn2d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =1
		if sum > min1 and s == 1:
			abs1 = abs(od-min1)
			abs2 = abs(d-min1)
			fn1d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =2
		if sum > max1 and s == 2:
			abs1 = abs(od-max1)
			abs2 = abs(d-max1)
			fn1u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =3
		if sum > max2 and s == 3:
			abs1 = abs(od-max2)
			abs2 = abs(d-max2)
			fn2u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =4
		ofn = fn
	if imax != 0 and imax != len(flar)-1:
		pu = pl2[imax+1]/sumar
		pd = pl2[imax-1]/sumar
		dpu = pmax-pu
		dpd = pmax-pd
	
		fd = (float(flar[imax-1])+fmax)/2.
		fu = (float(flar[imax+1])+fmax)/2.	
		am = (fu/dpu+fd/dpd)/(1./dpu+1./dpd)
	else:
		am = float(flar[imax])
	a1b = sumf/sum
	err1 = (fn1u-fn1d)/2.		

	sum = 0
	sumc = 0
	zero68 = 0
	zero95 = 0
	sumca = 0
	zero68a = 0
	zero95a = 0
	thmax = 0
	thmin = 100
	ml = 0
	sig1 = .682
	sig2 = .95
	min2 = (1.-sig2)/2.
	min1 = (1.-sig1)/2.
	max1 = min1+sig1
	max2 = min2+sig2
	pmax = 0
	s = 0
	#ofn = float(f[0].split()[0])
	sumf = 0
	for i in range(0,len(flat)):
		fn = flat[i]
		od = sum#/sumt
		pb = pl1[i]/sumat
		sum += pb
		sumf += fn*pb
		if pb > pmax:
			pmax = pb
			imax = i
			fmax = fn
		d = sum#/sumt	
		if d > min2 and s == 0:
			abs1 = abs(od-min2)
			abs2 = abs(d-min2)
			fn2d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =1
		if sum > min1 and s == 1:
			abs1 = abs(od-min1)
			abs2 = abs(d-min1)
			fn1d = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =2
		if sum > max1 and s == 2:
			abs1 = abs(od-max1)
			abs2 = abs(d-max1)
			fn1u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =3
		if sum > max2 and s == 3:
			abs1 = abs(od-max2)
			abs2 = abs(d-max2)
			fn2u = (ofn/abs1+fn/abs2)/(1./abs1+1./abs2)
			#print line, sum/sumt
			s =4
		ofn = fn
	if imax != 0 and imax != len(flat)-1:
		pu = pl1[imax+1]/sumat
		pd = pl1[imax-1]/sumat
		dpu = pmax-pu
		dpd = pmax-pd
	
		fd = (float(flat[imax-1])+fmax)/2.
		fu = (float(flat[imax+1])+fmax)/2.	
		am = (fu/dpu+fd/dpd)/(1./dpu+1./dpd)
	else:
		am = float(flat[imax])
	a2b = sumf/sum
	err2 = (fn1u-fn1d)/2.		


	#print sumf/sum,am
	corr = corr-a1b*a2b
	return a1b,err1,a2b,err2,chimin,corr,corr/(err1*err2)


def sigreg_c12(file,file2='',fac=1.,md='f'):
	#report the confidence region +/-1 for chi2
	dir = ''
	if md == 'l':
		f = file
	if md == 'f':
		f = np.loadtxt(file+'.dat').transpose()
	if file2 != '':
		f2 = np.loadtxt(file2+'.dat').transpose()
	chil = []
	chim = 1000
	
	fl = []
# 	for i in range(0,len(f)):
# 		a = float(f[i].split()[0])
# 		#if a > min and a < max:
# 		chiv = float(f[i].split()[-1])*fac
# 		if file2 != '':
# 			chiv = (chiv+float(f2[i].split()[-1])*fac)/2.
# 		chil.append((chiv,a))
# 		if chiv < chim:
# 			#better to fit a parabola to get these values
# 			chim = chiv	
# 			im = i
# 			am = a
	chill = f[1]
	if file2 != '':
		chill = (chill-min(f[1])+(f2[1]-min(f2[1])))/2.
	for i in range(0,len(chill)):
		chil.append((chill[i],f[0][i]))
		if chill[i] < chim:
			chim = chill[i]
			am = f[0][i]
			im = i
	#chim = min(chil)	
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

def sigreg_comb5(file,fac=1.,md='f'):
	#report the confidence region +/-1 for chi2
	chil = np.loadtxt(file+'0.dat').transpose()[1]
	al = np.loadtxt(file+'0.dat').transpose()[0]
	for i in range(1,5):
		chil += np.loadtxt(file+str(i)+'.dat').transpose()[1]
	chil = chil/5.
	#print(chil)
	chim = np.min(chil)
	im = np.argmin(chil)
	am = al[im]
	#print(len(chil),len(al))
	
	a1u = 2.
	a1d = 0
	a2u = 2.
	a2d = 0
	oa = 0
	ocd = 0
	s0 = 0
	s1 = 0
	for i in range(im+1,len(al)):
		chid = chil[i] - chim
		alph = al[i]
		if chid > 1. and s0 == 0:
			a1u = (alph/abs(chid-1.)+oa/abs(ocd-1.))/(1./abs(chid-1.)+1./abs(ocd-1.))
			s0 = 1
		if chid > 4. and s1 == 0:
			a2u = (alph/abs(chid-4.)+oa/abs(ocd-4.))/(1./abs(chid-4.)+1./abs(ocd-4.))
			s1 = 1
		#print(alph,chid)
		ocd = chid	
		oa = alph
	oa = 0
	ocd = 0
	s0 = 0
	s1 = 0
	for i in range(1,im):
		chid = chil[im-i] - chim
		alph = al[im-i]
		if chid > 1. and s0 == 0:
			a1d = (alph/abs(chid-1.)+oa/abs(ocd-1.))/(1./abs(chid-1.)+1./abs(ocd-1.))
			s0 = 1
		if chid > 4. and s1 == 0:
			a2d = (alph/abs(chid-4.)+oa/abs(ocd-4.))/(1./abs(chid-4.)+1./abs(ocd-4.))
			s1 = 1
		ocd = chid	
		oa = alph
	if a1u < a1d:
		a1u = 2.
		a1d = 0
	if a2u < a2d:
		a2u = 2.
		a2d = 0
	print(((a1d+a1u)/2.,(-a1d+a1u)/2.))		
	return am,a1d,a1u,a2d,a2u,chim	

def faccalc(nm,nb,nd):
	A = 2./(nm-nb-1.)/(nm-nb-4.)
	B = (nm-nb-2.)/(nm-nb-1.)/(nm-nb-4.)
	print(( A,B))
	mv = (1.+B*(nb-nd))/(1.+A+B*(nd+1.))
	md = (nm-1.)/(nm-nb-2.)*mv
	return mv,md		


def xibaoNS(sample,zmin,zmax,version='4',wm='fkp',bs=8,start=0,md='data',ns='',mom='0',covmd='QSO',rmin=35,rmax=180.,rmaxb=50.,m=1.,mb='',Bp=0.4,v='n',mockn='',angfac=0,damp='6.0',tempmd='',Nmock=1000,template='Challenge_matterpower',rec='',covv='',mumin=0,mumax=1,mupow=0,covfac=1.,npar=3):
	#does baofits, set mb='nobao' to do no BAO fit
	from baofit_pubtest import doxi_isolike
	from Cosmo import distance
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	if mupow != 0:
		muw += 'mup'+str(mupow)		

	outdir = ebossdir
	indir = ebossdir
	bsst = str(bs)+'st'+str(start)
	print(bsst)
	#if sample == 'ELG' or sample == 'LRGpCMASS':
	zw = ''
	if md == 'data':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		
		if sample == 'ELG' and zmin != 0.6:
			zw += 'zmin'+str(zmin)

		recf =''
		if rec == 'rec' or rec == '_rec':
			recf = '_rec'
		if version == '7':
			#dn = np.loadtxt(ebossdir+'2PCF_eBOSS_ELG_clustering_NGC_v7'+rec+'_z0.6z1.1.dat').transpose()
			#ds = np.loadtxt(ebossdir+'2PCF_eBOSS_ELG_clustering_SGC_v7'+rec+'_z0.6z1.1.dat').transpose()
			dn = np.loadtxt(ebossdir+'xi024ELG_v7NGC'+zw+'CZdata'+rec+muw+bsst+'.dat').transpose()
			ds = np.loadtxt(ebossdir+'xi024ELG_v7SGC'+zw+'CZdata'+rec+muw+bsst+'.dat').transpose()

		else:	
			dn = np.loadtxt(ebossdir+'xi024geboss'+sample+'_NGC'+version+recf+'_'+wz+wm+bsst+'.dat').transpose()
			ds = np.loadtxt(ebossdir+'xi024geboss'+sample+'_SGC'+version+recf+'_'+wz+wm+bsst+'.dat').transpose()
		print('number of rows in file:\n')
		print((len(dn[0])))
		ezver = '7'

	if md == 'dataJH':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		#d0N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[1][:maxind]
		#d2N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[2][:maxind]
		#d0S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[1][:maxind]
		#d2S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[2][:maxind]
		
		dn = np.loadtxt(dir+'dr16_qso_v'+version+'_N_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()
		ds = np.loadtxt(dir+'dr16_qso_v'+version+'_S_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()
		ezver = '7'


	if sample == 'ELGQPM':
		dn = np.loadtxt(ebossdir+'ELGmockxi_MV/qpm_mock_anymask_ELG_recon_specweights_NGC_'+mockn+'.mul').transpose()
		ds = np.loadtxt(ebossdir+'ELGmockxi_MV/qpm_mock_anymask_ELG_recon_specweights_SGC_'+mockn+'.mul').transpose()

	if sample == 'ELGEZ':
		dn = np.loadtxt(dirsci+'/EZmockELGv4/xi/xi024ELGNGCEZmock'+rec+mockn+muw+bsst+'.dat').transpose()
		ds = np.loadtxt(dirsci+'/EZmockELGv4/xi/xi024ELGSGCEZmock'+rec+mockn+muw+bsst+'.dat').transpose()
		indir = ''
		outdir = dirsci+'EZmockELGv4/BAOfits/'



	if md == 'EZmockave':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		ezver = version
		#dir = '/global/cscratch1/sd/ajross/ebossxi/'+ sample+'_v'+str(ezver)+'/'
		#outdir = '/global/cscratch1/sd/ajross/'
		#dirH = ''
		#xiave02NGCQSO_v5_EZ5st0.dat
		if sample == 'ELG' and zmin != 0.6:
			zw += 'zmin'+str(zmin)
		

		dn = np.loadtxt(dirH+'xiave0NGC'+zw+sample+ns+'_v'+str(ezver)+'_EZ'+rec+muw+bsst+'.dat').transpose()#[1][:maxind]
		ds = np.loadtxt(dirH+'xiave0SGC'+zw+sample+ns+'_v'+str(ezver)+'_EZ'+rec+muw+bsst+'.dat').transpose()#[1][:maxind]

	if md == 'ORmockave':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		ezver = version
		#dir = '/global/cscratch1/sd/ajross/ebossxi/'+ sample+'_v'+str(ezver)+'/'
		#outdir = '/global/cscratch1/sd/ajross/'
		#dirH = ''
		#xiave02NGCQSO_v5_EZ5st0.dat
		

		dn = np.loadtxt(dirH+'xiave024'+sample+'ORmockNGC'+rec+muw+bsst+'.dat').transpose()#[1][:maxind]
		ds = np.loadtxt(dirH+'xiave024'+sample+'ORmockSGC'+rec+muw+bsst+'.dat').transpose()#[1][:maxind]

	if md == 'EZmock':
		#dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		#dirout = '/Users/ashleyross/Dropbox/eboss/' 
		#dirH = '/Users/ashleyross/Dropbox/eboss/'
		ezver = version
		dir = '/global/cscratch1/sd/ajross/ebossxi/'+ sample+'_v'+str(ezver)+'/'
		outdir = '/global/cscratch1/sd/ajross/baofits/'
		dirH = ''
		#xiave02NGCQSO_v5_EZ5st0.dat
		dn = np.loadtxt(dir+'xi024'+sample+'_v'+str(ezver)+'NGCEZmock'+rec+str(mockn)+bsst+'.dat').transpose()
		ds = np.loadtxt(dir+'xi024'+sample+'_v'+str(ezver)+'SGCEZmock'+rec+str(mockn)+bsst+'.dat').transpose()

	if md == 'EZmockjh':
		#dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		#dirout = '/Users/ashleyross/Dropbox/eboss/' 
		#dirH = '/Users/ashleyross/Dropbox/eboss/'
		ezver = version
		#dir = '/global/cscratch1/sd/ajross/ebossxi/'+ sample+'_v'+str(ezver)+'/'
		dir = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys/'
		dirs = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys_obsolete/'
		outdir = '/global/cscratch1/sd/ajross/baofits/'
		dirH = ''
		zer = ''
		if mockn < 10:
			zer += '0'
		if mockn < 100:
			zer += '0'
		if mockn < 1000:
			zer += 	'0'
		mockstr = zer+str(mockn)	
		print((mockn,mockstr))
		mockn = str(mockn)
		dn = np.loadtxt(dir+'ez_qsov7_0.31_ngc_xi4m_wsys_ds5_'+mockstr+'.dat').transpose()
		ds = np.loadtxt(dirs+'ez_qsov7_0.31_sgc_xi4m_wsys_ds5_'+mockstr+'.dat').transpose()


	if md == 'mockave':
		try:
			indir = ''
			outdir = ''
			dn = np.loadtxt(indir+'xiave0NGC'+sample+'_EZ'+rec+'angfac'+str(angfac)+bsst+'.dat').transpose()
			ds = np.loadtxt(indir+'xiave0SGC'+sample+'_EZ'+rec+'angfac'+str(angfac)+bsst+'.dat').transpose()
		except:	
			outdir = ebossdir
			indir = ebossdir
			dn = np.loadtxt(indir+'xiave0NGC'+sample+'_EZ'+rec+'angfac'+str(angfac)+bsst+'.dat').transpose()
			ds = np.loadtxt(indir+'xiave0SGC'+sample+'_EZ'+rec+'angfac'+str(angfac)+bsst+'.dat').transpose()
		#if rec == 'rec':
		#	dn = np.loadtxt(ebossdir+'xiave_recon0NGCELG_EZ'+bsst+'.dat').transpose()
		#	ds = np.loadtxt(ebossdir+'xiave_recon0SGCELG_EZ'+bsst+'.dat').transpose()
	if md == 'qsozel':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		ezver = version

		from scipy.interpolate import interp1d
		nb = 200//bs
		dn = np.zeros((2,nb))
		zelf = np.loadtxt('/Users/ashleyross/Dropbox/eboss/2PCF_QSO_Zeldovich.dat').transpose()
		fz = interp1d(zelf[0], zelf[1], kind='cubic')
		for i in range(0,nb):
			r = bs/2.+bs*i
			xi = fz(r)
			dn[0][i] = r
			dn[1][i] = xi
		ds = dn
	rl = dn[0]
	#print( rl)
	if mumax == 1 and mumin == 0 and mupow == 0:

		mod = np.loadtxt('BAOtemplates/xi0'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
		modsmooth = np.loadtxt('BAOtemplates/xi0sm'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
	else:
		mod = np.loadtxt('BAOtemplates/xi'+muw+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
		modsmooth = np.loadtxt('BAOtemplates/xism'+muw+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
		

	if mb == 'nobao':		
		mod = modsmooth

			
	fn = 1.
	fs = 1.
	if sample == 'QSO':
		if version == 'v1.8':
			fn = 1./1.192
			fs = 822/857.
		if version == '4':
			fn = .89/1.82
			fs = .594/0.966	
		if version == '5' or version == 'test' or version == '5_1':	
			fn = .89/2.2
			fs = .594/1.24	
		fn = 1.
		fs = 1.
	print((fn,fs))	
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
	if covmd == 'QSO':
		covN = np.loadtxt(indir+'cov0EZmockv1.8ngc'+bsst+'.dat')*fn
		covS = np.loadtxt(indir+'cov0EZmockv1.8sgc'+bsst+'.dat')*fs
		#covN = np.loadtxt(dirH+'cov0NGCQSO_EZ5st0.dat')
		#covS= np.loadtxt(dirH+'cov0SGCQSO_EZ5st0.dat')
		
	if covmd == 'ELG':# and rec == '':
		covN = np.loadtxt(indir+'cov0NGC'+sample+'_EZ'+rec+'angfac'+str(angfac)+covv+bsst+'.dat')#*fn
		covS = np.loadtxt(indir+'cov0SGC'+sample+'_EZ'+rec+'angfac'+str(angfac)+covv+bsst+'.dat')#*fs

	if covmd == 'me':
		#covN = np.loadtxt(dirH+'cov0NGC'+zw+sample+'_v'+str(ezver)+'_EZ'+rec+muw+bsst+'.dat')
		#print(('open NGC cov '+dirH+'cov0NGC'+zw+sample+'_v'+str(ezver)+'_EZ'+rec+muw+bsst+'.dat'))
		#covS = np.loadtxt(dirH+'cov0SGC'+zw+sample+'_v'+str(ezver)+'_EZ'+rec+muw+bsst+'.dat')  
		covN = np.loadtxt(dirH+'cov0NGC'+sample+'_v'+str(ezver)+'_EZ'+rec+bsst+'.dat')
		covS = np.loadtxt(dirH+'cov0SGC'+sample+'_v'+str(ezver)+'_EZ'+rec+bsst+'.dat')  

	covS *= covfac
	covN *= covfac				
	covti = np.linalg.pinv(covN)+np.linalg.pinv(covS)
	covt = np.linalg.pinv(covti)
	
	dns = np.zeros(len(covN))
	if len(dn[1]) < len(dns):
		dns = np.zeros(len(dn[1]))
	#if len(dns) > len(dn[1]):
	#	dns = np.zeros(len(covN)//2)
	for i in range(0,len(dns)):
		x = (dn[1][i]/covN[i][i]+ds[1][i]/covS[i][i])/(1./covN[i][i]+1./covS[i][i])
		dns[i] = x
	#dns =np.dot(covt,  (np.dot(dn[1],np.linalg.pinv(covN))+np.dot(ds[1],np.linalg.pinv(covS))))
	chiln = doxi_isolike(dn[1],covN,mod,modsmooth,rl,bs=bs,rmin=rmin,rmax=rmax,rmaxb=rmaxb,v=v,wo=sample+'NGC'+version+rec+bsst+mb+'np'+str(npar),diro=outdir,Bp=Bp,Nmock=Nmock,npar=npar)

	wf = sample+md+mockn+version+rec+damp+mb+muw+bsst
	fo = open(outdir+'BAOxichilNGC'+wf+'.dat','w')
	for i in range(0,len(chiln)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chiln[i])+'\n')
	fo.close()
	an = sigreg_c12(outdir+'BAOxichilNGC'+wf)
	#print an
	print(((an[1]+an[2])/2.,(an[2]-an[1])/2.,min(chiln)))
	chils = doxi_isolike(ds[1],covS,mod,modsmooth,rl,bs=bs,rmin=rmin,rmax=rmax,rmaxb=rmaxb,v=v,wo=sample+'SGC'+version+rec+bsst+mb+'np'+str(npar),diro=outdir,Bp=Bp,Nmock=Nmock,npar=npar)
	fo = open(outdir+'BAOxichilSGC'+wf+'.dat','w')
	for i in range(0,len(chils)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chils[i])+'\n')
	fo.close()
	als = sigreg_c12(outdir+'BAOxichilSGC'+wf)
	#print als
	print(((als[1]+als[2])/2.,(als[2]-als[1])/2.,min(chils)))
	chilt = np.array(chiln)+np.array(chils)
	fo = open(outdir+'BAOxichilNScomb'+wf+'.dat','w')
	for i in range(0,len(chilt)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chilt[i])+'\n')
	fo.close()
	a = sigreg_c12(outdir+'BAOxichilNScomb'+wf)
	
	#print a
	print(((a[1]+a[2])/2.,(a[2]-a[1])/2.,a[-1]))
	chilns = doxi_isolike(dns,covt,mod,modsmooth,rl,bs=bs,rmin=rmin,rmax=rmax,rmaxb=rmaxb,v=v,wo=sample+'NScombf'+version+rec+bsst+mb+'np'+str(npar),diro=outdir,Bp=Bp,Nmock=Nmock,npar=npar)
	fo = open(outdir+'BAOxichilNScombf'+wf+'.dat','w')
	for i in range(0,len(chils)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chilns[i])+'\n')
	fo.close()
	alns = sigreg_c12(outdir+'BAOxichilNScombf'+wf)
	#print als
	print(((alns[1]+alns[2])/2.,(alns[2]-alns[1])/2.,min(chilns)))
	print(bsst)
	return True

def xibaoNSmu(sample,mu,zmin=0.6,zmax=1.1,version='4',wm='fkp',bs=8,start=0,md='data',mom='0',covmd='QSO',rmin=35,rmax=180.,rmaxb=50.,m=1.,mb='',Bp=0.4,v='n',mockn='',angfac=0,damp='6.0',tempmd='',Nmock=1000,template='Challenge_matterpower',rec='',covv='',mumin=0,mumax=1):
	#does baofits, set mb='nobao' to do no BAO fit
	from baofit_pubtest import doxi_isolike
	from Cosmo import distance
	wz = 'mz'+str(zmin)+'xz'+str(zmax)
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	outdir = ebossdir
	indir = ebossdir
	bsst = str(bs)+'st'+str(start)
	print(bsst)
	#if sample == 'ELG' or sample == 'LRGpCMASS':
	if md == 'data':
		recf =''
		if rec == 'rec':
			recf = '_rec'
		dn = np.loadtxt(ebossdir+'xi024geboss'+sample+'_NGC'+version+recf+'_'+wz+wm+bsst+'.dat').transpose()
		xin = dn[1]+P2(mu)*dn[2]+P4(mu)*dn[3]
		ds = np.loadtxt(ebossdir+'xi024geboss'+sample+'_SGC'+version+recf+'_'+wz+wm+bsst+'.dat').transpose()
		xis = ds[1]+P2(mu)*ds[2]+P4(mu)*ds[3]
	if sample == 'ELGEZ':
		dn = np.loadtxt(dirsci+'/EZmockELGv4/xi/xi024ELGNGCEZmock'+rec+mockn+muw+bsst+'.dat').transpose()
		xin = dn[1]+P2(mu)*dn[2]+P4(mu)*dn[3]
		ds = np.loadtxt(dirsci+'/EZmockELGv4/xi/xi024ELGSGCEZmock'+rec+mockn+muw+bsst+'.dat').transpose()
		xis = ds[1]+P2(mu)*ds[2]+P4(mu)*ds[3]
		indir = ''
		outdir = dirsci+'EZmockELGv4/BAOfits/'

	if md == 'mockave':
		#try:
		indir = ''
		outdir = ''
		dn = np.loadtxt(indir+'xiavemu'+str(mu)+'NGC'+sample+'_EZ'+rec+bsst+'.dat').transpose()
		xin = dn[1]
		ds = np.loadtxt(indir+'xiavemu'+str(mu)+'SGC'+sample+'_EZ'+rec+bsst+'.dat').transpose()
		xis = ds[1]
		#	ds = np.loadtxt(ebossdir+'xiave_recon0SGCELG_EZ'+bsst+'.dat').transpose()
	#rl = dn[0]
	##print rl

	mod0 = np.loadtxt('BAOtemplates/xi0'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
	mod2 = np.loadtxt('BAOtemplates/xi2'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
	mod4 = np.loadtxt('BAOtemplates/xi4'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
	mod = mod0 + P2(mu)*mod2 + P4(mu)*mod4

	modsmooth0 = np.loadtxt('BAOtemplates/xi0sm'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
	modsmooth2 = np.loadtxt('BAOtemplates/xi2sm'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
	modsmooth4 = np.loadtxt('BAOtemplates/xi4sm'+tempmd+'Challenge_matterpower'+damp+'.dat').transpose()[1]
	modsmooth = modsmooth0 + P2(mu)*modsmooth2 + P4(mu)*modsmooth4

	if mb == 'nobao':		
		mod = modsmooth

			
	fn = 1.
	fs = 1.
	rl = dn[0]
	if sample == 'QSO':
		if version == 'v1.8':
			fn = 1./1.192
			fs = 822/857.
		if version == '4':
			fn = .89/1.82
			fs = .594/0.966	
		if version == '5' or version == 'test' or version == '5_1':	
			fn = .89/2.2
			fs = .594/1.24	

	covN = np.loadtxt(indir+'covmu'+str(mu)+'NGC'+sample+'_EZ'+rec+covv+bsst+'.dat')#*fn
	covS = np.loadtxt(indir+'covmu'+str(mu)+'SGC'+sample+'_EZ'+rec+covv+bsst+'.dat')#*fs
					
	covti = np.linalg.pinv(covN)+np.linalg.pinv(covS)
	covt = np.linalg.pinv(covti)
	dns = np.zeros(len(covN))
	for i in range(0,len(dns)):
		x = (xin[i]/covN[i][i]+xis[i]/covS[i][i])/(1./covN[i][i]+1./covS[i][i])
		dns[i] = x
	chiln = doxi_isolike(xin,covN,mod,modsmooth,rl,rmin=rmin,rmax=rmax,rmaxb=rmaxb,v=v,wo=sample+'NGC'+version+rec+bsst+mb,diro=outdir,Bp=Bp,Nmock=Nmock)

	wf = sample+mockn+version+rec+damp+mb+'mu'+str(mu)+bsst
	fo = open(outdir+'BAOxichilNGC'+wf+'.dat','w')
	for i in range(0,len(chiln)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chiln[i])+'\n')
	fo.close()
	an = sigreg_c12(outdir+'BAOxichilNGC'+wf)
	#print an
	print(((an[1]+an[2])/2.,(an[2]-an[1])/2.,min(chiln)))
	chils = doxi_isolike(xis,covS,mod,modsmooth,rl,rmin=rmin,rmax=rmax,rmaxb=rmaxb,v=v,wo=sample+'SGC'+version+rec+bsst+mb,diro=outdir,Bp=Bp,Nmock=Nmock)
	fo = open(outdir+'BAOxichilSGC'+wf+'.dat','w')
	for i in range(0,len(chils)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chils[i])+'\n')
	fo.close()
	als = sigreg_c12(outdir+'BAOxichilSGC'+wf)
	#print als
	print(((als[1]+als[2])/2.,(als[2]-als[1])/2.,min(chils)))
	chilt = np.array(chiln)+np.array(chils)
	fo = open(outdir+'BAOxichilNScomb'+wf+'.dat','w')
	for i in range(0,len(chilt)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chilt[i])+'\n')
	fo.close()
	a = sigreg_c12(outdir+'BAOxichilNScomb'+wf)
	
	#print a
	print(((a[1]+a[2])/2.,(a[2]-a[1])/2.,a[-1]))
	chilns = doxi_isolike(dns,covt,mod,modsmooth,rl,rmin=rmin,rmax=rmax,rmaxb=rmaxb,v=v,wo=sample+'NScombf'+version+rec+bsst+mb,diro=outdir,Bp=Bp,Nmock=Nmock)
	fo = open(outdir+'BAOxichilNScombf'+wf+'.dat','w')
	for i in range(0,len(chils)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chilns[i])+'\n')
	fo.close()
	alns = sigreg_c12(outdir+'BAOxichilNScombf'+wf)
	#print als
	print(((alns[1]+alns[2])/2.,(alns[2]-alns[1])/2.,min(chilns)))
	print(bsst)
	return True


def xi2DBAOOR(md = 'mockORjh',samp='QSO',min=50,max=150,maxb=80.,bs=8,minz=.8,maxz=2.2,covf=1.,mockn=0,orver='mocks_blind_v1',orcos='mockb0',damp='',maxr=200,spat=0.003,spar=0.006,mina=.8,maxa=1.2,covmd='ORjh',Bp=0.4,Bt=0.4):
	from baofit_pub2D import Xism_arat_1C_an
	fout = samp+md+orver+orcos+covmd+'Bp'+str(Bp)+damp+'mn'+str(mockn)
	if md == 'mockORjh':
		covdir = '/global/project/projectdirs/eboss/octobers/cm_or/'+orver+'/'+orcos
		dir = '/global/project/projectdirs/eboss/octobers/twopcf/outerrim/blind/'+orver+'/'+orcos
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		zer = ''
		if mockn < 10:
			zer += '0'
		if mockn < 100:
			zer += '0'
		if mockn < 1000:
			zer += 	'0'
		mockstr = zer+str(mockn)	
		drl = np.loadtxt(dir+'/or_qso_0.26_xi4m_ds8_'+mockstr+'_mean.dat').transpose()[0]
		d0 = np.loadtxt(dir+'/or_qso_0.26_xi4m_ds8_'+mockstr+'_mean.dat').transpose()[1]
		d2 = np.loadtxt(dir+'/or_qso_0.26_xi4m_ds8_'+mockstr+'_mean.dat').transpose()[2]

	if md == 'mockORjh_notblind':
		covdir = '/global/project/projectdirs/eboss/octobers/cm_or/non-blind/mocks_v4/'+orver+'/'+orcos
		dir = '/global/project/projectdirs/eboss/octobers/twopcf/outerrim/non-blind/mocks_v4/'+orver+'/'+orcos
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		zer = ''
		if mockn < 10:
			zer += '0'
		if mockn < 100:
			zer += '0'
		if mockn < 1000:
			zer += 	'0'
		mockstr = zer+str(mockn)	
		drl = np.loadtxt(dir+'/or_qso_0.26_xi4m_ds8_'+mockstr+'_mean.dat').transpose()[0]
		d0 = np.loadtxt(dir+'/or_qso_0.26_xi4m_ds8_'+mockstr+'_mean.dat').transpose()[1]
		d2 = np.loadtxt(dir+'/or_qso_0.26_xi4m_ds8_'+mockstr+'_mean.dat').transpose()[2]

	dv = [] #empty list to become data vector
	dvb = [] #empty list to become data vector for setting bias prior

	rl = [] #empty list to become list of r values to evaluate model at	
	rlb  = [] #empty list to become list of r values to evaluate model at for bias prior
	mini = 0
	for i in range(0,len(d0)):
		r = drl[i]
		if r > min and r < max:
			dv.append(d0[i])
			rbc = r#.75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.) #correct for pairs should have slightly larger average pair distance than the bin center
			rl.append(rbc) 
			if mini == 0:
				mini = i #minimum index number to be used for covariance matrix index assignment
			if r < maxb:
				dvb.append(d0[i])
				rlb.append(rbc)
	for i in range(0,len(d2)):
		r = drl[i]
		if r > min and r < max:
			dv.append(d2[i])
			rbc = r#.75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.)
			rl.append(rbc)
			if r < maxb:
				dvb.append(d2[i])
				rlb.append(rbc)

	dv = np.array(dv)
	#print( len(dv))
	if covmd == 'ORjh':
		ca = np.loadtxt(covdir+'/or_qso_0.26_cov_gaussian_1_ds8.txt')[:40][:40]
		
		
		if len(ca) != len(d0)*2:
			print( 'MISMATCHED data and cov matrix!')
			print((len(ca),len(d0)))
		covm = np.zeros((len(dv),len(dv))) #will become covariance matrix to be used with data vector

		#need to cut it to correct size
		for i in range(0,len(ca)):
			if i < len(ca)/2:
				ri = drl[i]
				indi = i-mini
			else:
				#ri = (i-len(ca)/2)*bs+bs/2.+binc
				ri = drl[i-len(ca)//2]
				indi = len(dv)//2+i-mini-len(ca)//2	
			for j in range(0,len(ca)):		
				if j < len(ca)/2:
					#rj = j*bs+bs/2.+binc
					rj = drl[j]
					indj = j-mini
				else:
					#rj = (j-len(ca)/2)*bs+bs/2.+binc
					rj = drl[j-len(ca)//2]
					indj = len(dv)//2+j-mini-len(ca)//2
				if ri > min and ri < max and rj > min and rj < max:
					#print ri,rj,i,j,indi,indj
					indi = int(indi)
					indj = int(indj)
					covm[indi][indj] = ca[i][j]
		invc = np.linalg.pinv(covm) #the inverse covariance matrix to pass to the code
		covmb = np.zeros((len(dvb),len(dvb)))
		for i in range(0,len(ca)):
			if i < len(ca)/2:
				ri = drl[i]
				indi = i-mini
			else:
				ri = drl[i-len(ca)//2]
				indi = len(dvb)/2+i-mini-len(ca)/2	
			for j in range(0,len(ca)):
				if j < len(ca)/2:
					rj = drl[j-len(ca)//2]
					indj = j-mini
				else:
					rj = drl[j-len(ca)//2]
					indj = len(dvb)/2+j-mini-len(ca)/2
				if ri > min and ri < maxb and rj > min and rj < maxb:	
					#print(ri,rj,indi,indj)
					indi = int(indi)
					indj = int(indj)
					covmb[indi][indj] = ca[i][j]

		invcb = np.linalg.pinv(covmb)

	mod = 'pk_or'+damp+'.dat'		
	anst = Xism_arat_1C_an(dv,invc*covf,rl,mod,dvb,invcb*covf,rlb,bs=bs,dirout=dirout,fout=fout,spat=spat,spar=spar,mina=mina,maxa=maxa,Bp=Bp,Bt=Bt)
	print(anst)
	d = np.loadtxt(dirout+'2Dbaofits/arat'+fout+'1covchi.dat')
	



def xi2DBAONS(md='data',samp='LRGpCMASS',min=50.,max=150.,maxb=80.,bs=8,binc=0,minz=0.6,maxz=1.,covf=1.,ezver='7',ver='7_1',af='',wm='fkp',rec='',mockn=0,damp='',maxr=200,spat=0.003,spar=0.006,mina=.8,maxa=1.2,covm='me',Bp=0.4,Bt=0.4):
	from baofit_pub2D import Xism_arat_1C_an
	maxind = int(maxr/bs)
	print(maxind)
	bs = int(bs)
	if md == 'data':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		d0N = np.loadtxt(dir+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]

	if md == 'dataJH':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		#d0N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[1][:maxind]
		#d2N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[2][:maxind]
		#d0S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[1][:maxind]
		#d2S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[2][:maxind]
		#if ver == '7_1:'
		d0N = np.loadtxt(dir+'dr16_qso_v'+ver+'_N_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'dr16_qso_v'+ver+'_N_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[3][:maxind]
		d0S = np.loadtxt(dir+'dr16_qso_v'+ver+'_S_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'dr16_qso_v'+ver+'_S_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[3][:maxind]
		#if ver == '7_2':	
		#	d0N = np.loadtxt(dir+'dr16_qso_v7_1_N_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[1][:maxind]
		#	d2N = np.loadtxt(dir+'dr16_qso_v7_1_N_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[3][:maxind]
		#	d0S = np.loadtxt(dir+'dr16_qso_v7_1_S_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[1][:maxind]
		#	d2S = np.loadtxt(dir+'dr16_qso_v7_1_S_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()[3][:maxind]

	if md == 'dataN':
		#for on NERSC
		dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ezver)+'/'
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		d0N = np.loadtxt(dirH+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dirH+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dirH+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dirH+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
	if md == 'dataJB':
		dir = '/Users/ashleyross/DR16plots2/LRG/data/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		
		d0N = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[2][:maxind]

	if md == 'mockave':
		#dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		#dirout = '/Users/ashleyross/Dropbox/eboss/' 
		#dirH = '/Users/ashleyross/Dropbox/eboss/'
		dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ezver)+'/'
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		#xiave02NGCQSO_v5_EZ5st0.dat

		d0N = np.loadtxt(dirH+'xiave02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dirH+'xiave02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][maxind:]
		d0S = np.loadtxt(dirH+'xiave02SGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dirH+'xiave02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][maxind:]

	if md == 'mockaveJB':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		d0N = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[2][:maxind]

	if md == 'mock':
		dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ezver)+'/'
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		d0N = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'NGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'NGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'SGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'SGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]

	if md == 'mockjh':
		dir = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys/'
		dirs = '/global/project/projectdirs/eboss/octobers/twopcf/ezmock/qso_v7_sys_obsolete/'
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		zer = ''
		if mockn < 10:
			zer += '0'
		if mockn < 100:
			zer += '0'
		if mockn < 1000:
			zer += 	'0'
		mockstr = zer+str(mockn)	
		d0N = np.loadtxt(dir+'ez_qsov7_0.31_ngc_xi4m_wsys_ds5_'+mockstr+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'ez_qsov7_0.31_ngc_xi4m_wsys_ds5_'+mockstr+'.dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dirs+'ez_qsov7_0.31_sgc_xi4m_wsys_ds5_'+mockstr+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dirs+'ez_qsov7_0.31_sgc_xi4m_wsys_ds5_'+mockstr+'.dat').transpose()[2][:maxind]



	if covm == 'me':
		cN = np.loadtxt(dirH+'cov02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat')
		cS = np.loadtxt(dirH+'cov02SGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat')  
		
		if len(cN) != len(d0N)*2:
			print( 'MISMATCHED data and cov matrix!')
			print((len(cN),len(d0N)))
	dvN = [] #empty list to become data vector
	dvbN = [] #empty list to become data vector for setting bias prior
	dvS = [] #empty list to become data vector
	dvbS = [] #empty list to become data vector for setting bias prior

	rl = [] #empty list to become list of r values to evaluate model at	
	rlb  = [] #empty list to become list of r values to evaluate model at for bias prior
	mini = 0
	for i in range(0,len(d0N)):
		r = i*bs+bs/2.+binc
		if r > min and r < max:
			dvN.append(d0N[i])
			dvS.append(d0S[i])
			rbc = r#.75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.) #correct for pairs should have slightly larger average pair distance than the bin center
			rl.append(rbc) 
			if mini == 0:
				mini = i #minimum index number to be used for covariance matrix index assignment
			if r < maxb:
				dvbN.append(d0N[i])
				dvbS.append(d0S[i])
				rlb.append(rbc)
	for i in range(0,len(d2N)):
		r = i*bs+bs/2.+binc
		if r > min and r < max:
			dvN.append(d2N[i])
			dvS.append(d2S[i])
			rbc = r#.75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.)
			rl.append(rbc)
			if r < maxb:
				dvbN.append(d2N[i])
				dvbS.append(d2S[i])
				rlb.append(rbc)

	dvN = np.array(dvN)
	dvS = np.array(dvS)
	print(( len(dvN),len(dvS)))
	if covm == 'me':
		covmS = np.zeros((len(dvS),len(dvS))) #will become covariance matrix to be used with data vector
		covmN = np.zeros((len(dvS),len(dvS))) #will become covariance matrix to be used with data vector

		#need to cut it to correct size
		for i in range(0,len(cN)):
			if i < len(cN)/2:
				ri = i*bs+bs/2.+binc
				indi = i-mini
			else:
				ri = (i-len(cN)/2)*bs+bs/2.+binc
				indi = len(dvN)/2+i-mini-len(cN)/2	
			for j in range(0,len(cN)):		
				if j < len(cN)/2:
					rj = j*bs+bs/2.+binc
					indj = j-mini
				else:
					rj = (j-len(cN)/2)*bs+bs/2.+binc
					indj = len(dvN)/2+j-mini-len(cN)/2
				if ri > min and ri < max and rj > min and rj < max:
					#print ri,rj,i,j,indi,indj
					indi = int(indi)
					indj = int(indj)
					covmN[indi][indj] = cN[i][j]
					covmS[indi][indj] = cS[i][j]
		invcN = np.linalg.pinv(covmN) #the inverse covariance matrix to pass to the code
		invcS = np.linalg.pinv(covmS)
		covmbN = np.zeros((len(dvbN),len(dvbN)))
		covmbS = np.zeros((len(dvbN),len(dvbN)))
# 		for i in range(0,len(dvbN)):
# 			if i < len(dvbN)/2:
# 				indi = i
# 			else:
# 				indi = i-len(dvbN)/2+len(covmN)/2
# 			for j in range(0,len(dvbN)):
# 				if j < len(dvbN)/2:
# 					indj = j
# 				else:
# 					indj = j-len(dvbN)/2+len(covmN)/2
# 				indi = int(indi)
# 				indj = int(indj)
# 				covmbN[i][j] = covmN[indi][indj]
# 				covmbS[i][j] = covmS[indi][indj]
		for i in range(0,len(cN)):
			if i < len(cN)/2:
				ri = i*bs+bs/2.+binc
				indi = i-mini
			else:
				ri = (i-len(cN)/2)*bs+bs/2.+binc
				indi = len(dvbN)/2+i-mini-len(cN)/2	
			for j in range(0,len(cN)):
				if j < len(cN)/2:
					rj = j*bs+bs/2.+binc
					indj = j-mini
				else:
					rj = (j-len(cN)/2)*bs+bs/2.+binc
					indj = len(dvbN)/2+j-mini-len(cN)/2
				if ri > min and ri < maxb and rj > min and rj < maxb:	
					#print(ri,rj,indi,indj)
					indi = int(indi)
					indj = int(indj)
					covmbN[indi][indj] = cN[i][j]
					covmbS[indi][indj] = cS[i][j]

		invcbN = np.linalg.pinv(covmbN)
		invcbS = np.linalg.pinv(covmbS)
		invct = invcN+invcS
		invcbt = invcbN+invcbS
	if covm == 'JB':
		invct = mkinvcJB(min=min,max=max)
		if (maxb-min)/bs < 1:
			return 'ERROR, maxb needs to be larger to have any data to fit'
		invcbt = mkinvcJB(min=min,max=maxb)
	#if rec == 'rec':
	#	mod = 'Challenge_matterpower'+.dat' #BAO template used	post-recon
	#if rec == '':
	mod = 'Challenge_matterpower'+damp+'.dat'		
	fout = samp+rec+md+ver+str(bs)+covm+'Bp'+str(Bp)+damp
	if md == 'mock' or md == 'mockjh':
		fout += 'mn'+str(mockn)

	print((len(dvN),len(invct),len(dvbN),len(invcbt)))
	#ansN = Xism_arat_1C_an(dvN,invcN,rl,mod,dvbN,invcbN,rlb,bs=bs,fout=fout+'NGC')
	#print(ansN)
	#sigreg_2dme(
	#ansS = Xism_arat_1C_an(dvS,invcS,rl,mod,dvbS,invcbS,rlb,bs=bs,fout=fout+'SGC')
	#print(ansS)
	if covm == 'me':
		s2N = []
		s2S = []
		for i in range(0,len(dvN)):
			s2N.append(covmN[i][i])
			s2S.append(covmS[i][i])
		s2N = np.array(s2N)
		s2S = np.array(s2S)
		dvt = (dvN/s2N+dvS/s2S)/(1./s2N+1./s2S)
		dvbt = (dvbN/s2N[:len(dvbN)]+dvbS/s2S[:len(dvbN)])/(1./s2N[:len(dvbN)]+1./s2S[:len(dvbN)])	
	if covm == 'JB':
		#weighting by relative sums of fkp weights
		fsumS = 0.557
		fsumN = 1.21
		dvN = np.array(dvN)
		dvS = np.array(dvS)
		dvbN = np.array(dvbN)
		dvbS = np.array(dvbS)
		dvt = (dvN*fsumN+dvS*fsumS)/(fsumS+fsumN)
		dvbt = (dvbN*fsumN+dvbS*fsumS)/(fsumS+fsumN)
	anst = Xism_arat_1C_an(dvt,invct*covf,rl,mod,dvbt,invcbt*covf,rlb,bs=bs,dirout=dirout,fout=fout+'NScomb',spat=spat,spar=spar,mina=mina,maxa=maxa,Bp=Bp,Bt=Bt)
	print(anst)
	d = np.loadtxt(dirout+'2Dbaofits/arat'+fout+'NScomb1covchi.dat')
	#chi2f = np.loadtxt(dir+'2Dbaofits/arat'+fout+'NScomb1covchi.dat').transpose()
# 	chi2f = np.loadtxt(dir+'2Dbaofits/arat'+fout+'NScomb1covchigrid.dat').transpose()
# 	chi2f -= np.min(chi2f)
# 	ar = np.arange(.803,1.1931,.006)
# 	at = np.arange(.8015,1.19751,.003)
# 	#X,Y = np.meshgrid(chi2f[0],chi2f[1])
# 	#Z = chi2f[2].reshape(len(chi2f[0]),len(chi2f[0]))
# 	#plt.scatter(chi2f[0],chi2f[1],c=(chi2f[2]-np.min(chi2f[2])),vmax=1)
# 	if md == 'mockave':
# 		vmax=1
# 		levels = np.arange(0,vmax,vmax/10)
# 	if md == 'data':
# 		vmax = 25*2.3
# 		levels = [0,2.3,4.*2.3,9.*2.3,16*2.3]
# 	extent = (mina, maxa,mina,maxa)
# 	im = plt.imshow(chi2f,vmin=0,vmax=vmax,origin='lower',cmap = cm.PRGn,extent=extent,aspect='auto')
# 	#v = plt.axis()
# 	plt.colorbar(im)
# 	
# 	
# 	ct = plt.contour(chi2f, levels, colors='k', origin='lower',extent=extent)
# 	#plt.clabel(ct, inline=1, fontsize=12,fmt='%1.0f')
# 	#plt.axis(v)
# 	#ylim = plt.get(plt.gca(), 'ylim')
# 	#plt.setp(plt.gca(), ylim=ylim[::-1])
# 	plt.xlabel(r'$\alpha_{||}$')
# 	plt.ylabel(r'$\alpha_{\perp}$')
# 	plt.show()	
# 	if md == 'mockave':
# 		#plt.scatter(chi2f[0],chi2f[1],c=(chi2f[2]-np.min(chi2f[2])),vmax=1)
# 		plt.pcolormesh(ar,at,chi2f-np.min(chi2f),vmax=.1)
# 		plt.xlim(.98,1.02)
# 		plt.ylim(.98,1.02)
# 		plt.show()


# def mksubNERSC(start,N):
# 	#!/bin/bash
# 	#SBATCH --qos=shared
# 	#SBATCH --constraint=haswell
# 	#SBATCH --time=30
# 	#SBATCH --ntasks=1
# 	#SBATCH --mem=1GB
# 
# 	python eboss_DR16tools.py 1

def xi2DBAONandS(md='data',samp='LRGpCMASS',min=50.,max=150.,maxb=80.,bs=8,binc=0,minz=0.6,maxz=1.,covf=1.,mockn=0,ezver='5',ver='4',af='',wm='fkp',rec='',damp='',spat=0.003,spar=0.006,mina=.8,maxa=1.2,covm='JB',Bp=0.4):
	from baofit_pub2D import Xism_arat_1C_an
	maxind = int(200/bs)
	print(maxind)
	bs = int(bs)
	if md == 'data':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		d0N = np.loadtxt(dir+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]

	if md == 'dataJH':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		d0N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[2][:maxind]
# 		d0N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5.dat').transpose()[1][:maxind]
# 		d2N = np.loadtxt(dir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5.dat').transpose()[2][:maxind]
# 		d0S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5.dat').transpose()[1][:maxind]
# 		d2S = np.loadtxt(dir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5.dat').transpose()[2][:maxind]

	if md == 'dataN':
		#for on NERSC
		dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ezver)+'/'
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		d0N = np.loadtxt(dirH+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dirH+'xi024geboss'+samp+'_NGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dirH+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dirH+'xi024geboss'+samp+'_SGC'+ver+rd+'_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
	if md == 'dataJB':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		#rl = np.loadtxt(dir+'xi024geboss'+samp+'_SGC4_mz'+str(minz)+'xz'+str(maxz)+wm+str(bs)+'st'+str(binc)+'.dat').transpose()[0][:max]
		rd = ''
		if rec == 'rec':
			rd = '_rec'
		d0N = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'eBOSS_LRGpCMASS_clustering_COMB_v'+ver+'_rec_5mpc_shift0.mul').transpose()[2][:maxind]

	if md == 'mockave':
		#dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		#dirout = '/Users/ashleyross/Dropbox/eboss/' 
		#dirH = '/Users/ashleyross/Dropbox/eboss/'
		dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ezver)+'/'
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		#xiave02NGCQSO_v5_EZ5st0.dat

		d0N = np.loadtxt(dirH+'xiave02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dirH+'xiave02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][maxind:]
		d0S = np.loadtxt(dirH+'xiave02SGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dirH+'xiave02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat').transpose()[1][maxind:]

	if md == 'mockaveJB':
		dir = '/Users/ashleyross/Dropbox/eboss/' #change to wherever the data is
		dirout = '/Users/ashleyross/Dropbox/eboss/' 
		dirH = '/Users/ashleyross/Dropbox/eboss/'
		d0N = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.mul').transpose()[2][:maxind]

	if md == 'mock':
		dir = '/global/cscratch1/sd/ajross/ebossxi/'+ samp+'_v'+str(ezver)+'/'
		dirout = '/global/cscratch1/sd/ajross/'
		dirH = ''
		d0N = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'NGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2N = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'NGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]
		d0S = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'SGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[1][:maxind]
		d2S = np.loadtxt(dir+'xi024'+samp+'_v'+str(ezver)+'SGCEZmock'+rec+str(mockn)+str(bs)+'st'+str(binc)+'.dat').transpose()[2][:maxind]

	if covm == 'me':
		cN = np.loadtxt(dirH+'cov02NGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat')
		cS = np.loadtxt(dirH+'cov02SGC'+samp+'_v'+str(ezver)+'_EZ'+rec+af+str(bs)+'st'+str(binc)+'.dat')  
		
		if len(cN) != len(d0N)*2:
			print( 'MISMATCHED data and cov matrix!')
			print((len(cN),len(d0N)))
	dvN = [] #empty list to become data vector
	dvbN = [] #empty list to become data vector for setting bias prior
	dvS = [] #empty list to become data vector
	dvbS = [] #empty list to become data vector for setting bias prior

	rl = [] #empty list to become list of r values to evaluate model at	
	rlb  = [] #empty list to become list of r values to evaluate model at for bias prior
	mini = 0
	for i in range(0,len(d0N)):
		r = i*bs+bs/2.+binc
		if r > min and r < max:
			dvN.append(d0N[i])
			dvS.append(d0S[i])
			rbc = r#.75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.) #correct for pairs should have slightly larger average pair distance than the bin center
			rl.append(rbc) 
			if mini == 0:
				mini = i #minimum index number to be used for covariance matrix index assignment
			if r < maxb:
				dvbN.append(d0N[i])
				dvbS.append(d0S[i])
				rlb.append(rbc)
	for i in range(0,len(d2N)):
		r = i*bs+bs/2.+binc
		if r > min and r < max:
			dvN.append(d2N[i])
			dvS.append(d2S[i])
			rbc = r#.75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.)
			rl.append(rbc)
			if r < maxb:
				dvbN.append(d2N[i])
				dvbS.append(d2S[i])
				rlb.append(rbc)

	dvN = np.array(dvN)
	dvS = np.array(dvS)
	print(( len(dvN),len(dvS)))
	if covm == 'me':
		covmS = np.zeros((len(dvS),len(dvS))) #will become covariance matrix to be used with data vector
		covmN = np.zeros((len(dvS),len(dvS))) #will become covariance matrix to be used with data vector

		#need to cut it to correct size
		for i in range(0,len(cN)):
			if i < len(cN)/2:
				ri = i*bs+bs/2.+binc
				indi = i-mini
			else:
				ri = (i-len(cN)/2)*bs+bs/2.+binc
				indi = len(dvN)/2+i-mini-len(cN)/2	
			for j in range(0,len(cN)):		
				if j < len(cN)/2:
					rj = j*bs+bs/2.+binc
					indj = j-mini
				else:
					rj = (j-len(cN)/2)*bs+bs/2.+binc
					indj = len(dvN)/2+j-mini-len(cN)/2
				if ri > min and ri < max and rj > min and rj < max:
					#print ri,rj,i,j,indi,indj
					indi = int(indi)
					indj = int(indj)
					covmN[indi][indj] = cN[i][j]
					covmS[indi][indj] = cS[i][j]
		invcN = np.linalg.pinv(covmN) #the inverse covariance matrix to pass to the code
		invcS = np.linalg.pinv(covmS)
		covmbN = np.zeros((len(dvbN),len(dvbN)))
		covmbS = np.zeros((len(dvbN),len(dvbN)))
		for i in range(0,len(dvbN)):
			if i < len(dvbN)/2:
				indi = i
			else:
				indi = i-len(dvbN)/2+len(covmN)/2
			for j in range(0,len(dvbN)):
				if j < len(dvbN)/2:
					indj = j
				else:
					indj = j-len(dvbN)/2+len(covmN)/2
				indi = int(indi)
				indj = int(indj)
				covmbN[i][j] = covmN[indi][indj]
				covmbS[i][j] = covmS[indi][indj]
		invcbN = np.linalg.pinv(covmbN)
		invcbS = np.linalg.pinv(covmbS)
	if covm == 'JB':
		invct = mkinvcJB(min=min,max=max)
		if (maxb-min)/bs < 1:
			return 'ERROR, maxb needs to be larger to have any data to fit'
		invcbt = mkinvcJB(min=min,max=maxb)
	#if rec == 'rec':
	#	mod = 'Challenge_matterpower'+.dat' #BAO template used	post-recon
	#if rec == '':
	mod = 'Challenge_matterpower'+damp+'.dat'		
	fout = samp+rec+md+ver+str(bs)+covm+'Bp'+str(Bp)
	if md == 'mock':
		fout += 'mn'+str(mockn)

	print((len(dvN),len(invcN),len(dvbN),len(invcbN)))
	ansn = Xism_arat_1C_an(dvN,invcN*covf,rl,mod,dvbN,invcbN,rlb,bs=bs,dirout=dirout,fout=fout+'NGC',spat=spat,spar=spar,mina=mina,maxa=maxa,Bp=Bp,Bt=Bp)
	print(ansn)
	anss = Xism_arat_1C_an(dvS,invcS*covf,rl,mod,dvbS,invcbS,rlb,bs=bs,dirout=dirout,fout=fout+'SGC',spat=spat,spar=spar,mina=mina,maxa=maxa,Bp=Bp,Bt=Bp)
	print(anss)
	ds = np.loadtxt(dirout+'2Dbaofits/arat'+fout+'SGC1covchi.dat')
	dn = np.loadtxt(dirout+'2Dbaofits/arat'+fout+'NGC1covchi.dat')
	dt = ds+dn
	ans = sigreg_2dmea(dt)
	print(ans)
	return True


def mkinvcJB(min=50,max=150,bs=5,shift=0):
	dir = '/Users/ashleyross/DR16plots2/LRG/data/'
	d = np.loadtxt(dir+'EZmock_eBOSS_LRGpCMASS_COMB_v7_recon_average_'+str(bs)+'mpc.cov').transpose()
	siz = int((max-min)//bs*2)
	print(siz)
	cov = np.zeros((siz,siz))
	db = min//bs
	print(db)
	for i in range(0,len(d[0])):
		r1 = d[2][i]
		r2 = d[3][i]
		if r1 > min and r1 < max and r2 > min and r2 < max:
			ind1 = d[0][i]
			ind2 = d[1][i]
			if ind1 < 39:
				c1 = int(ind1)-	db
			if ind1 >= 39 and ind1 < 78:
				c1 = int(ind1)-db-39+int(siz//2)	
			if ind2 < 39:
				c2 = int(ind2)-	db
			if ind2 >= 39 and ind2 < 78:
				c2 = int(ind2)-db-39+int(siz//2)	
			#print(c1,c2)
			cov[int(c1)][int(c2)] = d[4][i]
	for i in range(0,len(cov)):
		for j in range(0,len(cov)):
			if cov[i][j] == 0:
				print(('0 cov element found at '+str(i)+' '+str(j)))
	icov = np.linalg.pinv(cov)
	return icov		
	
def faccalc(nm,nb,nd):
	A = 2./(nm-nb-1.)/(nm-nb-4.)
	B = (nm-nb-2.)/(nm-nb-1.)/(nm-nb-4.)
	print((A,B))
	mv = (1.+B*(nb-nd))/(1.+A+B*(nd+1.))
	md = (nm-1.)/(nm-nb-2.)*mv
	return mv,md		


def mksubfile_mocks(ind):
	
	fo = open('sub'+str(ind)+'.sh','w')
	fl = ''
# 	if ind < 1000:
# 		fl += '0'
# 	if ind < 100:
# 		fl += '0'
# 	if ind < 10:
# 		fl += '0'
	fl += str(ind)
	##print fl
	fo.write('#!/bin/bash\n')
	#fo.write('#$ -V -cwd\n')
	#fo.write('. /etc/profile.d/modules.sh \n')
	#fo.write('module add  apps/gcc/python/2.7.3 \n')
	#fo.write('/opt/gridware/apps/gcc/python/2.7.3/bin/python baofit.py '+file +' '+str(B)+' CMASS '+reg+' '+col+' '+tp +' '+str(bs)+' '+str(st)+'\n')
	fo.write('python eboss_DR16tools.py '+fl+'\n')
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


		
def plotallmocks(mom=0,reg='NGC',rec='',N=1000,bs=8,start=0,muw=''):
	from matplotlib import pyplot as plt
	dir = '/mnt/lustre/ashleyr/eboss/EZmockELGv4/'
	for i in range(1,1+N):
		d = np.loadtxt(dir+'xi/xi024ELG'+reg+'EZmock'+rec+str(i)+muw+str(bs)+'st'+str(start)+'.dat').transpose()
		plt.plot(d[0],d[0]**2.*d[mom/2+1],'k-')
	plt.show()
	return True

def putall2DBAO(start,N,samp='QSOmockjh45meBp0.4',damp='0.44.03.08.015.00',sigmax=0.2):
	dir = '/global/cscratch1/sd/ajross/2Dbaofits/'
	atl = []
	arl = []
	chil = []
	atpl = []
	arpl = []
	atel = []
	arel = []
	al =[]
	sl = []
	nt = 0
	ng = 0
	fo = open('BAO2D_EZmocks_eBOSS'+samp+damp+'DR16.dat','w')
	fo.write('#mock_number max_lik_arad max_like_atran arad atran sig_arad sig_atran chi2 alpha_iso sig_iso \n')
	for i in range(start,start+N):		
		an = sigreg_c12('/global/cscratch1/sd/ajross/baofits/BAOxichilNScombfQSOEZmockjh'+str(i)+'7'+damp+'5st0')
		#print(i,an)
		alph = (an[2]+an[1])/2.
		sig = (an[2]-an[1])/2.
		#if sig < sigmax and (an[0]-1.)/sig < 3.:
		#	try:
		f = np.loadtxt(dir+'arat'+samp+'mn'+str(i)+'NScomb1covchi.dat').transpose()
		chimin = np.min(f[2])
		ind = np.argmin(f[2])#np.where(f[2] == np.amin(f[2]))[0][0]
		ar = f[0][ind]
		at = f[1][ind]
	
		chil.append(chimin)
		atl.append(at)
		arl.append(ar)
		#al.append(an[0])
		#sl.append(sig)
		ans = sigreg_2dme(dir+'arat'+samp+damp+'mn'+str(i)+'NScomb1covchi')
		fo.write(str(i)+' '+str(ar)+' '+str(at)+' '+str(ans[0])+' '+str(ans[2])+' '+str(ans[1])+' '+str(ans[3])+' '+str(chimin)+' '+str(alph)+' '+str(sig)+'\n')
		arpl.append(ans[0])
		arel.append(ans[1])
		atpl.append(ans[2])
		atel.append(ans[3])
		if ans[1] < 0.04 and ans[3] < 0.025:
			ng += 1
		#a1b,err1,a2b,err2
		nt += 1
		print((i,nt))
			#except:
			#	print(str(i)+' bad')
	print(ng)		
	print((nt,np.mean(atl),np.std(atl),np.mean(arl),np.std(arl),np.median(atl),np.median(arl),np.mean(chil)))#,np.mean(al),np.mean(sl),np.std(al))
	print((nt,np.mean(atpl),np.std(atpl),np.mean(arpl),np.std(arpl),np.median(atpl),np.median(arpl),np.mean(atel),np.mean(arel)))
	return True

def putall2DBAO_ORblind(vb,N=100,start=0,damp='ORjhBp0.40.41.03.08.015.00',sigmax=0.2):
	dir = '/global/cscratch1/sd/ajross/2Dbaofits/'
	atl = []
	arl = []
	chil = []
	atpl = []
	arpl = []
	atel = []
	arel = []
	al =[]
	sl = []
	nt = 0
	ng = 0
	if vb < 10:
		orver = 'mocks_blind_v1'
	else:
		orver = 'mocks_blind_v2'
	orcos = 'mockb'+str(vb)	
	fo = open('BAO2D_ORblindmocks_eBOSS'+str(vb)+damp+'DR16.dat','w')
	fo.write('#mock_number max_lik_arad max_like_atran arad atran sig_arad sig_atran chi2 C_||,perp Corr_||,perp \n')
	for i in range(start,start+N):
		fl = dir+'aratQSOmockORjh'+orver+orcos+damp+'mn'+str(i)+'1covchi.dat'		
		f = np.loadtxt(fl).transpose()
		chimin = np.min(f[2])
		ind = np.argmin(f[2])#np.where(f[2] == np.amin(f[2]))[0][0]
		ar = f[0][ind]
		at = f[1][ind]
	
		chil.append(chimin)
		atl.append(at)
		arl.append(ar)
		ans = sigreg_2dEZ(fl)
		fo.write(str(i)+' '+str(ar)+' '+str(at)+' '+str(ans[0])+' '+str(ans[2])+' '+str(ans[1])+' '+str(ans[3])+' '+str(chimin)+' '+str(ans[-2])+' '+str(ans[-1])+'\n')
		arpl.append(ans[0])
		arel.append(ans[1])
		atpl.append(ans[2])
		atel.append(ans[3])
		if ans[1] < 0.04 and ans[3] < 0.025:
			ng += 1
		#a1b,err1,a2b,err2
		nt += 1
		print((i,nt))
			#except:
			#	print(str(i)+' bad')
	print(ng)		
	print((nt,np.mean(atl),np.std(atl),np.mean(arl),np.std(arl),np.median(atl),np.median(arl),np.mean(chil)))#,np.mean(al),np.mean(sl),np.std(al))
	print((nt,np.mean(atpl),np.std(atpl),np.mean(arpl),np.std(arpl),np.median(atpl),np.median(arpl),np.mean(atel),np.mean(arel)))
	return np.mean(atl),np.std(atl),np.mean(arl),np.std(arl)

def putall2DBAO_ORcat(vb,N=100,start=0,damp='ORjhBp0.40.41.03.08.015.00',sigmax=0.2,orver = 'catastrophic_redshifts'):
	dir = '/global/cscratch1/sd/ajross/2Dbaofits/'
	atl = []
	arl = []
	chil = []
	atpl = []
	arpl = []
	atel = []
	arel = []
	al =[]
	sl = []
	nt = 0
	ng = 0		
	fo = open('BAO2D_ORblindmocks_eBOSS'+str(vb)+orver+damp+'DR16.dat','w')
	fo.write('#mock_number max_lik_arad max_like_atran arad atran sig_arad sig_atran chi2 C_||,perp Corr_||,perp \n')
	for i in range(start,start+N):
		fl = dir+'aratQSOmockORjh_notblindmock'+str(vb)+orver+damp+'mn'+str(i)+'1covchi.dat'		
		f = np.loadtxt(fl).transpose()
		chimin = np.min(f[2])
		ind = np.argmin(f[2])#np.where(f[2] == np.amin(f[2]))[0][0]
		ar = f[0][ind]
		at = f[1][ind]
	
		chil.append(chimin)
		atl.append(at)
		arl.append(ar)
		ans = sigreg_2dEZ(fl)
		fo.write(str(i)+' '+str(ar)+' '+str(at)+' '+str(ans[0])+' '+str(ans[2])+' '+str(ans[1])+' '+str(ans[3])+' '+str(chimin)+' '+str(ans[-2])+' '+str(ans[-1])+'\n')
		arpl.append(ans[0])
		arel.append(ans[1])
		atpl.append(ans[2])
		atel.append(ans[3])
		if ans[1] < 0.04 and ans[3] < 0.025:
			ng += 1
		#a1b,err1,a2b,err2
		nt += 1
		print((i,nt))
			#except:
			#	print(str(i)+' bad')
	print(ng)		
	print((nt,np.mean(atl),np.std(atl),np.mean(arl),np.std(arl),np.median(atl),np.median(arl),np.mean(chil)))#,np.mean(al),np.mean(sl),np.std(al))
	print((nt,np.mean(atpl),np.std(atpl),np.mean(arpl),np.std(arpl),np.median(atpl),np.median(arpl),np.mean(atel),np.mean(arel)))
	return np.mean(atl),np.std(atl),np.mean(arl),np.std(arl)


def combORblind():
	fo = open('allORblindAJR.dat','w')
	fo.write('#using maximum likelihood position\n')
	fo.write('#bn <aperp> std(aperp) <arad> std(arad)\n')
	for i in range(0,25):
		a= putall2DBAO_ORblind(i)
		fo.write(str(i)+' '+str(a[0])+' '+str(a[1])+' '+str(a[2])+' '+str(a[3])+'\n')
		print('done with '+str(i))
	fo.close()
	return True

def combORcat(N=11,orver='catastrophic_redshifts',damp='ORjhBp0.40.44.03.08.015.00'):
	fo = open('allOR'+orver+damp+'AJR.dat','w')
	fo.write('#using maximum likelihood position\n')
	fo.write('#bn <aperp> std(aperp) <arad> std(arad)\n')
	for i in range(0,N):
		a= putall2DBAO_ORcat(i,orver=orver,damp=damp)
		fo.write(str(i)+' '+str(a[0])+' '+str(a[1])+' '+str(a[2])+' '+str(a[3])+'\n')
		print('done with '+str(i))
	fo.close()
	return True


def putall1DBAO(start=1,N=1000,samp='ELG',md='7_rec0.59304.07.015.01.08st0'):
	dir = '/global/cscratch1/sd/ajross/baofits/'
	ng = 0
	am = 0
	sm = 0
	ss = 0
	fo = open('BAO1D_EZmocks_'+samp+md+'.dat','w')
	fo.write('#mock_number a_ngc sigma_ngc chi2 a_sgc sigma_sgc chi2 a_sgcngcl sigma_sgcngcl chi2 a_sgcngcx sigma_sgcngcx chi2\n')
	regl = ['NGC','SGC','NScomb','NScombf']
	for i in range(start,start+N):		
		for reg in regl:
			an = sigreg_c12('/global/cscratch1/sd/ajross/baofits/BAOxichil'+reg+samp+str(i)+md)
		#print(i,an)
		#(a[2]-a[1])/2.
			sig = (an[2]-an[1])/2.
			al = (an[2]+an[1])/2.
			chi = an[-1]
			fo.write(str(al)+' '+str(sig)+' '+str(chi)+' ')
			if reg == 'NScombf':
				if sig < 0.2:
					ng += 1.
					am += al
					sm += sig
					ss += al*al
		fo.write('\n')	
		print(i)
	fo.close()
	am = am/ng
	sm = sm/ng
	ss = ss/ng	
	return ng,am,sm,sqrt(ss-am*am)

def putall1DBAO_comb5(start=1,N=1000,samp='ELG',md='7_rec0.59304.07.015.01.08st'):
	dir = '/global/cscratch1/sd/ajross/baofits/'
	ng = 0
	am = 0
	sm = 0
	ss = 0
	fo = open('BAO1D_EZmocks_'+samp+md+'comb5.dat','w')
	fo.write('#mock_number a_ngc sigma_ngc chi2 a_sgc sigma_sgc chi2 a_sgcngcl sigma_sgcngcl chi2 a_sgcngcx sigma_sgcngcx chi2\n')
	regl = ['NGC','SGC','NScomb','NScombf']
	for i in range(start,start+N):		
		for reg in regl:
			an = sigreg_comb5('/global/cscratch1/sd/ajross/baofits/BAOxichil'+reg+samp+str(i)+md)
		#print(i,an)
		#(a[2]-a[1])/2.
			sig = (an[2]-an[1])/2.
			al = (an[2]+an[1])/2.
			chi = an[-1]
			fo.write(str(al)+' '+str(sig)+' '+str(chi)+' ')
			if reg == 'NScombf':
				if sig < 0.2:
					ng += 1.
					am += al
					sm += sig
					ss += al*al
		fo.write('\n')	
		print(i)
	fo.close()
	am = am/ng
	sm = sm/ng
	ss = ss/ng	
	return ng,am,sm,sqrt(ss-am*am)
	
def BAO1Dstats(md,samp='ELG',reg='sgcngcl',sigmax=0.2):
	d = np.loadtxt(ebossdir+'BAO1D_EZmocks_'+samp+md+'.dat').transpose()
	if reg == 'ngc':
		acol = 0
	if reg == 'sgc':
		acol = 3
	if reg == 'sgcngcl':
		acol = 6
	if reg == 'sgcngcx':
		acol = 9
	scol = acol + 1
	ccol = acol + 2
	w = d[scol] < sigmax
	print('Nmocks,<alpha>,<sigma>,rms,<chi2>')
	print((len(d[0][w]),np.mean(d[acol][w]),np.mean(d[scol][w]),np.std(d[acol][w]),np.mean(d[ccol][w])))
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
	print (chin)
	lf = linfit(xl,yl,el)
	inl = np.array([1.,0])
	b,m = fmin(lf.chilin,inl)
	chilin = sum((yl-(m*xl+b))**2./el**2.)
	print (chilin)
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

def getCompstats(samp):
	if samp == 'LRG':
		ver = '7'
	if samp == 'QSO':
		ver = '7_1'
	fn = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_NGC_v'+ver+'.dat.fits')
	fs = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_SGC_v'+ver+'.dat.fits')
	fnr = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_NGC_v'+ver+'.ran.fits')
	fsr = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_SGC_v'+ver+'.ran.fits')
	wcn = fnr['COMP_BOSS'] <= 0.5
	wcs = fsr['COMP_BOSS'] <= 0.5
	areac = (len(fnr[wcn])+len(fsr[wcs]))/5000.
	print('area removed by completeness cut for '+samp+' is '+str(areac))
	wcn = (fnr['COMP_BOSS'] > 0.5) & (fnr['sector_SSR'] <= 0.5)
	wcs = (fsr['COMP_BOSS'] > 0.5) & (fsr['sector_SSR'] <= 0.5)
	areac = (len(fnr[wcn])+len(fsr[wcs]))/5000.
	print('area removed by SSR cut (after completeness) for '+samp+' is '+str(areac))
	wcn = (fnr['COMP_BOSS'] > 0.5) & (fnr['sector_SSR'] > 0.5)
	wcs = (fsr['COMP_BOSS'] > 0.5) & (fsr['sector_SSR'] > 0.5)
	arean = len(fnr[wcn])/5000.
	areas = len(fsr[wcs])/5000.
	areac = (arean+areas)
	print('area left '+samp+' is '+str(areac)+ ' and is '+ str(arean) +' in NGC and '+str(areas)+' in SGC')
	arean = sum(fnr[wcn]['COMP_BOSS'])/5000.
	areas = sum(fsr[wcs]['COMP_BOSS'])/5000.
	areac = (arean+areas)
	print('weighted area left '+samp+' is '+str(areac)+ ' and is '+ str(arean) +' in NGC and '+str(areas)+' in SGC')

	arean = len(fnr)/5000.
	areas = len(fsr)/5000.
	areac = (arean+areas)
	print('total area (post veto) for '+samp+' is '+str(areac)+ ' and is '+ str(arean) +' in NGC and '+str(areas)+' in SGC')
	
	wcn = fn['COMP_BOSS'] <= 0.5
	wcni = (fn['COMP_BOSS'] <= 0.5) & (fn['IMATCH'] == 1)
	wcs = fs['COMP_BOSS'] <= 0.5
	wcsi = (fs['COMP_BOSS'] <= 0.5) & (fs['IMATCH'] == 1)
	print('number of targets removed by completeness cut for '+samp+' is '+str(len(fn[wcn])+len(fs[wcs])))
	print('number of eboss redshifts removed by completeness cut for '+samp+' is '+str(len(fn[wcni])+len(fs[wcsi])))
	wcn = (fn['COMP_BOSS'] > 0.5) & (fn['sector_SSR'] <= 0.5)
	wcni = (fn['COMP_BOSS'] > 0.5) & (fn['IMATCH'] == 1) & (fn['sector_SSR'] <= 0.5)
	wcs = (fs['COMP_BOSS'] > 0.5) & (fs['sector_SSR'] <= 0.5)
	wcsi = (fs['COMP_BOSS'] > 0.5) & (fs['IMATCH'] == 1) & (fs['sector_SSR'] <= 0.5)
	print('number of targets removed by SSR cut (after completeness) for '+samp+' is '+str(len(fn[wcn])+len(fs[wcs])))
	print('number of eboss redshifts removed by SSR cut (after completeness) for '+samp+' is '+str(len(fn[wcni])+len(fs[wcsi])))

	wcn = (fn['COMP_BOSS'] > 0.5) & (fn['sector_SSR'] > 0.5)
	wcni = (fn['COMP_BOSS'] > 0.5) & (fn['IMATCH'] == 1) & (fn['sector_SSR'] > 0.5)
	wcs = (fs['COMP_BOSS'] > 0.5) & (fs['sector_SSR'] > 0.5)
	wcsi = (fs['COMP_BOSS'] > 0.5) & (fs['IMATCH'] == 1) & (fs['sector_SSR'] > 0.5)
	print('number of targets after completeness cuts for '+samp+' is '+str(len(fn[wcn])+len(fs[wcs])))
	print('number of eboss redshifts after completeness cuts for '+samp+' is '+str(len(fn[wcni])+len(fs[wcsi])))

def compLRG_BAORSD_RP():
	bf = np.loadtxt(ebossdir+'EZmock_eBOSS_LRGpCMASS_COMB_v7_recon_5mpc-rmin50.0-rmax150.0-bb-quad.txt').transpose()
	rf = np.loadtxt(ebossdir+'Ezmock_v7_RP_xi024.dat').transpose()
	w = (bf[3] == 1) & (bf[6] == 1)	
	aperpb = bf[1][w]
	aperpr = rf[2][w]
	eperpb = bf[2][w]
	eperpr = rf[3][w]
	plt.plot(aperpb,aperpr,'ko')
	datb = 1.026
	datr = 1.01
	detb = 0.0190
	detr = 0.022
	plt.plot(datb,datr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\perp}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	xl = [-1,2]
	yl = [-1,2]
	plt.plot(xl,yl,'r--')

	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPaperp.png')
	plt.show()
	plt.plot(eperpb,eperpr,'ko')
	plt.plot(detb,detr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\sigma\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{\perp}$')
	plt.xlim(0,0.05)
	plt.ylim(0,0.05)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPsigaperp.png')
	plt.show()
	
	aparb = bf[4][w]
	aparr = rf[0][w]
	eparb = bf[5][w]
	eparr = rf[1][w]
	print(np.mean(eperpb),np.mean(eperpr))
	print(np.std(aperpb),np.std(aperpr))
	print(np.mean(eparb),np.mean(eparr))
	print(np.std(aparb),np.std(aparr))
	plt.plot(aparb,aparr,'ko')
	darb = 0.956
	darr = 1.019
	derb = 0.027
	derr = 0.037

	plt.plot(darb,darr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{||}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPapar.png')
	plt.show()

	plt.plot(eparb,eparr,'ko')
	plt.plot(derb,derr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\sigma\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{||}$')
	plt.xlim(0,0.1)
	plt.ylim(0,0.1)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPsigaperp.png')
	plt.show()
	
	
	aisob = (aperpb**2.*aparb)**(1/3.)
	aisor = (aperpr**2.*aparr)**(1/3.)
	plt.plot(aisob,aisor,'ko')
	daisob = (datb**2.*darb)**(1/3.)
	daisor = (datr**2.*darr)**(1/3.)
	plt.plot(daisob,daisor,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(0.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$\xi$ BAO $\alpha_{\rm iso}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\rm iso}$')
	
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPasio.png')
	plt.show()
	epb = (aparb/aperpb)**(1/3.)
	epr = (aparr/aperpr)**(1/3.)
	plt.plot(epb,epr,'ko')
	depb = (darb/datb)**(1/3.)
	depr = (darr/datr)**(1/3.)
	plt.plot(depb,depr,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(0.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$\xi$ BAO $\epsilon$+1')
	plt.ylabel(r'$\xi$ RSD RP $\epsilon$+1')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPep.png')
	plt.show()
	print(np.std(aisob),np.std(aisor))
	print(np.std(epb),np.std(epr))
	ht = plt.hist(abs(aisob-aisor))
	yl = [0,1000]
	dda = abs(daisob-daisor)
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\alpha_{\rm iso}|$')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPalphhist.png')
	plt.show()


	ht = plt.hist(abs(epb-epr)/(eparb))
	yl = [0,1000]
	dda = abs(depb-depr)/derb
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\epsilon|$(BAO RP RSD)/$\sigma_{||,BAO}$')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPephist.png')
	plt.show()

def compLRG_BAORSD_HGM():
	bf = np.loadtxt(ebossdir+'EZmock_eBOSS_LRGpCMASS_COMB_v7_recon_5mpc-rmin50.0-rmax150.0-bb-quad.txt').transpose()
	rf = np.loadtxt(ebossdir+'EZmocksv7_HGM_P024_kmax015_015_FS_Anoiseprior_prerecon_fix.txt').transpose()
	w = (bf[3] == 1) & (bf[6] == 1)	 & (rf[1] != 0)
	aperpb = bf[1][w]
	aperpr = rf[2][w]
	eperpb = bf[2][w]
	eperpr = rf[3][w]
	plt.plot(aperpb,aperpr,'ko')
	datb = 1.026
	datr = 1.008
	detb = 0.0190
	detr = 0.026

	plt.plot(datb,datr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\perp}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	xl = [-1,2]
	yl = [-1,2]
	plt.plot(xl,yl,'r--')

	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_HGMaperp.png')
	plt.show()
	plt.plot(eperpb,eperpr,'ko')
	plt.plot(detb,detr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\sigma\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{\perp}$')
	plt.xlim(0,0.05)
	plt.ylim(0,0.05)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPsigaperp.png')
	plt.show()
	
	aparb = bf[4][w]
	aparr = rf[0][w]
	eparb = bf[5][w]
	eparr = rf[1][w]
	print(np.mean(eperpb),np.mean(eperpr))
	print(np.std(aperpb),np.std(aperpr))
	print(np.mean(eparb),np.mean(eparr))
	print(np.std(aparb),np.std(aparr))
	plt.plot(aparb,aparr,'ko')
	darb = 0.956
	darr = 0.994
	derb = 0.027
	derr = 0.036

	plt.plot(darb,darr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{||}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPapar.png')
	plt.show()

	plt.plot(eparb,eparr,'ko')
	plt.plot(derb,derr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\sigma\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{||}$')
	plt.xlim(0,0.1)
	plt.ylim(0,0.1)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPsigaperp.png')
	plt.show()
	
	
	aisob = (aperpb**2.*aparb)**(1/3.)
	aisor = (aperpr**2.*aparr)**(1/3.)
	plt.plot(aisob,aisor,'ko')
	daisob = (datb**2.*darb)**(1/3.)
	daisor = (datr**2.*darr)**(1/3.)
	plt.plot(daisob,daisor,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(0.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$\xi$ BAO $\alpha_{\rm iso}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\rm iso}$')
	
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPasio.png')
	plt.show()
	epb = (aparb/aperpb)**(1/3.)
	epr = (aparr/aperpr)**(1/3.)
	plt.plot(epb,epr,'ko')
	depb = (darb/datb)**(1/3.)
	depr = (darr/datr)**(1/3.)
	plt.plot(depb,depr,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(0.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$\xi$ BAO $\epsilon$+1')
	plt.ylabel(r'$\xi$ RSD RP $\epsilon$+1')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPep.png')
	plt.show()
	print(np.std(aisob),np.std(aisor))
	print(np.std(epb),np.std(epr))
	ht = plt.hist(abs(aisob-aisor))
	yl = [0,1000]
	dda = abs(daisob-daisor)
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\alpha_{\rm iso}|$')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPalphhist.png')
	plt.show()


	ht = plt.hist(abs(epb-epr)/(eparb))
	yl = [0,1000]
	dda = abs(depb-depr)/derb
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\epsilon|$(BAO RP RSD)/$\sigma_{||,BAO}$')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_RPephist.png')
	plt.show()


def compLRG_BAORSD_MV():	
	#rf = np.loadtxt(ebossdir+'Ezmock_v7_RP_xi024.dat').transpose()
	rf = np.loadtxt(ebossdir+'EZ_MVM_xi024_FS_fidcosmo_bugfree.txt').transpose()
	bf = np.loadtxt(ebossdir+'EZmock_eBOSS_LRGpCMASS_COMB_v7_recon_5mpc-rmin50.0-rmax150.0-bb-quad.txt')[:len(rf[0])].transpose()
	w = (bf[3] == 1) & (bf[6] == 1)	& (rf[0] != 0)
	aperpb = bf[1][w]
	aperpr = rf[1][w]
	eperpb = bf[2][w]
	#eperpr = rf[3][w]
	plt.plot(aperpb,aperpr,'ko')
	datb = 1.026
	datr = .991
	detb = 0.0190
	detr = 0.021
	plt.plot(datb,datr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD MVM $\alpha_{\perp}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	xl = [-1,2]
	yl = [-1,2]
	plt.plot(xl,yl,'r--')

	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVaperp.png')
	plt.show()
# 	plt.plot(eperpb,eperpr,'ko')
# 	plt.plot(detb,detr,'*',color='orange',markersize=14)
# 	plt.xlabel(r'$\xi$ BAO $\sigma\alpha_{\perp}$')
# 	plt.ylabel(r'$\xi$ RSD MV $\sigma\alpha_{\perp}$')
# 	plt.xlim(0,0.05)
# 	plt.ylim(0,0.05)
# 	plt.plot(xl,yl,'r--')
# 	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVsigaperp.png')
# 	plt.show()
	
	aparb = bf[4][w]
	aparr = rf[0][w]
	eparb = bf[5][w]
	eparr = rf[1][w]
	#print(np.mean(eperpb),np.mean(eperpr))
	print('rms BAOperp, RSDperp')
	print(np.std(aperpb),np.std(aperpr))
	#print(np.mean(eparb),np.mean(eparr))
	print('rms BAO||, RSD||')
	print(np.std(aparb),np.std(aparr))
	plt.plot(aparb,aparr,'ko')
	darb = 0.956
	darr = 1.011
	derb = 0.027
	derr = 0.029
	plt.plot(xl,yl,'r--')

	plt.plot(darb,darr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ BAO $\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD MVM $\alpha_{||}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVapar.png')
	plt.show()

# 	plt.plot(eparb,eparr,'ko')
# 	plt.plot(derb,derr,'*',color='orange',markersize=14)
# 	plt.xlabel(r'$\xi$ BAO $\sigma\alpha_{||}$')
# 	plt.ylabel(r'$\xi$ RSD MV $\sigma\alpha_{||}$')
# 	plt.xlim(0,0.1)
# 	plt.ylim(0,0.1)
# 	plt.plot(xl,yl,'r--')
# 	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVsigapar.png')
# 	plt.show()
	
	
	aisob = (aperpb**2.*aparb)**(1/3.)
	aisor = (aperpr**2.*aparr)**(1/3.)
	plt.plot(aisob,aisor,'ko')
	daisob = (datb**2.*darb)**(1/3.)
	daisor = (datr**2.*darr)**(1/3.)
	plt.plot(daisob,daisor,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$P$ BAO $\alpha_{\rm iso}$')
	plt.ylabel(r'$\xi$ RSD MVM $\alpha_{\rm iso}$')
	
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVaiso.png')
	plt.show()

	epb = (aparb/aperpb)**(1/3.)
	epr = (aparr/aperpr)**(1/3.)
	plt.plot(epb,epr,'ko')
	depb = (darb/datb)**(1/3.)
	depr = (darr/datr)**(1/3.)
	plt.plot(depb,depr,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$P$ BAO $\epsilon$+1')
	plt.ylabel(r'$\xi$ RSD MVM $\epsilon$+1')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVep.png')

	plt.show()

	ht = plt.hist(abs(aisob-aisor))
	yl = [0,1000]
	dda = abs(daisob-daisor)
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\alpha_{\rm iso}|$(BAO MVM RSD)')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVaisohist.png')
	plt.show()
	print(np.std(aisob),np.std(aisor))
	print(np.std(epb),np.std(epr))
	ht = plt.hist(abs(epb-epr)/(eparb))
	yl = [0,1000]
	dda = abs(depb-depr)/derb
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\epsilon|$(BAO MVM RSD)/$\sigma_{||,BAO}$')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSD_MVephist.png')
	plt.show()


def fixpkfile():
	f = open(ebossdir+'EZmocksv7_HGM_P024_kmax015_015_FS_Anoiseprior_prerecon.txt').readlines()
	fo = open(ebossdir+'EZmocksv7_HGM_P024_kmax015_015_FS_Anoiseprior_prerecon_fix.txt','w')
	for line in f:
		ln = line.split()
		if len(ln) == 8:
			fo.write(line)
		else:
			for i in range(0,8):
				fo.write(ln[i]+' ')
			fo.write('\n')		
	fo.close()	

def compLRGRSD_2xi():
	
	mf = np.loadtxt(ebossdir+'EZ_MVM_xi024_FS_fidcosmo_bugfree.txt').transpose()
	rf = np.loadtxt(ebossdir+'Ezmock_v7_RP_xi024.dat')[:len(mf[0])].transpose()

	w = (mf[1] != 0)# & (bf[3] != 0)	
	aperpb = mf[1][w]
	aperpr = rf[2][w]
	#eperpb = bf[4][w]
	eperpr = rf[3][w]
	plt.plot(aperpb,aperpr,'ko')
	datb = 0.991
	datr = 1.01
	detb = 0.021
	detr = 0.022

	plt.plot(datb,datr,'*',color='orange',markersize=14)
	plt.xlabel(r'$P$ RSD MVM $\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\perp}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	xl = [-1,2]
	yl = [-1,2]
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxiaperp.png')
	plt.show()

# 	plt.plot(eperpb,eperpr,'ko')
# 	plt.plot(detb,detr,'*',color='orange',markersize=14)
# 	plt.xlabel(r'$\xi$ RSD MVM $\sigma\alpha_{\perp}$')
# 	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{\perp}$')
# 	plt.xlim(0,0.05)
# 	plt.ylim(0,0.05)
# 	plt.plot(xl,yl,'r--')
# 	plt.show()
	
	aparb = mf[0][w]
	aparr = rf[0][w]
	#eparb = bf[2][w]
	eparr = rf[1][w]

	fb = mf[2][w]
	fr = rf[4][w]
	#eparb = bf[2][w]
	efr = rf[5][w]

	#print(np.mean(eperpb),np.mean(eperpr))
	print('aperp rms MVM,RP')
	print(np.std(aperpb),np.std(aperpr))
	#print(np.mean(eparb),np.mean(eparr))
	print('a|| rms MVM,RP')
	print(np.std(aparb),np.std(aparr))

	print('fsig8 rms MVM,RP')
	print(np.std(fb),np.std(fr))

	
	darb = 1.011
	darr = 1.019
	derb = 0.029
	derr = 0.037

	dfb = .467
	dfr = .431
	defb = 0.047
	defr = 0.044

	plt.plot(fb,fr,'ko')
	plt.plot(dfb,dfr,'*',color='orange',markersize=14)
	plt.xlabel(r'$P$ RSD MVM $f\sigma_8$')
	plt.ylabel(r'$\xi$ RSD RP $f\sigma_8$')
	plt.xlim(.3,.6)
	plt.ylim(.3,.6)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxifsig8.png')
	plt.show()


	plt.plot(aparb,aparr,'ko')
	plt.plot(darb,darr,'*',color='orange',markersize=14)
	plt.xlabel(r'$P$ RSD MVM $\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{||}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	plt.plot(xl,yl,'r--')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxiaper.png')
	plt.show()

# 	plt.plot(eparb,eparr,'ko')
# 	plt.plot(derb,derr,'*',color='orange',markersize=14)
# 	plt.xlabel(r'$P$ RSD HGM $\sigma\alpha_{||}$')
# 	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{||}$')
# 	plt.xlim(0,0.1)
# 	plt.ylim(0,0.1)
# 	plt.plot(xl,yl,'r--')
# 	plt.show()
	
	
	aisob = (aperpb**2.*aparb)**(1/3.)
	aisor = (aperpr**2.*aparr)**(1/3.)
	plt.plot(aisob,aisor,'ko')
	daisob = (datb**2.*darb)**(1/3.)
	daisor = (datr**2.*darr)**(1/3.)
	plt.plot(daisob,daisor,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$P$ RSD MVM $\alpha_{\rm iso}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\rm iso}$')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxiaiso.png')
	plt.show()

	epb = (aparb/aperpb)**(1/3.)
	epr = (aparr/aperpr)**(1/3.)
	plt.plot(epb,epr,'ko')
	depb = (darb/datb)**(1/3.)
	depr = (darr/datr)**(1/3.)
	plt.plot(depb,depr,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$P$ RSD MVM $\epsilon+1$')
	plt.ylabel(r'$\xi$ RSD RP $\epsilon+1$')
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxiep.png')
	plt.show()

	ht = plt.hist(abs(fb-fr)/efr)
	yl = [0,1000]
	dda = abs(dfb-dfr)/defr
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta f\sigma_8|$(MVM RP)/$\sigma$RP')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxihistfs8.png')
	plt.show()


	ht = plt.hist(abs(aisob-aisor))
	yl = [0,1000]
	dda = abs(daisob-daisor)
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\alpha_{\rm iso}|$(MVM RP)')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxihistaiso.png')
	plt.show()


	print(np.std(aisob),np.std(aisor))
	print(np.std(epb),np.std(epr))
	plt.hist(abs(epb-epr)/(eparr))
	hn = np.histogram(abs(epb-epr)/(eparr))
	yl = [0,1000]
	dda = abs(depb-depr)/sqrt(derr)
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\epsilon|$/$\sigma_{||,RP}$')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(hn[0])*1.1)
	plt.savefig(ebossdir+'LRGRSDcompplots/BAORSDxihistep.png')
	plt.show()


def compLRGRSD_pkxiRP():
	bf = np.loadtxt(ebossdir+'EZmocksv7_HGM_P024_kmax015_015_FS_Anoiseprior_prerecon_fix.txt').transpose()
	rf = np.loadtxt(ebossdir+'Ezmock_v7_RP_xi024.dat').transpose()
	w = (bf[1] != 0) & (bf[3] != 0)	
	aperpb = bf[3][w]
	aperpr = rf[2][w]
	eperpb = bf[4][w]
	eperpr = rf[3][w]
	plt.plot(aperpb,aperpr,'ko')
	datb = 1.008
	datr = 1.01
	detb = 0.026
	detr = 0.022
	plt.plot(datb,datr,'*',color='orange',markersize=14)
	plt.xlabel(r'$P$ RSD HGM $\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\perp}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	xl = [-1,2]
	yl = [-1,2]
	plt.plot(xl,yl,'r--')
	plt.show()
	plt.plot(eperpb,eperpr,'ko')
	plt.plot(detb,detr,'*',color='orange',markersize=14)
	plt.xlabel(r'$\xi$ RSD HGM $\sigma\alpha_{\perp}$')
	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{\perp}$')
	plt.xlim(0,0.05)
	plt.ylim(0,0.05)
	plt.plot(xl,yl,'r--')
	plt.show()
	
	aparb = bf[1][w]
	aparr = rf[0][w]
	eparb = bf[2][w]
	eparr = rf[1][w]
	print(np.mean(eperpb),np.mean(eperpr))
	print(np.std(aperpb),np.std(aperpr))
	print(np.mean(eparb),np.mean(eparr))
	print(np.std(aparb),np.std(aparr))
	plt.plot(aparb,aparr,'ko')
	darb = 0.994
	darr = 1.019
	derb = 0.036
	derr = 0.037

	plt.plot(darb,darr,'*',color='orange',markersize=14)
	plt.xlabel(r'$P$ RSD HGM $\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{||}$')
	plt.xlim(.8,1.2)
	plt.ylim(.8,1.2)
	plt.plot(xl,yl,'r--')
	plt.show()

	plt.plot(eparb,eparr,'ko')
	plt.plot(derb,derr,'*',color='orange',markersize=14)
	plt.xlabel(r'$P$ RSD HGM $\sigma\alpha_{||}$')
	plt.ylabel(r'$\xi$ RSD RP $\sigma\alpha_{||}$')
	plt.xlim(0,0.1)
	plt.ylim(0,0.1)
	plt.plot(xl,yl,'r--')
	plt.show()
	
	
	aisob = (aperpb**2.*aparb)**(1/3.)
	aisor = (aperpr**2.*aparr)**(1/3.)
	plt.plot(aisob,aisor,'ko')
	daisob = (datb**2.*darb)**(1/3.)
	daisor = (datr**2.*darr)**(1/3.)
	plt.plot(daisob,daisor,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$P$ RSD HGM $\alpha_{\rm iso}$')
	plt.ylabel(r'$\xi$ RSD RP $\alpha_{\rm iso}$')
	plt.show()

	epb = (aparb/aperpb)**(1/3.)
	epr = (aparr/aperpr)**(1/3.)
	plt.plot(epb,epr,'ko')
	depb = (darb/datb)**(1/3.)
	depr = (darr/datr)**(1/3.)
	plt.plot(depb,depr,'*',color='orange',markersize=14)
	plt.plot(xl,yl,'r--')
	plt.xlim(.9,1.1)
	plt.ylim(0.9,1.1)
	plt.xlabel(r'$P$ RSD HGM $\epsilon+1$')
	plt.ylabel(r'$\xi$ RSD RP $\epsilon+1$')
	plt.show()


	ht = plt.hist(abs(aisob-aisor))
	yl = [0,1000]
	dda = abs(daisob-daisor)
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\alpha_{\rm iso}|$(HGM RP)')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(ht[0])*1.1)
	plt.show()


	print(np.std(aisob),np.std(aisor))
	print(np.std(epb),np.std(epr))
	plt.hist(abs(epb-epr)/(eparr*eparb)**.5)
	hn = np.histogram(abs(epb-epr)/(eparr*eparb)**.5)
	yl = [0,1000]
	dda = abs(depb-depr)/sqrt(derr*derb)
	xl = [dda,dda]
	plt.plot(xl,yl,'r-')
	plt.xlabel(r'$|\Delta\epsilon|$/$\sqrt{\sigma_{||,HGM}\sigma_{||,RP}}$')
	plt.ylabel('number of EZ mocks')
	plt.ylim(0,np.max(hn[0])*1.1)
	plt.show()


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

def plotxirpp(reg,v='4',samp='ELG',bs=4,rmax=200,rtick=40,vmax=0.01,vmin=-0.01,zmin=.8,zmax=2.2,interp='spline16',subperp=False,wm='fkp',angfac=.75):
	d = np.loadtxt(ebossdir+'xirperprmugeboss'+samp+'_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+'.dat')#.transpose()
	xil = np.zeros((rmax/bs,rmax/bs))
	if subperp:
		ds = np.loadtxt(ebossdir+'xirprpgeboss'+samp+'_'+reg+v+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp1st0.dat').transpose()[1]
		for i in range(0,len(d)):
			d[i] -= angfac*ds[i]
# 		xip = np.zeros((rmax))
# 		for i in range(0,rmax):
# 			xiave = np.mean(d[i])
# 			xip[i] = xiave
# 		for i in range(0,rmax):
# 			for j in range(0,rmax):
# 				print(i,j,d[i][j])
# 				d[i][j] -= xip[i]
				#print(d[i][j],xip[i])			
		#print(xip)
	d = d.transpose()
	
	for i in range(0,rmax,bs):
		for j in range(0,rmax,bs):
			xi = 0
			nb = 0
			for bi in range(0,bs):
				for bj in range(0,bs):
					if xi != -1.:
						xi += d[i+bi][j+bj]
						nb += 1.
					#if subperp:
					#	xi -= xip[i+bi]
			xi = xi/nb#float(bs*bs)
			xil[i/bs][j/bs] = xi
	plt.imshow(xil,origin='lower',interpolation=interp,vmax=vmax,vmin=vmin)
	ax = plt.gca()
	ax.set_xticks(np.arange(0, rmax/bs, rtick/float(bs)))
	ax.set_xticklabels(np.arange(bs/2.,rmax+1,rtick))
	ax.set_yticks(np.arange(0, rmax/bs, rtick/float(bs)))
	ax.set_yticklabels(np.arange(bs/2.,rmax+1,rtick))
	plt.colorbar()
	plt.xlabel(r'$s_{\perp}h^{-1}$Mpc')
	plt.ylabel(r'$s_{||}h^{-1}$Mpc')
	if subperp:
		plt.title(reg+r' $\xi$'+'-'+str(angfac)+r'$w(r_{\perp}$)')
	if wm == 'fkp' and subperp == False:
		plt.title(reg+r' $\xi$')
	if wm == 'fkp.shuffled':
		plt.title(reg+r' $\xi_{\rm shuffled}$')		
	plt.show()
	sp = ''
	if subperp:
		sp = 'subperp'
	#plt.savefig('xirpp'+reg+sp+wm.strip('.')+'.png')
	plt.clf()
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

def plotxiELGcombNS(mom=0,minr=30,maxr=180,bs='8st0',v='test',rec='',zmin=.6,zmax=1.1,l1='',l2='',comb='',modplot=False):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcombNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .558
	d1 = np.loadtxt(ebossdir+'xi024gebossELG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossELG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .479
	dc = (d1*a1+d2*a2)/(a1+a2)
	rec = rec.strip('_')
	dmn = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGCELG_EZ'+rec+'angfac0'+bs+'.dat').transpose()
	dms = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGCELG_EZ'+rec+'angfac0'+bs+'.dat').transpose()
	errl = (1./(1./dmn[2]**2.+1./dms[2]**2.))**0.5
	sm = 0
	sx = 0
	for i in range(0,len(d1[0])):
		if d1[0][i] > minr and sm == 0:
			mini = i
			sm = 1
		if d1[0][i] > maxr and sx == 0:
			maxi = i
			sx = 1	
	
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
		plt.errorbar(dc[0][mini:maxi],dc[0][mini:maxi]**2.*dc[1][mini:maxi],dc[0][mini:maxi]**2.*errl[mini:maxi],fmt='o',color='steelblue')
		plt.plot(dc[0][mini:maxi],dc[0][mini:maxi]**2.*dc[1][mini:maxi],'ow',markersize=4,zorder=1000)
	if mom ==2:	
		plt.errorbar(dc[0][mini:maxi],dc[0][mini:maxi]*dc[2][mini:maxi],dc[0][mini:maxi]*errl[mini:maxi])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-20,30)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-1.5,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	if rec == 'rec':
		plt.title(r'post-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	else:
		plt.title(r'pre-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
# 	if l1 == '':
# 		if modplot:
# 			#if reg == 'ALL':
# 			plt.legend(labels=['model','NGC+SGC'])
# 			#if reg == 'NGC':
# 			#	plt.legend(labels=['model','eboss23+eboss25'])
# 			#if reg == 'SGC':
# 			#	plt.legend(labels=['model','eboss21+eboss22'])
# 		else:
# 			plt.legend(labels=[wm1,wm])
# 	else:
# 		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxiQSOcombNS(mom=0,minr=30,maxr=180,bs='5st0',v='test',rec='',zmin=.8,zmax=2.2,l1='',l2='',comb='',modplot=False):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSOcombNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .169
	d1 = np.loadtxt(ebossdir+'xi024gebossQSO_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossQSO_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .295
	dc = (d1*a1+d2*a2)/(a1+a2)
	rec = rec.strip('_')
	dmn = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGCQSO_EZ'+rec+bs+'.dat').transpose()
	dms = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGCQSO_EZ'+rec+bs+'.dat').transpose()
	errl = (1./(1./dmn[2]**2.+1./dms[2]**2.))**0.5
	sm = 0
	sx = 0
	for i in range(0,len(d1[0])):
		if d1[0][i] > minr and sm == 0:
			mini = i
			sm = 1
		if d1[0][i] > maxr and sx == 0:
			maxi = i
			sx = 1	
	
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
		plt.errorbar(dc[0][mini:maxi],dc[0][mini:maxi]**2.*dc[1][mini:maxi],dc[0][mini:maxi]**2.*errl[mini:maxi],fmt='s',color='forestgreen')
		plt.plot(dc[0][mini:maxi],dc[0][mini:maxi]**2.*dc[1][mini:maxi],'sw',markersize=4,zorder=1000)
	if mom ==2:	
		plt.errorbar(dc[0][mini:maxi],dc[0][mini:maxi]*dc[2][mini:maxi],dc[0][mini:maxi]*errl[mini:maxi])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-30,80)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-1.5,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	if rec == 'rec':
		plt.title(r'post-recon correlation function of ELGs '+' '+str(zmin)+' < z < '+str(zmax))
	else:
		plt.title(r'pre-recon correlation function of quasars '+' '+str(zmin)+' < z < '+str(zmax))
# 	if l1 == '':
# 		if modplot:
# 			#if reg == 'ALL':
# 			plt.legend(labels=['model','NGC+SGC'])
# 			#if reg == 'NGC':
# 			#	plt.legend(labels=['model','eboss23+eboss25'])
# 			#if reg == 'SGC':
# 			#	plt.legend(labels=['model','eboss21+eboss22'])
# 		else:
# 			plt.legend(labels=[wm1,wm])
# 	else:
# 		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True


def plotxi024ELGcombNS(mom=0,bs='8st0',v='5',rec='_rec',zmin=.6,zmax=1.1,l1='',l2='',comb='',modplot=True):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcombNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .558
	d1 = np.loadtxt(ebossdir+'xi024gebossELG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossELG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .479
	dc = (d1*a1+d2*a2)/(a1+a2)
	if rec == '_rec':
		#dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.59304.07.015.01.0.dat').transpose()
		dts = np.loadtxt(ebossdir+'xiave0SGCELG_EZrecangfac08st0.dat').transpose()
		dtn = np.loadtxt(ebossdir+'xiave0NGCELG_EZrecangfac08st0.dat').transpose()
		dt = (dts[1]*a1+dtn[1]*a2)/(a1+a2)
		et = (1./(1./dts[2]**2.+1./dtn[2]))**.5
	else:
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.563.06.010.015.00.dat').transpose()
		dts = np.loadtxt(ebossdir+'xiave0SGCELG_EZrecangfac08st0.dat').transpose()
		dtn = np.loadtxt(ebossdir+'xiave0NGCELG_EZrecangfac08st0.dat').transpose()
		dt = (dts[1]*a1+dtn[1]*a2)/(a1+a2)
		et = (1./(1./dts[2]**2.+1./dtn[2]))**.5
		
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
		plt.plot(dc[0],dc[0]*dc[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-25,45)
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


def plotxi024QSOcombNS(mom=0,bs='5st0',v='5',rec='',zmin=.8,zmax=2.2,l1='',l2='',comb='',modplot=True,modb=1.4):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSOcombNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a2 = .931935*6.
	d1 = np.loadtxt(ebossdir+'xi024gebossQSO_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossQSO_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a1 = .903293*4.
	dc = (d1*a1+d2*a2)/(a1+a2)
	dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.406.010.015.00.dat').transpose()
	
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*modb,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*modb,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		plt.plot(dc[0],dc[0]**2.*dc[1])
	if mom ==2:	
		plt.plot(dc[0],dc[0]*dc[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-50,100)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-4.,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'correlation function of quasars '+' '+str(zmin)+' < z < '+str(zmax))
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

def plotxi024LRGcombNS(mom=0,bs='5st0',v='5',rec='',zmin=.6,zmax=1.0,l1='',l2='',comb='',modplot=True,modb=2.2):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'LRGcombNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .1651
	d1 = np.loadtxt(ebossdir+'xi024gebossLRG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossLRG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .2508
	dc = (d1*a1+d2*a2)/(a1+a2)
	dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.406.010.015.00.dat').transpose()
	
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*modb,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*modb,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		plt.plot(dc[0],dc[0]**2.*dc[1])
	if mom ==2:	
		plt.plot(dc[0],dc[0]*dc[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-50,120)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-4.,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'correlation function of LRGs '+' '+str(zmin)+' < z < '+str(zmax))
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

def plotxi024LRG_jb(mom=0,min=30,max=180,bs='5st0',v='5',rec='rec',zmin=.6,zmax=1.0,l1='',l2='',comb='',modplot=False,modb=2.2):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'LRGcombNSjb'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dc = np.loadtxt(ebossdir+'eBOSS_LRGpCMASS_clustering_COMB_v5_rec_5mpc_shift0.mul').transpose()
	rl = []
	xil = []
	for i in range(0,len(dc[0])):
		r = dc[0][i]
		#print(r,min,max)
		if r > min and r < max:
			rl.append(r)
			if mom == 0:
				xil.append(dc[1][i])
			if mom == 2:
				xil.append(dc[2][i])
			#print(r,xil[i])		
	rl = np.array(rl)
	xil = np.array(xil)
	dcov = np.loadtxt(ebossdir+'eBOSS_LRGpCMASS_COMB_v5_recon_average_5mpc_shift0.cov').transpose()
	bs=5
	db = min//bs
	print(db)
	errl = []
	for i in range(0,len(dcov[0])):
		r1 = dcov[2][i]
		r2 = dcov[3][i]
		if r1 > min and r1 < max and r2 > min and r2 < max:
			ind1 = dcov[0][i]
			ind2 = dcov[1][i]
			#print(c1,c2)
			#cov[int(c1)][int(c2)] = d[4][i]
			if ind1 == ind2:
				if mom == 0:
					if ind1 < 39:
						errl.append(sqrt(dcov[4][i]))
				if mom == 2:
					if ind1 >= 39:
						errl.append(sqrt(dcov[4][i]))
	errl = np.array(errl)
	dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.406.010.015.00.dat').transpose()
	
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*xil*modb,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*xil*modb,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		#plt.plot(dc[0],dc[0]**2.*dc[1])
		plt.errorbar(rl,rl**2.*xil,rl**2.*errl,fmt='d',color='firebrick')
		plt.plot(rl,rl**2.*xil,'dw',markersize=4,zorder=1000)
	if mom ==2:	
		plt.errorbar(rl,rl*xil,rl**2.*errl,fmt='d',color='firebrick')

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-20,80)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-4.,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'post-rec correlation function of CMASS +LRGs '+' '+str(zmin)+' < z < '+str(zmax))
# 	if l1 == '':
# 		if modplot:
# 			#if reg == 'ALL':
# 			plt.legend(labels=['model','NGC+SGC'])
# 			#if reg == 'NGC':
# 			#	plt.legend(labels=['model','eboss23+eboss25'])
# 			#if reg == 'SGC':
# 			#	plt.legend(labels=['model','eboss21+eboss22'])
# 		else:
# 			plt.legend(labels=[wm1,wm])
# 	else:
# 		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True


def plotxi024LRGpCMASScombNS(mom=0,bs='5st0',v='5',rec='',zmin=.6,zmax=1.0,l1='',l2='',comb='',modplot=True,modb=2.2):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'LRGpCMASScombNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .5
	d1 = np.loadtxt(ebossdir+'xi024gebossLRGpCMASS_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossLRGpCMASS_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = 1.
	dc = (d1*a1+d2*a2)/(a1+a2)
	if rec == '':
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.406.010.015.00.dat').transpose()
	if rec == '_rec':
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.44.02.54.015.01.0.dat').transpose()
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*modb,'k:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*modb,'k:')	
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		plt.plot(dc[0],dc[0]**2.*dc[1])
	if mom ==2:	
		plt.plot(dc[0],dc[0]*dc[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-50,120)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-4.,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'correlation function of LRGs + CMASS'+' '+str(zmin)+' < z < '+str(zmax))
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

def plotxi024ELGcompNS(mom=0,bs='8st0',v='5',rec='',zmin=.6,zmax=1.1,l1='',l2='',modplot=False,comb=''):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcompNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = 1.2+3.1
	d1 = np.loadtxt(ebossdir+'xi024gebossELG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossELG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
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
		plt.plot(d2[0],d2[0]*d2[2])
		plt.plot(d1[0],d1[0]*d1[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-25,45)
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
			plt.legend(labels=['model','NGC','SGC'])
			#if reg == 'NGC':
			#	plt.legend(labels=['model','eboss23+eboss25'])
			#if reg == 'SGC':
			#	plt.legend(labels=['model','eboss21+eboss22'])
		else:
			plt.legend(labels=['NGC','SGC'])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxi024QSOcompNS(mom=0,bs='5st0',v='5',rec='',zmin=.8,zmax=2.2,l1='',l2='',modplot=False,comb=''):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'QSOcompNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .931935*6.
	d1 = np.loadtxt(ebossdir+'xi024gebossQSO_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossQSO_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .903293*4.
	dc = (d1*a1+d2*a2)/(a1+a2)
	dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.406.010.015.00.dat').transpose()
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
		plt.plot(d2[0],d2[0]*d2[2])
		plt.plot(d1[0],d1[0]*d1[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-50,100)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-4,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'correlation function of quasars '+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		if modplot:
			#if reg == 'ALL':
			plt.legend(labels=['model','NGC','SGC'])
			#if reg == 'NGC':
			#	plt.legend(labels=['model','eboss23+eboss25'])
			#if reg == 'SGC':
			#	plt.legend(labels=['model','eboss21+eboss22'])
		else:
			plt.legend(labels=['NGC','SGC'])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxi024LRGcompNS(mom=0,bs='5st0',v='5',rec='',zmin=.6,zmax=1.,l1='',l2='',modplot=False,comb=''):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'LRGcompNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .931935*6.
	d1 = np.loadtxt(ebossdir+'xi024gebossLRG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossLRG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .903293*4.
	dc = (d1*a1+d2*a2)/(a1+a2)
	dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.406.010.015.00.dat').transpose()
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
		plt.plot(d2[0],d2[0]*d2[2])
		plt.plot(d1[0],d1[0]*d1[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-50,120)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-4,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'correlation function of LRGs '+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		if modplot:
			#if reg == 'ALL':
			plt.legend(labels=['model','NGC','SGC'])
			#if reg == 'NGC':
			#	plt.legend(labels=['model','eboss23+eboss25'])
			#if reg == 'SGC':
			#	plt.legend(labels=['model','eboss21+eboss22'])
		else:
			plt.legend(labels=['NGC','SGC'])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxi024LRGpCMASScompNS(mom=0,bs='5st0',v='5',rec='',zmin=.6,zmax=1.,l1='',l2='',modplot=False,comb=''):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'LRGpCMASScompNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	a1 = .931935*6.
	d1 = np.loadtxt(ebossdir+'xi024gebossLRGpCMASS_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'xi024gebossLRGpCMASS_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+bs+'.dat').transpose()
	a2 = .903293*4.
	dc = (d1*a1+d2*a2)/(a1+a2)
	dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.406.010.015.00.dat').transpose()
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
		plt.plot(d2[0],d2[0]*d2[2])
		plt.plot(d1[0],d1[0]*d1[2])

	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-50,120)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-4,1)
	#plt.text(30,180,v,color='b')
	#plt.text(30,170,v2,color='r')
	plt.title(r'correlation function of LRGs + CMASS'+' '+str(zmin)+' < z < '+str(zmax))
	if l1 == '':
		if modplot:
			#if reg == 'ALL':
			plt.legend(labels=['model','NGC','SGC'])
			#if reg == 'NGC':
			#	plt.legend(labels=['model','eboss23+eboss25'])
			#if reg == 'SGC':
			#	plt.legend(labels=['model','eboss21+eboss22'])
		else:
			plt.legend(labels=['NGC','SGC'])
	else:
		plt.legend(labels=[l1,l2])
	pp.savefig()
	pp.close()
	return True

def plotxi_shuffsubcom(mom=0,reg='SGC',bs='8st0',v='4',rec='',zmin=.6,zmax=1.1,wm='',angfac=1.):
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages

	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+bs+'.dat').transpose()
	ds = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+'4_1'+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp.shuffled'+bs+'.dat').transpose()

	rl = d1[0]
	xil = d1[1]
	wp = np.loadtxt(ebossdir+'wrp024ELG'+reg+'_data'+bs+'.dat').transpose()
	ind = int(mom/2+1)
	xils = xil-angfac*wp[ind]
	pwr = 1.
	if mom == 0:
		pwr = 2.
	plt.plot(rl,rl**pwr*xils,'r-')
	plt.plot(rl,rl**pwr*ds[1],'b-')
	plt.plot(rl,rl**pwr*xil,'r--')
	plt.xlim(0,200)
	#plt.ylim(-75,75)
	xl = [130,136]
	tx =  138
	if mom == 0:
		xl = [12,18]
		tx =20
	ymin = plt.axis()[-2]
	ymax = plt.axis()[-1]

	yl = [.2*(ymax-ymin)+ymin,.2*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'r--')
	plt.text(tx,.19*(ymax-ymin)+ymin,'Standard')
	
	yl = [.15*(ymax-ymin)+ymin,.15*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'b-')
	plt.text(tx,.14*(ymax-ymin)+ymin,'Shuffled')
	
	yl = [.1*(ymax-ymin)+ymin,.1*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'r-')	
	plt.text(tx,.09*(ymax-ymin)+ymin,r'Standard -'+str(angfac)+r'$\xi(r_{\perp})$')
	plt.title(reg)
	plt.xlabel(r'$s$ $(h^{-1}$Mpc$)$')
	if pwr == 2:
		plt.ylabel(r'$s^2\xi_0$')
	if pwr == 1:
		if mom == 2:
			plt.ylabel(r'$s\xi_2$')
		if mom == 4:
			plt.ylabel(r'$s\xi_4$')	
	#plt.show()
	plt.savefig(ebossdir+'xi'+str(mom)+reg+'shuffsubcom.png')
	plt.clf()
	return True
	
	#plt.show()

def plotxi_submockcom(mom=0,reg='SGC',bs='8st0',mini=3,maxi=25,v='4',rec='',zmin=.6,zmax=1.1,wm='',angfac=1.,bf=1.):
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	recf = ''
	if rec == 'rec':
		recf = '_rec'
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+recf+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+bs+'.dat').transpose()
	ds = np.loadtxt(ebossdir+'xiave'+str(mom)+reg+'ELG_EZ'+rec+'angfac'+str(angfac)+str(bs)+'.dat').transpose()
	dsh = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+'4_1'+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp.shuffled'+bs+'.dat').transpose()
	mom = int(mom)
	#nb = len(ds[1])	
	rl = d1[0][mini:maxi]
	xil = d1[1]#[mini:maxi]
	wp = np.loadtxt(ebossdir+'wrp024ELG'+reg+'_data'+bs+'.dat').transpose()
	ind = int(mom/2+1)
	xils = xil-angfac*wp[ind]#[mini:maxi]
	xils = xils[mini:maxi]
	pwr = 1.
	if mom == 0:
		pwr = 2.
	#cov = np.loadtxt(ebossdir+'cov'+str(mom)+reg+'ELG_EZangfac'+str(angfac)+str(bs)+'.dat')[mini:maxi,mini:maxi]
	cov = np.loadtxt(ebossdir+'cov'+str(mom)+reg+'ELG_EZ'+rec+'angfac'+str(angfac)+str(bs)+'.dat')[mini:maxi,mini:maxi]#/2.
	diff = bf*ds[1][mini:maxi]-xils
	#diffsh = bf*dsh[1][mini:maxi]-ds[1][mini:maxi]
	facn = 1.
	facs = 1.

	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)*facs
	#chish = np.dot(np.dot(diffsh,np.linalg.pinv(cov)),diffsh)*facs
	print(chi)
	plt.plot(rl,rl**pwr*xils,'r-')
	#plt.plot(rl,rl**pwr*ds[1][mini:maxi],'b-')
	plt.errorbar(rl,rl**pwr*ds[1][mini:maxi]*bf,rl**pwr*ds[2][mini:maxi],fmt='b-')
	plt.plot(rl,rl**pwr*dsh[1][mini:maxi],'r--')
	plt.xlim(0,200)
	#plt.ylim(-75,75)
	xl = [130,136]
	tx =  138
	if mom == 0:
		xl = [12,18]
		tx =20
	ymin = plt.axis()[-2]
	ymax = plt.axis()[-1]

	yl = [.2*(ymax-ymin)+ymin,.2*(ymax-ymin)+ymin]
	#plt.plot(xl,yl,'r--')
	#plt.text(tx,.19*(ymax-ymin)+ymin,'Standard')
	
	yl = [.15*(ymax-ymin)+ymin,.15*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'b-')
	plt.text(tx,.14*(ymax-ymin)+ymin,'mocks-'+str(angfac)+r'$\xi(r_{\perp})$')
	
	yl = [.1*(ymax-ymin)+ymin,.1*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'r-')	
	plt.text(tx,.09*(ymax-ymin)+ymin,r'data -'+str(angfac)+r'$\xi(r_{\perp})$')
	plt.title(reg)
	plt.xlabel(r'$s$ $(h^{-1}$Mpc$)$')
	if pwr == 2:
		plt.ylabel(r'$s^2\xi_0$')
	if pwr == 1:
		if mom == 2:
			plt.ylabel(r'$s\xi_2$')
		if mom == 4:
			plt.ylabel(r'$s\xi_4$')	
	#plt.show()
	plt.savefig(ebossdir+'xi'+str(mom)+reg+rec+'mocksub'+str(angfac)+'com.png')
	plt.clf()
	return True

def plotxi_shuffmockcom(mom=0,reg='SGC',bs='8st0',mini=3,maxi=25,v='4',rec='',zmin=.6,zmax=1.1,wm='',angfac=.6,bf=1.):
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	recf = ''
	if rec == 'rec':
		recf = '_rec'
	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+recf+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+bs+'.dat').transpose()
	ds = np.loadtxt(ebossdir+'xiave'+str(mom)+reg+'ELG_EZshuffangfac0'+str(bs)+'.dat').transpose()
	dsh = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+'4_1'+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp.shuffled'+bs+'.dat').transpose()
	mom = int(mom)
	#nb = len(ds[1])	
	rl = d1[0][mini:maxi]
	xil = d1[1]#[mini:maxi]
	wp = np.loadtxt(ebossdir+'wrp024ELG'+reg+'_data'+bs+'.dat').transpose()
	ind = int(mom/2+1)
	xils = xil-angfac*wp[ind]#[mini:maxi]
	xils = xils[mini:maxi]
	pwr = 1.
	if mom == 0:
		pwr = 2.
	#cov = np.loadtxt(ebossdir+'cov'+str(mom)+reg+'ELG_EZangfac'+str(angfac)+str(bs)+'.dat')[mini:maxi,mini:maxi]
	cov = np.loadtxt(ebossdir+'cov'+str(mom)+reg+'ELG_EZshuffangfac0'+str(bs)+'.dat')[mini:maxi,mini:maxi]#/2.
	diff = bf*ds[1][mini:maxi]-dsh[1][mini:maxi]
	#diffsh = bf*dsh[1][mini:maxi]-ds[1][mini:maxi]
	facn = 1.
	facs = 1.

	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)*facs
	#chish = np.dot(np.dot(diffsh,np.linalg.pinv(cov)),diffsh)*facs
	print(chi)
	plt.plot(rl,rl**pwr*dsh[1][mini:maxi],'r-')
	#plt.plot(rl,rl**pwr*ds[1][mini:maxi],'b-')
	plt.errorbar(rl,rl**pwr*ds[1][mini:maxi]*bf,rl**pwr*ds[2][mini:maxi],fmt='b-')
	plt.plot(rl,rl**pwr*xils,'r--')
	plt.xlim(0,200)
	#plt.ylim(-75,75)
	xl = [130,136]
	tx =  138
	if mom == 0:
		xl = [12,18]
		tx =20
	ymin = plt.axis()[-2]
	ymax = plt.axis()[-1]

	yl = [.2*(ymax-ymin)+ymin,.2*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'r--')
	plt.text(tx,.19*(ymax-ymin)+ymin,'Standard data -'+str(angfac)+r'$\xi_{\perp}$')
	
	yl = [.15*(ymax-ymin)+ymin,.15*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'b-')
	plt.text(tx,.14*(ymax-ymin)+ymin,'shuffled mocks')
	
	yl = [.1*(ymax-ymin)+ymin,.1*(ymax-ymin)+ymin]
	plt.plot(xl,yl,'r-')	
	plt.text(tx,.09*(ymax-ymin)+ymin,r'shuffled data')
	plt.title(reg)
	plt.xlabel(r'$s$ $(h^{-1}$Mpc$)$')
	if pwr == 2:
		plt.ylabel(r'$s^2\xi_0$')
	if pwr == 1:
		if mom == 2:
			plt.ylabel(r'$s\xi_2$')
		if mom == 4:
			plt.ylabel(r'$s\xi_4$')	
	#plt.show()
	plt.savefig(ebossdir+'xi'+str(mom)+reg+rec+'mockshuffcom.png')
	plt.clf()
	return True


def plotxi_mockcomth(mom=0,sample='ELG',reg='SGC',bs='8st0',mini=3,maxi=25,damp='',v='4',rec='',angfac=0,bf=1.):
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	ds = np.loadtxt(ebossdir+'xiave'+str(mom)+reg+sample+'_EZ'+rec+'angfac'+str(angfac)+str(bs)+'.dat').transpose()
	dt = np.loadtxt('/Users/ashleyross/Github/LSSanalysis/BAOtemplates/xi'+str(mom)+'Challenge_matterpower'+damp+'.dat').transpose()
	
	#nb = len(ds[1])	
	rl = ds[0][mini:maxi]
	xil = ds[1]#[mini:maxi]
	#wp = np.loadtxt(ebossdir+'wrp024ELG'+reg+'_data'+bs+'.dat').transpose()
	ind = int(mom/2+1)
	xils = xil#-angfac*wp[ind]#[mini:maxi]
	xils = xils[mini:maxi]
	pwr = 1.
	if mom == 0:
		pwr = 2.
	plt.plot(rl,rl**pwr*xils,'r-')
	plt.plot(dt[0],dt[0]**pwr*dt[1]*bf,'b-')
	plt.xlim(0,200)
	#plt.ylim(-75,75)
	xl = [130,136]
	tx =  138
	if mom == 0:
		xl = [12,18]
		tx =20
	ymin = plt.axis()[-2]
	ymax = plt.axis()[-1]

	yl = [.2*(ymax-ymin)+ymin,.2*(ymax-ymin)+ymin]
	#plt.plot(xl,yl,'r--')
	plt.title(reg)
	plt.xlabel(r'$s$ $(h^{-1}$Mpc$)$')
	if pwr == 2:
		plt.ylabel(r'$s^2\xi_0$')
	if pwr == 1:
		if mom == 2:
			plt.ylabel(r'$s\xi_2$')
		if mom == 4:
			plt.ylabel(r'$s\xi_4$')	
	plt.show()
	#plt.savefig(ebossdir+'xi'+str(mom)+reg+rec+'mocksub'+str(angfac)+'com.png')
	plt.clf()
	return True


def plotxi_shuffsubcom_mubin(mom=0,reg='SGC',bs='8st0',v='4',rec='',zmin=.6,zmax=1.1,wm='',angfac=.5,mumin=0,mumax=1):
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	plt.clf()
	gf = ''
	if mumin != 0:
		gf += 'mum'+str(mumin)
	if mumax != 1.:
		gf += 'mux'+str(mumax)

	d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+gf+bs+'.dat').transpose()
	ds = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkpnosysshuff'+gf+bs+'.dat').transpose()

	rl = d1[0]
	xil = d1[1]
	wp = np.loadtxt(ebossdir+'wrp0'+gf+'ELG'+reg+'_data'+bs+'.dat').transpose()
	xils = xil-angfac*wp[mom/2+1]
	pwr = 1.
	if mom == 0:
		pwr = 2.
	plt.plot(rl,rl**pwr*xils,'r-')
	plt.plot(rl,rl**pwr*ds[1],'b-')
	plt.plot(rl,rl**pwr*xil,'r--')
	plt.xlim(0,200)
	plt.ylim(-75,75)
	xl = [12,18]
	yl = [-40,-40]
	plt.plot(xl,yl,'r--')
	plt.text(20,-42,'Standard')
	xl = [12,18]
	yl = [-48,-48]
	plt.plot(xl,yl,'b-')
	plt.text(20,-50,'Shuffled')
	xl = [12,18]
	yl = [-56,-56]
	plt.plot(xl,yl,'r-')	
	plt.text(20,-58,r'Standard -'+str(angfac)+r'$\xi(r_{\perp})$')
	plt.title(str(mumin)+r'$<\mu<$'+str(mumax))
	plt.xlabel(r'$s$ $(h^{-1}$Mpc$)$')
	plt.ylabel(r'$s^2\xi$')
	plt.show()
	#plt.savefig(ebossdir+'xishuffsubcom'+gf+'.png')
	return True
	


def plotxiELGcompth(mom=0,reg='SGC',bs='8st0',v='4',rec='',zmin=.6,zmax=1.1,thfac=.75,l1='',l2='',wm='',modplot=True,comb='',angfac=.5,md='subang'):
	#Plots comparison between NGC and SGC clustering and to theory for QSOs, no depth density correction
	from matplotlib import pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xi'+str(mom)+'ELGcompNS'+v+rec+'mz'+str(zmin)+'xz'+str(zmax)+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#d1 = np.loadtxt(ebossdir+'xi0gebosselg_'+chunk+v+'_mz'+str(zmin)+'xz'+str(zmax)+wm+bs+'.dat').transpose()
	if reg == 'NScomb':
		aves = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGCELG_EZv'+v+bs+'.dat').transpose()
		aven = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGCELG_EZv'+v+bs+'.dat').transpose()
		ave = aves
		ave[2] = (1./(1./aves[2]**2.+1./aven[2]**2.))**.5
		dn = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_NGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+bs+'.dat').transpose()
		ds = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_SGC'+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+bs+'.dat').transpose()
		rl = dn[0][:len(ave[0])]
		xiln = dn[1][:len(ave[0])]
		xils = ds[1][:len(ave[0])]
		
		if md == 'subang':
			print(('subtracting angular clustering with factor '+str(angfac)))
			wpn = np.loadtxt(ebossdir+'wrp024ELGNGC_data'+bs+'.dat').transpose()
			wps = np.loadtxt(ebossdir+'wrp024ELGSGC_data'+bs+'.dat').transpose()
			xilno = xiln
			xilso = xils
			xiln = xiln-angfac*wpn[mom/2+1][:len(ave[0])]
			xils = xils-angfac*wps[mom/2+1][:len(ave[0])]
		xilo = (xilso/aves[2]**2+xilno/aven[2]**2.)/(1./aves[2]**2.+1./aven[2]**2.)	
		xil = (xils/aves[2]**2+xiln/aven[2]**2.)/(1./aves[2]**2.+1./aven[2]**2.)
	else:
		if mom == 4:
			ave = np.loadtxt(ebossdir+'xiave0'+reg+'ELG_EZv'+v+bs+'.dat').transpose()*sqrt(9.) #error on hexadecapole is ~sqrt(9)x error on monopole
		else:	
			ave = np.loadtxt(ebossdir+'xiave'+str(mom)+reg+'ELG_EZv'+v+bs+'.dat').transpose()
		d1 = np.loadtxt(ebossdir+'xi'+str(mom)+'gebossELG_'+reg+comb+v+rec+'_mz'+str(zmin)+'xz'+str(zmax)+'fkp'+wm+bs+'.dat').transpose()
		rl = d1[0][:len(ave[0])]
		xil = d1[1][:len(ave[0])]

		if md == 'subang':
			print(('subtracting angular clustering with factor '+str(angfac)))
			wp = np.loadtxt(ebossdir+'wrp024ELG'+reg+'_data'+bs+'.dat').transpose()
			xilo = xil
			xil = xil-angfac*wp[mom/2+1][:len(ave[0])]
	if rec == '_rec':
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.593.02.04.015.01.0.dat').transpose()
	else:
		dt = np.loadtxt('BAOtemplates/xi'+str(mom)+'Challenge_matterpower0.563.04.07.015.00.dat').transpose()
	if wm == 'nosysshuff' or md == 'subang':
		dto = dt
		dt = dt - angfac*np.loadtxt('/Users/ashleyross/eBOSS/wrp024ELGSGC.dat').transpose()[mom/2+1]/.655**2.
	if modplot:
		if mom == 0:
			plt.plot(dt[0],dt[0]**2.*dt[1]*thfac,'k:')
			if md == 'subang':
				plt.plot(dt[0],dt[0]**2.*dto[1]*thfac,'r:')
		if mom == 2:
			plt.plot(dt[0],dt[0]*dt[1]*thfac,'k:')	
			if md == 'subang':
				plt.plot(dt[0],dt[0]*dto[1]*thfac,'r:')

		if mom == 4:
			plt.plot(dt[0],dt[0]*dt[1]*thfac,'k:')		
	#plt.plot(d1[0],d1[0]**2.*(d1[1]))
	#plt.plot(d2[0],d2[0]**2.*(d2[1]))
	if mom ==0:	
		
		plt.errorbar(rl,rl**2.*xil,rl**2.*ave[2],fmt='ko')
		if md == 'subang':
			plt.plot(rl,rl**2.*xilo,'r-')
	if mom ==2:	
		
		plt.errorbar(rl,rl*xil,rl*ave[2],fmt='ko')
		plt.plot(rl,rl*xilo,'r-')
	if mom ==4:	
		
		plt.errorbar(rl,rl*xil,rl*ave[2],fmt='ko')


	plt.xlim(10,200)
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylim(-25*thfac/.75,50*thfac/.75)
		plt.ylabel(r'$s^2\xi_0(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_2(s)$ ($h^{-1}$Mpc)',size=16)
		plt.ylim(-2.5,1.5)
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

def plotximuEZ(samp='ELG',covv='v7',rec='_rec',bs='5st0',mom=0):
	aves0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+rec+'mumax0.2'+bs+'.dat').transpose()
	aven0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+rec+'mumax0.2'+bs+'.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'k-')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'k--')
	aves0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.2mumax0.4'+bs+'.dat').transpose()
	aven0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.2mumax0.4'+bs+'.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'r-')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'r--')
	aves0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.4mumax0.6'+bs+'.dat').transpose()
	aven0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.4mumax0.6'+bs+'.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'b-')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'b--')
	aves0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.6mumax0.8'+bs+'.dat').transpose()
	aven0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.6mumax0.8'+bs+'.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'-',color='purple')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'--',color='purple')
	aves0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.8'+bs+'.dat').transpose()
	aven0 = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+rec+'mumin0.8'+bs+'.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'-',color='orange')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'--',color='orange')
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s^2\xi(\mu) $')
	plt.title('mean of post-recon EZ mocks in 0.2 mu bins')
	plt.show()
	tempdir = 'BAOtemplates/'
	aven0 = np.loadtxt(tempdir+'ximumax0.2Challenge_matterpower0.59304.07.015.01.0.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'k-')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'k--')
	aven0 = np.loadtxt(tempdir+'ximumin0.2mumax0.4Challenge_matterpower0.59304.07.015.01.0.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'r-')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'r--')
	aven0 = np.loadtxt(tempdir+'ximumin0.4mumax0.6Challenge_matterpower0.59304.07.015.01.0.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'b-')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'b--')
	aven0 = np.loadtxt(tempdir+'ximumin0.6mumax0.8Challenge_matterpower0.59304.07.015.01.0.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'-',color='purple')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'--',color='purple')
	aven0 = np.loadtxt(tempdir+'ximumin0.8Challenge_matterpower0.59304.07.015.01.0.dat').transpose()
	plt.plot(aven0[0],aven0[0]**2.*aven0[1],'-',color='orange')
	#plt.plot(aven0[0],aven0[0]**2.*aves0[1],'--',color='orange')
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s^2\xi(\mu) $')
	plt.title( 'post-recon BAO template in 0.2 mu bins')
	plt.show()

def plotxiNScompEZ(mom=0,samp='ELG',bs='5st0',v='7',mini=4,maxi=40,rec='',wm='fkp',mumin=0,mumax=1,bfs=1.,bfn=1.,tit='',covv='v7',zr='0.6xz1.1',mur='',angfac=0):
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
	pp = PdfPages(ebossdir+'xi'+str(mom)+samp+rec+'NScompEZ'+covv+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	#dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	if samp == 'ELG':
		ds = np.loadtxt(ebossdir+'xi024'+samp+'_'+covv+'SGCCZdata'+rec+bs+'.dat').transpose()[int(1+mom/2)]
		dn = np.loadtxt(ebossdir+'xi024'+samp+'_'+covv+'NGCCZdata'+rec+bs+'.dat').transpose()[int(1+mom/2)]
		colors = ['midnightblue','steelblue']
		ylim0 = [-60,50]
		ylim2 = [-2.,1.]
	if samp == 'QSO':
		ds = np.loadtxt(ebossdir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[int(1+mom/2)]
		dn = np.loadtxt(ebossdir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[int(1+mom/2)]
		colors = ['g','goldenrod']
		ylim0 = [-60,80]
		ylim2 = [-2.,1.]
		if maxi > 32:
			print('data xi only goes to 160 mpc/h!')
	if samp == 'LRGpCMASS':
		ds = np.loadtxt(ebossdir+'eBOSS_LRGpCMASS_clustering_SGC_v7_5mpc_shift0.mul').transpose()[int(1+mom/2)]
		dn = np.loadtxt(ebossdir+'eBOSS_LRGpCMASS_clustering_NGC_v7_5mpc_shift0.mul').transpose()[int(1+mom/2)]
		colors = ['firebrick','salmon']
		ylim0 = [-40,110]
		ylim2 = [-3.,1.]
		if maxi > 39:
			print('data xi only goes to 195 mpc/h!')

	aves = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+rec+bs+'.dat').transpose()
	aven = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+rec+bs+'.dat').transpose()
	covs = np.loadtxt(ebossdir+'cov'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+rec+bs+'.dat')[mini:maxi,mini:maxi]#*.857/1.4
	covn = np.loadtxt(ebossdir+'cov'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+rec+bs+'.dat')[mini:maxi,mini:maxi]#/2.
	#print(mini,maxi)
	#print(aves[1])
	diffs = bfs*aves[1][mini:maxi]*norm**2.-ds[mini:maxi]
	facn = 1.
	facs = 1.

	chis = np.dot(np.dot(diffs,np.linalg.pinv(covs)),diffs)*facs
	diffn = bfn*aven[1][mini:maxi]*norm**2.-dn[mini:maxi]
	chin = np.dot(np.dot(diffn,np.linalg.pinv(covn)),diffn)*facn
	diff = ds[mini:maxi]-dn[mini:maxi]
	cov = covn+covs
	chi = np.dot(np.dot(diff,np.linalg.pinv(cov)),diff)
	print(( chis,chin,chi))
	rl = aven[0][mini:maxi]
	if mom == 0:
		plt.plot(rl,rl**2.*aven[1][mini:maxi]*norm**2.,'--',color=colors[0])
		plt.plot(rl,rl**2.*aves[1][mini:maxi]*norm**2.,'--',color=colors[1])
		plt.errorbar(rl-.5,rl**2.*ds[mini:maxi],rl**2.*aves[2][mini:maxi]*norm**2.,fmt='s',color=colors[1],label=r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini))
		plt.errorbar(rl+.5,rl**2.*dn[mini:maxi],rl**2.*aven[2][mini:maxi]*norm**2.,fmt='o',markerfacecolor='w',color=colors[0],label=r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini))
		plt.legend(loc='lower left')
		plt.plot(rl+.5,rl**2.*dn[mini:maxi],'wo',markersize=4,zorder=100)
	if mom == 2:
		plt.plot(rl,rl*aven[1][mini:maxi]*norm**2.,'--',color=colors[0])
		plt.plot(rl,rl*aves[1][mini:maxi]*norm**2.,'--',color=colors[1])
		plt.errorbar(rl-.5,rl*ds[mini:maxi],rl*aves[2][mini:maxi]*norm**2.,fmt='s',color=colors[1])
		plt.errorbar(rl+.5,rl*dn[mini:maxi],rl*aven[2][mini:maxi]*norm**2.,fmt='o',color=colors[0])
		plt.plot(rl+.5,rl*dn[mini:maxi],'wo',markersize=4,zorder=100)
	
	plt.xlim(rl[0]-2.,rl[-1]+2.)
	if mom == 0:
		plt.ylim(ylim0[0],ylim0[1])
	else:
		plt.ylim(ylim2[0],ylim2[1])
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylabel(r'$s^2\xi_{'+str(mom)+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_{'+str(mom)+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
# 	if mom == 0:
# 		if samp != "LRGpCMASS":
# 			plt.text(30,-47,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color=colors[1])
# 			plt.text(30,-42,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color=colors[0])
# 		else:
# 			plt.text(30,-27,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color=colors[1])
# 			plt.text(30,-20,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color=colors[0])
# 			
# 	else:
# 		plt.text(30,.7,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color=colors[1])
# 		plt.text(30,.5,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color=colors[0])
	
	#plt.text(30,160,'Combined',color='k')
	if tit == '':
		if samp == 'LRGpCMASS':
			plt.title('DR16 LRGs+CMASS')
		else:	
			plt.title('DR16 '+samp+'s')
	else:
		plt.title(tit)		
	pp.savefig()
	pp.close()

	return True

def plotxiprepostcompEZ(mom=0,samp='ELG',bs='5st0',v='7',mini=4,maxi=40,wm='fkp',mumin=0,mumax=1,bfs=1.,bfn=1.,tit='NGC+SGC',covv='v7',zr='0.6xz1.1',mur='',angfac=0):
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
	pp = PdfPages(ebossdir+'xi'+str(mom)+samp+'prepostcompEZ'+covv+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	#ds = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_S'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	#dn = np.loadtxt(ebossdir+'xi'+mom+'gebossQSO_N'+v+'_mz0.8xz2.2'+wm+muwd+bs+'.dat').transpose()
	if samp == 'ELG':
		dsr = np.loadtxt(ebossdir+'xi024'+samp+'_'+covv+'SGCCZdata_rec'+bs+'.dat').transpose()[int(1+mom/2)]
		dnr = np.loadtxt(ebossdir+'xi024'+samp+'_'+covv+'NGCCZdata_rec'+bs+'.dat').transpose()[int(1+mom/2)]
		dsn = np.loadtxt(ebossdir+'xi024'+samp+'_'+covv+'SGCCZdata'+bs+'.dat').transpose()[int(1+mom/2)]
		dnn = np.loadtxt(ebossdir+'xi024'+samp+'_'+covv+'NGCCZdata'+bs+'.dat').transpose()[int(1+mom/2)]
		colors = ['dodgerblue','blue']
		msize = 7
		ylim0 = [-60,50]
		ylim2 = [-2.,1.]
	if samp == 'QSO':
		ds = np.loadtxt(ebossdir+'dr16_qso_obs_v7_sgc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[int(1+mom/2)]
		dn = np.loadtxt(ebossdir+'dr16_qso_obs_v7_ngc_0.31_z0.8x2.2_xi4m_dwfiber_rwfull_ds5[1].dat').transpose()[int(1+mom/2)]
		colors = ['g','goldenrod']
		ylim0 = [-60,80]
		ylim2 = [-2.,1.]
		if maxi > 32:
			print('data xi only goes to 160 mpc/h!')
	if samp == 'LRGpCMASS':
		ds = np.loadtxt(ebossdir+'eBOSS_LRGpCMASS_clustering_SGC_v7_5mpc_shift0.mul').transpose()[int(1+mom/2)]
		dn = np.loadtxt(ebossdir+'eBOSS_LRGpCMASS_clustering_NGC_v7_5mpc_shift0.mul').transpose()[int(1+mom/2)]
		colors = ['firebrick','salmon']
		ylim0 = [-40,110]
		ylim2 = [-3.,1.]
		if maxi > 39:
			print('data xi only goes to 195 mpc/h!')

	avesr = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ_rec'+bs+'.dat').transpose()
	avenr = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ_rec'+bs+'.dat').transpose()
	avetr = (avesr[1]/avesr[2]**2.+avenr[1]/avenr[2]**2.)/(1./avesr[2]**2.+1./avenr[2]**2.)
	etr = (1./(1./avesr[2]**2.+1./avenr[2]**2.))**.5

	avesn = np.loadtxt(ebossdir+'xiave'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+bs+'.dat').transpose()
	avenn = np.loadtxt(ebossdir+'xiave'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+bs+'.dat').transpose()
	avetn = (avesn[1]/avesn[2]**2.+avenn[1]/avenn[2]**2.)/(1./avesn[2]**2.+1./avenn[2]**2.)
	etn = (1./(1./avesn[2]**2.+1./avenn[2]**2.))**.5

	dtr = (dsr[mini:maxi]/avesr[2][mini:maxi]**2.+dnr[mini:maxi]/avenr[2][mini:maxi]**2.)/(1./avenr[2][mini:maxi]**2.+1./avesr[2][mini:maxi]**2.)
	dtn = (dsn[mini:maxi]/avesn[2][mini:maxi]**2.+dnn[mini:maxi]/avenn[2][mini:maxi]**2.)/(1./avenn[2][mini:maxi]**2.+1./avesn[2][mini:maxi]**2.)

	covsr = np.loadtxt(ebossdir+'cov'+str(mom)+'SGC'+samp+'_'+covv+'_EZ_rec'+bs+'.dat')[mini:maxi,mini:maxi]#*.857/1.4
	covnr = np.loadtxt(ebossdir+'cov'+str(mom)+'NGC'+samp+'_'+covv+'_EZ_rec'+bs+'.dat')[mini:maxi,mini:maxi]#/2.
	icovsr = np.linalg.pinv(covsr)
	icovnr = np.linalg.pinv(covnr)
	icovtr = icovsr+icovnr

	covsn = np.loadtxt(ebossdir+'cov'+str(mom)+'SGC'+samp+'_'+covv+'_EZ'+bs+'.dat')[mini:maxi,mini:maxi]#*.857/1.4
	covnn = np.loadtxt(ebossdir+'cov'+str(mom)+'NGC'+samp+'_'+covv+'_EZ'+bs+'.dat')[mini:maxi,mini:maxi]#/2.
	icovsn = np.linalg.pinv(covsn)
	icovnn = np.linalg.pinv(covnn)
	icovtn = icovsn+icovnn

	#print(mini,maxi)
	#print(aves[1])
	diffr = bfs*avetr[mini:maxi]*norm**2.-dtr#[mini:maxi]
	facn = 1.
	facs = 1.

	chir = np.dot(np.dot(diffr,icovtr),diffr)*facs
	diffn = bfn*avetn[mini:maxi]*norm**2.-dtn#[mini:maxi]
	chin = np.dot(np.dot(diffn,icovtn),diffn)*facn
	print(( chir,chin))
	rl = avenn[0][mini:maxi]
	if mom == 0:
		plt.plot(rl,rl**2.*avetn[mini:maxi]*norm**2.,'--',color=colors[1])
		plt.plot(rl,rl**2.*avetr[mini:maxi]*norm**2.,'--',color=colors[0])
		plt.errorbar(rl-.5,rl**2.*dtn,rl**2.*etn[mini:maxi]*norm**2.,fmt='v',color=colors[1],markersize=msize,label=r'pre-rec, $\chi^2=$'+str(chin)[:4]+'/'+str(maxi-mini))
		plt.errorbar(rl+.5,rl**2.*dtr,rl**2.*etr[mini:maxi]*norm**2.,fmt='^',color=colors[0],markeredgecolor='darkblue',markersize=msize,label='post-rec, $\chi^2=$'+str(chir)[:4]+'/'+str(maxi-mini))
		#plt.plot(rl+.5,rl**2.*dn[mini:maxi],'wo',markersize=4,zorder=100)
	if mom == 2:
		plt.plot(rl,rl*aven[1][mini:maxi]*norm**2.,'--',color=colors[0])
		plt.plot(rl,rl*aves[1][mini:maxi]*norm**2.,'--',color=colors[1])
		plt.errorbar(rl-.5,rl*ds[mini:maxi],rl*aves[2][mini:maxi]*norm**2.,fmt='s',color=colors[1])
		plt.errorbar(rl+.5,rl*dn[mini:maxi],rl*aven[2][mini:maxi]*norm**2.,fmt='o',color=colors[0])
		plt.plot(rl+.5,rl*dn[mini:maxi],'wo',markersize=4,zorder=100)
	
	plt.xlim(rl[0]-2.,rl[-1]+2.)
	if mom == 0:
		plt.ylim(ylim0[0],ylim0[1])
	else:
		plt.ylim(ylim2[0],ylim2[1])
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16)
	if mom == 0:
		plt.ylabel(r'$s^2\xi_{'+str(mom)+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	if mom == 2:
		plt.ylabel(r'$s\xi_{'+str(mom)+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
# 	if mom == 0:
# 		if samp != "LRGpCMASS":
# 			plt.text(30,-47,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color=colors[1])
# 			plt.text(30,-42,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color=colors[0])
# 		else:
# 			plt.text(30,-27,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color=colors[1])
# 			plt.text(30,-20,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color=colors[0])
# 			
# 	else:
# 		plt.text(30,.7,r'SGC, $\chi^2$/dof ='+str(chis)[:4]+'/'+str(maxi-mini),color=colors[1])
# 		plt.text(30,.5,r'NGC, $\chi^2$/dof ='+str(chin)[:4]+'/'+str(maxi-mini),color=colors[0])
	
	#plt.text(30,160,'Combined',color='k')
	plt.legend(loc='lower left')
	if tit == '':
		if samp == 'LRGpCMASS':
			plt.title('DR16 LRGs+CMASS')
		else:	
			plt.title('DR16 '+samp+'s')
	else:
		plt.title(tit)		
	pp.savefig()
	pp.close()

	return True


def plotxicompEZOR(mom=0,reg='SGC',samp='ELG',bs='5st0',v='7',mini=4,maxi=30,rec='_rec',wm='fkp',mumin=0,mumax=1,bfs=1.,bfn=1.,covv='v7',zr='0.6xz1.1',mur='',angfac=0):
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
	pp = PdfPages(ebossdir+'xi'+str(mom)+samp+rec+reg+'compEZOR'+covv+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()

	aveez = np.loadtxt(ebossdir+'xiave'+str(mom)+reg+samp+'_'+covv+'_EZ'+rec+bs+'.dat').transpose()
	aveor = np.loadtxt(ebossdir+'xiave024'+samp+'ORmock'+reg+rec+bs+'.dat').transpose()[1+mom//2]
	s, m0, e0, m2, e2, m4, e4, xi0, xi2, xi4 = np.load(ebossdir+'2PCF_EZmock_eBOSS_ELG_rec_comp.npz', allow_pickle=True)['data']

	facn = 1.
	facs = 1.
	
	colors = ['cadetblue','navy']

	rl = aveez[0][mini:maxi]
	if mom == 0:
		plt.plot(rl,rl**2.*aveez[1][mini:maxi]*norm**2.,'-',color=colors[0],label='EZmocks')
		plt.plot(rl*0.942,(rl*.942)**2.*aveor[mini:maxi]*norm**2.,'--',color=colors[1],label='OuterRim')
		plt.plot(rl,rl**2.*m0[1][mini:maxi]*norm**2.,':',color=colors[0],label='EZmocks, no sys')
	
	#plt.xlim(rl[0]-2.,rl[-1]+2.)
	#ylim0 = [-60,50]
	
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16,labelpad=-.75)
	plt.ylabel(r'$s^2\xi_{'+str(mom)+'}(s)$ ($h^{-2}$Mpc$^{2}$)',size=16)
	plt.title('Post-reconstruction')
	plt.legend()
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
	print((chins,maxi-mini,ds[0][mini],ds[0][maxi]))
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
	print((chis,chin,chi))
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

def plotxiELGNSbaofit(bs='5st0',v='4',a='',rec='',wm='fkp',mini=10,maxi=30,mom='0',covv='',angfac=0):
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
# 	if rec == '_rec' or rec == 'rec':
# 		covs = np.loadtxt(ebossdir+'cov_recon'+mom+'SGCELG_EZ'+bs+'.dat')
# 		covn = np.loadtxt(ebossdir+'cov_recon'+mom+'NGCELG_EZ'+bs+'.dat')
# 	if rec == '':	
	covs = np.loadtxt(ebossdir+'cov'+mom+'SGCELG_EZ'+rec+'angfac'+str(angfac)+bs+'.dat')
	covn = np.loadtxt(ebossdir+'cov'+mom+'NGCELG_EZ'+rec+'angfac'+str(angfac)+bs+'.dat')
	
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
	recf = ''
	if rec == 'rec':
		recf = '_rec'	
	dsw = np.loadtxt(ebossdir+'xi0gebossELG_SGC'+v+recf+'_mz0.6xz1.1'+wm+bs+'.dat').transpose()
	dsw_iso = dsw[1][mini:maxi]-dts[2]
	dnw = np.loadtxt(ebossdir+'xi0gebossELG_NGC'+v+recf+'_mz0.6xz1.1'+wm+bs+'.dat').transpose()
	dnw_iso = dnw[1][mini:maxi]-dtn[2]
	ddt_iso = (dnw_iso/etn**2.+dsw_iso/ets**2.)/(1./etn**2.+1./ets**2.)
	#print(ddt_iso)
	plt.errorbar(dnw[0][mini:maxi],ddt_iso*1.e3,et[mini:maxi]*1.e3,fmt='ko')
	plt.plot(dts[0],dtt_iso*1.e3,'k-')
	#plt.plot(dt[0],dt[0]**2.*dtn[1],'k--')
	plt.xlim(20,190)
	#plt.ylim(-35,80)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16,labelpad=-.75)
	plt.ylabel(r'$10^3$($\xi(s)-\xi_{no BAO}$)',size=16)
	#plt.text(30,71,r'$\alpha=1.044\pm0.042$',color='k',size=16)
	#plt.text(35,64,r'$\chi^2$/dof = 5.8/7',color='k',size=16)
	#plt.text(30,71,r'$\alpha=0.986\pm0.040$',color='k',size=16)
	covi = np.linalg.pinv(cov[mini:maxi,mini:maxi])
	diff = (ddt_iso-dtt_iso)
	chi2 = np.dot(diff,(np.dot(diff,covi)))
	print (chi2)
	plt.text(120,1.00,r'$\chi^2$/dof = '+str(round(chi2,1))+'/'+str(maxi-mini-5),color='k',size=16)
	#plt.title(r'BAO best-fit for v1.6 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True

def plotxiELGbaofit(bs='5st0',v='7',a='',rec='_rec',reg='NScombf',wm='fkp',minr=50,maxr=150,mom='0',covv='',angfac=0,zmin=0.6,zmax=1.1):
	#Plots comparison between QSO clustering and best-fit BAO theory
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiELGNSbaofit'+v+rec+reg+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dt = np.loadtxt(ebossdir+'ximodELG'+reg+v+rec+a+bs+'.dat').transpose()
	dt_iso = dt[1]-dt[2]
	zw = ''#'mz'+str(zmin)+'xz'+str(zmax)

	recf = ''
	if rec == 'rec':
		recf = '_rec'	
	dsw = np.loadtxt(ebossdir+'xi024ELG_v7SGC'+zw+'CZdata'+rec+bs+'.dat').transpose()
	for i in range(0,len(dsw[0])):
		if dsw[0][i] > minr:
			mini = i
			break
	for i in range(0,len(dsw[0])):
		if dsw[0][i] > maxr:
			maxi = i
			break
	rl = dsw[0][mini:maxi]
	print(rl)

	et = []
	if reg == 'NScombf':
		covn = np.loadtxt(ebossdir+'cov0NGC'+zw+'ELG_v'+str(v)+'_EZ'+rec+bs+'.dat')
		covs = np.loadtxt(ebossdir+'cov0SGC'+zw+'ELG_v'+str(v)+'_EZ'+rec+bs+'.dat')  
	
		covi = np.linalg.pinv(covn)+np.linalg.pinv(covs)
		cov = np.linalg.pinv(covi)
		
		ets = []
		etn = []
		for i in range(0,maxi):
			et.append(sqrt(cov[i][i]))
			etn.append(sqrt(covn[i][i]))
			ets.append(sqrt(covs[i][i]))
		etn = np.array(etn)[mini:maxi]
		ets = np.array(ets)[mini:maxi]
		et = np.array(et)[mini:maxi]			
		dnw = np.loadtxt(ebossdir+'xi024ELG_v7NGC'+zw+'CZdata'+rec+bs+'.dat').transpose()
		ds = dsw[1][mini:maxi]
		dn = dnw[1][mini:maxi]
		ddt = (dn/etn**2.+ds/ets**2.)/(1./et**2.)
	else:
		cov = np.loadtxt(ebossdir+'cov'+mom+reg+zw+'ELG_v'+str(v)+'_EZ'+rec+bs+'.dat')
		for i in range(0,maxi):
			et.append(sqrt(cov[i][i]))

		et = np.array(et)[mini:maxi]
		ddt = np.loadtxt(ebossdir+'xi024ELG_v7'+reg+zw+'CZdata'+rec+bs+'.dat').transpose()[1][mini:maxi]
	ddt_iso = ddt-dt[2]
	
	if reg == 'SGC':
		cs = 'steelblue'
		ps = 's'
		plt.text(130,2,'SGC only',color=cs)
	if reg == 'NScombf':
		cs = 'dodgerblue'
		ps =  '^'	
		#plt.text(130,2,'SGC+NGC',color=cs)
	plt.errorbar(rl,ddt_iso*1.e3,et*1.e3,fmt=ps,color=cs,markersize=7,markeredgecolor='darkblue')
	#if reg == 'NScombf':
	#	plt.plot(rl,ddt_iso*1.e3,'wo',markersize=4,zorder=1000)
	plt.plot(rl,dt_iso*1.e3,'-',color=cs)
	#plt.plot(dt[0],dt[0]**2.*dtn[1],'k--')
	plt.xlim(minr-10,maxr+10)
	#plt.ylim(-35,80)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16,labelpad=-.75)
	plt.ylabel(r'$10^3$($\xi(s)-\xi_{no BAO}$)',size=16)
	#plt.text(30,71,r'$\alpha=1.044\pm0.042$',color='k',size=16)
	#plt.text(35,64,r'$\chi^2$/dof = 5.8/7',color='k',size=16)
	#plt.text(30,71,r'$\alpha=0.986\pm0.040$',color='k',size=16)
	covi = np.linalg.pinv(cov[mini:maxi,mini:maxi])
	diff = (ddt_iso-dt_iso)
	chi2 = np.dot(diff,(np.dot(diff,covi)))
	print (chi2)
	#plt.text(120,1.00,r'$\chi^2$/dof = '+str(round(chi2,1))+'/'+str(maxi-mini-5),color='k',size=16)
	#plt.title(r'BAO best-fit for v1.6 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True

def plotxiQSObaofit(bs='5st0',v='7_1',ezv='7',a='',rec='',reg='NScombf',wm='fkp',minr=50,maxr=150,mom='0',covv=''):
	#Plots comparison between QSO clustering and best-fit BAO theory
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiQSONSbaofit'+v+rec+reg+wm+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	dt = np.loadtxt(ebossdir+'ximodQSO'+reg+v+rec+a+bs+'np3.dat').transpose()
	dt_iso = dt[1]-dt[2]
	zw = ''#'mz'+str(zmin)+'xz'+str(zmax)

	recf = ''
	if rec == 'rec':
		recf = '_rec'	

	dnw = np.loadtxt(ebossdir+'dr16_qso_v7_1_N_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()
	dsw = np.loadtxt(ebossdir+'dr16_qso_v7_1_S_0.31_z0.8x2.2_dwfiber_rwfull_covEZv7wsys_3mds5_ximulti_mean.dat').transpose()

	for i in range(0,len(dsw[0])):
		if dsw[0][i] > minr:
			mini = i
			break
	for i in range(0,len(dsw[0])):
		if dsw[0][i] > maxr:
			maxi = i
			break
	rl = dsw[0][mini:maxi]
	print(rl)

	et = []
	if reg == 'NScombf':
		covn = np.loadtxt(ebossdir+'cov0NGC'+zw+'QSO_v'+str(ezv)+'_EZ'+rec+bs+'.dat')
		covs = np.loadtxt(ebossdir+'cov0SGC'+zw+'QSO_v'+str(ezv)+'_EZ'+rec+bs+'.dat')  
	
		covi = np.linalg.pinv(covn)+np.linalg.pinv(covs)
		cov = np.linalg.pinv(covi)
		
		ets = []
		etn = []
		for i in range(0,maxi):
			et.append(sqrt(cov[i][i]))
			etn.append(sqrt(covn[i][i]))
			ets.append(sqrt(covs[i][i]))
		etn = np.array(etn)[mini:maxi]
		ets = np.array(ets)[mini:maxi]
		et = np.array(et)[mini:maxi]			
		ds = dsw[1][mini:maxi]
		dn = dnw[1][mini:maxi]
		ddt = (dn/etn**2.+ds/ets**2.)/(1./et**2.)
	else:
		return('NOT WRITTEN YET!!!')
		cov = np.loadtxt(ebossdir+'cov'+mom+reg+zw+'QSO_v'+str(v)+'_EZ'+rec+bs+'.dat')
		for i in range(0,maxi):
			et.append(sqrt(cov[i][i]))

		et = np.array(et)[mini:maxi]
		#ddt = np.loadtxt(ebossdir+'xi024ELG_v7'+reg+zw+'CZdata'+rec+bs+'.dat').transpose()[1][mini:maxi]
	ddt_iso = ddt-dt[2]
	
	if reg == 'SGC':
		cs = 'steelblue'
		ps = 's'
		plt.text(130,2,'SGC only',color=cs)
	if reg == 'NScombf':
		cs = 'orange'
		bc = 'gold'
		ps =  '^'	
		#plt.text(130,2,'SGC+NGC',color=cs)
	plt.errorbar(rl,ddt_iso*1.e3,et*1.e3,fmt=ps,color=cs,markersize=7,mfc=bc)
	#if reg == 'NScombf':
	#	plt.plot(rl,ddt_iso*1.e3,'wo',markersize=4,zorder=1000)
	plt.plot(rl,dt_iso*1.e3,'-',color=cs)
	#plt.plot(dt[0],dt[0]**2.*dtn[1],'k--')
	plt.xlim(minr-10,maxr+10)
	#plt.ylim(-35,80)
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)',size=16,labelpad=-.75)
	plt.ylabel(r'$10^3$($\xi(s)-\xi_{no BAO}$)',size=16)
	#plt.text(30,71,r'$\alpha=1.044\pm0.042$',color='k',size=16)
	#plt.text(35,64,r'$\chi^2$/dof = 5.8/7',color='k',size=16)
	#plt.text(30,71,r'$\alpha=0.986\pm0.040$',color='k',size=16)
	covi = np.linalg.pinv(cov[mini:maxi,mini:maxi])
	diff = (ddt_iso-dt_iso)
	chi2 = np.dot(diff,(np.dot(diff,covi)))
	print (chi2)
	#plt.text(120,1.00,r'$\chi^2$/dof = '+str(round(chi2,1))+'/'+str(maxi-mini-5),color='k',size=16)
	#plt.title(r'BAO best-fit for v1.6 eboss QSOs, 0.9 < z < 2.2')
	pp.savefig()
	pp.close()
	return True


def plotELGNSbaolike(v='7',p='3',Bp='0.59303.05.015.01.0',rec='_rec',bs='5st0',reg='NScomb'):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiELG'+reg+'baolik'+v+rec+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	
	db = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELGdata'+v+rec+Bp+bs+'.dat').transpose()
	a = db[0]
	ax = sigreg_c12(ebossdir+'BAOxichil'+reg+'ELGdata'+v+rec+Bp+bs)
	print (ax)
	alph = (ax[2]+ax[1])/2.
	err = (ax[2]-ax[1])/2.
	dnb = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELGdata'+v+rec+Bp+'nobao'+bs+'.dat').transpose()[1]
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
	#plt.text(.9,10,r'$\alpha=$'+str(round(alph,3))+r'$\pm$'+str(round(err,3)))
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	#plt.title(r'BAO likelihood for '+title)
	if reg == 'NScomb':
		lab = 'NGC + SGC chi2'
	if reg == 'NScombf':
		lab = 'NGC + SGC xi'
	if reg == 'NGC':
		lab = 'NGC only'
	if reg == 'SGC':
		lab = 'SGC only'

	#if rec == '_rec':
	#	plt.title(r'ELGs, $\xi$, post-recon '+lab)
	#if rec == '':
	#	plt.title(r'ELGs, $\xi$, pre-recon')
	pp.savefig()
	pp.close()
	return True

def plotQSObaolike(v='7_2',p='3',Bp='0.44.03.08.015.00',rec='',bs='5st0',reg='NScombf'):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiELG'+reg+'baolik'+v+rec+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	db = np.loadtxt(ebossdir+'BAOxichil'+reg+'QSOdataJH'+v+rec+Bp+bs+'.dat').transpose()
	a = db[0]
	ax = sigreg_c12(ebossdir+'BAOxichil'+reg+'QSOdataJH'+v+rec+Bp+bs)
	print (ax)
	alph = (ax[2]+ax[1])/2.
	err = (ax[2]-ax[1])/2.
	dnb = np.loadtxt(ebossdir+'BAOxichil'+reg+'QSOdataJH'+v+rec+Bp+'noBAO'+bs+'.dat').transpose()[1]
	chim = min(db[1])
	ol = np.ones((len(db[0])))
	plt.plot(db[0],db[1]-chim,'-',color='k',linewidth=2)
	plt.plot(db[0],dnb-chim,'--',color='k',linewidth=1)
	plt.plot(db[0],ol,'k:',linewidth=1)
	plt.text(0.825,1.1,r'$1\sigma$')
	#plt.plot(db[0],ol*4,'k:',linewidth=1)
	#plt.text(0.825,4.1,r'$2\sigma$')
	plt.plot(db[0],ol*9,'k:',linewidth=1)
	plt.text(0.825,9.1,r'$3\sigma$')
	plt.plot(db[0],ol*25,'k:',linewidth=1)
	plt.text(0.825,25.1,r'$5\sigma$')

	plt.ylim(0,50)
	plt.xlabel(r'$\alpha_{\rm iso,BAO}$',size=18)
	plt.ylabel(r'$\Delta\chi^2$',size=18)
	#plt.text(.9,10,r'$\alpha=$'+str(round(alph,3))+r'$\pm$'+str(round(err,3)))
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	#plt.title(r'BAO likelihood for '+title)
	if reg == 'NScomb':
		lab = 'NGC + SGC chi2'
	if reg == 'NScombf':
		lab = 'NGC + SGC xi'
	if reg == 'NGC':
		lab = 'NGC only'
	if reg == 'SGC':
		lab = 'SGC only'

	#if rec == '_rec':
	#	plt.title(r'ELGs, $\xi$, post-recon '+lab)
	#if rec == '':
	#	plt.title(r'ELGs, $\xi$, pre-recon')
	pp.savefig()
	pp.close()
	return True


def plotELGNSbaolikexiP(v='7',p='3',Bp='0.59303.05.015.01.0',rec='_rec',bs='5st0',reg='NScomb'):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiandPELG'+reg+'baolik'+v+rec+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	colors = ['dodgerblue','navy']
	db = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELGdata'+v+rec+Bp+bs+'.dat').transpose()
	a = db[0]
	ax = sigreg_c12(ebossdir+'BAOxichil'+reg+'ELGdata'+v+rec+Bp+bs)
	print (ax)
	alph = (ax[2]+ax[1])/2.
	err = (ax[2]-ax[1])/2.
	dnb = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELGdata'+v+rec+Bp+'nobao'+bs+'.dat').transpose()[1]
	chim = min(db[1])
	ol = np.ones((len(db[0])))
	plt.plot(db[0],db[1]-chim,'-',color=colors[0],linewidth=2,label=r'$\xi$')
	pkl = np.loadtxt(ebossdir+'sdss_DR16_ELG_BAO_Pk_baoiso_k-0.03-0.3_rec_analytic_sigmanl-4.dat').transpose()
	pklsm = np.loadtxt(ebossdir+'sdss_DR16_ELG_BAO_Pk_baoisosm_k-0.03-0.3_rec_analytic_sigmanl-4.dat').transpose()
	plt.plot(pkl[0],pkl[1]-np.min(pkl[1]),'-',color=colors[1],linewidth=1.5,label=r'$P(k)$')
	print((np.min(pkl[1]),chim))
	plt.legend()
	plt.plot(db[0],dnb-chim,'--',color=colors[0],linewidth=1)
	plt.plot(pkl[0],pklsm[1]-np.min(pkl[1]),'--',color=colors[1],linewidth=1)
	plt.plot(db[0],ol,'k:',linewidth=1)
	plt.text(0.825,1.1,r'$1\sigma$')
	plt.plot(db[0],ol*4,'k:',linewidth=1)
	plt.text(0.825,4.1,r'$2\sigma$')
	plt.plot(db[0],ol*9,'k:',linewidth=1)
	plt.text(0.825,9.1,r'$3\sigma$')
	plt.ylim(0,16)
	plt.xlim(.8,1.2)
	plt.xlabel(r'$\alpha_{\rm BAO}$',size=18)
	plt.ylabel(r'$\Delta\chi^2$',size=18)
	#plt.text(.9,10,r'$\alpha=$'+str(round(alph,3))+r'$\pm$'+str(round(err,3)))
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	#plt.title(r'BAO likelihood for '+title)
	if reg == 'NScomb':
		lab = 'NGC + SGC chi2'
	if reg == 'NScombf':
		lab = 'NGC + SGC xi'
	if reg == 'NGC':
		lab = 'NGC only'
	if reg == 'SGC':
		lab = 'SGC only'
	
	

	#if rec == '_rec':
	#	plt.title(r'ELGs, $\xi$, post-recon '+lab)
	#if rec == '':
	#	plt.title(r'ELGs, $\xi$, pre-recon')
	pp.savefig()
	pp.close()
	return True



def plotELGNSbaolikeNScompS(v='7',p='3',Bp='0.59304.07.015.01.0',rec='_rec',bs='5st0'):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiELGNScSbaolik'+v+rec+bs+'.pdf')
	plt.clf()
	plt.minorticks_on()
	
	db = np.loadtxt(ebossdir+'BAOxichilNScombELGdata'+v+rec+Bp+bs+'.dat').transpose()
	dbs = np.loadtxt(ebossdir+'BAOxichilSGCELGdata'+v+rec+Bp+bs+'.dat').transpose()
	#dbn = np.loadtxt(ebossdir+'BAOxichilNGCELGdata'+v+rec+Bp+bs+'.dat').transpose()
	a = db[0]
	#ax = sigreg_c12(ebossdir+'BAOxichil'+reg+'ELGdata'+v+rec+Bp+bs)
	dnb = np.loadtxt(ebossdir+'BAOxichilNScombELGdata'+v+rec+Bp+'nobao'+bs+'.dat').transpose()[1]
	dnbs = np.loadtxt(ebossdir+'BAOxichilSGCELGdata'+v+rec+Bp+'nobao'+bs+'.dat').transpose()[1]
	chim = min(db[1])
	chims = min(dbs[1])
	#chimn = min(dbn[1])
	ol = np.ones((len(db[0])))
	plt.plot(db[0],db[1]-chim,'-',color='k',linewidth=2.5,label='NGC+SGC')
	plt.plot(db[0],dbs[1]-chims,'-',color='steelblue',linewidth=1.5,label='SGC only')
	#plt.plot(db[0],dbn[1]-chimn,'-',color='midnightblue',linewidth=1.5,label='NGC only')
	plt.legend()
	plt.plot(db[0],dnb-chim,'--',color='k',linewidth=1.25)
	plt.plot(db[0],dnbs-chims,'--',color='steelblue',linewidth=.75)
	plt.plot(db[0],ol,'k:',linewidth=.5)
	plt.text(0.825,1.1,r'$1\sigma$')
	plt.plot(db[0],ol*4,'k:',linewidth=.5)
	plt.text(0.825,4.1,r'$2\sigma$')
	plt.plot(db[0],ol*9,'k:',linewidth=.5)
	plt.text(0.825,9.1,r'$3\sigma$')
	plt.ylim(0,16)
	plt.xlabel(r'$\alpha_{\rm BAO}$',size=18)
	plt.ylabel(r'$\Delta\chi^2$',size=18)
	#plt.text(.9,10,r'$\alpha=$'+str(round(alph,3))+r'$\pm$'+str(round(err,3)))
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	#plt.title(r'BAO likelihood for '+title)

	#if rec == '_rec':
	#	plt.title(r'ELGs, $\xi$, post-recon '+lab)
	#if rec == '':
	#	plt.title(r'ELGs, $\xi$, pre-recon')
	pp.savefig()
	pp.close()
	return True


def plotELGbaolikeprepost(v='4',reg='NScombf',bs='8st0',damppre='0.5933.06.010.015.00',damprec='0.59304.07.015.01.0'):
	#plot bao likelihood for QSOs
	from matplotlib import pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(ebossdir+'xiELGNSbaolikprepost.pdf')
	plt.clf()
	plt.minorticks_on()
	
	db = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELG'+v+'rec'+damprec+bs+'.dat').transpose()
	a = db[0]
	ax = sigreg_c12(db,md='l')
	print (ax)
	alph = (ax[2]+ax[1])/2.
	err = (ax[2]-ax[1])/2.
	dnb = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELG'+v+'rec'+damprec+'nobao'+bs+'.dat').transpose()[1]
	chim = min(db[1])
	ol = np.ones((len(db[0])))
	plt.plot(db[0],db[1]-chim,'-',color='k',linewidth=4)
	plt.plot(db[0],dnb-chim,'--',color='k',linewidth=3)

	db = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELG'+v+damppre+bs+'.dat').transpose()
	a = db[0]
	ax = sigreg_c12(db,md='l')
	print (ax)
	alph = (ax[2]+ax[1])/2.
	err = (ax[2]-ax[1])/2.
	dnb = np.loadtxt(ebossdir+'BAOxichil'+reg+'ELG'+v+damppre+'nobao'+bs+'.dat').transpose()[1]
	chim = min(db[1])
	plt.plot(db[0],db[1]-chim,'-',color='0.5',linewidth=2)
	plt.plot(db[0],dnb-chim,'--',color='0.5',linewidth=1)


	plt.plot(db[0],ol,'k:',linewidth=1)
	plt.text(0.825,1.1,r'$1\sigma$')
	plt.plot(db[0],ol*4,'k:',linewidth=1)
	plt.text(0.825,4.1,r'$2\sigma$')
	plt.plot(db[0],ol*9,'k:',linewidth=1)
	plt.text(0.825,9.1,r'$3\sigma$')
	plt.ylim(0,16)
	plt.xlabel(r'$\alpha_{\rm BAO}$',size=18)
	plt.ylabel(r'$\Delta\chi^2$',size=18)
	xl = [0.91,0.94]
	yl = [14.2,14.2]
	plt.plot(xl,yl,'k-',linewidth=4)
	yl = [13.2,13.2]
	plt.plot(xl,yl,'-',color='0.5',linewidth=2)
	plt.text(0.95,14,'post-recon',color='k')
	plt.text(0.95,13,'pre-recon',color='0.5')
	xl = [0.9,1.02]
	yl = [14.8,14.8]
	plt.plot(xl,yl,'k-',linewidth=1)
	xl = [0.9,1.02]
	yl = [12.5,12.5]
	plt.plot(xl,yl,'k-',linewidth=1)
	xl = [0.9,0.9]
	yl = [12.5,14.8]
	plt.plot(xl,yl,'k-',linewidth=1)
	xl = [1.02,1.02]
	plt.plot(xl,yl,'k-',linewidth=1)
	#plt.text(.9,10,r'$\alpha=$'+str(round(alph,3))+r'$\pm$'+str(round(err,3)))
	#plt.text(30,120,r'$\alpha=0.941\pm0.018$',color='k',size=16)
	#plt.title(r'BAO likelihood for '+title)
	#if rec == '_rec':
	#	plt.title(r'ELGs, $\xi$, post-recon')
	#if rec == '':
	#	plt.title(r'ELGs, $\xi$, pre-recon')
	pp.savefig()
	pp.close()
	return True

def plotbox(cnr,lw=1):
	plt.plot([cnr[0][0],cnr[0][0]],[cnr[1][0],cnr[1][1]],'k-')
	plt.plot([cnr[0][0],cnr[0][1]],[cnr[1][0],cnr[1][0]],'k-')
	plt.plot([cnr[0][1],cnr[0][1]],[cnr[1][0],cnr[1][1]],'k-')
	plt.plot([cnr[0][0],cnr[0][1]],[cnr[1][1],cnr[1][1]],'k-')

def Nzalleboss():
	bossf = fitsio.read('/Users/ashleyross/Dropbox/BOSS/galaxy_DR12v5_CMASSLOWZTOT_North.fits.gz')
	boss1 = (bossf['IMATCH'] == 1) & (bossf['Z'] > 0.43)
	lrgf = fitsio.read('/Users/ashleyross/Dropbox/eBOSS/eBOSS_LRG_full_ALLdata-vDR16.fits')
	lrg1 = (lrgf['IMATCH'] == 1) & (lrgf['Z'] > 0.01)  
	qsof = fitsio.read('/Users/ashleyross/Dropbox/eBOSS/eBOSS_QSO_full_ALL_v4.dat.fits')
	qsoc = (qsof['IMATCH'] == 1) & (qsof['Z'] > 0.01) & (qsof['Z'] < 4)
	qsol = (qsof['IMATCH'] == 2) & (qsof['Z'] > 0.01) & (qsof['Z'] < 4)
	elgf = fitsio.read('/Users/ashleyross/Dropbox/eBOSS/eBOSS_ELG_full_ALL_v4.dat_really_masked_cutted.fits')
	elg1 = elgf['IMATCH'] == 1
	plt.hist(bossf[boss1]['Z'],histtype='step',color='orange',bins=30,range=(0,3.5),log=True)
	plt.hist(lrgf[lrg1]['Z'],histtype='step',color='red',bins=30,range=(0,3.5),log=True)
	plt.hist(elgf[elg1]['Z'],histtype='step',color='blue',bins=30,range=(0,3.5),log=True)
	plt.hist(qsof[qsoc]['Z'],histtype='step',color='rebeccapurple',bins=30,range=(0,3.5),log=True)
	plt.hist(qsof[qsol]['Z'],histtype='step',color='mediumturquoise',bins=30,range=(0,3.5),log=True)
	plt.legend(['BOSS Galaxies','eBOSS LRGs','eBOSS ELGs','eBOSS Quasars','SDSS I/II/III Quasars'])
	plt.xlabel('Redshift')
	plt.ylabel('Number of redshifts')
	plt.ylim(100,1.5e6)
	plt.savefig('/Users/ashleyross/Dropbox/eBOSS/Nz4DR16unnorm.png')
	plt.show()

def plotnzQSO(smooth=4,sf=1.1):
	ds = np.loadtxt(ebossdir+'nbar_eBOSS_QSO_SGC_v7.dat').transpose()
	dn = np.loadtxt(ebossdir+'nbar_eBOSS_QSO_NGC_v7.dat').transpose()
	#ds4 = np.loadtxt(dir+'nbar-eboss_v1.84-QSO-S-eboss_v1.84.dat').transpose()
	#dn4 = np.loadtxt(dir+'nbar-eboss_v1.84-QSO-N-eboss_v1.84.dat').transpose()
	#pp = PdfPages('nzQSOv19f.pdf')
	plt.clf()
	plt.minorticks_on()
	zl = []
	nsl = []
	nnl = []
	nsl4 = []
	nnl4 = []
	Scol = 'g'
	Ncol = 'gold'
	colors = ['g','goldenrod']
	esl = []
	enl = []
	for i in range(0,len(ds[0]),smooth):
		z = 0
		ns = 0
		nn = 0
		ns4 = 0
		nn4 = 0
		vs = 0
		vn = 0
		nts = 0
		ntn = 0
		if i+smooth < len(ds[0]):
			for j in range(0,smooth):
				z += ds[0][i+j]
				ns += ds[3][i+j]
				nn += dn[3][i+j]
				#ns4 += ds4[3][i+j]
				#nn4 += dn4[3][i+j]
				nts += ds[-1][i+j]
				vs += ds[-2][i+j]
				ntn += dn[-1][i+j]
				vn += dn[-2][i+j]
			z = z/float(smooth)
			ns = ns/float(smooth)
			nn = nn/float(smooth)
			#fkpws = 1./(1.+1.e4*nts/vs)
			#fkpwn = 1./(1.+1.e4*ntn/vn)
			#print(z,fkpwn)
			es = sqrt(nts)/vs
			en = sqrt(ntn)/vn
			esl.append(es)
			enl.append(en)

			zl.append(z)
			nsl.append(ns)
			nnl.append(nn)
	zl = np.array(zl)
	nsl = np.array(nsl)
	nnl = np.array(nnl)
	esl = np.array(esl)
	enl = np.array(enl)
	plt.plot(zl,nsl*1e5,'-',color=colors[1],linewidth=2,label='SGC')
	plt.plot(zl,nsl*1e5*sf,':',color=colors[1],linewidth=2)
	plt.errorbar(zl,nsl*1e5,esl*1e5,fmt=',',color=colors[1],linewidth=2)
	plt.plot(zl,nnl*1e5,'--',color=colors[0],linewidth=2,label='NGC')
	plt.errorbar(zl,nnl*1e5,enl*1e5,fmt=',',color=colors[0],linewidth=2)
	#plt.plot(zl,nsl4*1e5,'--',color='paleturquoise',linewidth=2)
	#plt.plot(zl,nnl4*1e5,'--',color='lightcoral',linewidth=2)
	plt.xlim(0.5,2.5)
	xl = [.8,.8]
	yl = [0,10]
	plt.plot(xl,yl,'k:',linewidth=1)
	xl = [2.2,2.2]
	plt.plot(xl,yl,'k:',linewidth=1)
# 	cnr = [(1.8,2.13),(1.9,2.25)]
# 	plt.text(1.92,2.08,'NGC',fontsize=16,color=Ncol)
# 	plt.text(1.92,1.95,'SGC',fontsize=16,color=Scol)
# 	xl = [1.82,1.9]
# 	yl = [2.12,2.12]
# 	plt.plot(xl,yl,'-',color=Ncol,linewidth=2)
# 	yl = [1.99,1.99]
# 	plt.plot(xl,yl,'-',color=Scol,linewidth=2)
# 	plotbox(cnr)
	plt.legend()
	plt.ylim(0,2.3)
	plt.ylabel(r'$10^{5}$ $n$ ($h^{3}$Mpc$^{-3}$)',size=16)
	plt.xlabel('quasar redshift',size=16)

	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/nzquasarDR16werr.png')
	plt.show()
	#pp.close()
	w = (zl > 0.8) & (zl < 2.2)
	chi2 = sum((nsl[w]-nnl[w])**2./(esl[w]**2.+enl[w]**2.))
	
	chi2f = sum((nsl[w]*sf-nnl[w])**2./(esl[w]**2.+enl[w]**2.))
	print('chi2,chi2 sgc x'+str(sf)+', nbin')
	print(chi2,chi2f,len(nsl[w]))

	return True

def getnzMockLRG(reg='SGC',nbin = 55,zmin=0,zmax=1.1):
	from Cosmo import distance
	d = distance(.31,.69)
	vl = []
	sp = (zmax-zmin)/float(nbin)
	if reg == 'SGC':
		area = 1627.
	if reg == 'NGC':
		area = 2476.	
	zl = []	
	for i in range(0,nbin):
		zlo = i*sp
		zh = (i+1)*sp
		v = area/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		vl.append(v)
		zl.append(i*sp+sp/2.)
	vl = np.array(vl)
	zt = np.zeros(nbin)    
	for i in range(1,1001):
		zer = ''
		if i < 10:
			zer += '0'
		if i < 100:
			zer += 	'0'
		if i < 1000:
			zer += 	'0'
		#di = fitsio.read('/project/projectdirs/eboss/czhao/EZmock/LRG_v5/eBOSS_'+reg+'/EZmock_eBOSS_LRG_'+reg+'_v5_'+zer+str(i)+'.fits.gz')
		di = fitsio.read('/global/cscratch1/sd/jbautist/eboss/EZmocks/EZmock_LRG_v7.0_syst/eBOSS_LRG/EZmock_eBOSS_LRG_'+reg+'_v7_'+zer+str(i)+'.dat.fits ')
		zh = np.histogram(di['Z'],bins=nbin,range=(zmin,zmax),weights=di['WEIGHT_NOZ']*di['WEIGHT_CP'])
		zt += zh[0]
	zt = zt/1000./vl
	print(zt)
	#print(zl)
	#print(zh[1])
	cov = np.zeros((nbin,nbin))	
	for i in range(1,1001):
		zer = ''
		if i < 10:
			zer += '0'
		if i < 100:
			zer += 	'0'
		if i < 1000:
			zer += 	'0'
		#di = fitsio.read('/project/projectdirs/eboss/czhao/EZmock/LRG_v5/eBOSS_'+reg+'/EZmock_eBOSS_LRG_'+reg+'_v5_'+zer+str(i)+'.fits.gz')
		di = fitsio.read('/global/cscratch1/sd/jbautist/eboss/EZmocks/EZmock_LRG_v7.0_syst/eBOSS_LRG/EZmock_eBOSS_LRG_'+reg+'_v7_'+zer+str(i)+'.dat.fits ')

		zh = np.histogram(di['Z'],bins=nbin,range=(zmin,zmax),weights=di['WEIGHT_NOZ']*di['WEIGHT_CP'])
		zi = zh[0]/vl
		for j in range(0,nbin):
			for k in range(0,nbin):
				cov[j][k] += (zt[j]-zi[j])*(zt[k]-zi[k])
		print(i)
	cov = cov/999.
	fo = open('nzLRGEZ'+reg+'ave.dat','w')
	for i in range(0,nbin):
		fo.write(str(zl[i])+' '+str(zt[i])+' '+str(sqrt(cov[i][i]))+'\n')
	fo.close()
	fo = open('covnzLRGEZ'+reg+'.dat','w')
	for i in range(0,nbin):
		for j in range(0,nbin):
			fo.write(str(cov[i][j])+' ')
		fo.write('\n')				
	fo.close()
	return True

def getnzdataLRG(nbin = 30,zmin=.5,zmax=1.1):
	from Cosmo import distance
	d = distance(.31,.69)
	vl = []
	sp = (zmax-zmin)/float(nbin)
	#if reg == 'SGC':
	areas = 1627.
	#if reg == 'NGC':
	arean = 2476.	
	zl = []	
	for i in range(0,nbin):
		zlo = zmin+i*sp
		zh = zmin+(i+1)*sp
		v = 1./(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zlo)**3.)
		vl.append(v)
		zl.append(zmin+i*sp+sp/2.)
	vl = np.array(vl)
	dn = fitsio.read(ebossdir+'/eBOSS_LRG_full_NGC_v7_2.dat.fits')
	wn = dn['IMATCH'] == 1
	wnall = dn[wn]['WEIGHT_CP']*dn[wn]['WEIGHT_NOZ']*dn[wn]['WEIGHT_SYSTOT']
	zn = np.histogram(dn[wn]['Z'],bins=nbin,range=(zmin,zmax),weights=wnall)
	nzn = zn[0]/vl/arean
	ds = fitsio.read(ebossdir+'/eBOSS_LRG_full_SGC_v7_2.dat.fits')
	ws = ds['IMATCH'] == 1
	wsall = ds[ws]['WEIGHT_CP']*ds[ws]['WEIGHT_NOZ']*ds[ws]['WEIGHT_SYSTOT']
	zs = np.histogram(ds[ws]['Z'],bins=nbin,range=(zmin,zmax),weights=wsall)
	nzs = zs[0]/vl/areas
	plt.plot(zl,nzn,zl,nzs)
	plt.show()


def getnzdataLRG_syssplit(sys,reg,nbin = 20,zmin=.6,zmax=1.,fkpp=1.2):
	from healpix import radec2thphi
	import healpy as hp
	#from Cosmo import distance
	if sys == 'star':
		res = 256
		ms = np.loadtxt(ebossdir+'allstars17.519.9Healpixall256.dat')
		print(len(ms))
		nest = True
	else:
		res = 512
		nest = False
		ms  = fitsio.read(ebossdir+'SDSSimageprop_Nside512.fits')[sys]
	f = fitsio.read(ebossdir+'eBOSS_LRG_clustering_'+reg+'_v7_2.dat.fits')
	thphi = radec2thphi(f['RA'],f['DEC'])	
	pix = hp.ang2pix(res,thphi[0],thphi[1],nest=nest)
	sysg = np.zeros(len(f))
	nb = 0
	for i in range(0,len(f)):
		px = pix[i]
		sysv = ms[int(px)]
		if sysv*0 == 0:
			sysg[i] = sysv
		else:
			nb += 1	
	medsys = np.median(sysg)
	print(medsys,nb)
	#plt.hist(sysg)
	#plt.show()
	fr = fitsio.read(ebossdir+'eBOSS_LRG_clustering_'+reg+'_v7_2.ran.fits')
	thphi = radec2thphi(fr['RA'],fr['DEC'])	
	pixr = hp.ang2pix(res,thphi[0],thphi[1],nest=nest)
	sysr = np.zeros(len(fr))
	nb = 0
	for i in range(0,len(fr)):
		px = pixr[i]
		sysv = ms[int(px)]
		if sysv*0 == 0:
			sysr[i] = sysv
		else:
			nb += 1
	print(nb)	
	wr = sysr > medsys
	norml = sum(fr[~wr]['WEIGHT_SYSTOT'])
	normh = sum(fr[wr]['WEIGHT_SYSTOT'])
	print(norml,normh)
	wg = sysg > medsys
	print(sum(f[wg]['WEIGHT_CP']*f[wg]['WEIGHT_NOZ']*f[wg]['WEIGHT_SYSTOT'])/normh,sum(f[~wg]['WEIGHT_CP']*f[~wg]['WEIGHT_NOZ']*f[~wg]['WEIGHT_SYSTOT'])/norml)

	nh = np.histogram(f[wg]['Z'],bins=nbin,weights=f[wg]['WEIGHT_CP']*f[wg]['WEIGHT_NOZ']*f[wg]['WEIGHT_SYSTOT'])
	nl = np.histogram(f[~wg]['Z'],bins=nbin,weights=f[~wg]['WEIGHT_CP']*f[~wg]['WEIGHT_NOZ']*f[~wg]['WEIGHT_SYSTOT'])
	fkph = np.histogram(f[wg]['Z'],bins=nbin,weights=f[wg]['WEIGHT_CP']*f[wg]['WEIGHT_NOZ']*f[wg]['WEIGHT_SYSTOT']*f[wg]['WEIGHT_FKP'])
	fkpl = np.histogram(f[~wg]['Z'],bins=nbin,weights=f[~wg]['WEIGHT_CP']*f[~wg]['WEIGHT_NOZ']*f[~wg]['WEIGHT_SYSTOT']*f[~wg]['WEIGHT_FKP'])
	fkph = fkph[0]/nh[0]
	fkpl = fkpl[0]/nl[0]
	eh = np.sqrt(nh[0])/fkph**fkpp
	el = np.sqrt(nl[0])/fkpl**fkpp
	zl = []
	for i in range(0,nbin):
		zl.append((nh[1][i]+nh[1][i+1])/2.)
	nh = nh[0]/normh
	eh = eh/normh
	nl = nl[0]/norml
	el = el/norml

	plt.errorbar(zl,nh,eh,fmt='r-,',label='high '+sys.split('_')[0])
	plt.errorbar(zl,nl,el,fmt='b-,',label='low '+sys.split('_')[0])
	plt.legend()
	plt.show()
	chi2 = np.sum((nh-nl)**2./(eh**2.+el**2.))
	print(chi2)
	return True


def plotnzLRG(smooth=2,fkpp=1.2):
	ds = np.loadtxt(ebossdir+'nbar_eBOSS_LRG_SGC_v7.dat').transpose()
	dn = np.loadtxt(ebossdir+'nbar_eBOSS_LRG_NGC_v7.dat').transpose()
	#ds4 = np.loadtxt(dir+'nbar-eboss_v1.84-QSO-S-eboss_v1.84.dat').transpose()
	#dn4 = np.loadtxt(dir+'nbar-eboss_v1.84-QSO-N-eboss_v1.84.dat').transpose()
	#pp = PdfPages('nzQSOv19f.pdf')
	plt.clf()
	plt.minorticks_on()
	zl = []
	nsl = []
	nnl = []
	nsl4 = []
	nnl4 = []
	eln = []
	els = []
	colors = ['firebrick','salmon']
	ntots = 0
	for i in range(0,len(ds[0]),smooth):
		z = 0
		ns = 0
		nn = 0
		ns4 = 0
		nn4 = 0
		nts = 0
		ntn = 0
		vs = 0
		vn = 0
		
		if i+smooth < len(ds[0]):
			for j in range(0,smooth):
				z += ds[0][i+j]
				ns += ds[3][i+j]
				nn += dn[3][i+j]
				nts += ds[-1][i+j]
				vs += ds[-2][i+j]
				ntn += dn[-1][i+j]
				vn += dn[-2][i+j]
				ntots += ds[-1][i+j]
			
				#ns4 += ds4[3][i+j]
				#nn4 += dn4[3][i+j]
			z = z/float(smooth)
			ns = ns/float(smooth)
			nn = nn/float(smooth)
			fkpws = 1./(1.+1.e4*nts/vs)
			fkpwn = 1./(1.+1.e4*ntn/vn)
			print(z,fkpwn)
			es = sqrt(nts)/vs/fkpws**fkpp
			en = sqrt(ntn)/vn/fkpwn**fkpp
			#print(ns,nts/vs,es,es/(nts/vs),nts,sqrt(nts)/nts)
			#ns4 = ns4/float(smooth)
			#nn4 = nn4/float(smooth)
			zl.append(z)
			nsl.append(nts/vs)
			nnl.append(ntn/vn)
			eln.append(en)
			els.append(es)
			#nsl4.append(ns4)
			#nnl4.append(nn4)
	print('total in SGC '+str(ntots))
	zl = np.array(zl)
	nsl = np.array(nsl)
	nnl = np.array(nnl)
	eln = np.array(eln)
	els = np.array(els)
	#nsl4 = np.array(nsl4)
	#nnl4 = np.array(nnl4)
	elsm = np.loadtxt(ebossdir+'nzLRGEZSGCave.dat').transpose()[2]
	elnm = np.loadtxt(ebossdir+'nzLRGEZNGCave.dat').transpose()[2]
	plt.plot(zl,nsl*1e4,'-',color=colors[1],linewidth=2,label='SGC')
	plt.errorbar(zl[:len(elsm)],nsl[:len(elsm)]*1e4,elsm*1e4,fmt=',',color=colors[1])
	plt.plot(zl,nnl*1e4,'--',color=colors[0],linewidth=2,label='NGC')
	plt.errorbar(zl[:len(elsm)],nnl[:len(elsm)]*1e4,elnm*1e4,fmt=',',color=colors[0])
	#plt.plot(zl,nsl4*1e5,'--',color='paleturquoise',linewidth=2)
	#plt.plot(zl,nnl4*1e5,'--',color='lightcoral',linewidth=2)
	plt.xlim(0.5,1.1)
	xl = [.6,.6]
	yl = [0,10]
	plt.plot(xl,yl,'k:',linewidth=1)
	xl = [1.,1.0]
	plt.plot(xl,yl,'k:',linewidth=1)
# 	cnr = [(1.8,2.13),(1.9,2.25)]
# 	plt.text(1.92,2.08,'NGC',fontsize=16,color=Ncol)
# 	plt.text(1.92,1.95,'SGC',fontsize=16,color=Scol)
# 	xl = [1.82,1.9]
# 	yl = [2.12,2.12]
# 	plt.plot(xl,yl,'-',color=Ncol,linewidth=2)
# 	yl = [1.99,1.99]
# 	plt.plot(xl,yl,'-',color=Scol,linewidth=2)
# 	plotbox(cnr)
	plt.legend()
	plt.ylim(0,1.3)
	plt.ylabel(r'$10^{4}$ $n$ ($h^{3}$Mpc$^{-3}$)',size=16)
	plt.xlabel('LRG redshift',size=16)

	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/nzLRGDR16werr.png')
	plt.show()
	#pp.close()
	covs = np.loadtxt(ebossdir+'covnzLRGEZSGC.dat')[30:50,30:50]
	covn = np.loadtxt(ebossdir+'covnzLRGEZNGC.dat')[30:50,30:50]
	covt = covs+covn
	icov = np.linalg.inv(covt)
	dns = (nsl-nnl)[30:50]
	chi2 = np.dot(dns,(np.dot(dns,icov)))
	print(chi2)
	chi2d = 0
	for i in range(0,20):
		chi2d += dns[i]**2./covt[i][i]
	print(chi2d)
	plt.plot(zl[:len(elsm)],elsm,label='mocks')
	plt.plot(zl,els,label='Poisson with FKW factor')
	plt.xlim(0.5,1.1)
	plt.ylim(0,1e-5)
	plt.legend()
	plt.show()	
	plt.plot(zl[:len(elsm)],elnm,label='mocks')
	plt.plot(zl,eln,label='Poisson with FKW factor')
	plt.xlim(0.5,1.1)
	plt.ylim(0,1e-5)
	plt.legend()
	plt.show()	

	return True

def plotnzLRGpCMASS(smooth=2):
	ds = np.loadtxt(ebossdir+'nbar_eBOSS_LRG_SGC_v7.dat').transpose()
	dn = np.loadtxt(ebossdir+'nbar_eBOSS_LRG_NGC_v7.dat').transpose()
	dsc = np.loadtxt(ebossdir+'nbar_DR12v5_CMASS_South_om0p31_Pfkp10000.dat').transpose()
	dnc = np.loadtxt(ebossdir+'nbar_DR12v5_CMASS_North_om0p31_Pfkp10000.dat').transpose()
	#ds4 = np.loadtxt(dir+'nbar-eboss_v1.84-QSO-S-eboss_v1.84.dat').transpose()
	#dn4 = np.loadtxt(dir+'nbar-eboss_v1.84-QSO-N-eboss_v1.84.dat').transpose()
	#pp = PdfPages('nzQSOv19f.pdf')
	plt.clf()
	plt.minorticks_on()
	zl = []
	nsl = []
	nnl = []
	nscl = []
	nncl = []
	colors = ['darkred','r','orange']
	for i in range(0,len(ds[0]),smooth):
		z = 0
		ns = 0
		nn = 0
		ns4 = 0
		nn4 = 0
		if i+smooth < len(ds[0]):
			for j in range(0,smooth):
				z += ds[0][i+j]
				ns += ds[3][i+j]
				nn += dn[3][i+j]
				#ns4 += ds4[3][i+j]
				#nn4 += dn4[3][i+j]
			z = z/float(smooth)
			ns = ns/float(smooth)
			nn = nn/float(smooth)
			#ns4 = ns4/float(smooth)
			#nn4 = nn4/float(smooth)
			zl.append(z)
			nsl.append(ns)
			nnl.append(nn)
			#nsl4.append(ns4)
			#nnl4.append(nn4)
	for i in range(0,len(dsc[0]),smooth*2):
		z = 0
		nsc = 0
		nnc = 0
		ns4 = 0
		nn4 = 0
		if i+smooth < len(dsc[0]):
			for j in range(0,smooth*2):
				#z += ds[0][i+j]
				nsc += dsc[3][i+j]
				nnc += dnc[3][i+j]
				#ns4 += ds4[3][i+j]
				#nn4 += dn4[3][i+j]
			#z = z/float(smooth)
			nsc = nsc/float(smooth*2)
			nnc = nnc/float(smooth*2)
			#ns4 = ns4/float(smooth)
			#nn4 = nn4/float(smooth)
			#zl.append(z)
			nscl.append(nsc)
			nncl.append(nnc)
			#nsl4.append(ns4)
			#nnl4.append(nn4)

	zl = np.array(zl)
	nsl = np.array(nsl)
	nnl = np.array(nnl)
	nscl = np.array(nscl)
	nncl = np.array(nncl)
	print(len(zl),len(nsl),len(nscl))
	#nsl4 = np.array(nsl4)
	#nnl4 = np.array(nnl4)
	#plt.plot(zl,nsl*1e4,'-',color=colors[1],linewidth=2,label='SGC')
	plt.plot(zl,nnl*1e4,'--',color=colors[1],linewidth=2,label='eBOSS LRGs')
	plt.plot(zl[:len(nncl)],nncl*1e4,':',color=colors[2],linewidth=2,label='BOSS CMASS')
	plt.plot(zl[:len(nncl)],(nnl[:len(nncl)]+nncl)*1e4,'-',color=colors[0],linewidth=2,label='CMASS+eBOSS LRGs')
	#plt.plot(zl,nsl4*1e5,'--',color='paleturquoise',linewidth=2)
	#plt.plot(zl,nnl4*1e5,'--',color='lightcoral',linewidth=2)
	plt.xlim(0.5,1.1)
	xl = [.6,.6]
	yl = [0,10]
	plt.plot(xl,yl,'k:',linewidth=1)
	xl = [1.,1.0]
	plt.plot(xl,yl,'k:',linewidth=1)
# 	cnr = [(1.8,2.13),(1.9,2.25)]
# 	plt.text(1.92,2.08,'NGC',fontsize=16,color=Ncol)
# 	plt.text(1.92,1.95,'SGC',fontsize=16,color=Scol)
# 	xl = [1.82,1.9]
# 	yl = [2.12,2.12]
# 	plt.plot(xl,yl,'-',color=Ncol,linewidth=2)
# 	yl = [1.99,1.99]
# 	plt.plot(xl,yl,'-',color=Scol,linewidth=2)
# 	plotbox(cnr)
	plt.legend()
	plt.ylim(0,3.)
	plt.ylabel(r'$10^{4}$ $n$ ($h^{3}$Mpc$^{-3}$)',size=16)
	plt.xlabel('redshift',size=16)

	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/nzLRGpCMASSDR16.png')
	#pp.close()
	return True

def mkcombNbar_LRGcmass(reg):
	nl = np.loadtxt(ebossdir+'nbar_eBOSS_LRG_'+reg+'_v7_2.dat').transpose()
	if reg == 'NGC':
		cr = 'North'
	else:
		cr = 'South'
	nc = np.loadtxt(ebossdir+'nbar_DR12v5_CMASS_'+cr+'_om0p31_Pfkp10000.dat').transpose()
	w = nl[0] < 1.
	ntl = []
	for i in range(0,len(nl[0][w])):
		vcc = nc[-2][i*2]+nc[-2][i*2+1]
		ncc = nc[-1][i*2]+nc[-1][i*2+1]
		nbc = ncc/vcc
		nt = nl[3][i] +nbc
		ntl.append(nt)
	fo = open(ebossdir+'nbar_eBOSS_LRGpCMASS_'+reg+'_v7_2.dat','w')
	for i in range(0,len(ntl)):
		fo.write(str(nl[0][i])+' '+str(nl[1][i])+' '+str(nl[2][i])+' '+str(ntl[i])+'\n')	
	fo.close()

def convert_ELGnbar():
	areas = [117.5,240,237.5,131.9]
	nn = np.loadtxt(ebossdir+'nbar_eBOSS_ELG_v7.dat').transpose()
	nl = []
	nh = []
	nm = []
	nzs = []
	nzn = []
	for i in range(1,len(nn[0])-1,2):
		nl.append((i-1)/2*0.01)
		nh.append((i-1)/2*0.01+0.01)
		nm.append((i-1)/2*0.01+0.005)
		nzsb = ((areas[0]*nn[1][i]+areas[1]*nn[3][i])/(areas[0]+areas[1])+(areas[0]*nn[1][i+1]+areas[1]*nn[3][i+1])/(areas[0]+areas[1]))/2.
		nzs.append(nzsb)
		nznb = ((areas[2]*nn[5][i]+areas[3]*nn[7][i])/(areas[2]+areas[3])+(areas[2]*nn[5][i+1]+areas[3]*nn[7][i+1])/(areas[2]+areas[3]))/2.
		nzn.append(nznb)
	fon = open(ebossdir+'nbar_eBOSS_ELG_NGC_v7.dat','w')	
	fos = open(ebossdir+'nbar_eBOSS_ELG_SGC_v7.dat','w')
	for i in range(0,len(nzn)):
		fon.write(str(nm[i])+' '+str(nl[i])+' '+str(nh[i])+' '+str(nzn[i])+'\n')
		fos.write(str(nm[i])+' '+str(nl[i])+' '+str(nh[i])+' '+str(nzs[i])+'\n')
	fon.close()
	fos.close()
	
def getveff(f1,f2,P1,P2,zmin,zmax,area,omega_matter=.31):
	from Cosmo import distance
	d = distance(omega_matter,1.-omega_matter)
	d1 = np.loadtxt(ebossdir+'nbar_eBOSS_'+f1+'.dat').transpose()
	d2 = np.loadtxt(ebossdir+'nbar_eBOSS_'+f2+'.dat').transpose()
	w = (d1[0]>zmin) & (d1[0]<zmax)
	w2 = (d2[0]>zmin) & (d2[0]<zmax)
	if len(d1[0][w]) != len(d2[0][w2]):
		return "ERROR, will not work if zbinning is different"
	vtot = 0
	veff = 0
	for i in range(0,len(d1[0][w])):
		vi = area/(360.*360./pi)*4.*pi/3.*(d.dc(d1[2][w][i])**3.-d.dc(d1[1][w][i])**3.)
		vtot += vi
		veff += vi*(d1[3][w][i]*P1)/(1.+d1[3][w][i]*P1)*(d2[3][w2][i]*P2)/(1.+d2[3][w2][i]*P2)
	return veff,vtot	

def getveff1(f1,P1,zmin,zmax,area,omega_matter=.31):
	from Cosmo import distance
	d = distance(omega_matter,1.-omega_matter)
	d1 = np.loadtxt(ebossdir+'nbar'+f1+'.dat').transpose()
	w = (d1[0]>zmin) & (d1[0]<zmax)
	vtot = 0
	veff = 0
	for i in range(0,len(d1[0][w])):
		vi = area/(360.*360./pi)*4.*pi/3.*(d.dc(d1[2][w][i])**3.-d.dc(d1[1][w][i])**3.)
		vtot += vi
		veff += vi*(d1[3][w][i]*P1)/(1.+d1[3][w][i]*P1)**2.
	return veff,vtot	


def plotfootsizevscomp(sp=0.02):
	fqn = fitsio.read(ebossdir+'eboss_QSO_full_NGC_v7_1.ran.fits')
	fqs = fitsio.read(ebossdir+'eboss_QSO_full_sGC_v7_1.ran.fits')
	fln = fitsio.read(ebossdir+'eboss_LRG_full_NGC_v7.ran.fits')
	fls = fitsio.read(ebossdir+'eboss_LRG_full_sGC_v7.ran.fits')
	cmp = 0
	cmpl = []
	aql = []
	all = []
	while cmp < 1:
		wqn = fqn['COMP_BOSS'] >= cmp
		nqn = len(fqn[wqn])
		wqs = fqs['COMP_BOSS'] >= cmp
		nqs = len(fqs[wqs])
		aq = (nqn+nqs)/5000.
		wln = fln['COMP_BOSS'] >= cmp
		nln = len(fln[wln])
		wls = fls['COMP_BOSS'] >= cmp
		nls = len(fls[wls])
		al = (nln+nls)/5000.
		print(cmp,aq,al)
		cmpl.append(cmp)
		aql.append(aq)
		all.append(al)
		cmp += sp
	plt.plot(cmpl,aql,'-',color='forestgreen',label='quasars')
	plt.plot(cmpl,all,'--',color='firebrick',label = 'LRGs')
	plt.legend()
	xl = [.5,.5]
	yl = [0,10000]
	plt.plot(xl,yl,'k:')
	plt.ylim(2500,6000)	
	plt.xlabel('Completeness cut')
	plt.ylabel('area (square degrees)')
	plt.savefig(ebossdir+'areavscomp.png')
	plt.show()
	#plt.show()
	#pp.savefig()
	#pp.close()
	return True

def plotfootsizevsssr(sp=0.02):
	fqn = fitsio.read(ebossdir+'eboss_QSO_full_NGC_v7_1.ran.fits')
	fqs = fitsio.read(ebossdir+'eboss_QSO_full_sGC_v7_1.ran.fits')
	fln = fitsio.read(ebossdir+'eboss_LRG_full_NGC_v7.ran.fits')
	fls = fitsio.read(ebossdir+'eboss_LRG_full_sGC_v7.ran.fits')
	cmp = 0
	cmpl = []
	aql = []
	all = []
	while cmp < 1:
		wqn = fqn['sector_SSR'] >= cmp
		nqn = len(fqn[wqn])
		wqs = fqs['sector_SSR'] >= cmp
		nqs = len(fqs[wqs])
		aq = (nqn+nqs)/5000.
		wln = fln['sector_SSR'] >= cmp
		nln = len(fln[wln])
		wls = fls['sector_SSR'] >= cmp
		nls = len(fls[wls])
		al = (nln+nls)/5000.
		print(cmp,aq,al)
		cmpl.append(cmp)
		aql.append(aq)
		all.append(al)
		cmp += sp
	plt.plot(cmpl,aql,'-',color='forestgreen',label='quasars')
	plt.plot(cmpl,all,'--',color='firebrick',label = 'LRGs')
	plt.legend()
	xl = [.5,.5]
	yl = [0,10000]
	plt.plot(xl,yl,'k:')
	plt.ylim(2500,6000)	
	plt.xlabel('Redshift Success cut')
	plt.ylabel('area (square degrees)')
	plt.savefig(ebossdir+'areavsssr.png')
	plt.show()
	#plt.show()
	#pp.savefig()
	#pp.close()
	return True


def getebpcstats(ver='7_2'):
	samp = 'LRGpCMASS'
	f = fitsio.read(ebossdir+'eBOSS_'+samp+'_clustering_NGC_v'+ver+'.dat.fits')
	fs = fitsio.read(ebossdir+'eBOSS_'+samp+'_clustering_SGC_v'+ver+'.dat.fits')
	#statistics for eboss galaxies
	w = f['ISCMASS'] == 0
	ws = fs['ISCMASS'] == 0
	print('total number of eBOSS redshifts, SGC, NGC, combined:')
	print(len(fs[ws]),len(f[w]),len(f[w])+len(fs[ws]))
	print('effective number of eBOSS galaxies, SGC, NGC, combined:')
	nn = np.sum(f[w]['WEIGHT_CP']*f[w]['WEIGHT_NOZ'])
	ns = np.sum(fs[ws]['WEIGHT_CP']*fs[ws]['WEIGHT_NOZ'])
	print(ns,nn,ns+nn)
	w = f['ISCMASS'] == 1
	ws = fs['ISCMASS'] == 1
	print('total number of CMASS redshifts, SGC, NGC, combined:')
	print(len(fs[ws]),len(f[w]),len(f[w])+len(fs[ws]))
	print('effective number of CMASS galaxies, SGC, NGC, combined:')
	nn = np.sum(f[w]['WEIGHT_CP']+f[w]['WEIGHT_NOZ']-1)
	ns = np.sum(fs[ws]['WEIGHT_CP']+fs[ws]['WEIGHT_NOZ']-1)
	print(ns,nn,ns+nn)

	w &= f['IN_EBOSS_FOOT'] == 1
	ws &= fs['IN_EBOSS_FOOT'] == 1
	print('total number of CMASS redshifts in eboss footprint, SGC, NGC, combined:')
	print(len(fs[ws]),len(f[w]),len(f[w])+len(fs[ws]))
	print('effective number of CMASS galaxies in eboss footprint, SGC, NGC, combined:')
	nn = np.sum(f[w]['WEIGHT_CP']+f[w]['WEIGHT_NOZ']-1)
	ns = np.sum(fs[ws]['WEIGHT_CP']+fs[ws]['WEIGHT_NOZ']-1)
	print(ns,nn,ns+nn)

	w = f['ISCMASS'] == 1
	ws = fs['ISCMASS'] == 1
	w &= f['IN_EBOSS_FOOT'] == 0
	ws &= fs['IN_EBOSS_FOOT'] == 0
	print('total number of CMASS redshifts not in eboss footprint, SGC, NGC, combined:')
	print(len(fs[ws]),len(f[w]),len(f[w])+len(fs[ws]))
	print('effective number of CMASS galaxies not eboss footprint, SGC, NGC, combined:')
	nn = np.sum(f[w]['WEIGHT_CP']+f[w]['WEIGHT_NOZ']-1)
	ns = np.sum(fs[ws]['WEIGHT_CP']+fs[ws]['WEIGHT_NOZ']-1)
	print(ns,nn,ns+nn)

	

def getstats(samp,ver='7_2'):
	if samp == 'QSO':
		zmin = 0.8
		zmax = 2.2
	if samp == 	'LRG':
		zmin = 0.6
		zmax = 1.
	f = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_NGC_v'+ver+'.dat.fits')
	fs = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_SGC_v'+ver+'.dat.fits')
	w = f['IMATCH'] == 1
	ws = fs['IMATCH'] == 1
	print('NGC eboss redshifts: '+str(len(f[w])))
	print('SGC eboss redshifts: '+str(len(fs[ws])))
	w = f['IMATCH'] == 2
	ws = fs['IMATCH'] == 2
	print('NGC legacy redshifts: '+str(len(f[w])))
	print('SGC legacy redshifts: '+str(len(fs[ws])))
	w = (f['IMATCH'] == 1) | (f['IMATCH'] == 2)
	ws = (fs['IMATCH'] == 1) | (fs['IMATCH'] == 2)
	print('NGC Neff: '+str(np.sum(f[w]['WEIGHT_NOZ']*f[w]['WEIGHT_CP'])))
	print('SGC Neff: '+str(np.sum(fs[ws]['WEIGHT_NOZ']*fs[ws]['WEIGHT_CP'])))
	w &= (f['COMP_BOSS'] > 0.5) & (f['sector_SSR'] > 0.5)
	ws &= (fs['COMP_BOSS'] > 0.5) & (fs['sector_SSR'] > 0.5)
	print('NGC Neff after completeness cuts: '+str(np.sum(f[w]['WEIGHT_NOZ']*f[w]['WEIGHT_CP'])))
	print('SGC Neff after completeness cuts: '+str(np.sum(fs[ws]['WEIGHT_NOZ']*fs[ws]['WEIGHT_CP'])))
	print('NGC Nztot after completeness cuts: '+str(len(f[w])))
	print('SGC Nztot after completeness cuts: '+str(len(fs[ws])))

	w &= (f['Z'] > zmin) & (f['Z'] < zmax)
	ws &= (fs['Z'] > zmin) & (fs['Z'] < zmax)
	print('NGC Neff after redshift cuts: '+str(np.sum(f[w]['WEIGHT_NOZ']*f[w]['WEIGHT_CP'])))
	print('SGC Neff after redshift cuts: '+str(np.sum(fs[ws]['WEIGHT_NOZ']*fs[ws]['WEIGHT_CP'])))
	print('NGC Nztot after redshift cuts: '+str(len(f[w])))
	print('SGC Nztot after redshift cuts: '+str(len(fs[ws])))

	w &= f['IMATCH'] == 2
	ws &= fs['IMATCH'] == 2
	print('NGC Neff legacy after redshift cuts: '+str(np.sum(f[w]['WEIGHT_NOZ']*f[w]['WEIGHT_CP'])))
	print('SGC Neff legacy after redshift cuts: '+str(np.sum(fs[ws]['WEIGHT_NOZ']*fs[ws]['WEIGHT_CP'])))
	print('NGC Nztot legacy after redshift cuts: '+str(len(f[w])))
	print('SGC Nztot legacy after redshift cuts: '+str(len(fs[ws])))


	w = f['IMATCH'] == 4
	ws = fs['IMATCH'] == 4
	print('NGC stars: '+str(len(f[w])))
	print('SGC stars: '+str(len(fs[ws])))
	w = f['IMATCH'] == 7
	ws = fs['IMATCH'] == 7
	print('NGC z fail: '+str(len(f[w])))
	print('SGC z fail: '+str(len(fs[ws])))
	w = f['IMATCH'] == 3
	ws = fs['IMATCH'] == 3
	print('NGC close pairs: '+str(len(f[w])))
	print('SGC close pairs: '+str(len(fs[ws])))
	w = f['IMATCH'] == 9
	ws = fs['IMATCH'] == 9
	print('NGC bad type: '+str(len(f[w])))
	print('SGC bad type: '+str(len(fs[ws])))
	w = f['IMATCH'] == 13
	ws = fs['IMATCH'] == 13
	print('NGC legacy stars: '+str(len(f[w])))
	print('SGC legacy stars: '+str(len(fs[ws])))



def catzqso(ver='7_2',zmin=.8,zmax=2.2):
	samp = 'QSO'
	f = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_NGC_v'+ver+'.dat.fits')
	fs = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_SGC_v'+ver+'.dat.fits')
	wl = (f['IMATCH'] == 2) & (f['RANDOM_SELECT'] == 1)
	print(len(f[wl]))
	w = ((f['IMATCH'] == 1) | (f['IMATCH'] == 2)) & (f['Z'] > zmin) & (f['Z'] < zmax) & (f['RANDOM_SELECT'] == 1)
	ws = ((fs['IMATCH'] == 1) | (fs['IMATCH'] == 2)) & (fs['Z'] > zmin) & (fs['Z'] < zmax) & (fs['RANDOM_SELECT'] == 1)
	print(len(f[w]),len(fs[ws]))
	wnq = w & (f['IS_QSO_10K'] != 1)
	wnqs = ws & (fs['IS_QSO_10K'] != 1)
	wbz = w & (abs(f['Z']-f['Z_10K'])>0.01*(1+f['Z']))
	wbzs = ws & (abs(fs['Z']-fs['Z_10K'])>0.01*(1+fs['Z']))
	wbt = wnq | wbz
	wbts = wnqs | wbzs
	print(len(f[wnq]),len(fs[wnqs]))
	print(len(f[wbt]),len(fs[wbts]))
	print('total inspections '+str(len(f[w])+len(fs[ws])))
	print('total catastrophic failures '+str(len(f[wbt])+len(fs[wbts])))
	


def zgoodvsxy(samp,ver,nbin=20,ylim=(0.9,1.05)):
	if samp == 'QSO':
		pt = '^'
		pc = 'forestgreen'
		frac = 0.05
		ylab = r'quasar spectra:  $N_{\rm good}/N_{\rm total}$ '
	if samp == 'LRG':
		pt = 'o'
		pc = 'firebrick'
		frac = 0.02
		ylab = r'LRG spectra:  $N_{\rm good}/N_{\rm total}$ '
	f = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_NGC_v'+ver+'.dat.fits')
	fs = fitsio.read(ebossdir+'eBOSS_'+samp+'_full_SGC_v'+ver+'.dat.fits')
	#print(f.dtype.names)
	wgood = (f['IMATCH'] == 1) | (f['IMATCH'] == 4) | (f['IMATCH'] == 9)
	wgoods = (fs['IMATCH'] == 1) | (fs['IMATCH'] == 4) | (fs['IMATCH'] == 9)	
	wall = (f['IMATCH'] == 1) | (f['IMATCH'] == 4) | (f['IMATCH'] == 9) | (f['IMATCH'] == 7)
	walls = (fs['IMATCH'] == 1) | (fs['IMATCH'] == 4) | (fs['IMATCH'] == 9) | (fs['IMATCH'] == 7)
	
	hg = np.histogram(f[wgood]['XFOCAL'],bins=nbin)
	ha = np.histogram(f[wall]['XFOCAL'],bins=hg[1])	
	hgw = np.histogram(f[wgood]['XFOCAL'],bins=hg[1],weights=f[wgood]['WEIGHT_NOZ'])
	
	#print(f.dtype.names)
	hgs = np.histogram(fs[wgoods]['XFOCAL'],bins=hg[1])
	has = np.histogram(fs[walls]['XFOCAL'],bins=hg[1])	
	hgws = np.histogram(fs[wgoods]['XFOCAL'],bins=hg[1],weights=fs[wgoods]['WEIGHT_NOZ'])

	xl = np.zeros(nbin)
	for i in range(0,len(xl)):
		xl[i] = (hg[1][i]+hg[1][i+1])/2.
	
	#plt.plot(xl,(hg[0]+hgs[0])/(ha[0]+has[0]))
	nt = (ha[0]+has[0])
	err = (nt*frac)**.5/nt
	plt.errorbar(xl,(hgw[0]+hgws[0])/nt,err,fmt=pt,color=pc,label=r'with $w_{\rm noz}$')
	plt.plot(xl,(hg[0]+hgs[0])/(ha[0]+has[0]),'--',color='gray',label=r'without $w_{\rm noz}$')
	plt.legend()
	plt.plot(xl,np.ones(len(xl)),'k:')
	plt.ylim(ylim[0],ylim[1])
	plt.xlabel('focal plane x position (mm)')
	plt.ylabel(ylab)
	
	chi2x = np.sum(((hgw[0]+hgws[0])/nt-1.)**2./err**2.)
	print('chi2 for x focal null is '+str(chi2x))

	#plt.plot(xl,(hgw[0]+hgws[0])/(ha[0]+has[0]))
	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/n'+samp+'vsxfoc.png')

	plt.show()

	hgy = np.histogram(f[wgood]['YFOCAL'],bins=nbin)
	hay = np.histogram(f[wall]['YFOCAL'],bins=hgy[1])	
	hgwy = np.histogram(f[wgood]['YFOCAL'],bins=hgy[1],weights=f[wgood]['WEIGHT_NOZ'])
	
	#print(f.dtype.names)
	hgsy = np.histogram(fs[wgoods]['YFOCAL'],bins=hgy[1])
	hasy = np.histogram(fs[walls]['YFOCAL'],bins=hgy[1])	
	hgwsy = np.histogram(fs[wgoods]['YFOCAL'],bins=hgy[1],weights=fs[wgoods]['WEIGHT_NOZ'])

	yl = np.zeros(nbin)
	for i in range(0,len(yl)):
		yl[i] = (hgy[1][i]+hgy[1][i+1])/2.
	
	#plt.plot(xl,(hg[0]+hgs[0])/(ha[0]+has[0]))
	nty = (hay[0]+hasy[0])
	erry = (nty*frac)**.5/nty
	plt.errorbar(yl,(hgwy[0]+hgwsy[0])/nty,erry,fmt=pt,color=pc,label=r'with $w_{\rm noz}$')
	plt.plot(yl,(hgy[0]+hgsy[0])/(hay[0]+hasy[0]),'--',color='gray',label=r'without $w_{\rm noz}$')
	plt.legend()
	plt.plot(yl,np.ones(len(yl)),'k:')
	plt.ylim(ylim[0],ylim[1])
	plt.xlabel('focal plane y position (mm)')
	plt.ylabel(ylab)
	chi2y = np.sum(((hgwy[0]+hgwsy[0])/nty-1.)**2./erry**2.)
	print('chi2 for y focal null is '+str(chi2y))

	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/n'+samp+'vsyfoc.png')

	#plt.plot(xl,(hgw[0]+hgws[0])/(ha[0]+has[0]))

	plt.show()
	

def zsysplotsQSO(ver='7_1',frac=0.05,ylim=(0.9,1.05)):
	pt = '^'
	pc = 'forestgreen'
	
	nvsn = np.loadtxt(ebossdir+'nvsSPECSN2_IQSONGC'+ver+'.dat').transpose()
	nvss = np.loadtxt(ebossdir+'nvsSPECSN2_IQSOSGC'+ver+'.dat').transpose()
	nwt = nvsn[1]+nvss[1]
	nnwt = nvsn[3]+nvss[3]
	nt = nvsn[2]+nvss[2]
	err = (nt*nnwt/nt*(1.-nnwt/nt))**.5/nt
	normw = sum(nt)/sum(nwt)
	plt.errorbar(nvsn[0],nwt/nt,err,fmt=pt,color=pc,label=r'with $w_{\rm noz}$')
	chi2 = sum((1.-normw*nwt/nt)**2./err**2.)	
	norm = sum(nt)/sum(nnwt)
	chi2null = sum((1.-norm*nnwt/nt)**2./err**2.)	
	print(chi2,chi2null)
	plt.plot(nvsn[0],nnwt/nt,'--',color='gray',label=r'without $w_{\rm noz}$')
	plt.legend()
	plt.plot(nvsn[0],np.ones(len(nvsn[0])),'k:')
	plt.xlabel(r'Spectrograph $(S/N)^2$ in the i-band')
	plt.ylabel(r'quasar spectra:  $N_{\rm good}/N_{\rm total}$ ')
	plt.ylim(ylim[0],ylim[1])
	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/MNRAS_revision/nQSOvsSNI.png')
	plt.show()
	
	nvsn = np.loadtxt(ebossdir+'nvsIDQSONGC'+ver+'.dat').transpose()
	nvss = np.loadtxt(ebossdir+'nvsIDQSOSGC'+ver+'.dat').transpose()
	nwt = nvsn[1]+nvss[1]
	nnwt = nvsn[3]+nvss[3]
	nt = nvsn[2]+nvss[2]
	#err = (nt*frac)**.5/nt
	err = (nt*nnwt/nt*(1.-nnwt/nt))**.5/nt
	normw = sum(nt)/sum(nwt)
	plt.errorbar(nvsn[0],nwt/nt,err,fmt=pt,color=pc,label=r'with $w_{\rm noz}$')
	chi2 = sum((1.-normw*nwt/nt)**2./err**2.)	
	norm = sum(nt)/sum(nnwt)
	chi2null = sum((1.-norm*nnwt/nt)**2./err**2.)	
	print(chi2,chi2null)
	plt.plot(nvsn[0],nnwt/nt,'--',color='gray',label=r'without $w_{\rm noz}$')
	plt.plot(nvsn[0],np.ones(len(nvsn[0])),'k:')
	plt.xlabel('Fiber ID')
	plt.ylabel(r'quasar spectra:  $N_{\rm good}/N_{\rm total}$ ')
	plt.ylim(ylim[0],ylim[1])
	xl = [499,499]
	yl = [0,2]
	plt.plot(xl,yl,'k:')
	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/MNRAS_revision/nQSOvsID.png')
	plt.show()
	return True

def zsysplotsLRG(ver='7_2',frac=0.02,ylim=(0.9,1.05)):
	pt = 'o'
	pc = 'firebrick'
	
	nvsn = np.loadtxt(ebossdir+'nvsSPECSN2_ILRGNGC'+ver+'.dat').transpose()
	nvss = np.loadtxt(ebossdir+'nvsSPECSN2_ILRGSGC'+ver+'.dat').transpose()
	nwt = nvsn[1]+nvss[1]
	nnwt = nvsn[3]+nvss[3]
	nt = nvsn[2]+nvss[2]
	#err = (nt*frac)**.5/nt
	err = (nt*nnwt/nt*(1.-nnwt/nt))**.5/nt
	w = (err == 0)
	err[w] = 1.
	normw = sum(nt)/sum(nwt)
	plt.errorbar(nvsn[0][~w],nwt[~w]/nt[~w],err[~w],fmt=pt,color=pc,label=r'with $w_{\rm noz}$')
	chi2 = sum((1.-normw*nwt/nt)**2./err**2.)	
	norm = sum(nt)/sum(nnwt)
	chi2null = sum((1.-normw*norm*nnwt/nt)**2./err**2.)	
	print(chi2,chi2null)
	plt.plot(nvsn[0][~w],nnwt[~w]/nt[~w],'--',color='gray',label=r'without $w_{\rm noz}$')
	plt.legend()
	plt.plot(nvsn[0],np.ones(len(nvsn[0])),'k:')
	plt.xlabel(r'Spectrograph $(S/N)^2$ in the i-band')
	plt.ylabel(r'LRG spectra:  $N_{\rm good}/N_{\rm total}$')
	plt.ylim(ylim[0],ylim[1])
	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/MNRAS_revision/nLRGvsSNI.png')
	
	plt.show()
	
	nvsn = np.loadtxt(ebossdir+'nvsIDLRGNGC'+ver+'.dat').transpose()
	nvss = np.loadtxt(ebossdir+'nvsIDLRGSGC'+ver+'.dat').transpose()
	nwt = nvsn[1]+nvss[1]
	nnwt = nvsn[3]+nvss[3]
	nt = nvsn[2]+nvss[2]
	#err = (nt*frac)**.5/nt
	err = (nt*nnwt/nt*(1.-nnwt/nt))**.5/nt
	normw = sum(nt)/sum(nwt)
	plt.errorbar(nvsn[0],nwt/nt,err,fmt=pt,color=pc,label=r'with $w_{\rm noz}$')
	chi2 = sum((1.-normw*nwt/nt)**2./err**2.)	
	norm = sum(nt)/sum(nnwt)
	chi2null = sum((1.-norm*nnwt/nt)**2./err**2.)	
	print(chi2,chi2null)
	plt.plot(nvsn[0],nnwt/nt,'--',color='gray',label=r'without $w_{\rm noz}$')
	plt.plot(nvsn[0],np.ones(len(nvsn[0])),'k:')
	plt.xlabel('Fiber ID')
	plt.ylabel(r'LRG spectra:  $N_{\rm good}/N_{\rm total}$')
	plt.ylim(ylim[0],ylim[1])
	xl = [499,499]
	yl = [0,2]
	plt.plot(xl,yl,'k:')

	plt.savefig('/Users/ashleyross/eBOSSCatalogPaper/MNRAS_revision/nLRGvsID.png')
	plt.show()
	return True
	
	

def sysplotsQSO6pan(ver='7_1'):
	from numpy import loadtxt as load
	from numpy import ones,zeros
	#pp = PdfPages(ebossdir+'nQSOdr16vsys.pdf')

	fig = plt.figure(figsize=(8.5,7))
	lskn = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2nosysvSKY_I.dat').transpose()
	lsks = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2nosysvSKY_I.dat').transpose()
	lskt = (lskn[-2]/lskn[-1]**2.+lsks[-2]/lsks[-1]**2.)/(1./lskn[-1]**2.+1./lsks[-1]**2.)
	lske = (1./(1./lskn[-1]**2.+1./lsks[-1]**2.))**.5
	chi = 0
	for i in range(0,len(lskt)):
		chi += (lskt[i]-1.)**2./lske[i]**2.
	print(chi)	
	cskn = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2vSKY_I.dat').transpose()
	csks = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2vSKY_I.dat').transpose()
	cskt = (cskn[-2]/cskn[-1]**2.+csks[-2]/csks[-1]**2.)/(1./cskn[-1]**2.+1./csks[-1]**2.)
	cske = (1./(1./cskn[-1]**2.+1./csks[-1]**2.))**.5
	chi = 0
	for i in range(0,len(cskt)):
		chi += (cskt[i]-1.)**2./cske[i]**2.
	print(chi)	

	leen = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2nosysvEBV.dat').transpose()
	les = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2nosysvEBV.dat').transpose()
	cen = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2vEBV.dat').transpose()
	ces = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2vEBV.dat').transpose()
	let = (leen[1]/leen[2]**2.+les[1]/les[2]**2.)/(1./leen[2]**2.+1./les[2]**2.)
	lee = (1./(1./leen[2]**2.+1./les[2]**2.))**.5
	chi = 0
	for i in range(0,len(let)):
		chi += (let[i]-1.)**2./lee[i]**2.
		#print chi

	print(chi)
	
	cet = (cen[1]/cen[2]**2.+ces[1]/ces[2]**2.)/(1./cen[2]**2.+1./ces[2]**2.)
	cee = (1./(1./cen[2]**2.+1./ces[2]**2.))**.5
	chi = 0
	#chim = 0
	for i in range(0,len(cet)):
		chi += (cet[i]-1.)**2./cee[i]**2.
		#print chi
	print(chi)

	lan = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2nosysvAIRMASS.dat').transpose()
	las = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2nosysvAIRMASS.dat').transpose()
	can = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2vAIRMASS.dat').transpose()
	cas = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2vAIRMASS.dat').transpose()
		
	lat = (lan[1]/lan[2]**2.+las[1]/las[2]**2.)/(1./lan[2]**2.+1./las[2]**2.)
	lae = (1./(1./lan[2]**2.+1./las[2]**2.))**.5
	chi = 0
	for i in range(0,len(lat)):
		chi += (lat[i]-1.)**2./lae[i]**2.
	print(chi)	
	cat = (can[1]/can[2]**2.+cas[1]/cas[2]**2.)/(1./can[2]**2.+1./cas[2]**2.)
	cae = (1./(1./can[2]**2.+1./cas[2]**2.))**.5
	chi = 0
	for i in range(0,len(cat)):
		chi += (cat[i]-1.)**2./cae[i]**2.
	print(chi)	

	lsen = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2nosysvPSF_I.dat').transpose()
	lses = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2nosysvPSF_I.dat').transpose()
	csen = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2vPSF_I.dat').transpose()
	cses = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2vPSF_I.dat').transpose()
		
	lset = (lsen[1]/lsen[2]**2.+lses[1]/lses[2]**2.)/(1./lsen[2]**2.+1./lses[2]**2.)
	lsee = (1./(1./lsen[2]**2.+1./lses[2]**2.))**.5
	chi = 0
	for i in range(0,len(lset)):
		chi += (lset[i]-1.)**2./lsee[i]**2.
	print(chi)	
	cset = (csen[1]/csen[2]**2.+cses[1]/cses[2]**2.)/(1./csen[2]**2.+1./cses[2]**2.)
	csee = (1./(1./csen[2]**2.+1./cses[2]**2.))**.5
	chi = 0
	for i in range(0,len(cat)):
		chi += (cset[i]-1.)**2./csee[i]**2.
	print(chi)	

	lstn = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2nosysvstar.dat').transpose()
	lsts = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2nosysvstar.dat').transpose()
	cstn = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2vstar.dat').transpose()
	csts = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2vstar.dat').transpose()
		
	lstt = (lstn[1]/lstn[2]**2.+lsts[1]/lsts[2]**2.)/(1./lstn[2]**2.+1./lsts[2]**2.)
	lste = (1./(1./lstn[2]**2.+1./lsts[2]**2.))**.5
	chi = 0
	for i in range(0,len(lstt)):
		chi += (lstt[i]-1.)**2./lste[i]**2.
	print(chi)	
	cstt = (cstn[1]/cstn[2]**2.+csts[1]/csts[2]**2.)/(1./cstn[2]**2.+1./csts[2]**2.)
	cste = (1./(1./cstn[2]**2.+1./csts[2]**2.))**.5
	chi = 0
	for i in range(0,len(cstt)):
		chi += (cstt[i]-1.)**2./cste[i]**2.
	print(chi)	

	ldn = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2nosysvDEPTH_G.dat').transpose()
	lds = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2nosysvDEPTH_G.dat').transpose()
	cdn = load(ebossdir+'ngebossQSO_NGC'+ver+'_mz0.8xz2.2vDEPTH_G.dat').transpose()
	cds = load(ebossdir+'ngebossQSO_SGC'+ver+'_mz0.8xz2.2vDEPTH_G.dat').transpose()
		
	ldt = (ldn[1]/ldn[2]**2.+lds[1]/lds[2]**2.)/(1./ldn[2]**2.+1./lds[2]**2.)
	lde = (1./(1./ldn[2]**2.+1./lds[2]**2.))**.5
	chi = 0
	for i in range(0,len(ldt)):
		chi += (ldt[i]-1.)**2./lde[i]**2.
	print(chi)	
	cdt = (cdn[1]/cdn[2]**2.+cds[1]/cds[2]**2.)/(1./cdn[2]**2.+1./cds[2]**2.)
	cde = (1./(1./cdn[2]**2.+1./cds[2]**2.))**.5
	chi = 0
	for i in range(0,len(cdt)):
		chi += (cdt[i]-1.)**2./cde[i]**2.
	print(chi)	


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
	
	pt = '^'
	pc = 'forestgreen'
	lc = 'k'
	ax.xaxis.set_ticks(np.arange(5,20,5.))
	ax.set_xlabel('i-band sky background (nm)',size=16)
	ax.plot(lskn[0],ols,':',color='gray')
	ax.plot(lskn[0],lskt,'--',color=lc) 	
	ax.errorbar(cskn[0]-.07,cskt,cske,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor=pc)
	cnr = [(5,14.5),(.78,.88)]
	plotbox(cnr)
	xl = [5.3,6.9]
	yl = [.81,.81]
	ax.plot(xl,yl,'--',color=lc)
	ax.text(7.2,.8,'raw',color=lc,size=16)
	xl = [5.9]
	yl = [.845]
	ax.plot(xl,yl,pt,color=pc,markeredgecolor=pc)
	ax.text(7.2,.835,'corrected',color=pc,size=16)
	ax2 = fig.add_subplot(2,3,2,sharey=ax)
	ax2.minorticks_on()
# 	for tick in ax2.xaxis.get_major_ticks():
# 		tick.label.set_fontsize(8)
	ax2.plot(leen[0],ols,'k:')
	ax2.plot(leen[0]+.001,let,'--',color=lc) 
	ax2.errorbar(cen[0]-.001,cet,cee,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor=pc)
	ax2.set_xlabel('E[B-V] (mag)',size=16)
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
	ax3.set_xlabel('airmass',size=16)
	ax3.plot(lan[0],ols,'k:')
	ax3.plot(lan[0]+.007,lat,'--',color=lc)
	ax3.errorbar(can[0]-.007,cat,cae,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor=pc)
	ax3.set_xlim(1.01,2)
	ax3.xaxis.set_ticks(np.arange(1.2,2.,.2))
	ax3.set_ylim(.751,1.2)
	for ylabel_i in ax3.axes.get_yticklabels():
		ylabel_i.set_visible(False)

	ax4 = fig.add_subplot(2,3,4)
	ax4.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax4.set_xlabel(r'i-band seeing ($^{\prime\prime}$)',size=16)
	ax4.plot(lsen[0],ols,'k:')
	ax4.plot(lsen[0]+.007,lset,'--',color=lc)
	ax4.errorbar(csen[0]-.007,cset,csee,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor=pc)
	ax4.set_xlim(0.7,1.99)
	ax4.set_ylim(.751,1.2)
	ax4.text(.25,1.6,r'quasar number density / average number density',size=16,rotation='vertical')
	ax5 = fig.add_subplot(2,3,5,sharey=ax4)
	ax5.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax5.set_xlabel(r'$N_{\rm star}/{\rm deg}^2$',size=16)
	ax5.plot(lstn[0]/0.0524558,ols,'k:')
	ax5.plot(lstn[0]/0.0524558,lstt,'--',color=lc)
	ax5.errorbar(cstn[0]/0.0524558,cstt,cste,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor=pc)
	ax5.set_xlim(31/0.0524558,299/0.0524558)
	ax5.set_ylim(.751,1.2)
	for ylabel_5 in ax5.axes.get_yticklabels():
		ylabel_5.set_visible(False)

	ax6 = fig.add_subplot(2,3,6,sharey=ax4)
	ax6.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax6.set_xlabel(r'$g$-band depth (mag)',size=16)
	ax6.plot(ldn[0],ols,'k:')
	ax6.plot(ldn[0]+.007,ldt,'--',color=lc)
	ax6.errorbar(cdn[0]-.007,cdt,cde,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor=pc)
	#ax6.set_xlim(21.81,22.9)
	ax6.set_xlim(22.3,23.7)
	ax6.set_ylim(.751,1.2)
	for ylabel_6 in ax6.axes.get_yticklabels():
		ylabel_6.set_visible(False)


	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=.3)
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

	fig.savefig('/Users/ashleyross/eBOSSCatalogPaper/nQSOv'+ver+'vsys.png')
	#pp.close()
	return True

def sysplotsLRG6pan(ver='7'):
	from numpy import loadtxt as load
	from numpy import ones,zeros
	#pp = PdfPages(ebossdir+'nQSOdr16vsys.pdf')

	fig = plt.figure(figsize=(8.5,7))
	lskn = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0nosysvSKY_I.dat').transpose()
	lsks = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0nosysvSKY_I.dat').transpose()
	lskt = (lskn[-2]/lskn[-1]**2.+lsks[-2]/lsks[-1]**2.)/(1./lskn[-1]**2.+1./lsks[-1]**2.)
	lske = (1./(1./lskn[-1]**2.+1./lsks[-1]**2.))**.5
	chi = 0
	for i in range(0,len(lskt)):
		chi += (lskt[i]-1.)**2./lske[i]**2.
	print(chi)	
	cskn = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0vSKY_I.dat').transpose()
	csks = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0vSKY_I.dat').transpose()
	cskt = (cskn[-2]/cskn[-1]**2.+csks[-2]/csks[-1]**2.)/(1./cskn[-1]**2.+1./csks[-1]**2.)
	cske = (1./(1./cskn[-1]**2.+1./csks[-1]**2.))**.5
	chi = 0
	for i in range(0,len(cskt)):
		chi += (cskt[i]-1.)**2./cske[i]**2.
	print(chi)	

	leen = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0nosysvEBV.dat').transpose()
	les = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0nosysvEBV.dat').transpose()
	cen = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0vEBV.dat').transpose()
	ces = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0vEBV.dat').transpose()
	let = (leen[1]/leen[2]**2.+les[1]/les[2]**2.)/(1./leen[2]**2.+1./les[2]**2.)
	lee = (1./(1./leen[2]**2.+1./les[2]**2.))**.5
	chi = 0
	for i in range(0,len(let)):
		chi += (let[i]-1.)**2./lee[i]**2.
		#print chi

	print(chi)
	
	cet = (cen[1]/cen[2]**2.+ces[1]/ces[2]**2.)/(1./cen[2]**2.+1./ces[2]**2.)
	cee = (1./(1./cen[2]**2.+1./ces[2]**2.))**.5
	chi = 0
	#chim = 0
	for i in range(0,len(cet)):
		chi += (cet[i]-1.)**2./cee[i]**2.
		#print chi
	print(chi)

	lan = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0nosysvAIRMASS.dat').transpose()
	las = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0nosysvAIRMASS.dat').transpose()
	can = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0vAIRMASS.dat').transpose()
	cas = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0vAIRMASS.dat').transpose()
		
	lat = (lan[1]/lan[2]**2.+las[1]/las[2]**2.)/(1./lan[2]**2.+1./las[2]**2.)
	lae = (1./(1./lan[2]**2.+1./las[2]**2.))**.5
	chi = 0
	for i in range(0,len(lat)):
		chi += (lat[i]-1.)**2./lae[i]**2.
	print(chi)	
	cat = (can[1]/can[2]**2.+cas[1]/cas[2]**2.)/(1./can[2]**2.+1./cas[2]**2.)
	cae = (1./(1./can[2]**2.+1./cas[2]**2.))**.5
	chi = 0
	for i in range(0,len(cat)):
		chi += (cat[i]-1.)**2./cae[i]**2.
	print(chi)	

	lsen = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0nosysvPSF_I.dat').transpose()
	lses = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0nosysvPSF_I.dat').transpose()
	csen = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0vPSF_I.dat').transpose()
	cses = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0vPSF_I.dat').transpose()
		
	lset = (lsen[1]/lsen[2]**2.+lses[1]/lses[2]**2.)/(1./lsen[2]**2.+1./lses[2]**2.)
	lsee = (1./(1./lsen[2]**2.+1./lses[2]**2.))**.5
	chi = 0
	for i in range(0,len(lset)):
		chi += (lset[i]-1.)**2./lsee[i]**2.
	print(chi)	
	cset = (csen[1]/csen[2]**2.+cses[1]/cses[2]**2.)/(1./csen[2]**2.+1./cses[2]**2.)
	csee = (1./(1./csen[2]**2.+1./cses[2]**2.))**.5
	chi = 0
	for i in range(0,len(cat)):
		chi += (cset[i]-1.)**2./csee[i]**2.
	print(chi)	

	lstn = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0nosysvstar.dat').transpose()
	lsts = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0nosysvstar.dat').transpose()
	cstn = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0vstar.dat').transpose()
	csts = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0vstar.dat').transpose()
		
	lstt = (lstn[1]/lstn[2]**2.+lsts[1]/lsts[2]**2.)/(1./lstn[2]**2.+1./lsts[2]**2.)
	lste = (1./(1./lstn[2]**2.+1./lsts[2]**2.))**.5
	chi = 0
	for i in range(0,len(lstt)):
		chi += (lstt[i]-1.)**2./lste[i]**2.
	print(chi)	
	cstt = (cstn[1]/cstn[2]**2.+csts[1]/csts[2]**2.)/(1./cstn[2]**2.+1./csts[2]**2.)
	cste = (1./(1./cstn[2]**2.+1./csts[2]**2.))**.5
	chi = 0
	for i in range(0,len(cstt)):
		chi += (cstt[i]-1.)**2./cste[i]**2.
	print(chi)	

	ldn = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0nosysvDEPTH_G.dat').transpose()
	lds = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0nosysvDEPTH_G.dat').transpose()
	cdn = load(ebossdir+'ngebossLRG_NGC'+ver+'_mz0.6xz1.0vDEPTH_G.dat').transpose()
	cds = load(ebossdir+'ngebossLRG_SGC'+ver+'_mz0.6xz1.0vDEPTH_G.dat').transpose()
		
	ldt = (ldn[1]/ldn[2]**2.+lds[1]/lds[2]**2.)/(1./ldn[2]**2.+1./lds[2]**2.)
	lde = (1./(1./ldn[2]**2.+1./lds[2]**2.))**.5
	chi = 0
	for i in range(0,len(ldt)):
		chi += (ldt[i]-1.)**2./lde[i]**2.
	print(chi)	
	cdt = (cdn[1]/cdn[2]**2.+cds[1]/cds[2]**2.)/(1./cdn[2]**2.+1./cds[2]**2.)
	cde = (1./(1./cdn[2]**2.+1./cds[2]**2.))**.5
	chi = 0
	for i in range(0,len(cdt)):
		chi += (cdt[i]-1.)**2./cde[i]**2.
	print(chi)	


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
	pt = 'o'
	pc = 'firebrick'
	lc = 'k'
	ax.xaxis.set_ticks(np.arange(5,20,5.))
	ax.set_xlabel('i-band sky background (nm)',size=16)
	ax.plot(lskn[0],ols,':',color='gray')
	ax.plot(lskn[0],lskt,'--',color=lc) 	
	ax.errorbar(cskn[0]-.07,cskt,cske,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor='k')
	cnr = [(5,14.5),(.78,.88)]
	plotbox(cnr)
	xl = [5.3,6.9]
	yl = [.81,.81]
	ax.plot(xl,yl,'--',color=lc)
	ax.text(7.2,.8,'raw',color=lc,size=16)
	xl = [5.9]
	yl = [.845]
	ax.plot(xl,yl,pt,color=pc,markeredgecolor='k')
	ax.text(7.2,.835,'corrected',color=pc,size=16)
	ax2 = fig.add_subplot(2,3,2,sharey=ax)
	ax2.minorticks_on()
# 	for tick in ax2.xaxis.get_major_ticks():
# 		tick.label.set_fontsize(8)
	ax2.plot(leen[0],ols,'k:')
	ax2.plot(leen[0]+.001,let,'--',color=lc) 
	ax2.errorbar(cen[0]-.001,cet,cee,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor='k')
	ax2.set_xlabel('E[B-V] (mag)',size=16)
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
	ax3.set_xlabel('airmass',size=16)
	ax3.plot(lan[0],ols,'k:')
	ax3.plot(lan[0]+.007,lat,'--',color=lc)
	ax3.errorbar(can[0]-.007,cat,cae,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor='k')
	ax3.set_xlim(1.01,2)
	ax3.xaxis.set_ticks(np.arange(1.2,2.,.2))
	ax3.set_ylim(.751,1.2)
	for ylabel_i in ax3.axes.get_yticklabels():
		ylabel_i.set_visible(False)

	ax4 = fig.add_subplot(2,3,4)
	ax4.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax4.set_xlabel(r'i-band seeing ($^{\prime\prime}$)',size=16)
	ax4.plot(lsen[0],ols,'k:')
	ax4.plot(lsen[0]+.007,lset,'--',color=lc)
	ax4.errorbar(csen[0]-.007,cset,csee,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor='k')
	ax4.set_xlim(0.7,1.99)
	ax4.set_ylim(.751,1.2)
	ax4.text(.25,1.6,r'LRG number density / average number density',size=16,rotation='vertical')
	ax5 = fig.add_subplot(2,3,5,sharey=ax4)
	ax5.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax5.set_xlabel(r'$N_{\rm star}/{\rm deg}^2$',size=16)
	ax5.plot(lstn[0]/0.0524558,ols,'k:')
	ax5.plot(lstn[0]/0.0524558,lstt,'--',color=lc)
	ax5.errorbar(cstn[0]/0.0524558,cstt,cste,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor='k')
	ax5.set_xlim(31/0.0524558,299/0.0524558)
	ax5.set_ylim(.751,1.2)
	for ylabel_5 in ax5.axes.get_yticklabels():
		ylabel_5.set_visible(False)

	ax6 = fig.add_subplot(2,3,6,sharey=ax4)
	ax6.minorticks_on()
	#for tick in ax3.xaxis.get_major_ticks():
	#	tick.label.set_fontsize(8)
	ax6.set_xlabel(r'$g$-band depth (mag)',size=16)
	ax6.plot(ldn[0],ols,'k:')
	ax6.plot(ldn[0]+.007,ldt,'--',color=lc)
	ax6.errorbar(cdn[0]-.007,cdt,cde,fmt=pt,color=pc,markersize=7,elinewidth=2,markeredgecolor='k')
	#ax6.set_xlim(21.81,22.9)
	ax6.set_xlim(22.3,23.7)
	ax6.set_ylim(.751,1.2)
	for ylabel_6 in ax6.axes.get_yticklabels():
		ylabel_6.set_visible(False)


	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=.3)
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

	fig.savefig('/Users/ashleyross/eBOSSCatalogPaper/nLRGv'+ver+'vsys.png')
	#pp.close()
	return True


def plotObsdata(reg='NGC',mjdmin=56659,ver='7',elgper=0.05):
	#cdir = '/uufs/chpc.utah.edu/common/home/sdss02/ebosswork/eboss/sandbox/lss/catalogs/input/'
	cmgaln = fitsio.read(ebossdir+'galaxy_DR12v5_CMASS_North.fits.gz')
	plt.plot(cmgaln['RA'],cmgaln['DEC'],',',color='gray',zorder=0)
	cmgals = fitsio.read(ebossdir+'galaxy_DR12v5_CMASS_South.fits.gz')
	ral = cmgals['RA']
	ral[ral > 300] -=360
	plt.plot(ral,cmgals['DEC'],',',color='gray',zorder=1)

	#lrgf = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'+ver+'/eBOSS_LRG_full_noveto_ALL_v'+ver+'.dat.fits')
	lrgf = fitsio.read(ebossdir+'eBOSS_LRG_full_noveto_ALL_v'+ver+'.dat.fits')
	wlobs = ((lrgf['IMATCH'] == 1) | (lrgf['IMATCH'] == 4) | (lrgf['IMATCH'] == 7) | (lrgf['IMATCH'] == 9) ) & (lrgf['MJD'] > 56659)
	chunkl = ['boss214','boss217','eboss1','eboss16','eboss17','eboss18','eboss2','eboss20','eboss24','eboss26','eboss27','eboss3','eboss4','eboss5','eboss6','eboss7','eboss9']
	wlobs &= np.in1d(lrgf['CHUNKspec'],chunkl)
	print((len(lrgf[wlobs])))
	wlnobs = (lrgf['IMATCH'] == 0) & (lrgf['MJD'] > 56659)
	ral = lrgf['RA']
	#if reg == 'SGC':
	ral[ral > 300] -=360
	plt.plot(ral[wlnobs],lrgf[wlnobs]['DEC'],'k,')
	#qsof = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'+ver+'/eBOSS_QSO_full_noveto_ALL_v'+ver+'.dat.fits')
	qsof = fitsio.read(ebossdir+'eBOSS_QSO_full_noveto_ALL_v'+ver+'.dat.fits')
	wqobs = ((qsof['IMATCH'] == 1) | (qsof['IMATCH'] == 4) | (qsof['IMATCH'] == 7) | (qsof['IMATCH'] == 9) ) & (qsof['MJD'] > 56659)
	chunkl = ['boss214','boss217','eboss1','eboss16','eboss17','eboss18','eboss2','eboss20','eboss24','eboss26','eboss27','eboss3','eboss4','eboss5','eboss6','eboss7','eboss9']
	wqobs &= np.in1d(qsof['CHUNK_DR16Q'],chunkl)
	print((len(qsof[wqobs])))
	wqnobs = (qsof['IMATCH'] == 0) & (qsof['MJD'] > 56659)
	raq = qsof['RA']
	#if reg == 'SGC':
	raq[raq > 300] -=360
	plt.plot(raq[wqnobs],qsof[wqnobs]['DEC'],'k,')
	plt.plot(ral[wlobs],lrgf[wlobs]['DEC'],'r,')	
	plt.plot(raq[wqobs],qsof[wqobs]['DEC'],'y,')
	#elgf = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'+ver+'/eBOSS_ELG_full_ALL_v'+ver+'.dat.fits')
	#elgf = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/7/input/ELG/ELG.v5_13_0.rrv3bis.all.fits')
	elgf = fitsio.read(ebossdir+'ELG.v5_13_0.rrv3bis.all.fits')
	rae = elgf['ra']
	#if reg == 'SGC':
	rae[rae > 300] -=360
	#weobs = ((elgf['IMATCH'] == 1) | (elgf['IMATCH'] == 4) | (elgf['IMATCH'] == 7) | (elgf['IMATCH'] == 9) ) 
	weobs = elgf['hasfiber'] == 1
	print((len(elgf[weobs])))
	randl = np.random.rand(len(elgf))
	weobs &= (randl < elgper)
	print((len(elgf[weobs])))
	plt.plot(rae[weobs],elgf[weobs]['dec'],'b,')
	plt.ylim(-12,75)
	plt.xlabel('right ascension (J2000)')
	plt.ylabel('declination (J2000)')
	plt.savefig('DR16obstarg.png')
	plt.show()
	return True

def plotObsdata_sindec(reg='NGC',mjdmin=56659,ver='7',elgper=0.2):
	#cdir = '/uufs/chpc.utah.edu/common/home/sdss02/ebosswork/eboss/sandbox/lss/catalogs/input/'
	cmgaln = fitsio.read(ebossdir+'galaxy_DR12v5_CMASS_North.fits.gz')
	
	
	decl = np.sin(cmgaln['DEC']*np.pi/180.)
	plt.plot(cmgaln['RA'],decl,',',color='gray',zorder=0)
	cmgals = fitsio.read(ebossdir+'galaxy_DR12v5_CMASS_South.fits.gz')
	ral = cmgals['RA']
	ral[ral > 300] -=360
	decl = np.sin(cmgals['DEC']*np.pi/180.)
	plt.plot(ral,decl,',',color='gray',zorder=1)

	#lrgf = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'+ver+'/eBOSS_LRG_full_noveto_ALL_v'+ver+'.dat.fits')
	lrgf = fitsio.read(ebossdir+'eBOSS_LRG_full_noveto_ALL_v'+ver+'.dat.fits')
	wlobs = ((lrgf['IMATCH'] == 1) | (lrgf['IMATCH'] == 4) | (lrgf['IMATCH'] == 7) | (lrgf['IMATCH'] == 9) ) & (lrgf['MJD'] > 56659)
	chunkl = [b'boss214',b'boss217',b'eboss1',b'eboss16',b'eboss17',b'eboss18',b'eboss2',b'eboss20',b'eboss24',b'eboss26',b'eboss27',b'eboss3',b'eboss4',b'eboss5',b'eboss6',b'eboss7',b'eboss9']
	#chunkl = np.unique(lrgf['CHUNKspec'])
	print(lrgf['CHUNKspec'].dtype)
	lrgf['CHUNKspec'] = lrgf['CHUNKspec'].astype('U')
	print(np.unique(lrgf['CHUNKspec']))
	wlobs &= np.in1d(lrgf['CHUNKspec'],chunkl)
	print((len(lrgf[wlobs])))
	wlnobs = (lrgf['IMATCH'] == 0) & (lrgf['MJD'] > 56659)
	ral = lrgf['RA']
	#if reg == 'SGC':
	ral[ral > 300] -=360
	decl = np.sin(lrgf[wlnobs]['DEC']*np.pi/180.)
	plt.plot(ral[wlnobs],decl,'k,')
	#qsof = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'+ver+'/eBOSS_QSO_full_noveto_ALL_v'+ver+'.dat.fits')
	qsof = fitsio.read(ebossdir+'eBOSS_QSO_full_noveto_ALL_v'+ver+'.dat.fits')
	wqobs = ((qsof['IMATCH'] == 1) | (qsof['IMATCH'] == 4) | (qsof['IMATCH'] == 7) | (qsof['IMATCH'] == 9) ) & (qsof['MJD'] > 56659)
	chunkl = [b'boss214',b'boss217',b'eboss1',b'eboss16',b'eboss17',b'eboss18',b'eboss2',b'eboss20',b'eboss24',b'eboss26',b'eboss27',b'eboss3',b'eboss4',b'eboss5',b'eboss6',b'eboss7',b'eboss9']
	print(np.unique(qsof['CHUNK_DR16Q']))
	wqobs &= np.in1d(qsof['CHUNK_DR16Q'],chunkl)
	print((len(qsof[wqobs])))
	wqnobs = (qsof['IMATCH'] == 0) & (qsof['MJD'] > 56659)
	raq = qsof['RA']
	#if reg == 'SGC':
	raq[raq > 300] -=360
	decl = np.sin(qsof[wqnobs]['DEC']*np.pi/180.)
	plt.plot(raq[wqnobs],decl,'k,')
	decl = np.sin(lrgf[wlobs]['DEC']*np.pi/180.)
	plt.plot(ral[wlobs],decl,'r,')	
	decl = np.sin(qsof[wqobs]['DEC']*np.pi/180.)
	plt.plot(raq[wqobs],decl,'y,')
	#elgf = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/'+ver+'/eBOSS_ELG_full_ALL_v'+ver+'.dat.fits')
	#elgf = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/sandbox/lss/catalogs/versions/7/input/ELG/ELG.v5_13_0.rrv3bis.all.fits')
	elgf = fitsio.read(ebossdir+'ELG.v5_13_0.rrv3bis.all.fits')
	rae = elgf['ra']
	#if reg == 'SGC':
	rae[rae > 300] -=360
	#weobs = ((elgf['IMATCH'] == 1) | (elgf['IMATCH'] == 4) | (elgf['IMATCH'] == 7) | (elgf['IMATCH'] == 9) ) 
	weobs = elgf['hasfiber'] == 1
	print((len(elgf[weobs])))
	randl = np.random.rand(len(elgf))
	weobs &= (randl < elgper)
	print((len(elgf[weobs])))
	decl = np.sin(elgf[weobs]['dec']*np.pi/180.)
	plt.plot(rae[weobs],decl,'b,')
	plt.ylim(np.sin(-12*np.pi/180.),np.sin(75*np.pi/180.))
	plt.xlim(-60,300)
	plt.xlabel('right ascension J2000 (degrees)',size=16,labelpad=-.75)
	plt.ylabel('declination J2000 (degrees)',size=16)
	plt.yticks([])
	latl = [np.sin(60*pi/180.),np.sin(40*pi/180.),np.sin(20*pi/180.),0]
	xl = [-60,350]
	yl = [latl[0],latl[0]]
	plt.text(-58,np.sin(62*pi/180.),r'60$^o$')
	plt.plot(xl,yl,'k:')
	yl = [latl[1],latl[1]]
	plt.plot(xl,yl,'k:')
	plt.text(-58,np.sin(41.5*pi/180.),r'40$^o$')
	yl = [latl[2],latl[2]]
	plt.plot(xl,yl,'k:')
	plt.text(-58,np.sin(21*pi/180.),r'20$^o$')
	yl = [latl[3],latl[3]]
	plt.plot(xl,yl,'k:')
	plt.text(-58,np.sin(.5*pi/180.),r'0$^o$')

	plt.savefig('DR16obstargsindec.png',dpi=300)
	plt.show()
	return True

def plotfoot_sindec(infile,outfile,cmap=cm.coolwarm,tit=''):
	from matplotlib import pyplot as plt
	from matplotlib import rc

	from matplotlib.backends.backend_pdf import PdfPages
	import fitsio   
	from mpl_toolkits.basemap import Basemap
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	plt.rc('font', family='serif', size=16)
	plt.clf()
	fig = plt.figure(figsize=(12,6.7))
	ax = fig.add_subplot(111)
	for file in infile:
		f = fitsio.read(file)
		ral = f['RA']
		ral[ral > 300] -=360
		decl = np.sin(f['DEC']*np.pi/180.)
		map = plt.scatter(ral,decl,c=f['COMP_BOSS'],s=.1,cmap=cmap,vmin=.5,lw=0)
	#cbar = plt.colorbar()
	#cbar.set_label(r'$n_g$ (arcmin$^{-2}$)', rotation=270,fontsize=14)
	plt.xlabel('right ascension J2000 (degrees)',fontsize=16)
	plt.ylabel('declination J2000 (degrees)',fontsize=16)

	plt.yticks([])
	latl = [np.sin(60*pi/180.),np.sin(40*pi/180.),np.sin(20*pi/180.),0]
	xl = [-60,350]
	yl = [latl[0],latl[0]]
	plt.text(-58,np.sin(61.5*pi/180.),r'60$^o$')
	plt.plot(xl,yl,'k:')
	yl = [latl[1],latl[1]]
	plt.plot(xl,yl,'k:')
	plt.text(-58,np.sin(41.*pi/180.),r'40$^o$')
	yl = [latl[2],latl[2]]
	plt.plot(xl,yl,'k:')
	plt.text(-58,np.sin(21*pi/180.),r'20$^o$')
	yl = [latl[3],latl[3]]
	plt.plot(xl,yl,'k:')
	plt.text(-58,np.sin(.5*pi/180.),r'0$^o$')
	plt.xlim(-60,300)
	plt.ylim(np.sin(-12*np.pi/180.),np.sin(70*np.pi/180.))
	plt.gca().yaxis.set_label_coords(-.05,.5)
	plt.gca().xaxis.set_label_coords(.5,-.05)
	axins=inset_axes(plt.gca(),width="40%",height="3%",loc=4,borderpad=1.)
	#axins=inset_axes(plt.gca(),width="40%",height="3%",position=(.6,.2))
	plt.colorbar(map,cax=axins,orientation='horizontal',ticks=[.5, .6, .7, .8, .9, 1.])
	plt.title(tit)

	plt.savefig(outfile)
	#plt.show()
	#pp.savefig()
	#pp.close()
	return True


def getchunkstats():	
	chunkl = ['boss214','boss217','eboss1','eboss16','eboss17','eboss18','eboss2','eboss20','eboss21','eboss22','eboss23','eboss24','eboss25','eboss26','eboss27','eboss3','eboss4','eboss5','eboss6','eboss7','eboss9']
	chunklp = ['boss214       ','boss217       ','eboss1        ','eboss16       ','eboss17       ','eboss18       ','eboss2        ','eboss20       ','eboss21       ','eboss22       ','eboss23       ','eboss24       ','eboss25       ','eboss26       ','eboss27       ','eboss3        ','eboss4        ','eboss5        ','eboss6        ','eboss7        ','eboss9        ']

	fp = fitsio.read('/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/v5_13_0/platelist.fits')
	wp = fp['PLATEQUALITY'] == 'good'
	fp = fp[wp]
	for i in range(0,len(chunkl)):
		chunk = chunkl[i]
		if chunk == 'boss214' or chunk == 'boss217':
			dir = '/uufs/chpc.utah.edu/common/home/sdss05/software/svn.sdss.org/repo/eboss/ebosstilelist/trunk/outputs/eboss0'
		else:
			dir = '/uufs/chpc.utah.edu/common/home/sdss05/software/svn.sdss.org/repo/eboss/ebosstilelist/trunk/outputs/'+chunk
		f = fitsio.read(dir+'/final-'+chunk+'.fits')
		ntiletot = np.unique(f['TILE'])
		wt = ntiletot > 0
		
		chunkp = chunklp[i]
		wc = fp['CHUNK'] == chunkp
		ntileobs = np.unique(fp[wc]['TILEID'])
		print((chunk,len(ntiletot[wt]),len(ntileobs)))
	return True
	
		
			
		

def BAOsigreccompELG(cp='w',cd='orange'):
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	d = np.loadtxt('/Users/ashleyross/Dropbox/eBOSS/BAO1D_EZmocks_ELG70.5933.06.010.015.005st0.dat').transpose() #pre recon
	#dx = np.loadtxt('/Users/ashleyross/Dropbox/BAO-DES/Halogen/V2.0.0_DNF/3D_likelihoods/ximuwSM_sigt8sigr8_z0.6-1.0_bs8_rmin30/BAOchi2fitresultsz0.6-1.0_bs8_s88_rmin30.dat').transpose()
	fig = plt.figure(figsize=(7,5))
	dp = np.loadtxt('/Users/ashleyross/Dropbox/eBOSS/BAO1D_EZmocks_ELG7_rec0.59304.07.015.01.05st0.dat').transpose() #post recon

	plt.clf()
	plt.minorticks_on()
	pp = PdfPages('/Users/ashleyross/Dropbox/eboss/ELGBAOprepostsigDR16.pdf')
	#plt.plot((dx[2]-dx[1])/2.,dw[2],'o',color=cp,markeredgecolor='k')
	plt.plot(d[1],dp[1],'o',color=cp,markeredgecolor='b')
	#plt.plot(dezxi[0],dezpk[0],'o',color=EZcol,markeredgecolor='k')
	#plt.plot(dq[1],dq[2],'^',color=QPMcol,markeredgecolor='k')
	xl = [0,2]
	yl = [0,2]
	plt.plot(xl,yl,'--k')
	xl = [0.061]
	yl = [0.036]
	plt.plot(xl,yl,'*',color=cd,markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	#xl = [0.036]
	#yl = [0.037]
	#plt.plot(xl,yl,'*',color='y',markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	#xl = [0.050]
	#yl = [0.047]
	#plt.plot(xl,yl,'*',color='teal',markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	
	plt.xlim( 0.0, 0.1 )
	plt.ylim(0.0,0.1)
	plt.xlabel (r' $\sigma(\alpha)$ pre recon', fontsize=18)
	plt.ylabel (r' $\sigma(\alpha)$ post recon', fontsize=18)
	plt.title('DR16 ELGs BAO')
	pp.savefig()
	pp.close()
	return True

def BAOsigQSOraratcomp(cp='w',cd='orange'):
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	
	d = np.loadtxt('/Users/ashleyross/Dropbox/eboss/BAO2D_EZmocks_eBOSSQSODR16.dat').transpose() #pre recon
	#dx = np.loadtxt('/Users/ashleyross/Dropbox/BAO-DES/Halogen/V2.0.0_DNF/3D_likelihoods/ximuwSM_sigt8sigr8_z0.6-1.0_bs8_rmin30/BAOchi2fitresultsz0.6-1.0_bs8_s88_rmin30.dat').transpose()
	fig = plt.figure(figsize=(7,5))

	plt.clf()
	plt.minorticks_on()
	pp = PdfPages('/Users/ashleyross/Dropbox/eboss/BAOQSOarat_sysmock.pdf')
	#plt.plot((dx[2]-dx[1])/2.,dw[2],'o',color=cp,markeredgecolor='k')
	plt.plot(d[-3],d[-2],'o',color=cp,markeredgecolor='k')
	#plt.plot(dezxi[0],dezpk[0],'o',color=EZcol,markeredgecolor='k')
	#plt.plot(dq[1],dq[2],'^',color=QPMcol,markeredgecolor='k')
	xl = [0,2]
	yl = [0,2]
	#plt.plot(xl,yl,'--k')
	xl = [0.042]
	yl = [0.025]
	plt.plot(xl,yl,'*',color=cd,markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	#xl = [0.036]
	#yl = [0.037]
	#plt.plot(xl,yl,'*',color='y',markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	#xl = [0.050]
	#yl = [0.047]
	#plt.plot(xl,yl,'*',color='teal',markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	
	plt.xlim( 0.0, 0.1 )
	plt.ylim(0.0,0.1)
	plt.xlabel (r' $\sigma(\alpha_{||})$ ', fontsize=18)
	plt.ylabel (r' $\sigma(\alpha_{\perp})$ ', fontsize=18)
	plt.title('DR16 quasar BAO')
	pp.savefig()
	pp.close()
	return True


def BAOsigSGC(cp='w',cd='orange'):
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	d = np.loadtxt('/Users/ashleyross/Dropbox/eboss/BAOfitsELGEZSGC4rec0.59304.07.015.01.08st0.dat').transpose() #pre recon
	#dx = np.loadtxt('/Users/ashleyross/Dropbox/BAO-DES/Halogen/V2.0.0_DNF/3D_likelihoods/ximuwSM_sigt8sigr8_z0.6-1.0_bs8_rmin30/BAOchi2fitresultsz0.6-1.0_bs8_s88_rmin30.dat').transpose()
	fig = plt.figure(figsize=(7,5))
	dp = np.loadtxt('/Users/ashleyross/Dropbox/eboss/BAOfitsELGEZNScombf4rec0.59304.07.015.01.08st0.dat').transpose() #post recon

	plt.clf()
	plt.minorticks_on()
	pp = PdfPages('/Users/ashleyross/Dropbox/ELGBAO/BAOsgcfullsig.pdf')
	#plt.plot((dx[2]-dx[1])/2.,dw[2],'o',color=cp,markeredgecolor='k')
	plt.plot(d[1],dp[1],'o',color=cp,markeredgecolor='k')
	#plt.plot(dezxi[0],dezpk[0],'o',color=EZcol,markeredgecolor='k')
	#plt.plot(dq[1],dq[2],'^',color=QPMcol,markeredgecolor='k')
	xl = [0,2]
	yl = [0,2]
	plt.plot(xl,yl,'--k')
	xl = [0.021]
	yl = [0.024]
	plt.plot(xl,yl,'*',color=cd,markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	#xl = [0.036]
	#yl = [0.037]
	#plt.plot(xl,yl,'*',color='y',markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	#xl = [0.050]
	#yl = [0.047]
	#plt.plot(xl,yl,'*',color='teal',markersize=20,markeredgecolor='k',markeredgewidth=1.5)
	
	plt.xlim( 0.0, 0.1 )
	plt.ylim(0.0,0.1)
	plt.xlabel (r' $\sigma(\alpha)$ SGC only', fontsize=18)
	plt.ylabel (r' $\sigma(\alpha)$ SGC + NGC', fontsize=18)
	pp.savefig()
	pp.close()
	return True

def plotmockGauss(mm=1.,sfac=1.):
	import matplotlib.mlab as mlab
	from matplotlib.backends.backend_pdf import PdfPages
	dez = np.loadtxt('/Users/ashleyross/Dropbox/eboss/BAOfitsELGEZNScombf4rec0.59304.07.015.01.08st0.dat').transpose()
	plt.clf()
	plt.minorticks_on()
	pp = PdfPages('/Users/ashleyross/Dropbox/ELGBAO/BAOcompgauss.pdf')
	al = []
	sl = []
	n = 0
	for i in range(0,len(dez[0])):
		if dez[1][i] < .1:
			n += 1.
			al.append(dez[0][i])
			sl.append(dez[1][i])
	print( n)		
	ma = sum(al)/n
	dl = (np.array(al)-mm)/(np.array(sl)*sfac)
	#print max(dl),min(dl),n,ma,sum(sl)/n
	plt.hist(dl,25,normed=1,color='steelblue',histtype='stepfilled')#,histtype='step')
	plt.hist(dl,25,normed=1,color='k',histtype='step')
	x = np.linspace(-4,4,100)
	y = mlab.normpdf( x, 0, 1)
	from scipy.stats import kstest
	print(( kstest(dl,'norm')))
	plt.plot(x,y,'k--')
	plt.xlabel (r'$(\alpha-\langle \alpha\rangle)/\sigma(\alpha)$', fontsize=18)
	plt.ylabel ('relative number of mocks', fontsize=18)
	plt.ylim(0,.55)
	plt.text(-4,.5,r'K-S $D_n$: 0.022 (p-value 0.71)')
	plt.gca().xaxis.set_label_coords(.5,-.075)
	pp.savefig()
	pp.close()
	return True

def plotLRGBAO_Pk0():
	d = np.loadtxt(ebossdir+'LRGpk0BAODR16.txt').transpose()
	kl = d[0]
	dat = (d[1]-d[4])/d[4]
	err = d[2]/d[4]
	mod = (d[3]-d[4])/d[4]
	plt.errorbar(kl,dat,err,fmt='o',mfc='w',ms=8,elinewidth=1.75,capsize=4,mew=1.5,color='k')
	plt.plot(kl,mod,'k--')
	plt.xlabel(r'$k$ ($h$Mpc$^{-1}$)')
	plt.ylabel(r'$(P_0 - P_{\rm 0, smooth})/P_{\rm 0,smooth})$',labelpad=0)
	plt.savefig(ebossdir+'LRGPkBAO_mono.png')
	
	plt.show()

def plotLRGBAO_Pk2():
	d = np.loadtxt(ebossdir+'LRGpk2BAODR16.txt').transpose()
	d0 = np.loadtxt(ebossdir+'LRGpk0BAODR16.txt').transpose()
	kl = d[0]
	dat = (d[1]-d[4])/d0[4]
	err = d[2]/d0[4]
	mod = (d[3]-d[4])/d0[4]
	plt.errorbar(kl,dat,err,fmt='o',mfc='w',ms=8,elinewidth=1.75,capsize=4,mew=1.5,color='k')
	plt.plot(kl,mod,'k--')
	plt.xlabel(r'$k$ ($h$Mpc$^{-1}$)')
	plt.ylabel(r'$(P_2 - P_{\rm 2, smooth})/P_{\rm 0,smooth})$',labelpad=0)
	plt.savefig(ebossdir+'LRGPkBAO_quad.png')
	plt.show()
	
def plotLRGBAO_Pkring():
	d2 = np.loadtxt(ebossdir+'LRGpk2BAODR16.txt').transpose()
	d0 = np.loadtxt(ebossdir+'LRGpk0BAODR16.txt').transpose()
	kl = d0[0]
	dat0 = (d0[1]-d0[4])/d0[4]
	dat2 = (d2[1]-d2[4])/d0[4]
	datm = np.zeros((2*(len(kl)+2),2*(len(kl)+2)))
	for i in range(0,2*(len(kl)+2)):
		#if i < len(kl):
		#	kp = kl[-1-i]
		#else:
		#	kp = kl[i-len(kl)]
		kp = -.3+.005+.01*i
		for j in range(0,2*(len(kl)+2)):
			#if j < len(kl):
			#	kr = kl[-1-j]
			#else:	
			#	kr = kl[j-len(kl)]
			kr = -.3+.005+.01*j
			kt = sqrt(kp**2.+kr**2.)
			if kt < 0.3 and kt > 0.02:
				mu = kr/kt
				kind = int((kt-0.02)/0.01)
				pt = dat0[kind]+dat2[kind]*P2(mu)
				datm[i][j] = pt
	#err = d[2]/d0[4]
	#mod = (d[3]-d[4])/d0[4]
	#plt.errorbar(kl,dat,err,fmt='ko')
	#plt.plot(kl,mod,'k--')
	plt.imshow(datm,origin='lower',interpolation='bicubic',cmap='gray')
	locs,labels = plt.xticks()
	locs = np.array([4.5 ,14.5,  29.5,  44.5,  54.5  ])
	labs = []
	for loc in locs:
		kval = abs(-.3+.005+.01*loc)
		labs.append(str(kval.round(3)))
		
	plt.xticks(locs,labs)
	plt.yticks(locs,labs)
	print(locs)
	plt.xlim(0,len(datm)-1)
	plt.ylim(0,len(datm)-1)	
	plt.xlabel(r'$k_{\perp}$ ($h$Mpc$^{-1}$)')
	plt.ylabel(r'$k_{||}$ ($h$Mpc$^{-1}$)')
	plt.savefig(ebossdir+'LRGpkBAOring.png')
	plt.show()
	


if __name__ == '__main__':
	#do bao fit on EZ mocks
	import sys
	#fl = ''
	#ind = int(sys.argv[1])
	#LRGpCMASS post-rec 2D
	xi2DBAONS(md='dataJB',bs=5,damp='0.44.02.54.015.01.0',rec='rec',ver='7',covm = 'JB',min=50,max=150,Bp=0.4,af='angfac0')#,spat=0.001,spar=0.001,mina=.95,maxa=1.05)
	#pre-rec
	#xi2DBAONS(md='data',bs=5,damp='0.406.010.015.00',rec='',ver='5_1')
	
	#quasars

	#quasar mockave on NERSC
	#xibaoNS('QSO',0.8,2.2,md='qsomockave',version='5',wm='',bs=5,mom='024',covmd='QSO',tempmd='iso',rmin=50,rmax=150)
	#xibaoNS('QSO',0.8,2.2,md='qsomockave',version='5',wm='',bs=5,mom='024',covmd='me',tempmd='',damp='0.42.04915.00',rmin=50,rmax=150)
	#covf = 1.
	#test result for mockave on my computer and assuming zel'dovic
	#xibaoNS('QSO',0.8,2.2,md='qsozel',version='5',wm='',bs=5,mom='024',covmd='me',tempmd='',damp='0.42.03615.00',rmin=50,rmax=150)

	#xibaoNS('QSO',0.8,2.2,md='EZmockave',version='7',wm='',bs=5,mom='024',covmd='me',tempmd='',damp='0.44.03.08.015.00',rmin=50,rmax=150)
	#xibaoNS('QSO',0.8,2.2,md='EZmockave',version='7',wm='',bs=5,mom='024',covmd='QSO',tempmd='iso',damp='6.0',rmin=50,rmax=150)

	#print(covf)
	#
	#xi2DBAONandS(md='mockave',covf=covf,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.42.04915.00',rec='',ezver=5,min=50,max=150,Bp=0.4,covm='me',spat=0.002,spar=0.004)
	#print(covf)
	
	#data
	#xibaoNS('QSO',0.8,2.2,version='7_1',wm='fkp',bs=5,md='dataJH',covmd='me',tempmd='',damp='0.44.03.08.015.00',rmin=50,rmax=150)
	#xibaoNS('QSO',0.8,2.2,version='7_2',wm='fkp',bs=5,md='dataJH',covmd='me',tempmd='',damp='0.44.03.08.015.00',rmin=50,rmax=150)
	#xibaoNS('QSO',0.8,2.2,version='7_2',wm='fkp',bs=5,md='dataJH',covmd='me',tempmd='',damp='0.44.03.08.015.00',rmin=50,rmax=150,mb='nobao')
	#2D LRG
	#xi2DBAONS(md='data',samp='LRGpCMASS',min=50.,max=150.,maxb=80.,bs=8,binc=0,minz=0.6,maxz=1.,ver='4',af='',wm='fkp',rec='',damp='',spat=0.003,spar=0.006,mina=.8,maxa=1.2,covm='JB',Bp=0.4)
	#2D quasars, with varying damping terms
	#fiducial from Richard, Sigma_|| = 8, Sigma_perp = 3, Sigma_s = 4
#	print('Sigma_|| = 8, Sigma_perp = 3, Sigma_s = 4')
#	xi2DBAONS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.44.03.08.015.00',rec='',ver='7_2',ezver='7',min=50,max=150,Bp=0.4,covm='me')
# 	#xi2DBAONandS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.42.04915.00',rec='',ver='7',ezver='7',min=50,max=150,Bp=0.4,covm='me')
# 	print('Sigma_|| = 11, Sigma_perp = 3, Sigma_s = 4')
# 	xi2DBAONS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.44.03.011.015.00',rec='',ver='7',ezver='7',min=50,max=150,Bp=0.4,covm='me')
# 	print('Sigma_|| = 5, Sigma_perp = 3, Sigma_s = 4')
# 	xi2DBAONS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.44.03.05.015.00',rec='',ver='7',ezver='7',min=50,max=150,Bp=0.4,covm='me')
# 	print('Sigma_|| = 8, Sigma_perp = 2, Sigma_s = 4')
# 	xi2DBAONS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.44.02.08.015.00',rec='',ver='7',ezver='7',min=50,max=150,Bp=0.4,covm='me')
# 	print('Sigma_|| = 8, Sigma_perp = 4, Sigma_s = 4')
# 	xi2DBAONS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.44.04.08.015.00',rec='',ver='7',ezver='7',min=50,max=150,Bp=0.4,covm='me')
# 	print('Sigma_|| = 8, Sigma_perp = 3, Sigma_s = 2')
# 	xi2DBAONS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.42.03.08.015.00',rec='',ver='7',ezver='7',min=50,max=150,Bp=0.4,covm='me')
# 	print('Sigma_|| = 8, Sigma_perp = 3, Sigma_s = 6')
# 	xi2DBAONS(md='dataJH',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.46.03.08.015.00',rec='',ver='7',ezver='7',min=50,max=150,Bp=0.4,covm='me')

	#quasar mocks on NERSC
	#print('Sigma_|| = 8, Sigma_perp = 3, Sigma_s = 4')
	#xi2DBAONS(md='mockave',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.44.03.08.015.00',rec='',ver='7',ezver='7',min=50,max=150,maxr=170,Bp=0.1,Bt=0.4,covm='me')
	
	#xi2DBAONS(md='mock',covf=1.,samp='QSO',minz=0.8,maxz=2.2,bs=5,damp='0.44.04.08.015.00',rec='',mockn=int(sys.argv[1]),ezver=5,min=50,max=150,Bp=0.4,covm='me')
	#xibaoNS('QSO',0.8,2.2,md='qsomock',version='5',wm='',bs=5,mom='024',tempmd='',mockn=int(sys.argv[1]),covmd='me',damp='0.42.04915.00',rmin=50,rmax=150)
	#make all xi files pre and post reconstruction
	#calcwrp_mockELGEZ(ind,reg='SGC',rec='')
	#calcwrp_mockELGEZ(ind,reg='NGC',rec='')
	#calcwrp_mockELGEZ(ind,reg='SGC',rec='rec')
	#calcwrp_mockELGEZ(ind,reg='NGC',rec='rec')
	#calcxi_mockEZ(ind,reg='SGC',samp='LRG',pcmass=True,bs=5,rec='rec')
	#calcxi_mockEZ(ind,reg='NGC',samp='LRG',pcmass=True,bs=5,rec='rec')
#	calcxi_mockELGEZ(ind,reg='SGC',bs=8,start=0,rec='')
#	calcxi_mockELGEZ(ind,reg='NGC',bs=8,start=0,rec='')
#	calcxi_mockELGEZ(ind,reg='SGC',bs=8,start=0,rec='shuff')
#	calcxi_mockELGEZ(ind,reg='NGC',bs=8,start=0,rec='shuff')
 	
# 	calcxi_mockELGEZ(ind,reg='SGC',bs=10,start=0,rec='')
# 	calcxi_mockELGEZ(ind,reg='NGC',bs=10,start=0,rec='')
#calc mock xi for different start
#  	calcxi_mockELGEZ(ind,reg='SGC',bs=8,start=2,rec='rec')
#  	calcxi_mockELGEZ(ind,reg='NGC',bs=8,start=2,rec='rec')
#  	calcxi_mockELGEZ(ind,reg='SGC',bs=8,start=4,rec='rec')
#  	calcxi_mockELGEZ(ind,reg='NGC',bs=8,start=4,rec='rec')
#  	calcxi_mockELGEZ(ind,reg='SGC',bs=8,start=6,rec='rec')
#  	calcxi_mockELGEZ(ind,reg='NGC',bs=8,start=6,rec='rec')

# 	calcxi_mockELGEZ(ind,reg='SGC',bs=6,start=0,rec='rec')
# 	calcxi_mockELGEZ(ind,reg='NGC',bs=6,start=0,rec='rec')

#ELG data
#prerec
#	print('PRE RECON data, rmin = 50')
#	xibaoNS('ELG',.6,1.1,'7',rec='',covv='',md='data',covmd='me',damp='0.5933.06.010.015.00',bs=5,rmin=50,rmax=150,rmaxb=58)
#	print('PRE RECON data, rmin = 30')
#	xibaoNS('ELG',.6,1.1,'7',rec='',covv='',md='data',covmd='me',damp='0.5933.06.010.015.00',bs=5,rmin=30,rmax=150,rmaxb=50)
#	xibaoNS('ELG',.6,1.1,'5_2_badphotAD',rec='',covv='',md='data',covmd='ELG',damp='0.5933.06.010.015.00',bs=8,rmin=30,rmax=150,rmaxb=50,mb='nobao')
#postrec
# 	xibaoNS('ELG',.6,1.1,'4',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('4 done')
# 	xibaoNS('ELG',.6,1.1,'4_a',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('4_a done')
# 	xibaoNS('ELG',.6,1.1,'5_a',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('5_a done')
# 	xibaoNS('ELG',.6,1.1,'5_b',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('5_b done')
# 	xibaoNS('ELG',.6,1.1,'5_c',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('5_c done')
# 	xibaoNS('ELG',.6,1.1,'5_d',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('5_d done')
# 	xibaoNS('ELG',.6,1.1,'5_e',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('5_e done')
#	xibaoNS('ELG',.6,1.1,'5_1',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
#	print('5_1 done')
#	xibaoNS('ELG',.6,1.1,'5_1AD',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
#	print('5_1AD done')
#	xibaoNS('ELG',.6,1.1,'5_1_wAD',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
#	print('5_1_wAD done')
#	xibaoNS('ELG',.6,1.1,'5_2AD',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
#	print('5_2AD done')
# 	print('POST RECON data, rmin = 50, bs = 8')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58)

# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin 0.6,damp 47')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=55,start=0)
# 	print('POST RECON data, rmin = 50, bs = 5, start=0 zmin 0.7, damp 47')
# 	xibaoNS('ELG',.7,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=55,start=0)

# 	print('POST RECON data, rmin = 50, bs = 5, start=0,mup=1')
# 	xibaoNS('ELG',.6,1.1,'7',mupow=1,rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=55,start=0)
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.7, damp = 35')
# 	xibaoNS('ELG',.7,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0)
#	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35')
#	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0)
#	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35, no bias prior')
#	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0,Bp=100)
#	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35, no An')
#	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0,npar=0)
#	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35, 4th order')
#	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0,npar=3)

# 	# 	print('POST RECON data, rmin = 50, bs = 5, start=1, zmin=0.6, damp = 35')
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35,start 1')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=1)
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35,start 2')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=2)
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35, start 3')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=3)
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35, start 4')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=4)
# 
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 47, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=1, zmin=0.6, damp = 35, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=1,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=1, zmin=0.6, damp = 47, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=1,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=2, zmin=0.6, damp = 35, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=2,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=2, zmin=0.6, damp = 47, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=2,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=3, zmin=0.6, damp = 35, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=3,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=3, zmin=0.6, damp = 47, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=3,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 35, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=4,mb='nobao')
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 47, nobao')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=4,mb='nobao')
# 
# 
# 	print('POST RECON data, rmin = 50, bs = 5, start=0, zmin=0.6, damp = 47')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=0)
# 	print('POST RECON data, rmin = 50, bs = 5, start=1')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=1)
# 	print('POST RECON data, rmin = 50, bs = 5, start=2')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=2)
# 	print('POST RECON data, rmin = 50, bs = 5, start=3')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=3)
# 	print('POST RECON data, rmin = 50, bs = 5, start=4')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=58,start=4)

# 	print('POST RECON data, rmin = 30')
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',md='data',covmd='me',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)

#	xibaoNS('ELG',.6,1.1,'5_2_badphotAD',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50,mb='nobao')
#	print('5_2badphotAD done')

# 	xibaoNS('ELG',.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('5 done')
# 	xibaoNS('ELG',.6,1.1,'5_1',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=30,rmax=150,rmaxb=50)
# 	print('5_1 done')
# 	xibaoNS('ELG',.6,1.1,'5_1',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,start=2,rmin=30,rmax=150,rmaxb=50)
# 	print('5_1 2 done')
# 	xibaoNS('ELG',.6,1.1,'5_1',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,start=4,rmin=30,rmax=150,rmaxb=50)
# 	print('5_1 4 done')
# 	xibaoNS('ELG',.6,1.1,'5_1',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,start=6,rmin=30,rmax=150,rmaxb=50)
# 	print('5_1 6 done')

	#xibaoNS('ELG',.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
#as function of mu
# 	print('mu = 0.05')
# 	xibaoNSmu('ELG',0.05,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.15')
# 	xibaoNSmu('ELG',0.15,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.25')
# 	xibaoNSmu('ELG',0.25,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.35')
# 	xibaoNSmu('ELG',0.35,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.45')
# 	xibaoNSmu('ELG',0.45,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.55')
# 	xibaoNSmu('ELG',0.55,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.65')
# 	xibaoNSmu('ELG',0.65,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.75')
# 	xibaoNSmu('ELG',0.75,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.85')
# 	xibaoNSmu('ELG',0.85,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')
# 	print('mu = 0.95')
# 	xibaoNSmu('ELG',0.95,.6,1.1,'5',rec='rec',covv='',md='data',covmd='ELG',damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58,mb='')

#ELG mocks
# for i in range(int(sys.argv[1]),1001):
# 	calcxi_mockEZ(i,ver='7',rec='_rec')
# 	calcxi_mockEZ(i,ver='7',rec='')
# 	calcxi_mockEZ(i,ver='7',rec='',reg='NGC')
# 	calcxi_mockEZ(i,ver='7',rec='_rec',reg='NGC')
# 	print(i)

# mockave rec as function of mu:
# 	bs = 5
# 	print('mock ave POST RECON, rmin =50, mu 0,0.2, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mumax=0.2,mumin=0,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=.1,Bp=100)
# 	print('mock ave POST RECON, rmin =50, mu 0.2,0.4, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mumax=0.4,mumin=0.2,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=.1,Bp=100)
# 	print('mock ave POST RECON, rmin =50, mu 0.4,0.6, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mumax=0.6,mumin=0.4,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=.1,Bp=100)
# 	print('mock ave POST RECON, rmin =50, mu 0.6,0.8, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mumax=0.8,mumin=0.6,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=.1,Bp=100)
# 	print('mock ave POST RECON, rmin =50, mu 0.8,1, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mumax=1,mumin=0.8,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=.1,Bp=100)
#mockave mupow 1
# 	bs = 5
# 	print('mock ave POST RECON, rmin =50, mupow = 1, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mupow=1,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=0.1)
# 	bs = 5
# 	print('mock ave POST RECON, rmin =50, mupow = 0, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mupow=0,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=0.1)

# mockave rec and no rec:
# 	bs = 5
# 	print('mockave POST RECON, rmin =50, zmin 0.6, damp 58 bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59305.08.515.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)

# 	print('mockave POST RECON, rmin =50, zmin 0.6, damp 47 bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)
# 	print('mockave POST RECON, rmin =50, zmin 0.6, damp 35 bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59303.05.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)

# 	bs = 5
# 	print('mockave POST RECON, rmin =50, zmin 0.6, damp 58 bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',ns='_nosys',rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59305.08.515.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)
# 
# 	print('mockave POST RECON, rmin =50, zmin 0.6, damp 47 bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',ns='_nosys',rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)
# 	print('mockave POST RECON, rmin =50, zmin 0.6, damp 35 bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',ns='_nosys',rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59303.05.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)


# 	print('mockave POST RECON, rmin =50, zmin 0.7, damp 47 bin size'+str(bs))
# 	xibaoNS('ELG',.7,1.1,'7',rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)

# 	print('mock ave POST RECON, rmin =50, mumax=0.7, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',mumax=0.7,rec='_rec',covv='',covmd='me',md = 'EZmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=.1)
# 	print('PRE RECON, rmin =50, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',rec='',covv='',covmd='me',md = 'EZmockave',damp='0.5933.06.010.015.00',bs=bs,rmin=50,rmax=150,rmaxb=58)
# 	print('PRE RECON, rmin =30, bin size'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',rec='',covv='',covmd='me',md = 'EZmockave',damp='0.5933.06.010.015.00',bs=bs,rmin=30,rmax=150,rmaxb=50)

	bs = 5
# 	print('OR mockave POST RECON, rmin =50, bin size damp 35'+str(bs))
# 	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',covmd='me',md = 'ORmockave',damp='0.59303.05.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)
#	print('OR mockave POST RECON, rmin =50, bin size damp 47'+str(bs))
#	xibaoNS('ELG',.6,1.1,'7',rec='_rec',covv='',covmd='me',md = 'ORmockave',damp='0.59304.07.015.01.0',bs=bs,rmin=50,rmax=150,rmaxb=55,covfac=1.)


#pre-recon BAO fit	
	#xibaoNS('ELGEZ',.6,1.1,'4',rec='',covv='',covmd='ELG',mockn=str(ind),damp='0.5933.06.010.015.00',bs=5,rmin=50,rmax=150,rmaxb=55)
# 	sind = int(sys.argv[1])
# 	max = int(sys.argv[2])
# 	for ind in range(sind,max):
# 		xibaoNS('ELG',.6,1.1,'7',rec='',md='EZmock',covv='',covmd='me',mockn=str(ind),damp='0.5933.06.010.015.00',bs=8,rmin=50,rmax=150,rmaxb=58)
# 		xibaoNS('ELG',.6,1.1,'7',rec='_rec',md='EZmock',covv='',covmd='me',mockn=str(ind),damp='0.59304.07.015.01.0',bs=8,start=0,rmin=49.9,rmax=150,rmaxb=58)
# 		print('MOCK '+str(ind)+ 'DONE')
#different starts
#	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59304.07.015.01.0',bs=8,start=2,rmin=49.9,rmax=150,rmaxb=58)
#	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59304.07.015.01.0',bs=8,start=4,rmin=49.9,rmax=150,rmaxb=58)
#	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59304.07.015.01.0',bs=8,start=6,rmin=49.9,rmax=150,rmaxb=58)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59304.07.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59303.05.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=55)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59303.05.015.01.0',bs=10,rmin=50,rmax=150,rmaxb=60)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59303.05.015.01.0',bs=6,rmin=50,rmax=150,rmaxb=56)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.59303.05.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.593.02.04.015.01.0',bs=5,rmin=50,rmax=150,rmaxb=55)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.593.02.04.015.01.0',bs=10,rmin=50,rmax=150,rmaxb=60)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.593.02.04.015.01.0',bs=6,rmin=50,rmax=150,rmaxb=56)
# 	xibaoNS('ELGEZ',.6,1.1,'4',rec='rec',covv='',covmd='ELG',mockn=str(ind),damp='0.593.02.04.015.01.0',bs=8,rmin=50,rmax=150,rmaxb=58)
	
