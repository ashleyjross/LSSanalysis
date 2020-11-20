import healpy as hp
from healpix import radec2thphi,thphi2radec
import fitsio
import numpy as np
from math import *
from matplotlib import pyplot as plt
from optimize import fmin
from astropy.table import Table

bossdir = '/Users/ashleyross/Dropbox/BOSS/'
app = '.fits.gz'
ic = 3

def luptm(nmag,bnd):
	#convert fluxes to magnitudes
	b = []
	b.append(1.4*10**-10.)
	b.append(.9*10**-10.)
	b.append(1.2*10**-10.)
	b.append(1.8*10**-10.)
	b.append(7.4*10**-10.)
	return -2.5/log(10.)*(asinh((nmag/10.**9.)/(2.*b[bnd]))+log(b[bnd]))


def mkpixl_fits(file,res,min=0,max=1,c=999,wm='cp',zmin=.35,zmax=1.,md='',app='fits.gz'):
    npix = 12*res*res
    f = fitsio.read(bossdir+file+app) #written for mksample format, is the double fits something I did?
    #if wm == 'stm' or wm == 'stmsee' or wm == 'stmseem':
     #  fstw = open(file+'-stweights.dat')
    if wm == 'seem' or wm == 'stseem' or wm == 'stmseem':
        fseew = open(file+'-seeweights.dat')    
    pixl = np.zeros(npix)
    n = 0
    g = 0
    #assumes file is already masked properly and has equatorial ra,dec coordinates in degrees
    for i in range (0,len(f)):
        k = 0
        #if wm == 'stm' or wm == 'stmsee' or wm == 'stmseem':
        #   wstm = float(fstw.readline())
        if wm == 'seem' or wm == 'stseem' or wm == 'stmseem':
            wsee = float(fseew.readline())  
        if c != 999:
            if c == 'ifib':
                t = f[i]['FIBER2FLUX'][ic] 
                fibi = luptm(t,ic)-f[i]['EXTINCTION'][ic] 
                if fibi < min or fibi > max:
                    k = 1
        if md == 'lze':
            zpsf = luptm(f[i]['PSFFLUX'][4],4)-f[i]['EXTINCTION'][4]
            ipsf = luptm(f[i]['PSFFLUX'][3],3)-f[i]['EXTINCTION'][3]
            zmod = luptm(f[i]['MODELFLUX'][4],4)-f[i]['EXTINCTION'][4]
            imod = luptm(f[i]['MODELFLUX'][3],3)-f[i]['EXTINCTION'][3]
            rmod = luptm(f[i]['MODELFLUX'][2],2)-f[i]['EXTINCTION'][2]
            fracpsfr = f[i]['FRACPSF'][2]
            frdev = f[i]['DEVFLUX'][2]
            frexp = f[i]['EXPFLUX'][2]
            cmodfluxr = fracpsfr*frdev+frexp*(1.-fracpsfr)
            rcmod = luptm(cmodfluxr,2)-f[i]['EXTINCTION'][2]
            gmod = luptm(f[i]['MODELFLUX'][1],1)-f[i]['EXTINCTION'][1]
            cp = 0.7*(gmod-rmod)+1.2*(rmod-imod-.18)
            if zpsf - zmod < 9.125 -0.46*zmod or ipsf-imod < 0.2+.02*(20.0-imod):
                k = 1
            if rcmod > 19.5 or rcmod > 13.4 + cp/0.3:
                k = 1   

        z = f[i]['Z']
        if z < zmin or z > zmax:
            k = 1       
        if k == 0:
            ra,dec = f[i]['RA'],f[i]['DEC']
            w = 1. #if there are weights to be used in the file, put them here
            if wm == 'cp':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_FKP']
            if wm == 'st':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_STAR']*f[i]['WEIGHT_FKP']
            if wm == 'stg':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_STARGAIA']*f[i]['WEIGHT_FKP']
            if wm == 'stm':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*wstm*f[i]['WEIGHT_FKP']
            if wm == 'stmsee':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*wstm*f[i]['WEIGHT_SEEING']*f[i]['WEIGHT_FKP']
            if wm == 'see':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']*f[i]['WEIGHT_FKP']
            if wm == 'stsee':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*f[i]['WEIGHT_SEEING']*f[i]['WEIGHT_STAR']*f[i]['WEIGHT_FKP']
            if wm == 'seem':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*wsee*f[i]['WEIGHT_FKP']           
            if wm == 'stmseem':
                w = (f[i]['WEIGHT_NOZ']+f[i]['WEIGHT_CP']-1.)*wsee*wstm*f[i]['WEIGHT_FKP']      

            th,phi = radec2thphi(ra,dec)
            p = hp.ang2pix(res,th,phi)
            pixl[p] += 1.*w
            g += 1.
        else:
            n += 1.
    print(g,n)          
    return pixl

def mkranpixfile(reg,res=256):
    npo = 12*res**2
    ml = np.zeros(npo)
    ranf = fitsio.read(bossdir+'random0_DR12v5_CMASS_'+reg+'.fits.gz')
    th,phi = radec2thphi(ranf['RA'],ranf['DEC'])
    pixr = hp.ang2pix(res,th,phi)
    for pix in pixr:
        ml[pix] += 1
    fo = bossdir+'nran_'+reg+str(res)+'.dat'
    np.savetxt(fo, ml)
    return True



def ngvstar_ifib(file,ifibm,ifibx,reg='South',map='sdss',sysmin=50,sysmax=350,res=256,wm='cp',zmin=.35,zmax=1.,md='',app='.fits.gz'):
    #file is string for mksample fits file without 
    #hardcoded columns for relevant quantities, probably a better way to write this...
    stl = []
    wstl = []
    errl = []
    if map == 'sdss':
        fsys = open('maps/allstars17.519.9Healpixall'+str(res)+'.dat').readlines()
    if map == 'gaia':
        fsys = fitsio.read('maps/Gaia.dr2.bGT10.12g17.hp256.fits')['hpstardens']

    pixl = mkpixl_fits(file,res,min=ifibm,max=ifibx,c='ifib',wm=wm,zmin=zmin,zmax=zmax,md=md,app=app)
    
    fr = bossdir+'nran_'+reg+str(res)+'.dat'
    ml = np.loadtxt(fr)
    npo = 12*res**2

    
    

    
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
            nt += ml[i]
            nbt += pixl[i]
            bins = int((sysv-sysmin)*sysm)
            if bins >= 0 and bins < nsysbin:
                binnbs[bins] += pixl[i]
                binns[bins] += ml[i]
            else:
                bs += pixl[i] #count numbers outside of sysmin/sysmax
                bsr +=  ml[i]
        else:
            n0 += ml[i] #count numbers inside bad pixels in sys map
            ng0 += pixl[i]
                    
    print ('total number of randoms(or pixels)/objects '+str(nt)+'/'+str(nbt))
    print ('number of randoms(or pixels)/objects where sys = 0 '+str(n0)+'/'+str(ng0))
    print ('number of randoms(or pixels)/objects outside tested range '+str(bsr)+'/'+str(bs)  )       
    ave = nbt/nt
    print ('average number of objects per random (or pixel) is '+ str(ave))
    zw = ''
    if zmin != 0.35:
        zw += 'mz'+str(zmin)
    if zmax != 1.:
        zw += 'xz'+str(zmax)    
    fs = open(bossdir+'n'+file+str(res)+zw+'ifib'+str(ifibm)+str(ifibx)+str(wm)+'vnst'+map+'.dat','w')
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
    plt.errorbar(xl,yl,el,fmt='ko')
    plt.title(reg +' '+map)
    plt.show()
    return xl,yl,el

def ifiblinfits_gaia():
	#combines NGC+SGC for five ifib bins used in previous papers
	
	rn = ngvstar_ifib('galaxy_DR12v5_CMASS_North',20,20.3,reg='North',map='gaia',sysmin=150,sysmax=2000)
	rs = ngvstar_ifib('galaxy_DR12v5_CMASS_South',20,20.3,reg='South',map='gaia',sysmin=150,sysmax=2000)
	y0t = []
	e0t = []
	for i in range(0,len(rn[1])):
		yt = (rn[1][i]/rn[2][i]**2.+rs[1][i]/rs[2][i]**2.)/(1./rn[2][i]**2.+1./rs[2][i]**2.)
		et = sqrt(1./(1./rn[2][i]**2.+1./rs[2][i]**2.))
		y0t.append(yt)
		e0t.append(et)
	lf = linfit(rn[0],y0t,e0t)
	inl = np.array([1.,0])
	b0,m0 = fmin(lf.chilin,inl)
	print (y0t[4],rn[1][4],e0t[4],rn[2][4]) #check that answers make sense
	print (b0,m0)
	#return True
	rn = ngvstar_ifib('galaxy_DR12v5_CMASS_North',20.3,20.6,reg='North',map='gaia',sysmin=150,sysmax=2000)
	rs = ngvstar_ifib('galaxy_DR12v5_CMASS_South',20.3,20.6,reg='South',map='gaia',sysmin=150,sysmax=2000)
	y1t = []
	e1t = []
	for i in range(0,len(rn[1])):
		yt = (rn[1][i]/rn[2][i]**2.+rs[1][i]/rs[2][i]**2.)/(1./rn[2][i]**2.+1./rs[2][i]**2.)
		et = sqrt(1./(1./rn[2][i]**2.+1./rs[2][i]**2.))
		y1t.append(yt)
		e1t.append(et)
	lf = linfit(rn[0],y1t,e1t)
	inl = np.array([1.,0])
	b1,m1 = fmin(lf.chilin,inl)
	print (b1,m1)
	rn = ngvstar_ifib('galaxy_DR12v5_CMASS_North',20.6,20.9,reg='North',map='gaia',sysmin=150,sysmax=2000)
	rs = ngvstar_ifib('galaxy_DR12v5_CMASS_South',20.6,20.9,reg='South',map='gaia',sysmin=150,sysmax=2000)
	y2t = []
	e2t = []
	for i in range(0,len(rn[1])):
		yt = (rn[1][i]/rn[2][i]**2.+rs[1][i]/rs[2][i]**2.)/(1./rn[2][i]**2.+1./rs[2][i]**2.)
		et = sqrt(1./(1./rn[2][i]**2.+1./rs[2][i]**2.))
		y2t.append(yt)
		e2t.append(et)
	lf = linfit(rn[0],y2t,e2t)
	inl = np.array([1.,0])
	b2,m2 = fmin(lf.chilin,inl)
	print (b2,m2)
	rn = ngvstar_ifib('galaxy_DR12v5_CMASS_North',20.9,21.2,reg='North',map='gaia',sysmin=150,sysmax=2000)
	rs = ngvstar_ifib('galaxy_DR12v5_CMASS_South',20.9,21.2,reg='South',map='gaia',sysmin=150,sysmax=2000)
	y3t = []
	e3t = []
	for i in range(0,len(rn[1])):
		yt = (rn[1][i]/rn[2][i]**2.+rs[1][i]/rs[2][i]**2.)/(1./rn[2][i]**2.+1./rs[2][i]**2.)
		et = sqrt(1./(1./rn[2][i]**2.+1./rs[2][i]**2.))
		y3t.append(yt)
		e3t.append(et)
	lf = linfit(rn[0],y3t,e3t)
	inl = np.array([1.,0])
	b3,m3 = fmin(lf.chilin,inl)
	print (b3,m3)
	rn = ngvstar_ifib('galaxy_DR12v5_CMASS_North',21.2,30,reg='North',map='gaia',sysmin=150,sysmax=2000)
	rs = ngvstar_ifib('galaxy_DR12v5_CMASS_South',21.2,30,reg='South',map='gaia',sysmin=150,sysmax=2000)
	y4t = []
	e4t = []
	for i in range(0,len(rn[1])):
		yt = (rn[1][i]/rn[2][i]**2.+rs[1][i]/rs[2][i]**2.)/(1./rn[2][i]**2.+1./rs[2][i]**2.)
		et = sqrt(1./(1./rn[2][i]**2.+1./rs[2][i]**2.))
		y4t.append(yt)
		e4t.append(et)
	lf = linfit(rn[0],y4t,e4t)
	inl = np.array([1.,0])
	b4,m4 = fmin(lf.chilin,inl)
	print (b4,m4)
	fo = open(bossdir+'nstlinfits256gaia.dat','w')
	fo.write(str(b0)+' '+str(m0)+'\n')
	fo.write(str(b1)+' '+str(m1)+'\n')
	fo.write(str(b2)+' '+str(m2)+'\n')
	fo.write(str(b3)+' '+str(m3)+'\n')
	fo.write(str(b4)+' '+str(m4)+'\n')
	fo.close()
	return True

def assignstweights_gaia(NS='North',rc=0,dc=1,res=256):
	#npix = 12*res*res
	f = Table.read(bossdir+'galaxy_DR12v5_CMASS_'+NS+app)
	
	lfits = open(bossdir+'nstlinfits256gaia.dat').readlines()
	fsys = fitsio.read('maps/Gaia.dr2.bGT10.12g17.hp256.fits')['hpstardens']
	wts = np.ones(len(f))
	#fo = open(v+'-'+NS+'-stweights.dat','w')
	for i in range(0,len(f)):
		t = f[i]['FIBER2FLUX'][ic] 
		fibi = luptm(t,ic)-f[i]['EXTINCTION'][ic] 
		if fibi < 20.45:
			b0 = float(lfits[0].split()[0])
			m0 = float(lfits[0].split()[1])
			b1 = float(lfits[1].split()[0])
			m1 = float(lfits[1].split()[1])
			x0 = 20.15
			x1 = 20.45
		if fibi >= 20.45 and fibi < 20.75:
			b0 = float(lfits[1].split()[0])
			m0 = float(lfits[1].split()[1])
			b1 = float(lfits[2].split()[0])
			m1 = float(lfits[2].split()[1])
			x0 = 20.45
			x1 = 20.75
		if fibi >= 20.75 and fibi < 21.05:
			b0 = float(lfits[2].split()[0])
			m0 = float(lfits[2].split()[1])
			b1 = float(lfits[3].split()[0])
			m1 = float(lfits[3].split()[1])
			x0 = 20.75
			x1 = 21.05
		if fibi >= 21.05:
			b0 = float(lfits[3].split()[0])
			m0 = float(lfits[3].split()[1])
			b1 = float(lfits[4].split()[0])
			m1 = float(lfits[4].split()[1])
			x0 = 21.05
			x1 = 21.35
		mm = (m1-m0)/(x1-x0)
		bm = m0-mm*x0
		mb = (b1-b0)/(x1-x0)
		bb = b0-mb*x0
		mf = mm*fibi+bm
		bf = mb*fibi+bb		
		ra,dec = f[i][rc],f[i][dc]
		th,phi = radec2thphi(ra,dec)
		p = hp.ang2pix(res,th,phi)
		nst = float(fsys[p])
		wst = 1./(mf*nst+bf)
		wts[i] = wst
	f['WEIGHT_STARGAIA'] = wts
	f.write(bossdir+'galaxy_DR12v5_CMASS_gaia'+NS+'.fits',format='fits',overwrite=True)
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
