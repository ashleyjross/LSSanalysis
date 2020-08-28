import fitsio
import numpy as np
#from desitarget.io import read_targets_in_hp, read_targets_in_box, read_targets_in_cap
import astropy.io.fits as fits
import glob
import os
import healpy as hp
from matplotlib import pyplot as plt
import sys
#from math import *

bmzls = b'N' #if in desi environment

elgandlrgbits = [1,5,6,7,8,9,11,12,13]

R_G=3.214 # http://legacysurvey.org/dr8/catalogs/#galactic-extinction-coefficients
R_R=2.165
R_Z=1.211


def mask(fl):
	keep = (fl['NOBS_G']>0) & (fl['NOBS_R']>0) & (fl['NOBS_Z']>0)
	print(len(fl[keep]))
	for bit in elgandlrgbits:
		keep &= ((fl['MASKBITS'] & 2**bit)==0)
	print(len(fl[keep]))
	flk = fl[keep]
	return flk

#for healpix coordinates
def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.		
    
def splitcat(cat):
    NN = cat['PHOTSYS'] == bmzls
    d1 = (cat['PHOTSYS'] != bmzls) & (cat['RA'] < 300) & (cat['RA'] > 100) & (cat['DEC'] > -20)
    d2 = (d1==0) & (NN ==0) & (cat['DEC'] > -30)
    return cat[NN],cat[d1],cat[d2]

def splitcat_ind(cat):
    NN = cat['PHOTSYS'] == bmzls
    d1 = (cat['PHOTSYS'] != bmzls) & (cat['RA'] < 300) & (cat['RA'] > 100) & (cat['DEC'] > -20)
    d2 = (d1==0) & (NN ==0) & (cat['DEC'] > -30)
    return NN,d1,d2    

def splitcatS_ind(cat):
    #NN = cat['PHOTSYS'] == bmzls
    d1 = (cat['RA'] < 300) & (cat['RA'] > 100) & (cat['DEC'] > -20) & (cat['DEC'] < 32.375)
    d2 = (d1==0) & (cat['DEC'] > -30) & (cat['DEC'] < 32.375)
    return d1,d2

def plotvshp_compmc(r1,delg,dlrg,sys,rng,gdzm=20,ebvm=0.15,useMCeff=True,correctstar=True,title='',effac=1.,south=True,fout=''):
    w = hpq['GALDEPTH_Z'] > gdzm
    w &= hpq['EBV'] < ebvm
    if useMCeff:
        w &= mcl > 0
    if sys != 'gdc' and sys != 'rdc' and sys != 'zdc':
        sm = hpq[w][sys]
        xlab = sys
    else:
        if sys == 'gdc':
            print('g band depth, extinction corrected')
            sm = hpq[w]['GALDEPTH_G']*10.**(-0.4*R_G*hpq[w]['EBV'])
            xlab = 'g band depth, extinction corrected'
        if sys == 'rdc':
            sm = hpq[w]['GALDEPTH_R']*10.**(-0.4*R_R*hpq[w]['EBV'])
            xlab = 'r band depth, extinction corrected'
        if sys == 'zdc':
            sm = hpq[w]['GALDEPTH_Z']*10.**(-0.4*R_Z*hpq[w]['EBV'])
            xlab = 'z band depth, extinction corrected'
    ds = np.ones(len(delg))
    print(len(ds),len(delg),len(w),len(sm))
    hdnoc = np.histogram(sm,weights=dlrg[w],range=rng)
    #print(hd1)
    hr1 = np.histogram(sm,weights=r1[w],bins=hdnoc[1],range=rng)
    xl = []
    for i in range(0,len(hr1[0])):
        xl.append((hr1[1][i]+hr1[1][i+1])/2.)

    nlrg = hdnoc[0]/hr1[0]/(sum(dlrg[w])/sum(r1[w]))
    elrg = np.sqrt(hdnoc[0])/hr1[0]/(sum(dlrg[w])/sum(r1[w]))
    plt.errorbar(xl,nlrg,elrg,fmt='ko',label='LRG')
    fo = open('nlrg'+fout+'vsEBV_ebf'+ebf+'_Rv'+rv+'.dat','w')
    for i in range(0,len(xl)):
    	fo.write(str(xl[i])+' '+str(nlrg[i])+' '+str(elrg[i])+'\n')
    fo.close()	
    dmcse = mcl**effac
    hd1 = np.histogram(sm,weights=delg[w]*ds[w]/dmcse[w],bins=hdnoc[1],range=rng)
    nelg = hd1[0]/hr1[0]/(sum(delg[w]*ds[w]/dmcse[w])/sum(r1[w]))
    eelg = np.sqrt(hd1[0])/hr1[0]/(sum(delg[w]*ds[w]/dmcse[w])/sum(r1[w]))
    fo = open('nelgMC'+fout+'vsEBV_ebf'+ebf+'_Rv'+rv+'.dat','w')
    for i in range(0,len(xl)):
    	fo.write(str(xl[i])+' '+str(nelg[i])+' '+str(eelg[i])+'\n')
    fo.close()	

    plt.plot(xl,nelg,'r-',label='ELG+MC')
    plt.plot(xl,np.ones(len(xl)),'k:',label='null')
    plt.legend()#(['raw','with stellar density weights','+sed ext MC','just sed MC','old MC','null']))
    plt.ylabel('relative density')
    plt.xlabel(xlab)
    plt.ylim(0.7,1.3)
    plt.title(title)
    plt.show()    

def plotvsEBVvar():
	ebfl = ['0.9_Rv3.1','1.1_Rv3.1','1_Rv2.6','1_Rv3.6']
	#rvl = ['3.1']
	regl = ['DECaLS_SGC','DECaLS_NGC','BMzLS']
	typel = ['elgMC','lrg']
	for type in typel:
		for reg in regl:
			plt.clf()
			plt.title(type+" "+reg)
			plt.ylim(0.7,1.3)
			plt.ylabel('relative density')
			plt.xlabel('SFD E(B-V)')
			for ebfac in ebfl:
				#for r in rvl:
				delg = np.loadtxt('n'+type+reg+'vsEBV_ebf'+ebfac+'.dat').transpose()
				plt.errorbar(delg[0],delg[1],delg[2],label='EBVx'+ebfac)
			plt.plot(delg[0],np.ones(len(delg[0])),'k:')
			plt.legend()
			plt.savefig('n'+type+reg+'vsEBV.png')
			plt.show()		
	

if __name__ == '__main__':
	ebf = sys.argv[1]
	rv = sys.argv[2]
	es = 'ebvfac'+str(ebf)+'Rv'+str(rv)

	flrg = fitsio.read('/global/cscratch1/sd/ajross/combbricks/mysweeps/LRGdr8_south'+es+'.fits')
	lrgs = mask(flrg)
	lsdnl,lsdsl = splitcatS_ind(lrgs)
	flrgn = fitsio.read('/global/cscratch1/sd/ajross/combbricks/mysweeps/LRGdr8_north'+es+'.fits')
	lrgn = mask(flrgn)

	felg = fitsio.read('/global/cscratch1/sd/ajross/combbricks/mysweeps/ELGdr8_south'+es+'.fits')
	elgs = mask(felg)
	esdnl,esdsl = splitcatS_ind(elgs)
	felgn = fitsio.read('/global/cscratch1/sd/ajross/combbricks/mysweeps/ELGdr8_north'+es+'.fits')
	elgn = mask(felgn)

	#full random file is available, easy to read some limited number; take 1.5x ELG to start with
	lelg = len(felg)
	rall = fitsio.read('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits',rows=np.arange(int(1.5*lelg)))
	print('check this is correct:')
	print(np.unique(rall['PHOTSYS']),bmzls)

	rm = mask(rall)
	rbml,rdnl,rdsl = splitcat_ind(rm)


	#Some information is in pixelized map
	#get nside and nest from header
	pixfn      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/pixweight/pixweight-dr8-0.31.1.fits'
	hdr        = fits.getheader(pixfn,1)
	nside,nest = hdr['HPXNSIDE'],hdr['HPXNEST']
	hpq = fitsio.read(pixfn)

	def puthealpix(fl):
		dth,dphi = radec2thphi(fl['RA'],fl['DEC'])
		dpix = hp.ang2pix(nside,dth,dphi,nest)
		allpix = np.zeros(12*nside*nside)
		for pix in dpix:
			allpix[pix] += 1.
		return allpix

	def assignpix(fl):
		dth,dphi = radec2thphi(fl['RA'],fl['DEC'])
		dpix = hp.ang2pix(nside,dth,dphi,nest)
		return dpix

	#espix = assignpix(elgs)
	#enpix = assignpix(elgn)
	#lspix = assignpix(lrgs)
	#lnpix = assignpix(lrgn)
	#rpix = assignpix(rm)


	#get MC efficiency map
	mydir = '/project/projectdirs/desi/users/ajross'
	mcf = fitsio.read(mydir+'/ELGMCeffHSCHPsedext.fits') #new with improved snr cut
	mmc = np.mean(mcf['EFF'])
	mcl = np.ones(12*nside*nside)
	for i in range(0,len(mcf)):
		pix = mcf['HPXPIXEL'][i]
		mcl[pix] = mcf['EFF'][i]/mmc
	
	#put into full sky maps (probably not necessary but easier to keep straight down the line)
	pixlrbm = puthealpix(rm[rbml])
	pixlrdn = puthealpix(rm[rdnl])
	pixlrds = puthealpix(rm[rdsl])
	w = lrgn['DEC'] >  32.375
	pixllbm = puthealpix(lrgn[w])
	pixlldn = puthealpix(lrgs[lsdnl])
	pixllds = puthealpix(lrgs[lsdsl])
	w = elgn['DEC'] >  32.375
	pixlgbm = puthealpix(elgn[w])
	pixlgdn = puthealpix(elgs[esdnl])
	pixlgds = puthealpix(elgs[esdsl])

	title = 'DECaLS South, EBV factor '+str(ebf)+' Rv ='+str(rv)
	effac=2.
	slp = 0
	b = 1.
	ws = 1./(slp*hpq['STARDENS']+b)

	plotvshp_compmc(pixlrds,pixlgds,pixllds,'EBV',(0,0.15),title=title,effac=effac,fout='DECaLS_SGC')

	title = 'DECaLS North, EBV factor '+str(ebf)+' Rv ='+str(rv)
	effac=2.
	slp = 0
	b = 1.
	ws = 1./(slp*hpq['STARDENS']+b)

	plotvshp_compmc(pixlrdn,pixlgdn,pixlldn,'EBV',(0,0.15),title=title,effac=effac,fout='DECaLS_NGC')

	title = 'BMzLS, EBV factor '+str(ebf)+' Rv ='+str(rv)
	effac=1.
	slp = 0
	b = 1.
	ws = 1./(slp*hpq['STARDENS']+b)

	plotvshp_compmc(pixlrbm,pixlgbm,pixllbm,'EBV',(0,0.15),title=title,effac=effac,fout='BMzLS')


