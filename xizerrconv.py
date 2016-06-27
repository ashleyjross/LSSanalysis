from numpy import loadtxt
from math import *

def P2(mu):
	return .5*(3.*mu**2.-1.)
	
def P4(mu):
	return .125*(35.*mu**4.-30.*mu**2.+3.)


def mkxifile_zerrconv(z=0.8,sigz=0.029,sp=1.,bias=1.8,rmin=10.,rmax=300,rsd='',muww='muw',a='',v='y',gam=-1.7,file='MICE_matterpower',dir='',mun=0,beta=0.4,sfog=0,amc=0.0,sigt=6.,sigr=10.,mult=1.,sigs=15.,mumin=0,mumax=1):
	#Santi used zspec=0.45 to 1.2
	from random import gauss
	dir = '/Users/ashleyross/DR12/'
	if file == 'MICE_matterpower':
		dir = '/Users/ashleyross/DESY1/'
	f0 = loadtxt(dir+'xi0'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
	f2 = loadtxt(dir+'xi2'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
	f4 = loadtxt(dir+'xi4'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
	f0sm = loadtxt(dir+'xi0sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
	f2sm = loadtxt(dir+'xi2sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
	f4sm = loadtxt(dir+'xi4sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
	spf = 1.
	r = rmin
	from Cosmo import distance
	if file == 'Challenge_matterpower':
		d = distance(.31,.69)
	if file == 'MICE_matterpower':
		d = distance(.25,.75)
	ff = d.omz(z)**.557
	betad = ff/bias
	betaf = betad/beta	
	zmin = z-sigz*(1.+z)*5.*sqrt(2.)
	zmax = z+sigz*(1.+z)*5.*sqrt(2.)
	#zmin = z-.15
	#zmax = z+.2
	nz = 4000
	dz = (zmax-zmin)/float(nz)
	d0 = d.dc(z)
	dzl = []
	rzl = []
	wl = []
	for i in range(0,nz):
		zb = zmin + dz/2.+dz*i
		#print zb
		dzl.append(d.dc(zb)-d0)
		rzl.append(d.dc(zb)/d0)
		#rzl.append(1.)
		wl.append(1./(sqrt(2.)*sigz*(1.+z)*sqrt(2.*pi))*exp(-.5*((zb-z)/(sqrt(2.)*sigz*(1.+z)))**2.))
	print sum(wl)*dz
	sumw = sum(wl)
	#print wl
	#fmuw = loadtxt('/Users/ashleyross/DESY1/FdaHvsmu_z0.61.0_zerr0.03_10e3n2.6_b1.5.dat').transpose()
	fmuw = loadtxt('/Users/ashleyross/DESY1/FmuobsdaHvsmu_z0.61.0_zerr0.03_10e3n2.6_b1.5.dat').transpose()
	#muwt = sum(fmuw[1])
	muwt = 0
	#for i in range(0, len(fmuw[1])):
	#	if i > int(mumin*100) and i < int(mumax*100):
	for i in range(int(100*mumin),int(100*mumax)):	
		if muww == 'muw':
			muwt += fmuw[1][i]
		else:
			muwt += 1.
				
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	fo = open('xizconv'+muww+file+muw+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigz)+rsd+'sp'+str(sp)+'.dat','w')
	nzs = 10000
	while r < rmax:
		sumxi = 0
		sumxin = 0
		mmin = int(100*mumin)
		mmax = int(100*mumax)
		for m in range(mmin,mmax):
			mu = .005+0.01*m
			summ = 0
			summn = 0
			rt = sqrt(1.-mu**2.)*r
			rr = mu*r
			for i in range(0,nz):
				#dzt = gauss(0,sigz*(1.+z))-gauss(0,sigz*(1.+z))
				#indz = int((z+dzt-zmin)/dz)
				#if indz >= 0 and indz < nz: 
				rrp = rr + dzl[i]
				#if abs(rrp) < rr:
				#	print rrp,rr,r
				#rtp = rt*rzl[i]
				rtp = rt
				rp = sqrt(rrp**2.+rtp**2.)
				mup = abs(rrp/rp)
				if rp >= 10. and rp < 299:
					indd = int((rp-10.)/spf)
					indu = indd + 1
					fac = (rp-10.)/spf-indd
					if fac > 1 or fac < 0:
						print fac,rp,spf,indd
					xi0 = (f0[1][indu]*fac+(1.-fac)*f0[1][indd])*(1.+2/3.*betad+.2*betad**2.)/(1.+2/3.*beta+.2*beta**2.)
					xi2 = (f2[1][indu]*fac+(1.-fac)*f2[1][indd])*(4/3.*betad+4/7.*betad**2.)/(4/3.*beta+4/7.*beta**2.)
					xi4 = (f4[1][indu]*fac+(1.-fac)*f4[1][indd])*(betad/beta)**2.					
					xi0n = f0sm[1][indu]*fac+(1.-fac)*f0sm[1][indd]*(1.+2/3.*betad+.2*betad**2.)/(1.+2/3.*beta+.2*beta**2.)
					xi2n = (f2sm[1][indu]*fac+(1.-fac)*f2sm[1][indd])*(4/3.*betad+4/7.*betad**2.)/(4/3.*beta+4/7.*beta**2.)
					xi4n = (f4sm[1][indu]*fac+(1.-fac)*f4sm[1][indd])*(betad/beta)**2.					
				else:
					if rp < 10.:
						xi0,xi2,xi4 = wmod3(rp,mup,f0[1][0],10.,gam=gam)
						xi0n,xi2n,xi4n = wmod3(rp,mup,f0sm[1][0],10.,gam=gam)
						#xi2 = f2[1][0]
						#xi4 = f4[1][0]
						#xi2n = f2sm[1][0]
						#xi4n = f4sm[1][0]

					if rp >= 300:
						xi0 = f0[1][-1]*(1.+2/3.*betad+.2*betad**2.)/(1.+2/3.*beta+.2*beta**2.)
						xi2 = f2[1][-1]*(4/3.*betad+4/7.*betad**2.)/(4/3.*beta+4/7.*beta**2.)
						xi4 = f4[1][-1]*(betad/beta)**2.					
						xi0n = f0sm[1][-1]*(1.+2/3.*betad+.2*betad**2.)/(1.+2/3.*beta+.2*beta**2.)
						xi2n = f2sm[1][-1]*(4/3.*betad+4/7.*betad**2.)/(4/3.*beta+4/7.*beta**2.)
						xi4n = f4sm[1][-1]*(betad/beta)**2.					
					
						xi0,xi2,xi4 = 0,0,0	
						xi0n,xi2n,xi4n = 0,0,0
				if rsd == 'norsd':
					xi2,xi4,xi2n,xi4n = 0,0,0,0
				xip = xi0+P2(mup)*xi2+P4(mup)*xi4
				xipn = xi0n+P2(mup)*xi2n+P4(mup)*xi4n
				#if xip > 10:
				#	print xip,rp,indd,indu,f4[1][indu],f4[1][indd],P4(mup),mup
				summ += xip*wl[i]
				summn += xipn*wl[i]
			xi = summ/sumw#/float(nzs)#
			xin = summn/sumw#/float(nzs)#
			if muww == 'muw':
				sumxi += xi*fmuw[1][m]
				sumxin += xin*fmuw[1][m]
			else:
				sumxi += xi
				sumxin += xin
		sumxi = sumxi/muwt
		sumxin = sumxin/muwt		
		fo.write(str(r)+' '+str(sumxi)+' '+str(sumxin)+'\n')		
		print r,sumxi,sumxin
		r += sp	 
	fo.close()
# 	from matplotlib import pyplot as plt
# 	from numpy import loadtxt as load
# 	d = load('xizconvmuw'+file+muw+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigz)+'sp'+str(sp)+'.dat').transpose()
# 	dd = load('xiavemuw'+muw+'SM.dat').transpose()
# 	plt.plot(d[0],d[0]**2.*d[1]*1.3,'k-',dd[0],dd[0]**2.*dd[1],'r-')
# 	plt.show()
# 	plt.plot(d[0]/.97,d[0]**2.*d[1]*1.3,'k-',dd[0],dd[0]**2.*dd[1],'r-')
# 	plt.show()		
	return True
	
def wmod3(r,mu,mf,r0,sp=1.,gam=-1.7):
	xi10 = mf
	xi0 = xi10*(r/r0)**(gam)
	xi2 = xi0-3.*xi0/(3.+gam)
	xi4 = xi0+2.5*3.*xi0/(3.+gam)-3.5*5.*xi0/(5.+gam)
	return xi0,xi2,xi4
	
