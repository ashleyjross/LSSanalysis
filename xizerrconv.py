from numpy import loadtxt
from math import *
from EH import simulate
from scipy.special import jn

def sph_jn(l,x):
	return sqrt(pi/(2.*x))*jn(l+.5,x)

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

def mkxifile_zerrconvc(z=0.8,sigz=0.029,sp=1.,bias=1.8,rmin=10.,rmax=300,rsd='',muww='',a='',v='y',gam=-1.7,file='MICE_matterpower',mun=0,beta=0.4,sfog=0,sigt=6.,sigr=10.,sigs=15.,mumin=0,mumax=1):
	#Santi used zspec=0.45 to 1.2
	from random import gauss
	zc = zerrconv(file=file,mun=mun,beta=beta,sfog=sfog,sigt=sigt,sigr=sigr,sigs=sigs,gam=gam)
	wl,dzl,rzl = zc.calcwl(z,sigz)
	spf = 1.
	r = rmin
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)

	fo = open('xizconvc'+muww+file+muw+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigz)+rsd+'sp'+str(sp)+'.dat','w')
	nzs = 10000
	while r < rmax:
		
		xi,xin = zc.calcxi_zerrconv(r,wl,dzl,rzl,mumin=mumin,mumax=mumax)		
		fo.write(str(r)+' '+str(xi)+' '+str(xin)+'\n')		
		print r,xi,xin
		r += sp	 
	fo.close()
	return True

def mkxifile_zerrconvc_combzerr(sp=1.,bias=1.8,rmin=10.,rmax=300,rsd='',muww='',a='',v='y',gam=-1.7,file='MICE_matterpower',mun=0,beta=0.4,sfog=0,sigt=6.,sigr=10.,sigs=15.,mumin=0,mumax=1):
	#Santi used zspec=0.45 to 1.2
	from random import gauss
	spf = 1.
	r = rmin
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)

	fo = open('xizconvc'+muww+file+muw+str(beta)+str(sfog)+str(sigt)+str(sigr)+'combzerr'+rsd+'sp'+str(sp)+'.dat','w')
	
	zc = zerrconv(file=file,mun=mun,beta=beta,sfog=sfog,sigt=sigt,sigr=sigr,sigs=sigs,gam=gam)
	zl = [0.7,0.85,.95]
	sigl = [0.029,0.035,0.049]
	wxil = [0.5,.25,.25]
	wll= []
	dzll = []
	rzll = []
	for i in range(0,len(zl)):
		wl,dzl,rzl = zc.calcwl(zl[i],sigl[i])
		wll.append(wl)
		dzll.append(dzl)
		rzll.append(rzl)
	while r < rmax:
		xis = 0
		xisn = 0
		for i in range(0,len(zl)):
			wl,dzl,rzl = wll[i],dzll[i],rzll[i]		
			xi,xin = zc.calcxi_zerrconv(r,wl,dzl,rzl,mumin=mumin,mumax=mumax)
			xis += wxil[i]*xi
			xisn += wxil[i]*xin		
		fo.write(str(r)+' '+str(xis)+' '+str(xisn)+'\n')		
		print r,xis,xisn
		r += sp	 
	fo.close()
	return True

def mkxifile_zerrconvc_combzsigl(sp=1.,bias=1.8,rmin=10.,rmax=300,rsd='',muww='',a='',v='y',gam=-1.7,file='MICE_matterpower',mun=0,beta=0.4,sfog=0,sigt=6.,sigr=10.,sigs=15.,mumin=0,mumax=1):
	#Santi used zspec=0.45 to 1.2
	from random import gauss
	spf = 1.
	r = rmin
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)

	fo = open('xizconvc'+muww+file+muw+str(beta)+str(sfog)+str(sigt)+str(sigr)+'combzsigl'+rsd+'sp'+str(sp)+'.dat','w')
	
	zc = zerrconv(file=file,mun=mun,beta=beta,sfog=sfog,sigt=sigt,sigr=sigr,sigs=sigs,gam=gam)
	sigl = [0.031,0.029,0.028,0.029,0.033,0.038,0.044,0.052]
	#sigl = [0.029,0.029,0.035,0.049]
	wl,dzl,rzl = zc.calcwlave(.8,sigl)
	while r < rmax:
		xis = 0
		xisn = 0
		xi,xin = zc.calcxi_zerrconv(r,wl,dzl,rzl,mumin=mumin,mumax=mumax)
		fo.write(str(r)+' '+str(xi)+' '+str(xin)+'\n')		
		print r,xi,xin
		r += sp	 
	fo.close()
	return True


def compmocks(b=1.4,beta=.4):
	from matplotlib import pyplot as plt
	d = loadtxt('/Users/ashleyross/DESY1/xiSMlampavemu00.25.dat').transpose()
	dth = loadtxt('xizconvcMICE_matterpowermumax0.2'+str(beta)+'06.010.0combzsiglsp5.0.dat').transpose()
	plt.plot(d[0],d[0]**2.*d[1],'ko')
	plt.plot(dth[0],dth[0]**2.*dth[1]*b,'k-')
	plt.show()
	return True


class zerrconv:
	def __init__(self,file='MICE_matterpower',dir='',bias=1.8,mun=0,beta=0.4,sfog=0,amc=0.0,sigt=6.,sigr=10.,mult=1.,sigs=15.,gam=-1.7):
		dir = 'BAOtemplates/'
		self.f0 = loadtxt(dir+'xi0'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
		self.f2 = loadtxt(dir+'xi2'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
		self.f4 = loadtxt(dir+'xi4'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
		self.f0sm = loadtxt(dir+'xi0sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
		self.f2sm = loadtxt(dir+'xi2sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
		self.f4sm = loadtxt(dir+'xi4sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat').transpose()
		from Cosmo import distance
		if file == 'Challenge_matterpower':
			d = distance(.31,.69)
		if file == 'MICE_matterpower':
			d = distance(.25,.75)
		self.d = d
		self.beta = beta
		self.bias = bias
		self.gam = gam
		

	def calcwl(self,z,sigz):
		self.ff = self.d.omz(z)**.557
		self.betad = self.ff/self.bias
		#self.betad = self.beta
		self.betaf = self.betad/self.beta	

		zmin = z-sigz*(1.+z)*5.*sqrt(2.)
		zmax = z+sigz*(1.+z)*5.*sqrt(2.)
		nz = 4000
		dz = (zmax-zmin)/float(nz)
		d0 = self.d.dc(z)
		dzl = []
		rzl = []
		wl = []
		for i in range(0,nz):
			zb = zmin + dz/2.+dz*i
			dzl.append(self.d.dc(zb)-d0)
			rzl.append(self.d.dc(zb)/d0)
			wl.append(1./(sqrt(2.)*sigz*(1.+z)*sqrt(2.*pi))*exp(-.5*((zb-z)/(sqrt(2.)*sigz*(1.+z)))**2.))
		print sum(wl)*dz
		return wl,dzl,rzl

	def calcwlave(self,z,sigzl):
		#wl becomes average of the inputs from sigzl
		self.ff = self.d.omz(z)**.557
		self.betad = self.ff/self.bias
		self.betaf = self.betad/self.beta	

		zmin = z-max(sigzl)*(1.+z)*5.*sqrt(2.)
		zmax = z+max(sigzl)*(1.+z)*5.*sqrt(2.)
		nz = 4000
		dz = (zmax-zmin)/float(nz)
		d0 = self.d.dc(z)
		dzl = []
		rzl = []
		wl = []
		for i in range(0,nz):
			zb = zmin + dz/2.+dz*i
			dzl.append(self.d.dc(zb)-d0)
			rzl.append(self.d.dc(zb)/d0)
			w = 0
			for j in range(0,len(sigzl)):
				sigz = sigzl[j]
				w += 1./(sqrt(2.)*sigz*(1.+z)*sqrt(2.*pi))*exp(-.5*((zb-z)/(sqrt(2.)*sigz*(1.+z)))**2.)
			w = w/float(len(sigzl))	
			wl.append(w)
		print sum(wl)*dz
		return wl,dzl,rzl


	def calcxi_zerrconv(self,r,wl,dzl,rzl,sp=1.,rmin=10.,rmax=300,rsd='',mumin=0,mumax=1):
		#Santi used zspec=0.45 to 1.2
		#calculate just for one r
		bias = self.bias
		from random import gauss
		spf = 1.
		
		sumw = sum(wl)
		nz = len(wl)
		muwt = 0
		for i in range(int(100*mumin),int(100*mumax)):	
			muwt += 1.
		
		muw = ''
		if mumin != 0:
			muw += 'mumin'+str(mumin)
		if mumax != 1:
			muw += 'mumax'+str(mumax)
		nzs = 10000
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
				rrp = rr + dzl[i]
				rtp = rt
				rp = sqrt(rrp**2.+rtp**2.)
				mup = abs(rrp/rp)
				if rp >= 10. and rp < 299:
					indd = int((rp-10.)/spf)
					indu = indd + 1
					fac = (rp-10.)/spf-indd
					if fac > 1 or fac < 0:
						print fac,rp,spf,indd
					xi0 = (self.f0[1][indu]*fac+(1.-fac)*self.f0[1][indd])*(1.+2/3.*self.betad+.2*self.betad**2.)/(1.+2/3.*self.beta+.2*self.beta**2.)
					xi2 = (self.f2[1][indu]*fac+(1.-fac)*self.f2[1][indd])*(4/3.*self.betad+4/7.*self.betad**2.)/(4/3.*self.beta+4/7.*self.beta**2.)
					xi4 = (self.f4[1][indu]*fac+(1.-fac)*self.f4[1][indd])*(self.betad/self.beta)**2.					
					xi0n = self.f0sm[1][indu]*fac+(1.-fac)*self.f0sm[1][indd]*(1.+2/3.*self.betad+.2*self.betad**2.)/(1.+2/3.*self.beta+.2*self.beta**2.)
					xi2n = (self.f2sm[1][indu]*fac+(1.-fac)*self.f2sm[1][indd])*(4/3.*self.betad+4/7.*self.betad**2.)/(4/3.*self.beta+4/7.*self.beta**2.)
					xi4n = (self.f4sm[1][indu]*fac+(1.-fac)*self.f4sm[1][indd])*(self.betad/self.beta)**2.					
				else:
					if rp < 10.:
						xi0,xi2,xi4 = wmod3(rp,mup,self.f0[1][0],10.,gam=self.gam)
						xi0n,xi2n,xi4n = wmod3(rp,mup,self.f0sm[1][0],10.,gam=self.gam)
						#xi2 = f2[1][0]
						#xi4 = f4[1][0]
						#xi2n = f2sm[1][0]
						#xi4n = f4sm[1][0]

					if rp >= 300:#input files don't go beyond 300, just using maximum value
						xi0 = self.f0[1][-1]*(1.+2/3.*self.betad+.2*self.betad**2.)/(1.+2/3.*self.beta+.2*self.beta**2.)
						xi2 = self.f2[1][-1]*(4/3.*self.betad+4/7.*self.betad**2.)/(4/3.*self.beta+4/7.*self.beta**2.)
						xi4 = self.f4[1][-1]*(self.betad/self.beta)**2.					
						xi0n = self.f0sm[1][-1]*(1.+2/3.*self.betad+.2*self.betad**2.)/(1.+2/3.*self.beta+.2*self.beta**2.)
						xi2n = self.f2sm[1][-1]*(4/3.*self.betad+4/7.*self.betad**2.)/(4/3.*self.beta+4/7.*self.beta**2.)
						xi4n = self.f4sm[1][-1]*(self.betad/self.beta)**2.						
		
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
			sumxi += xi
			sumxin += xin
		sumxi = sumxi/muwt
		sumxin = sumxin/muwt			
		return sumxi,sumxin

	
def wmod3(r,mu,mf,r0,sp=1.,gam=-1.7):
	xi10 = mf
	xi0 = xi10*(r/r0)**(gam)
	xi2 = xi0-3.*xi0/(3.+gam)
	xi4 = xi0+2.5*3.*xi0/(3.+gam)-3.5*5.*xi0/(5.+gam)
	return xi0,xi2,xi4

def mkxifile_3dewig(sp=1.,a='',v='y',file='Challenge_matterpower',dir='',mun=0,beta=0.4,sfog=0,amc=0.0,sigt=6.,sigr=10.,mult=1.,sigs=15.):
	dir = 'BAOtemplates/'
	f0 = open(dir+'xi0'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat','w')
	f2 = open(dir+'xi2'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat','w')
	f4 = open(dir+'xi4'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat','w')
	f0mc = open(dir+'xi0sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat','w')
	f2mc = open(dir+'xi2sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat','w')
	f4mc = open(dir+'xi4sm'+file+str(beta)+str(sfog)+str(sigt)+str(sigr)+str(sigs)+str(mun)+'.dat','w')
	r = 10.

	while r < 300:
		xid = xi3elldfile_dewig(r,file,dir,beta=beta,sfog=sfog,sigt=sigt,sigr=sigr,sigs=sigs,mun=mun)
		f0.write(str(r)+' '+str(xid[0])+'\n')
		f2.write(str(r)+' '+str(xid[1])+'\n')
		f4.write(str(r)+' '+str(xid[2])+'\n')
		f0mc.write(str(r)+' '+str(xid[3])+'\n')
		f2mc.write(str(r)+' '+str(xid[4])+'\n')
		f4mc.write(str(r)+' '+str(xid[5])+'\n')
		r += sp
		if v == 'y':
			print r
	f0.close()
	f2.close()
	f4.close()
	f0mc.close()
	f2mc.close()
	f4mc.close()
	return True

def xi3elldfile_dewig(r,file='Challenge_matterpower',dir='',beta=0.4,sigt=3.0,sigr=3.0,sfog=3.5,max=51,mun=1.,sigs=15.,ns=.95,pw='n'):
	from scipy.integrate import quad
	#f = open('/Users/ashleyr/BOSS/spec/camb_MWcos.dat')
	#print dir
	mult = 1.
	dir = 'powerspectra/'
	if file=='Challenge_matterpower':
		om = 0.31
		lam = 0.69
		h = .676
		nindex = .963
		ombhh = .022
	if file == 'MICE_matterpower':
		om = 0.25
		lam = .75
		h = .7
		ombhh = .044*0.7*.7	
		nindex = .949

	f = open(dir+file+'.dat').readlines()
	if pw == 'y':
		fo = open('P02'+file+'beta'+str(beta)+'sigs'+str(sfog)+'sigxy'+str(sigt)+'sigz'+str(sigr)+'Sk'+str(sigs)+'.dat','w')
		fo.write('# k P0 P2 P4 Psmooth\n')

	s = simulate(omega=om,lamda=lam,h=h,nindex=nindex,ombhh=ombhh)
	pl2 = []
	pl4 = []
	beta0 = 0.4
	for i in range(0,100):
		pl2.append(P2(i/100.+.005))
	for i in range(0,100):
		pl4.append(P4(i/100.+.005))
	mul = []
	anipolyl = []
	for i in range(0,100):
		mu = i/100.+.005
		mul.append(mu)		
		anipoly = (1.+beta*mu**2.)**2.
		anipoly2 = (1.+beta*mu**2.)**2.*pl2[i]
		anipoly4 = (1.+beta*mu**2.)**2.*pl4[i]
		anipolyl.append((anipoly,anipoly2,anipoly4))
	k0 = float(f[0].split()[0])
	k1 = float(f[1].split()[0])
	ldk = log(k1)-log(k0)
	#print ldk
	xi0 = 0
	xi2 = 0
	xi4 = 0
	xism0 = 0
	xism2 = 0
	xism4 = 0
	#print max
	dmk = r/10.
	k = float(f[0].split()[0])
	if file == 'camb_Nacc':
		norm = 1.
	else:
		norm = float(f[0].split()[1])/s.Psmooth(k,0)*mult
	#print norm
	b = 2.
	ff =beta*b
	for i in range(0,len(f)-1):
		k = float(f[i].split()[0])
		pk = float(f[i].split()[1])*mult
		pk0 = 0
		pk2 = 0
		pk4 = 0
		pksm0 = 0
		pksm2 = 0
		pksm4 = 0
		pksm = s.Psmooth(k,0)*norm
		dpk = pk-pksm
		for m in range(0,100):
			#mu = (1.-mul[m])			
			mu = mul[m]

			if mun == 'n':
				mu = (1.-mul[m])
			if sfog > 0:
				F = 1./(1.+k**2.*mu**2.*sfog**2./2.)**2.
			else:
				F = (1.+k**2.*sfog**2./2.)**2./(1.+k**2.*(1.-mu)**2.*sfog**2./2.)**2.
			if mun == 'b':
				mus2 = mu**2.
				F = (1.+beta*mus2*(1.-exp(-0.5*(k*sigs)**2.)))**2.*1./(1.+k**2.*mu**2.*(sfog)**2./2.)**2.
			#C *doesn't include damping*
			S = mun*exp(-0.5*(k*sigs)**2.)
			C = (1.+beta*mu*mu*(1.-S))*1./(1.+k**2.*mu**2.*(sfog)**2./2.)
			sigv2 = (1-mu**2.)*sigt**2./4.+mu**2.*sigr**2./4.
			damp = exp(-1.*k*k*sigv2)	
			#damp1 = exp(-0.25*k**2.*mu**2.*sigr**2.*(1.+ff)**2.)
			#damp2 = exp(-0.25*k**2.*(1.-mu**2.)*sigt**2.)
			#damp3 = 1./b*S*exp(-0.5*k**2.*sigt*sigr)
			#damp4 = (1.-exp(-0.25*k**2.*mu**2.*(2.*ff+ff**2.)*sigt*sigr))
			#damp = damp1*damp2+damp3*damp4 	
			
			anipoly = anipolyl[m]
			pk0 += C**2.*(dpk*damp**2.+pksm)
			pk2 += (dpk*damp**2.+pksm)*pl2[m]*C**2.
			pk4 += (dpk*damp**2.+pksm)*pl4[m]*C**2.
			pksm0 += pksm*C**2.
			pksm2 += pksm*pl2[m]*C**2.
			pksm4 += pksm*pl4[m]*C**2.
		pk0 = pk0/100.
		pk2 = 5.*pk2/100.
		pk4 = 9.*pk4/100.
		pksm0 = pksm0/100.
		pksm2 = 5.*pksm2/100.
		pksm4 = 9.*pksm4/100.
		if pw == 'y':
			fo.write(str(k)+' '+str(pk0)+' '+str(pk2)+' '+str(pk4)+' '+str(pksm0)+' '+str(pksm2)+' '+str(pksm4)+' '+str(pk)+' '+str(pksm)+'\n')
		dk = ldk*k
		xi0 += dk*k*k*pk0*sph_jn(0,k*r)*exp(-1.*dmk*k**2.)
		xi2 += dk*k*k*pk2*sph_jn(2,k*r)*exp(-1.*dmk*k**2.)
		xi4 += dk*k*k*pk4*sph_jn(4,k*r)*exp(-1.*dmk*k**2.)
		xism0 += dk*k*k*pksm0*sph_jn(0,k*r)*exp(-1.*dmk*k**2.)
		xism2 += dk*k*k*pksm2*sph_jn(2,k*r)*exp(-1.*dmk*k**2.)
		xism4 += dk*k*k*pksm4*sph_jn(4,k*r)*exp(-1.*dmk*k**2.)
	if pw == 'y':
		fo.close()
		from matplotlib import pyplot as plt
		from numpy import loadtxt as load
		d = load('P02'+file+'beta'+str(beta)+'sigs'+str(sfog)+'sigxy'+str(sigt)+'sigz'+str(sigr)+'Sk'+str(sigs)+'.dat').transpose()
		plt.xlim(0,.3)
		plt.plot(d[0],d[-2]/d[-1])
		plt.show()
	return xi0/(2.*pi*pi),-1.*xi2/(2.*pi*pi),xi4/(2.*pi*pi),xism0/(2.*pi*pi),-1.*xism2/(2.*pi*pi),xism4/(2.*pi*pi)
	
