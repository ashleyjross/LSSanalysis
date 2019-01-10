#from math import *
from cmath import *
from Cosmo import *
from romberg import rom
from numpy import array
from numpy.linalg import *
from numpy import loadtxt as load
#from pklin import *
#from simulate import *
from bessel import *
from scipy.special import jn
from legendre import *
#from nz import *


def sph_jn(l,x):
	return sqrt(pi/(2.*x))*jn(l+.5,x)

def diserr(z1,z2,ang,dz,dz2):
	from Cosmo import distance
	angrad = ang*pi/180.
	d = distance(.275,.725)
	cz1 = d.dc(z1)
	cz2 = d.dc(z2)
	cra = cos(angrad)
	dsq = cz1**2.+cz2**2.-2.*cz1*cz2*cra
	dis = sqrt(dsq)
	dd1 = (d.dc(z1+dz)-d.dc(z1-dz))/2. 
	#print cz1,dd1
	dd2 = (d.dc(z2+dz2)-d.dc(z2-dz2))/2.
	ddsq = (2.*dd1*cz1-2.*dd1*cz2*cra)**2.+(2.*dd2*cz2-2.*dd2*cz1*cra)**2.
	dd = sqrt(ddsq)
	derr = .5*dd/dis
	print dd1,dd2,ddsq,dis
	return derr

def thbao(z1,z2,rbao=105.):
	from Cosmo import distance
	d = distance(.275,.725)
	cz1 = d.dc(z1)
	cz2 = d.dc(z2)
	cra = rbao**2./(cz1**2.+cz2**2.-2.*cz1*cz2)
	ang = acos(cra)*180./pi
	return ang


def mkACmeanpc():
	thl = []
	ml = []
	mlth = []
	stdl = []
	for i in range(0,30):
		ml.append(0)
		mlth.append(0)
		stdl.append(0)
	for i in range(1,83):
		if i != 42:
			f = open('ACpc0/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-0.03').readlines()
			for j in range(0,len(f)):
				ml[j] += float(f[j].split()[2])/81.
				mlth[j] += float(f[j].split()[1])/81.
	for i in range(1,83):
		if i != 42:
			f = open('ACpc0/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-0.03').readlines()
			for j in range(0,len(f)):
				stdl[j] += (ml[j]-float(f[j].split()[2]))**2./81.
	fo = open('ACpcmean0.05sigz0.045.dat','w')
	f = open('ACpc0/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-1_zm-0.50_dz-0.03').readlines()
	for i in range(0,30):
		fo.write(f[i].split()[0]+' '+str(ml[i])+' '+str(sqrt(stdl[i]))+' '+str(mlth[i])+'\n')
	fo.close()
class p_proj:
	def __init__(self,dNdZfile='nzTH0.40.60.05',zmax=2.,zmin=.01,zmed=.5):
		self.d = distance()
		self.zmed = zmed
		self.dndz,minmin,self.maxmax,self.step = zbinnerb(dNdZfile,zmax,zmin)
		
		
	def Wfunc(self,z):
		#from cmath import *
		d = distance()
		rz = abs(d.dc(z)-d.dc(self.zmed))
		#rz = d.dc(z)
		return self.nz(z)*exp(-1.*sqrt(-1.)*self.kz*rz)

	def Wsquare(self,kz):
		self.kz = kz
		ans = 0
		for i in range(0,len(self.dndz)):
			z = i/1000.
			ans += self.Wfunc(z)*.001
		#return (2997.9*ans)**2.
		return abs(ans)**2.
		#return (rom(0.01,.5,self.Wfunc)+rom(.5,1.,self.Wfunc))**2.

	def W(self,kz):
		self.kz = kz
		ans = 0
		for i in range(0,len(self.dndz)):
			z = i/1000.
			ans += self.Wfunc(z)*.001
		#return (2997.9*ans)**2.
		return abs(ans)


	def pprojfunc(self,z):
		dis = self.d.dc(z)
		p = self.s.P(self.k/dis,z)
		return self.nz(z)**2.*2997.9/self.d.Hz(z)*p*self.k**3./(2.*pi**2.)*self.d.dc(z)
		
	def pprojlimb(self,k):
		self.s = simulate()
		self.k = k
		ans = 0
		for i in range(0,1000):
			z = i/1000.+.0005
			ans += self.pprojfunc(z)*.001
		return ans*pi/k


	def nz(self,z):
		
		ind = int(z*1000)
		if ind < len(self.dndz) and ind >0:
			return self.dndz[ind]
		else:
			return 0

def tophatGaussconv(zmin,zmax,sig):
	from random import random,gauss
	zphot = (zmax-zmin)*random()+zmin
	ztrue = gauss(zphot,sig)
	return ztrue

def mkClsigfile_limb(dz,tz=.2,zsig=.05,nz=1000000.,zmed=.5,nbar=.0003,fsky=1/8.,dl=10,lmax=300):
	from Cosmo import distance
	from numpy import zeros
	n = tz/dz
	print n
	n = int(n)
	print n
	zls = zeros((n,2000),'f')
	d = distance(.25,.75)
	zmin = zmed-tz/2.
	zmax = zmed+tz/2.
	zb = []
	for i in range(0,n+1):
		zb.append(zmin+i*dz)
	nl = []
	for i in range(0,n):
		z1 = zb[i]
		z2 = zb[i+1]
		dmin = d.dc(z1)
		dmax = d.dc(z2)
		vol = 4/3.*pi*(dmax**3.-dmin**3.)
		nn = vol*nbar
		np = nn/(4.*pi)
		nl.append(np)

		for j in range(0,nz):
			z = tophatGaussconv(z1,z2,zsig)
			zind = int(1000*z)
			zls[i][zind] += 1.
	print 'zl done'
	
	fo = open('clsigTHG'+str(tz)+str(dz)+'.dat','w')
	print nl	
	l = dl
	while l < lmax:
		cl = 0
		sigcl = 0
		for i in range(0,n):
			cli = cl_limb_zl(l,zls[i])
			#print cli
			cl += cli
			sigcl += 2./(2.*l+1)*(cli+1./nl[i])**2.
		p=1
		while p < n:
			for i in range(0,n):
				if i+p >= n:
					break
				#cli = cl_limbcross_zl(l,zls[i],zls[i+p])
				#cl += cli
				sigcl += 2./(2.*l+1)*(cli)**2.
			p+=1
		fo.write(str(l)+' '+str(cl)+' '+str(sigcl)+'\n')
		l += dl
	fo.close()
	return True
	
def mkClsigfilepz_limb(zz,wf,tz=.2,zsig=.05,nz=1000000.,zmed=.5,nbar=.0003,fsky=1/8.,dl=10,lmax=300):
	from Cosmo import distance
	from numpy import zeros
	from random import random
	n = tz/zsig
	print n
	zls = zeros((2000),'f')
	d = distance(.25,.75)
	zmin = zmed-tz/2.
	zmax = zmed+tz/2.
	zsg1 = zz-zsig/2.
	zsg2 = zz+zsig/2.
	nl = []
	dmin = d.dc(zmin)
	dmax = d.dc(zmax)
	vol = 4/3.*pi*(dmax**3.-dmin**3.)
	nn = vol*nbar
	np = nn/(4.*pi)
	gco = 1./(zsig*sqrt(2.*pi))
	for i in range(0,nz):
		z = (zmax-zmin)*random()+zmin
		sig = abs((z-zz)/zsig)
		w = exp(-.5*sig**2.*wf)
		zind = int(1000*z)
		zls[zind] += w
	print 'zl done'
	sum = 0
	for i in range(0,len(zls)):
		sum += zls[i]#*.001
	print sum
	neff = nn/nz*sum
	print neff
	fo = open('clsigPZ'+str(zz)+str(zsig)+str(wf)+'.dat','w')
	l = dl
	#minn = (zmin-zz)/zsig
	#maxn = (zmax-zz)/zsig
	#frn = 1.5*sqrt(pi)*erf(maxn)+4*exp(-1.*maxn**2.)*(-.5*maxn**3.-.75*maxn)-(1.5*sqrt(pi)*erf(minn)+4*exp(-1.*minn**2.)*(-.5*minn**3.-.75*minn))
	#frd = sqrt(pi)*erf(maxn)-2.*exp(-1.*maxn**2.)*maxn-(sqrt(pi)*erf(minn)-2.*exp(-1.*minn**2.)*minn)
	#fr = (1.-exp(-.5*n**2.))/(sqrt(pi/2.)*erf(n/sqrt(2.)))
	npf = neff/(4.*pi)
	#fr = frn/frd**2.
	#print fr
	while l < lmax:
		cli = cl_limb_zl(l,zls)
		sigcl = 2./(2.*l+1)*(cli+1./npf)**2.
		fo.write(str(l)+' '+str(cli)+' '+str(sigcl)+'\n')
		l += dl
	fo.close()
	return True

def mkClsigfilepz_series_limb(dz,wf,tz=.2,zsig=.05,nz=1000000.,zmed=.5,nbar=.0003,fsky=1/8.,dl=10,lmax=300):
	from Cosmo import distance
	from numpy import zeros
	from random import random
	n = tz/dz
	print n
	n = int(n)
	print n
	zls = zeros((n,2000),'f')
	d = distance(.25,.75)
	zmin = zmed-tz/2.
	zmax = zmed+tz/2.
	zb = []
	for i in range(0,n):
		zb.append(zmin+i*dz+dz/2.)

	nl = []
	dmin = d.dc(zmin)
	dmax = d.dc(zmax)
	vol = 4/3.*pi*(dmax**3.-dmin**3.)
	nn = vol*nbar
	np = nn/(4.*pi)
	gco = 1./(zsig*sqrt(2.*pi))
	neffl = []
	for j in range(0,n):
		zz = zb[j]
		for i in range(0,nz):
			z = (zmax-zmin)*random()+zmin
			sig = abs((z-zz)/zsig)
			w = exp(-.5*sig**2.*wf)
			zind = int(1000*z)
			zls[j][zind] += w
		sum = 0
		for i in range(0,len(zls[j])):
			sum += zls[j][i]#*.001
		print sum
		neff = nn/nz*sum
		neffl.append(neff/(4.*pi))
	print neffl

	print 'zl done'
	fo = open('clsigPZs'+str(dz)+str(zsig)+str(wf)+'.dat','w')
	l = dl
	while l < lmax:
		cl = 0
		sigcl = 0
		for i in range(0,n):
			cli = cl_limb_zl(l,zls[i])
			#print cli
			cl += cli
			sigcl += 2./(2.*l+1)*(cli+1./neffl[i])**2.
		p=1
		while p < n:
			for i in range(0,n):
				if i+p >= n:
					break
				#cli = cl_limbcross_zl(l,zls[i],zls[i+p])
				#cl += cli
				sigcl += 2./(2.*l+1)*(cli)**2.
			p+=1
		fo.write(str(l)+' '+str(cl)+' '+str(sigcl)+'\n')
		l += dl


	while l < lmax:
		cli = cl_limb_zl(l,zls)
		sigcl = 2./(2.*l+1)*(cli+1./npf)**2.
		fo.write(str(l)+' '+str(cli)+' '+str(sigcl)+'\n')
		l += dl
	fo.close()
	return True
	
	

def mkradistPL(gam=-1.8,np=1000000.):
	from random import random
	rl = []
	oldr = random()
	dmin = .01
	pmin = 1.
	for i in range(0,np):
		prob = random()
		d = dmin*exp(log(prob)/gam)
		cf = random()
		if cf < .5:
			r = oldr +d
		else:
			r = oldr - d
		rl.append(r)
		oldr = r
	return rl
			
def mkGradfile(errmin=.025,errmax=.075):
	from random import random, gauss
	rl = mkradistPL()
	fo = open('radistGauss.dat','w')
	for i in range(0,len(rl)):
		err = (errmax-errmin)*random()+errmin
		nr = gauss(rl[i],err)
		fo.write(str(nr)+' '+str(err)+' '+str(rl[i])+'\n')
	fo.close()
	

def mkWsqfile(nzfile,zm):
	p = p_proj(dNdZfile=nzfile,zmed=zm)
	kmin = .00001
	mult = 1.01
	kmax = 10000
	k = kmin
	fo = open('Wsq'+nzfile+'.dat','w')
	while k < kmax:
		w = p.Wsquare(k)
		fo.write(str(k)+' '+str(w.real)+'\n')
		k = k*mult
	fo.close()

def mkWfile(nzfile,zm):
	p = p_proj(dNdZfile=nzfile,zmed=zm)
	kmin = .00001
	mult = 1.01
	kmax = 10000
	k = kmin
	fo = open('W'+nzfile+'.dat','w')
	while k < kmax:
		w = p.W(k)
		fo.write(str(k)+' '+str(w.real)+'\n')
		k = k*mult
	fo.close()

def Perr(k,zm,ngal,v,dk=.001):
	s = simulate(.25,.75)
	nbar = ngal/v
	kvol = 4.*pi*k**2*dk
	sig = sqrt((2.*pi)**3./(kvol*v))*(1.+1./(nbar*s.P(k,0)))*s.P(k,0)
	#sig = sqrt((2.*pi)**3./(kvol*v))*s.P(k,0)
	return sig

def nr2veff(file,pk):
	f = open(file).readlines()
	#p = 10000.
	sum = 0
	#fo = open(file.strip('nr.dat')+'veff.dat','w')
	for i in range(0,len(f)):
		n = float(f[i].split()[1])
		v = float(f[i].split()[0])
		veff = v*(n*pk/(1.+n*pk))**2.
		#fo.write(f[i].split()[2]+' '+str(veff)+'\n')
		sum += veff
	return sum


def xicov(file,r1,r2,dk=.03,kmin=.001,kmax=10.,bias=1.4):
	#from Cosmo import *
	#d = distance(.25,.75)
	s = simulate(.25,.75)
	#v = 4.*pi/3.*d.dc(zm)**3.*0.177675
	
	#nbar = ngal/v
	#print v,nbar
	#veff = zm
	k = kmin
	err2 = 0
	xi = 0
	while k < kmax:
		kvol = 4.*pi*k**2*dk*k
		pk = s.P(k,0)
		#sig = sqrt((2.*pi)**3./(kvol*v))*(pk+1./nbar)
		sig = sqrt((2.*pi)**3./(kvol*nr2veff('nzDR7spece2c'+file+'nr.dat',pk)))*(pk*bias**2.)
		err2 += (k*k*dk*sig)**2.*sin(r1*k)/r1*sin(r2*k)/r2
		xi += k*k*dk*pk*sin(r1*k)/r1*bias**2.
		k += dk*k
		#print k
	return xi/(2.*pi**2),err2/(2.*pi**2)/(2.*pi**2)
	#return err2/(2.*pi**2)

def xicovl(r1,r2,pkl,vfl,kl,dk=.03,bias=1.):
	err2 = 0
	xi = 0
	for i in range(0,len(kl)):
		k = kl[i]
		#kvol = 4.*pi*k**2*dk*k
		pk = pkl[i]
		#sig = sqrt((2.*pi)**3./(kvol*v))*(pk+1./nbar)
		#sig = sqrt((2.*pi)**3./(kvol*vfl[i]))*(pkl[i]*bias**2.)
		#err2 += (k*k*dk*sig)**2.*sin(r1*k)/r1*sin(r2*k)/r2
		err2 += k*dk*(pkl[i]*bias**2.)**2.*sin(r1*k)/r1*sin(r2*k)/r2/vfl[i]
		xi += k*k*dk*pk*sin(r1*k)/r1
	#return xi/(2.*pi**2),err2/(2.*pi**2)
	#return err2/(2.*pi**2)/(2.*pi**2)
	return err2/(2.*pi**2)


def mkxicovfile(file,kmin=.001,kmax=10.,bias=1.4,dk=.03):
	s = simulate(.25,.75)
	fo = open('covxi'+file+'a.dat','w')
	r0 = 1.1059869434
	mult = 1.35285110279/r0
	rl = []
	r = r0
	kl = []
	vfl = []
	pkl = []
	k = kmin
	while k < kmax:
		kl.append(k)
		pk = s.P(k,0)
		pkl.append(pk)
		vfl.append(nr2veff('nzDR7spece2c'+file+'nr.dat',pk))
		k += k*dk
	print 'lists done'
	while r < 210.:
		rl.append(r)
		r = r*mult
	for i in range(0,len(rl)):
		for j in range(0,len(rl)):
			c = xicovl(rl[i],rl[j],pkl,vfl,kl,bias=1.4,dk=.03)
			if j != len(rl)-1:
				fo.write(str(c)+' ')
			else:
				fo.write(str(c)+'\n')
	fo.close()


def cl(l,zf='THzsig0.030.450.55nzDES',zmed=.5,md='norm',bias=2.,kmax=.2,om=.279):
	from simulate import wp, simulate
	sum = 0
	s = wp(om=om,lam=1.-om)
	s.Dzf = 1.
	s.Ps = s.Pspline()	
	sim = simulate(omega=om,lamda=1.-om,h=.7,h0=1.,ombaryon=.046,sig8=.8)

	f = open('phiL/phi'+str(l)+str(om)+zf+'.dat').readlines()
	if md == 'RS':
		frs = open('phiL/phiRS'+str(l)+str(om)+zf+str(bias)+'.dat').readlines()
	for i in range(0,len(f)):
		ln = f[i].split()
		k = float(ln[0])
		if k > kmax:
			break
		dk = .03*k
		#damp = exp(-1.*k**2.*sigv**2.)
		damp = 1.
		rs = 0
		if md == 'RS':
			rs += float(frs[i].split()[1])
		sum += dk*k*k*s.Psplev(k)*damp*(float(ln[1])+rs)**2.
		#print k, sum
	return 2/pi*sum

def cl_fnl(l,zf='THzsig0.030.450.55nzDES',zmed=.5,md='norm',bias=2.,kmax=.2,om=.279):
	from simulate import wp, simulate
	sum = 0
	sum1 = 0
	sum2 = 0
	pm = 0
	try:
		fp = open('pk4fnl'+str(om)+'.dat').readlines()
	except:
		fp = open('pk4fnl'+str(om)+'.dat','w')
		pm = 1
		s = wp(om=om,lam=1.-om)
		s.Dzf = 1.
		s.Ps = s.Pspline()	
		
	sim = simulate(omega=om,lamda=1.-om,h=.7,h0=1.,ombaryon=.046,sig8=.8)
	
	#s = simulate(.274,.726)
	f = open('phiL/phi'+str(l)+str(om)+zf+'.dat').readlines()
	fD1 = open('phiL/phi'+str(l)+str(om)+zf+'D1.dat').readlines()
	if md == 'RS':
		frs = open('phiL/phiRS'+str(l)+zf+str(bias)+'.dat').readlines()
	#from Cosmo import *
	#d = distance(.25,.75)
	
	#sigv = 7.*d.D(zmed)
	for i in range(0,len(f)):
		ln = f[i].split()
		psiD1 = float(fD1[i].split()[1])
		k = float(ln[0])
		if pm == 0:
			kp = float(fp[i].split()[0])
			if kp != k:
				return 'mismatched k values!'
			p = float(fp[i].split()[1])
		else:
			p = s.Psplev(k)	
			fp.write(str(k)+' '+str(p)+'\n')
		if k > kmax:
			break
		dk = .03*k
		#damp = exp(-1.*k**2.*sigv**2.)
		damp = 1.
		rs = 0
		if md == 'RS':
			rs += float(frs[i].split()[1])
		
		sum += dk*k*k*p*damp*(float(ln[1])+rs)**2.
		#if k < .1:
		tk = sim.T(k)
		sum1 += dk*p*damp*(float(ln[1])*psiD1)/tk
		sum2 += dk*p*damp*(psiD1)**2./k/k/tk**2.
		#print k, sum
	if pm == 1:
		fp.close()
	return 2/pi*sum,2/pi*sum1,2/pi*sum2


def cl_cross(l,zf1,zf2='THzsig0.030.450.55nzDES',zmed=.5,md='norm',beta=0.7,nf='n'):
	sum = 0
	s = simulate(.25,.75)
	if nf == 'n':
		f = open('phiL/phi'+str(l)+zf1+'.dat').readlines()
		f2 = open('phiL/phi'+str(l)+zf2+'.dat').readlines()
	else:
		f = open('phiL/phi'+str(l)+zf1+'nf.dat').readlines()
		f2 = open('phiL/phi'+str(l)+zf2+'nf.dat').readlines()
	
	if md == 'RS':
		frs = open('phiL/phiRS'+str(l)+zf1+str(beta)+'.dat').readlines()
		frs2 = open('phiL/phiRS'+str(l)+zf2+str(beta)+'.dat').readlines()
	d = distance(.25,.75)
	
	sigv = 7.*d.D(zmed)
	slp = (float(f[1].split()[1])*float(f2[1].split()[1])-float(f[0].split()[1])*float(f2[0].split()[1]))/.00003
	#sum = (float(f[0].split()[1])*float(f2[0].split()[1])-.0005*slp)*.001*.0005**2.*s.P(.0005,0)
	#phisq = (float(f[0].split()[1])*float(f2[0].split()[1])-.0005*slp)
	#print sum,slp,phisq,s.P(.0005,0),s.P(.001,0)
	for i in range(0,len(f)):
		ln = f[i].split()
		ln2 = f2[i].split()
		k = float(ln[0])
		#if k*sigv > 5.:
		#if k > .03:
		#	break
		#dk = float(f[i+1].split()[0])
		if i == 0:
			dk = k
		else:
			if i < len(f) - 1:
				dk = (float(f[i+1].split()[0])-float(f[i-1].split()[0]))/2.
			else:
				dk = .03*k
		damp = exp(-1.*k**2.*sigv**2.)
		#damp = 1.
		rs = 0
		rs2 = 0
		if md == 'RS':
			rs += float(frs[i].split()[1])
			rs2 += float(frs2[i].split()[1])
		phisq = (float(ln[1])+rs)*(float(ln2[1])+rs2)
		sum += dk*k*k*s.P(k,0)*damp*phisq
		oldk = k
		#print k, sum,phisq
	return 2/pi*sum


def cl_pc(l,zf='zsig0.030.450.55',md='norm'):
	sum = 0
	s = simulate(.25,.75)
	f = open('phiL/phisqpc'+str(l)+zf+'.dat').readlines()
			
	sigv = 7.
	for i in range(0,len(f)):
		ln = f[i].split()
		k = float(ln[0])
		#if k*sigv > 5.:
		#if k > .03:
		#	break
		#dk = float(f[i+1].split()[0])-k
		dk = .03*k
		#damp = exp(-1.*k**2.*sigv**2.)
		damp = 1.
		rs = 0
		if md == 'RS':
			rs += float(frs[i].split()[1])
		sum += dk*k*k*s.P(k,0)*damp*(float(ln[1])+rs)
		#print k, sum
	return 2/pi*sum


def cl_limb(l,zf='TH0.40.60.05',zmed =1.06,ddz=.01,om=.279):
	sum = 0
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	nzf = open('nz'+zf+'.dat').readlines()
	dl = []
	rl = []
	Hl = []
	sumz = 0
	#sigv = 7.*d.D(zmed)
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom'+str(om)+'.dat').readlines()
	except:
		d = distance(om,1.-om)

		hf = open('hlistom'+str(om)+'.dat','w')
		hm = 1
	try:
		df = open('dlistom'+str(om)+'.dat').readlines()
	except:
		df = open('dlistom'+str(om)+'.dat','w')
		dm = 1
	try:
		rf = open('rlistom'+str(om)+'.dat').readlines()
	except:
		rf = open('rlistom'+str(om)+'.dat','w')
		rm = 1
	for i in range(0,len(nzf)):
		z = i*ddz
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if dm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if dm == 1:
			dz = d.Hz(z)
			Hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			Hl.append(float(hf[i]))

		sumz += float(nzf[i].split()[1])*ddz
	norm = 1./sumz
	pkl = []
	s = wp(om=om,lam=1.-om)
	s.Dzf = 1.
	s.Ps = s.Pspline()	
	sim = simulate(omega=om,lamda=1.-om,h=.7,h0=1.,ombaryon=.046,sig8=.8)
	try:
		pklf = open('pklfiles/pklfileom'+str(om)+str(l)+'.dat')
		for line in pklf:
			pkl.append(float(line))
	except:
		pklf = open('pklfiles/pklfileom'+str(om)+str(l)+'.dat','w')
		
		for i in range(1,len(nzf)):
			k = (l+.5)/rl[i]
			pk = s.Psplev(k)
			pkl.append(pk)
			pklf.write(str(pk)+'\n')
	for i in range(1,len(nzf)):
		k = (l+.5)/rl[i]
		#damp = exp(-1.*k**2.*sigv**2.)
		damp = 1.
		sum += ddz*(float(nzf[i].split()[1])*norm)**2.*dl[i]**2.*pkl[i-1]*damp*Hl[i]/rl[i]**2.	
	return sum/2997.9

def cl_limb_fnl(l,zf='TH0.40.60.05',zmed =1.06,ddz=.01,om=.279):
	sum = 0
	sum1 = 0
	sum2 = 0
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	nzf = open('nz'+zf+'.dat').readlines()
	dl = []
	rl = []
	Hl = []
	sumz = 0
	#sigv = 7.*d.D(zmed)
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom'+str(om)+'.dat').readlines()
	except:
		d = distance(om,1.-om)

		hf = open('hlistom'+str(om)+'.dat','w')
		hm = 1
	try:
		df = open('dlistom'+str(om)+'.dat').readlines()
	except:
		df = open('dlistom'+str(om)+'.dat','w')
		dm = 1
	try:
		rf = open('rlistom'+str(om)+'.dat').readlines()
	except:
		rf = open('rlistom'+str(om)+'.dat','w')
		rm = 1
	for i in range(0,len(nzf)):
		z = i*ddz
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if dm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if dm == 1:
			dz = d.Hz(z)
			Hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			Hl.append(float(hf[i]))

		sumz += float(nzf[i].split()[1])*ddz
	norm = 1./sumz
	pkl = []
	sim = simulate(omega=om,lamda=1.-om,h=.7,h0=1.,ombaryon=.046,sig8=.8)
	try:
		pklf = open('pklfiles/pklfileom'+str(om)+str(l)+'.dat')
		for line in pklf:
			pkl.append(float(line))
	except:
		pklf = open('pklfiles/pklfileom'+str(om)+str(l)+'.dat','w')
		s = wp(om=om,lam=1.-om)
		s.Dzf = 1.
		s.Ps = s.Pspline()	
		
		for i in range(1,len(nzf)):
			k = (l+.5)/rl[i]
			pk = s.Psplev(k)
			pkl.append(pk)
			pklf.write(str(pk)+'\n')
	for i in range(1,len(nzf)):
		k = (l+.5)/rl[i]
		#damp = exp(-1.*k**2.*sigv**2.)
		damp = 1.
		sum += ddz*(float(nzf[i].split()[1])*norm)**2.*dl[i]**2.*pkl[i-1]*damp*Hl[i]/rl[i]**2.
		if k < .1:
			sum1 += ddz*(float(nzf[i].split()[1])*norm)**2.*pkl[i-1]*damp*Hl[i]/rl[i]**2./k/k/sim.T(k)*dl[i]
			sum2 += ddz*(float(nzf[i].split()[1])*norm)**2.*pkl[i-1]*damp*Hl[i]/rl[i]**2./k/k/k/k/sim.T(k)**2.
	return sum/2997.9,sum1/2997.9,sum2/2997.9


def cl_limbcross_zl(l,zl1,zl2):
	sum = 0
	d = distance(.25,.75)
	dl = []
	rl = []
	Hl = []
	sumz = 0
	sumz2 = 0
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom.25.dat').readlines()
	except:
		d = distance(.25,.75)

		hf = open('hlistom.25.dat','w')
		hm = 1
	try:
		df = open('dlistom.25.dat').readlines()
	except:
		df = open('dlistom.25.dat','w')
		dm = 1
	try:
		rf = open('rlistom.25.dat').readlines()
	except:
		rf = open('rlistom.25.dat','w')
		rm = 1
	for i in range(0,len(zl1)):
		z = i*.001
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if dm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if dm == 1:
			dz = d.Hz(z)
			Hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			Hl.append(float(hf[i]))

		sumz += float(zl1[i])*.001
		sumz2 += float(zl2[i])*.001
	pkl = []
	try:
		pklf = open('pklfiles/pklfileom.25'+str(l)+'.dat')
		for line in pklf:
			pkl.append(float(line))
	except:
		pklf = open('pklfiles/pklfileom.25'+str(l)+'.dat','w')
		s = simulate(.25,.75)
		for i in range(1,len(zl1)):
			k = (l+.5)/rl[i]
			pk = s.P(k,0)
			pkl.append(pk)
			pklf.write(str(pk)+'\n')
	norm = 1./sumz
	norm2 = 1./sumz2
	for i in range(1,len(zl1)):
		k = (l+.5)/rl[i]
		damp = 1.
		sum += .001*(float(zl1[i])*norm)*(float(zl2[i])*norm2)*dl[i]**2.*pkl[i-1]*damp*Hl[i]/rl[i]**2.	
	return sum/d.c


def cl_limb_zl(l,zl):
	sum = 0
	dl = []
	rl = []
	Hl = []
	sumz = 0
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom.25.dat').readlines()
	except:
		d = distance(.25,.75)

		hf = open('hlistom.25.dat','w')
		hm = 1
	try:
		df = open('dlistom.25.dat').readlines()
	except:
		df = open('dlistom.25.dat','w')
		dm = 1
	try:
		rf = open('rlistom.25.dat').readlines()
	except:
		rf = open('rlistom.25.dat','w')
		rm = 1
	for i in range(0,len(zl)):
		z = i*.001
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if dm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if dm == 1:
			dz = d.Hz(z)
			Hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			Hl.append(float(hf[i]))

		sumz += float(zl[i])*.001
	norm = 1./sumz
	pkl = []
	try:
		pklf = open('pklfiles/pklfileom.25'+str(l)+'.dat')
		for line in pklf:
			pkl.append(float(line))
	except:
		pklf = open('pklfiles/pklfileom.25'+str(l)+'.dat','w')
		s = simulate(.25,.75)
		for i in range(1,len(zl)):
			k = (l+.5)/rl[i]
			pk = s.P(k,0)
			pkl.append(pk)
			pklf.write(str(pk)+'\n')
	for i in range(1,len(zl)):
		k = (l+.5)/rl[i]
		damp = 1.
		sum += .001*(float(zl[i])*norm)**2.*dl[i]**2.*pkl[i-1]*damp*Hl[i]/rl[i]**2.	
	return sum/2997.9


def cl_limbcross(l,zf1,zf2='TH0.40.60.05',zmed =1.06):
	sum = 0
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	nzf = open('nz'+zf1+'.dat').readlines()
	nzf2 = open('nz'+zf2+'.dat').readlines()
	d = distance(.25,.75)
	dl = []
	rl = []
	Hl = []
	sumz = 0
	sumz2 = 0
#	sigv = 7.*d.D(zmed)
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom.25.dat').readlines()
	except:
		d = distance(.25,.75)

		hf = open('hlistom.25.dat','w')
		hm = 1
	try:
		df = open('dlistom.25.dat').readlines()
	except:
		df = open('dlistom.25.dat','w')
		dm = 1
	try:
		rf = open('rlistom.25.dat').readlines()
	except:
		rf = open('rlistom.25.dat','w')
		rm = 1
	for i in range(0,len(nzf)):
		z = i*.001
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if dm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if dm == 1:
			dz = d.Hz(z)
			Hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			Hl.append(float(hf[i]))

		sumz += float(nzf[i].split()[1])*.001
		sumz2 += float(nzf2[i].split()[1])*.001
	pkl = []
	try:
		pklf = open('pklfiles/pklfileom.25'+str(l)+'.dat')
		for line in pklf:
			pkl.append(float(line))
	except:
		pklf = open('pklfiles/pklfileom.25'+str(l)+'.dat','w')
		s = simulate(.25,.75)
		for i in range(1,len(nzf)):
			k = (l+.5)/rl[i]
			pk = s.P(k,0)
			pkl.append(pk)
			pklf.write(str(pk)+'\n')
	norm = 1./sumz
	norm2 = 1./sumz2
	for i in range(1,len(nzf)):
		k = (l+.5)/rl[i]
		#damp = exp(-1.*k**2.*sigv**2.)
		damp = 1.
		sum += .001*(float(nzf[i].split()[1])*norm)*(float(nzf2[i].split()[1])*norm2)*dl[i]**2.*pkl[i-1]*damp*Hl[i]/rl[i]**2.	
	return sum/d.c



def mkphi(l,zf='TH0.40.60.05',dz=.01,om=.279):
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	nzf = open('nz'+zf+'.dat').readlines()
	d = distance(om,1.-om)
	fo = open('phiL/phi'+str(l)+str(om)+zf+'.dat','w')
	k = .0001
	dl = []
	rl = []
	hl = []
	sumz = 0
	for i in range(1,len(nzf)):
		sumz += float(nzf[i].split()[1])*dz
	print sumz
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom'+str(om)+'.dat').readlines()
	except:
		hf = open('hlistom'+str(om)+'.dat','w')
		hm = 1
	try:
		df = open('dlistom'+str(om)+'.dat').readlines()
	except:
		df = open('dlistom'+str(om)+'.dat','w')
		dm = 1
	try:
		rf = open('rlistom'+str(om)+'.dat').readlines()
	except:
		rf = open('rlistom'+str(om)+'.dat','w')
		rm = 1
	for i in range(0,len(nzf)):	
		z = i*dz
		#print z
		if dm == 1:
			ddz = d.D(z)
			dl.append(ddz)
			df.write(str(ddz)+'\n')
		else:
			dl.append(float(df[i]))
		if rm == 1:
			ddz = d.dc(z)
			rl.append(ddz)
			rf.write(str(ddz)+'\n')
		else:
			rl.append(float(rf[i]))
		if hm == 1:
			ddz = d.Hz(z)
			hl.append(ddz)
			hf.write(str(ddz)+'\n')
		else:
			hl.append(float(hf[i]))
	if dm == 1:
		df.close()
	if hm == 1:
		hf.close()
	if rm == 1:
		rf.close()
	cfac = d.c
	norm = 1./sumz
	#print norm
	nzint = 0
	for i in range(0,len(nzf)):
		zf = norm*float(nzf[i].split()[1])
		nzint += zf**2.*dz
	#print nzint
	sphjm = 0
	try:
		fsph = open('sphjfiles/sphjn'+str(l)+str(om)+'n4mk1.03.dat')
	except:
		fsph = open('sphjfiles/sphjn'+str(l)+str(om)+'n4mk1.03.dat','w')
		sphjm = 1
	while k < 10:
		sum = 0
		for i in range(1,len(nzf)):
			zf = norm*float(nzf[i].split()[1])
			#for j in range(0,10):
			#	ind = i*10+j
			if sphjm == 1:
				#print i,len(rl)
				sphjl = sph_jn(l,k*abs(rl[i]))
				
				fsph.write(str(sphjl)+'\n')
			else:
				sphjl = float(fsph.readline())
			sum += (zf*dl[i]*sphjl)*dz#**2.*.001*hl[i]/cfac#*.1
		fo.write(str(k)+' '+str(sum)+'\n')
		k = k*1.03
		#print k
	if sphjm == 1:
		fsph.close()
	
	fo.close()

def mkphi_D1(l,zf='TH0.40.60.05',dz=.01,om=.279):
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	nzf = open('nz'+zf+'.dat').readlines()
	d = distance(om,1.-om)
	fo = open('phiL/phi'+str(l)+str(om)+zf+'D1.dat','w')
	k = .0001
	dl = []
	rl = []
	hl = []
	sumz = 0
	for i in range(1,len(nzf)):
		sumz += float(nzf[i].split()[1])*dz
	print sumz
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom'+str(om)+'.dat').readlines()
	except:
		hf = open('hlistom'+str(om)+'.dat','w')
		hm = 1
	try:
		df = open('dlistom'+str(om)+'.dat').readlines()
	except:
		df = open('dlistom'+str(om)+'.dat','w')
		dm = 1
	try:
		rf = open('rlistom'+str(om)+'.dat').readlines()
	except:
		rf = open('rlistom'+str(om)+'.dat','w')
		rm = 1
	for i in range(0,len(nzf)):	
		z = i*dz
		#print z
		if dm == 1:
			ddz = d.D(z)
			dl.append(ddz)
			df.write(str(ddz)+'\n')
		else:
			dl.append(float(df[i]))
		if rm == 1:
			ddz = d.dc(z)
			rl.append(ddz)
			rf.write(str(ddz)+'\n')
		else:
			rl.append(float(rf[i]))
		if hm == 1:
			ddz = d.Hz(z)
			hl.append(ddz)
			hf.write(str(ddz)+'\n')
		else:
			hl.append(float(hf[i]))
	if dm == 1:
		df.close()
	if hm == 1:
		hf.close()
	if rm == 1:
		rf.close()
	cfac = d.c
	norm = 1./sumz
	#print norm
	nzint = 0
	for i in range(0,len(nzf)):
		zf = norm*float(nzf[i].split()[1])
		nzint += zf**2.*dz
	#print nzint
	sphjm = 0
	try:
		fsph = open('sphjfiles/sphjn'+str(l)+str(om)+'n4mk1.03.dat')
	except:
		fsph = open('sphjfiles/sphjn'+str(l)+str(om)+'n4mk1.03.dat','w')
		sphjm = 1
	while k < 10:
		sum = 0
		for i in range(1,len(nzf)):
			zf = norm*float(nzf[i].split()[1])
			#for j in range(0,10):
			#	ind = i*10+j
			if sphjm == 1:
				#print i,len(rl)
				sphjl = sph_jn(l,k*abs(rl[i]))
				
				fsph.write(str(sphjl)+'\n')
			else:
				sphjl = float(fsph.readline())
			sum += (zf*sphjl)*dz#**2.*.001*hl[i]/cfac#*.1
		fo.write(str(k)+' '+str(sum)+'\n')
		k = k*1.03
		#print k
	if sphjm == 1:
		fsph.close()
	
	fo.close()


def mkphil(l,z):
	d = distance(.25,.75)
	phi = []
	k = .001
	dl = []
	rl = []
	hl = []
	sumz = 0
	for i in range(1,zl):
		sumz += float(nzf[i].split()[1])*.001
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom.25.dat').readlines()
	except:
		hf = open('hlistom.25.dat','w')
		hm = 1
	try:
		df = open('dlistom.25.dat').readlines()
	except:
		df = open('dlistom.25.dat','w')
		dm = 1
	try:
		rf = open('rlistom.25.dat').readlines()
	except:
		rf = open('rlistom.25.dat','w')
		rm = 1
	for i in range(0,2000):
		z = i*.001
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if rm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if hm == 1:
			dz = d.Hz(z)
			hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			hl.append(float(hf[i]))
	cfac = d.c
	norm = 1./sumz
	nzint = 0
	for i in range(0,zl):
		zf = norm*zl[i]
		nzint += zf**2.*.001
	sphjm = 0
	try:
		fsph = open('sphjfiles/sphjn'+str(l)+'mk1.03.dat')
	except:
		fsph = open('sphjfiles/sphjn'+str(l)+'mk1.03.dat','w')
		sphjm = 1
	while k < 10:
		sum = 0
		for i in range(1,len(nzf)):
			zf = norm*float(nzf[i].split()[1])
			if sphjm == 1:
				sphjl = sph_jn(l,k*abs(rl[i]))
				fsph.write(str(sphjl)+'\n')
			else:
				sphjl = float(fsph.readline())
			sum += (zf*dl[i]*sphjl)*.001
		phil.append((k,sum))
		k = k*1.03
	return phil


def mkphinf(l,zf='TH0.40.60.05',sk=.001,dk=1.03):
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	nzf = open('nz'+zf+'.dat').readlines()
	d = distance(.25,.75)
	fo = open('phiL/phi'+str(l)+zf+'nf.dat','w')
	k = sk
	dl = []
	rl = []
	hl = []
	sumz = 0
	for i in range(1,len(nzf)):
		sumz += float(nzf[i].split()[1])*.001
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom.25.dat').readlines()
	except:
		hf = open('hlistom.25.dat','w')
		hm = 1
	try:
		df = open('dlistom.25.dat').readlines()
	except:
		df = open('dlistom.25.dat','w')
		dm = 1
	try:
		rf = open('rlistom.25.dat').readlines()
	except:
		rf = open('rlistom.25.dat','w')
		rm = 1
	for i in range(0,2000):
		z = i*.001
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if rm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if hm == 1:
			dz = d.Hz(z)
			hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			hl.append(float(hf[i]))
	cfac = d.c
	norm = 1./sumz
	#print norm
	nzint = 0
	for i in range(0,len(nzf)):
		zf = norm*float(nzf[i].split()[1])
		nzint += zf**2.*.001
	#print nzint
	sphjm = 1
	while k < 10:
		sum = 0
		for i in range(1,len(nzf)):
			zf = norm*float(nzf[i].split()[1])
			#for j in range(0,10):
			#	ind = i*10+j
			if sphjm == 1:
				sphjl = sph_jn(l,k*abs(rl[i]))
			sum += (zf*dl[i]*sphjl)*.001#**2.*.001*hl[i]/cfac#*.1
		fo.write(str(k)+' '+str(sum)+'\n')
		k = k*dk
	fo.close()


def pcdist(z1,z2,sigz=.03,binmin=.45,binmax=.55,dz=.01):
	from time import time
	t0 = time()
	zerr1 = sigz*(1.+z1)
	zerr1sq = zerr1**2.
	zerr2 = sigz*(1.+z2)
	zerr2sq = zerr2**2.
	sum = 0
	mult = 1./(2.*pi*zerr1*zerr2)
	#dz = .001
	indm = int(1/dz)
	binminind = int(binmin*1./dz)
	binmaxind = int(binmax*1./dz)
	zl1 = []
	zl2 = []
	f1m = 0
	try:
		f1 = open('zdist/zdist'+str(z1)+str(sigz)+str(dz)+'.dat').readlines()
	except:
		f1m = 1
		f1 = open('zdist/zdist'+str(z1)+str(sigz)+str(dz)+'.dat','w')
	if z1 != z2:
		f2m = 0
		try:
			f2 = open('zdist/zdist'+str(z2)+str(sigz)+str(dz)+'.dat').readlines()
		except:
			f2m = 1
			f2 = open('zdist/zdist'+str(z2)+str(sigz)+str(dz)+'.dat','w')
	for i in range(0,indm):
		if f1m == 1:
			p =exp(-1.*(z1-(i*dz+dz/2.))**2./(2.*zerr1sq))
			zl1.append(p)
			f1.write(str(p)+'\n')
		else:
			p = float(f1[i])
			zl1.append(p)
		if z1 != z2:
			if f2m == 1:
				p =exp(-1.*(z2-(i*dz+dz/2.))**2./(2.*zerr2sq))
				zl2.append(p)
				f2.write(str(p)+'\n')
			else:
				zl2.append(float(f2[i]))
		else:
			zl2.append(p)
	#print time()-t0
# 	thm = 0
# 	thl = []
# 	#f1.close()
# 	#f2.close()
# 	try:
# 		f = open('thetaf'+str(binmin)+str(binmax)+'.dat').readlines()
# 	except:
# 		thm = 1
# 		f = open('thetaf'+str(binmin)+str(binmax)+'.dat','w')
# 	if thm == 0:
# 		for i in range(0,1000):
# 			for j in range(0,1000):
# 				thl.append(int(f[j*1000+i]))
# 	else:
# 		for i in range(0,1000):
# 			for j in range(0,1000):
# 				binave = (i+j)/2
# 				if binave > binminind and binave < binmaxind:
# 					thl.append(1)
# 					f.write('1 \n')
# 				else:
# 					f.write('0 \n')
# 					thl.append(0)
	#f.close()
				
	#print time()-t0			
	for i in range(0,indm):
		for j in range(0,indm):
			binave = (i+j)/2
			if binave > binminind and binave < binmaxind:
				sum += zl1[i]*zl2[j]
			#sum += zl1[i]*zl2[j]*thl[j*1000+i]
	#print time()-t0
	return sum*mult*dz**2.
			
def mkpcdistfile(binmin,binmax,sigz=0.03):
	f = open('pcdist'+str(binmin)+str(binmax)+str(sigz)+'.dat','w')
	for i in range(0,100):
		print i
		for j in range(0,100):
			#print j
			f.write(str(pcdist(i/100.,j/100.,binmin=binmin,binmax=binmax))+'\n')
	f.close()

def mkphipc(l,zmin,zmax,zcmin=.4,zcmax=.6):
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	#pcf = open('pczdpc'+zf+'nzDES.dat').readlines()
	#dzf = open('zdpc'+zf+'nzDES.dat').readlines()
	from time import time
	t0 = time()
	nzf = open('nzDES.dat').readlines()
	norm1 = 0
	norm2 = 0
	normf = 0
	nzpc1 = []
	nzpc2 = []
	for i in range(0,len(nzf)):
		nzpc1.append(0)
		nzpc2.append(0)
	zcminind = int(1000*zcmin)
	zcmaxind = int(1000*zcmax)
	minind = int(1000*zmin)
	maxind = int(1000*zmax)
	for i in range(zcminind,zcmaxind):
		normf += float(nzf[i].split()[1])*.001
	zl = []
	for i in range(0,len(nzf)):
		zl.append(float(nzf[i].split()[1])/normf)
		
	print zcminind,maxind
	for i in range(zcminind,zcmaxind):
		sum = 0
		#zf = float(nzf[i].split()[1])/normf
		jmin = 2*minind-i
		if jmin < 0:
			jmin = 0
		jmax = 2*maxind-i
		if jmax < 0:
			jmax = 0
		for j in range(jmin,jmax):
			sum += zl[i]*.001*zl[j]
		nzpc1[i] = sum
	fo = open('nzpctemp.dat','w')
	for i in range(0,len(nzpc1)):
		fo.write(str(i/1000.)+' '+str(nzpc1[i])+'\n')
	print minind,zcmaxind
	for i in range(minind,zcmaxind):
		sum = 0
		zf = float(nzf[i].split()[1])/normf
		jmin = 2*minind-i
		if jmin < 0:
			jmin = 0
		jmax = 2*maxind-i
		if jmax < 0:
			jmax = 0
		for j in range(jmin,jmax):
			sum += zf*.001*float(nzf[j].split()[1])/normf
		nzpc2[i] = sum
	normpc1 = 0
	normpc2 = 0
	for i in range(zcminind,zcmaxind):
		normpc1 += .001*nzpc1[i]
		normpc2 += .001*nzpc2[i]
	t = 0
	ta = 0
	tb = 0
	pcdf = open('pcdist0.450.550.03.dat').readlines()
	
	for i in range(zcminind,zcmaxind):
		#print i
		for j in range(zcminind,zcmaxind):
			ta += sqrt(zl[i]*zl[j])*.001
			#if (i+j)/2 > minind and (i+j)/2 < maxind:
			pcdfind = (i/10)*100 +j/10
			#if i == 991:
			#	print i,j,i*10,j/10, pcdfind
			t += sqrt(zl[i]*zl[j])*.001*float(pcdf[pcdfind])
	fac = ta/t			
	print t,ta,normpc1,normpc2,normf
# 	for i in range(0,1000*zmax):
# 		if i > zcminind and i < zcmaxind:
# 			norm1 += float(nzf[i].split()[1])*.001
# 	for i in range(1000*zmin,len(nzf)):
# 		if i > zcminind and i < zcmaxind:
# 			norm2 += float(nzf[i].split()[1])*.001
	d = distance(.25,.75)
	fo = open('phiL/phisqpc'+str(l)+str(zmin)+str(zmax)+'zsig0.03DES.dat','w')
	k = .001
	dl = []
	rl = []
	hl = []
	dm = 0
	rm = 0
	hm = 0
	try:
		hf = open('hlistom.25.dat').readlines()
	except:
		hf = open('hlistom.25.dat','w')
		hm = 1
	try:
		df = open('dlistom.25.dat').readlines()
	except:
		df = open('dlistom.25.dat','w')
		dm = 1
	try:
		rf = open('rlistom.25.dat').readlines()
	except:
		rf = open('rlistom.25.dat','w')
		rm = 1
	for i in range(0,2000):
		z = i*.001
		if dm == 1:
			dz = d.D(z)
			dl.append(dz)
			df.write(str(dz)+'\n')
		else:
			dl.append(float(df[i]))
		if rm == 1:
			dz = d.dc(z)
			rl.append(dz)
			rf.write(str(dz)+'\n')
		else:
			rl.append(float(rf[i]))
		if hm == 1:
			dz = d.Hz(z)
			hl.append(dz)
			hf.write(str(dz)+'\n')
		else:
			hl.append(float(hf[i]))
	cfac = d.c
	nzint = 0
	phikl = []
	kl = []
	sphl  =[]
	sphjm = 0
	try:
		fsph = open('sphjfiles/sphjn'+str(l)+'mk1.03.dat')
	except:
		fsph = open('sphjfiles/sphjn'+str(l)+'mk1.03.dat','w')
		sphjm = 1
	for i in range(0,len(nzf)):
		sphl.append(0)
	print fac
	while k < 10:
		#print k
		sum = 0
		for i in range(1,len(nzf)):
			if sphjm == 1:
				sphjl = sph_jn(l,k*abs(rl[i]))
				fsph.write(str(sphjl)+'\n')
			else:
				sphjl = float(fsph.readline())
			#print i,len(nzf),len(sphl)
			sphl[i] = sphjl
		for i in range(zcminind,zcmaxind):
# 			jmin = 2*minind-i
# 			if jmin < 0:
# 				jmin = 0
# 			jmax = 2*maxind-i
# 			if jmax < 0:
# 				jmax = 0
# 			for j in range(jmin,jmax):			#print i
			for j in range(zcminind,zcmaxind):
			#	if (i+j)/2 >= minind and (i+j)/2 <= maxind:
					#zl1 = nzf[i].split()
					#zl2 = nzf[j].split()
					#zf1 = float(zl1[1])/normpc1/normf
					#zf2 = float(zl2[1])/normpc1/normf
				zf1 = zl[i]#/normpc1
				zf2 = zl[j]#/normpc1
				pcdfind = (i/10)*100 +j/10
				sum += (zf1*zf2*dl[i]*dl[j]*sphl[i]*sphl[j])*.001**2.*float(pcdf[pcdfind])
				#zf1 = nzpc1[i]/normpc1
				#zf2 = nzpc1[j]/normpc1
				#zf2 = nzpc2[j]/normpc2
# 				if i !=j:
# 					sum += 1.*(zf1*zf2*dl[i]*dl[j]*sphl[i]*sphl[j])*.001**2.*pcdist(i/1000.,j/1000.)
# 				else:
# 					sum += (zf1*zf2*dl[i]*dl[j]*sphl[i]*sphl[j])*.001**2.*pcdist(i/1000.,j/1000.)
		fo.write(str(k)+' '+str(sum*fac)+'\n')
		k = k*1.03
	fo.close()

def nzsqpc(zf='zsig0.030.450.55'):
	#nzf = open('nz'+zf+'nzDES.dat').readlines()
	pcf = open('pczdpc'+zf+'nzDES.dat').readlines()
	dzf = open('zdpc'+zf+'nzDES.dat').readlines()
	d = distance(.25,.75)
	fo = open('nzsqpc'+zf+'.dat','w')
	sumz = 0
	pcsum = 0
	for i in range(1,len(dzf)):
		sumz += float(dzf[i].split()[1])*.001
		if i < 200:
			pcsum += float(pcf[i].split()[1])*.01
	
	norm = 1./sumz
	normpc = 1./pcsum
	print norm,normpc
	nzl = []
	for i in range(0,2000):
		nzl.append(0)
	sum = 0
	
	for j in range(1,len(pcf)):
		#print j,pcf[j]
		pcl = pcf[j].split()
		z = float(pcl[0])
		ind = int(z*1000)
		pczf = float(pcl[1])
		if pczf > 0:
			pcfac = float(pcl[1])*.01
			for i in range(2,len(dzf),2):
				zf = norm*float(dzf[i].split()[1])
				if zf > 0:
					ind1 = ind-i/2
					ind2 = ind+i/2
					if ind1 > 0 and ind1 < 1995:
						nzl[ind1] += pcfac*zf*.001
					if ind2 >0 and ind2 < 1995:
						nzl[ind2] += pcfac*zf*.001
				#for j in range(0,10):
				#	ind = i*10+j
				#**2.*.001*hl[i]/cfac#*.1
	sum = 0
	for i in range(0,2000):
		sum += nzl[i]
	for i in range(0,2000):
		fo.write(str(i/1000.)+' '+str(nzl[i])+'\n')
	fo.close()

def mkphiRS(l,zf='TH0.40.60.05',zmed=.5,bias=2.,dz=.01,om=.279):
	nzf = open('nz'+zf+'.dat').readlines()
	d = distance(om,1.-om)
	beta = d.omz(zmed)**.557/bias
	print beta
	fo = open('phiL/phiRS'+str(l)+str(om)+zf+str(bias)+'.dat','w')
	k = .001
	dl = []
	rl = []
	hl = []
	sumz = 0
	for i in range(1,len(nzf)):
		sumz += float(nzf[i].split()[1])*dz

	for i in range(0,len(nzf)):
		z = i*dz
		dl.append(d.D(z))
		rl.append(d.dc(z))
		hl.append(d.Hz(z))
	cfac = d.c
	norm = 1./sumz
	print norm
	nzint = 0
	for i in range(0,len(nzf)):
		zf = norm*float(nzf[i].split()[1])
		nzint += zf**2.*dz
	print nzint
	ld = l - 2
	lu = l + 2
	
	lf = 0
	try:
		fsph = open('sphjfiles/sphjn'+str(l)+str(om)+'mk1.03.dat')
	except:
		fsph = open('sphjfiles/sphjn'+str(l)+str(om)+'mk1.03.dat','w')
		lf = 1
	lfd = 0
	try:
		fsphd = open('sphjfiles/sphjn'+str(ld)+str(om)+'mk1.03.dat')
	except:
		fsphd = open('sphjfiles/sphjn'+str(ld)+str(om)+'mk1.03.dat','w')
		lfd = 1
	lfu = 0
	try:
		fsphu = open('sphjfiles/sphjn'+str(lu)+str(om)+'mk1.03.dat')
	except:
		fsphu = open('sphjfiles/sphjn'+str(lu)+str(om)+'mk1.03.dat','w')
		lfu = 1
	l = float(l)
	while k < 10:
		sum = 0
		for i in range(1,len(nzf)):
			zf = norm*float(nzf[i].split()[1])
			#for j in range(0,10):
			#	ind = i*10+j
			lfac1 = (2.*l**2.+2.*l-1.)/((2.*l+3.)*(2.*l-1.))
			lfac2 = l*(l-1.)/((2.*l-1.)*(2.*l+1.))
			lfac3 = (l+1.)*(l+2.)/((2.*l+1.)*(2.*l+3.))
			if lf == 0:
				jlfac1 = lfac1*float(fsph.readline())
			else:
				sph = sph_jn(l,k*abs(rl[i]))
				jlfac1 = lfac1*sph
				fsph.write(str(sph)+'\n')
			if lfd == 0:
				jlfac2 = lfac2*float(fsphd.readline())
			else:
				sphd = sph_jn(ld,k*abs(rl[i]))
				jlfac2 = lfac2*sphd
				fsphd.write(str(sphd)+'\n')
			if lfu == 0:
				jlfac3 = lfac3*float(fsphu.readline())
			else:
				sphu = sph_jn(lu,k*abs(rl[i]))
				jlfac3 = lfac3*sphu
				fsphu.write(str(sphu)+'\n')
			#jlfac2 = lfac2*sph_jn(ld,k*abs(rl[i]))
			#jlfac3 = lfac3*sph_jn(lu,k*abs(rl[i]))
			sum += (zf*dl[i]*(jlfac1-jlfac2-jlfac3))*.001#**2.*.001*hl[i]/cfac#*.1
		fo.write(str(k)+' '+str(sum*beta)+'\n')
		k = k*1.03
	fo.close()

def om2cl(omfile,clguess='cltemp.dat'):

	f = open(clguess).readlines()	
	clg = []
	for i in range(0,len(f)):
		if i >= 200:
			break
		clg.append(float(f[i].split()[1]))
	fom = open(omfile).readlines()
	thetal = []
	oml = []
	for line in fom:
		thetal.append(float(line.split()[0])*pi/180.)
		oml.append(float(line.split()[1]))
	#for j in range(0,len(thetal)):
		for l in range(0,len(clg)):
			cl = oml[j]
			for i in range(0,len(clg)):
				if l != i:
					cl -= (2.*i+1.)/(4.*pi)*legendre(i,cos(thetal[j]))*clg[i]
			cl = cl*(4.*pi/(2.*l+1.))/legendre(l,cos(thetal[j]))
			#print l,cl
			clg[l] = cl
	return clg

def almfunc(omf,l):
	om = open(omf).readlines()
	bins = abs(float(om[0].split()[0])-float(om[1].split()[0]))*pi/180.
	sum = 0
	for j in range(0,10000):
		cth = -1. + .0002*j + .0001
		theta = acos(cth)
		if theta < .43*pi/180.:
			bin = 0
		else:
			bin = int((theta -.43*pi/180.)/bins)
		if bin >= len(om):
			bin = len(om)-1
			omb = 0
		else:
			ln = om[bin].split()
			omb = float(ln[1])
		if cth < .999:
			sum += omb*legendre(l,cth)*.0002
		else:
			for i in range(0,20):
				cth = cth + .00001
				if cth >= 1.:
				#	print j,i
					break
				theta = acos(cth)
				if theta < .5*pi/180. and theta > 5e-5:
					omb = 0.0107217869515*(theta/(.43*pi/180.))**(-.8)
				#print omb, legendre(l,cth)
				sum += omb*legendre(l,cth)*.00001
		#print bin, sum, omb
# 		if i == len(om) -1:
# 			nj = int((pi-theta)/bins)
# 			for j in range(1,nj):
# 				theta = theta + bins
# 				sum += omb*legendre(l,cos(theta))*sin(theta)*bins
	return sum*2.*pi
		
		

def cl2om_fnl(theta,zf='TH0.40.60.05',maxl=1000,md = 'n'):
	om = 0
	om1 = 0
	om2 = 0
	rad = theta*pi/180.
	f = open('cl'+zf+'_fnl.dat').readlines()
	oldpl = 1.
	cr = 0
	for l in range(1,maxl):
		ln = f[l].split()
		c = float(ln[1])
		c1 = float(ln[2])
		c2 = float(ln[3])
		pl = legendre(l,cos(rad))
		if pl < 0 and oldpl >0:
			cr += 1
			#print l
			#if cr == 7:
			#	break
			omp = om
			omp1 = om1
			omp2 = om2
		om += c*pl*(2.*l+1.)/(4.*pi)
		#if l < 60:
		om1 += c1*pl*(2.*l+1.)/(4.*pi)
		om2 += c2*pl*(2.*l+1.)/(4.*pi)
		oldpl = pl
		if l == 500:
			c5 = c
		if l == 999:
			c9 = c

	if l == 999:
		plin = (log(c9)-log(c5))/(log(999)-log(500))
		for l in range(1000,3000):
			pl = legendre(l,cos(rad))
			if pl < 0 and oldpl >0:
				cr += 1
				#print l
				#if cr == 7:
				#	break
				omp = om
			om += l**plin*c9/(999**plin)*pl*(2.*l+1.)/(4.*pi)
			oldpl = pl

		#print l,om
	
	return omp,omp1,omp2

def cl2om(theta,zf='TH0.40.60.05',maxl=1000,md = 'n',ll='n'):
	om = 0
	omrs = 0
	rad = theta*pi/180.
	f = open('cl'+zf+'nodamp.dat').readlines()
	if md == 'rs':
		try:
			frs = open('clrs'+zf+'.dat').readlines()
		except:
			print 'no z-dist file'
			md = 'n'
	oldpl = 1.
	cr = 0
	for l in range(0,maxl):
		#if l < 30:
		#	c = cl(l,zf)
		#else:
		#	c = cl_limb(l,zf)
		c = float(f[l].split()[1])
		if ll == 'y' and l < 5:
			c = c*4.
		if ll == 'y' and l > 5 and l < 10:
			c = c*1.5
		if md == 'n':
			crs = c
		else:
			if l < 30:
				crs = float(frs[l].split()[1])
			else:
				crs = c
		if l == 500:
			c1 = c
		if l == 999:
			c2 = c
		pl = legendre(l,cos(rad))
		if pl < 0 and oldpl >0:
			cr += 1
			#print l
			#if cr == 7:
			#	break
			omp = om
		om += c*pl*(2.*l+1.)/(4.*pi)
		omrs += crs*pl*(2.*l+1.)/(4.*pi)
		oldpl = pl
		#print l,om
	
	if l == 999:
		plin = (log(c2)-log(c1))/(log(999)-log(500))
		for l in range(1000,3000):
			pl = legendre(l,cos(rad))
			if pl < 0 and oldpl >0:
				cr += 1
				#print l
				#if cr == 7:
				#	break
				omp = om
			om += l**plin*c2/(999**plin)*pl*(2.*l+1.)/(4.*pi)
			omrs += l**plin*c2/(999**plin)*pl*(2.*l+1.)/(4.*pi)
			oldpl = pl
	if md == 'rs':		
		return omrs,om
	else:
		return omp


def cl2om_camb(theta,minl=2,maxl=2000,crx=5,nm='',md=''):
	om = 0
	omrs = 0
	rad = theta*pi/180.
	if md == 'me':
		f= open('cl'+nm+'.dat').readlines()
	else:
		f = open('/Users/ashleyross/CAMB_sources/counts_test_scalCls'+nm+'.dat').readlines()
		maxl = int(f[-1].split()[0])
	cr = 0
	oldpl = 1.
	omp = 0
	for l in range(2,maxl):
		if l < minl:
			c = float(f[0].split()[-2])*2.*pi/(1.*2*(2+1))
		else:
			c = float(f[l-minl].split()[-2])*2.*pi/(1.*l*(l+1.))
		pl = legendre(l,cos(rad))
		om += c*pl*(2.*l+1.)/(4.*pi)#*exp(-1.*(0.1*l/180.)**2.)
		if pl < 0 and oldpl >0:
		#if pl/oldpl < 0:
			cr += 1
			#print l
			omp = om
			if cr == crx:
				break
		
		oldpl = pl
	#omp = om	
# 	cu = float(f[-1].split()[-2])
# 	for l in range(1000,maxl):
# 		pl = legendre(l,cos(rad))
# 		c = cu*2.*pi/(1.*l*(l+1.))
# 		om += c*pl*(2.*l+1.)/(4.*pi)#*exp(-1.*(0.1*l/180.)**2.)
# 		if pl < 0 and oldpl >0:
# 			cr += 1
# 			omp = om
# 		
# 		oldpl = pl
	return omp

def cl2om_tom(theta,minl=2,maxl=1000,crx=5,nm='',md=''):
	om = 0
	omrs = 0
	rad = theta*pi/180.
	
	f = open('Cltomqso_100.dat').readlines()
	cr = 0
	for l in range(2,maxl):
		if l < minl:
			c = float(f[0].split()[-2])*2.*pi/(1.*2*(2+1))
		else:
			c = float(f[l-minl].split()[1])*2.*pi/(1.*l*(l+1.))
		pl = legendre(l,cos(rad))
		om += c*pl*(2.*l+1.)/(4.*pi)
		if pl < 0 and oldpl >0:
			cr += 1
			#print l
			omp = om
			#if cr == crx:
			#	break
		
		oldpl = pl
	#print '\n'
	return omp


def mkw_camb(nm,th0=19.5,sp=1.,nbin=20):

	fo = open('w'+nm+str(sp)+'.dat','w')
	for i in range(0,nbin):
		thi = th0-sp*i
		w = cl2om_camb(thi,crx=11,nm=nm)
		fo.write(str(thi)+' '+str(w)+'\n')
	fo.close()
	return True

def mkw_camblg(nm,z=.74,md='',om=.25,ob=.0459):
	from simulate import wp, simulate
	#ol = 1.-om
	#s = simulate(omega=om,lamda=ol,h=.7,h0=1.,ombaryon=ob,sig8=.8)
	#pw = wp(om=om,lam=ol,omb=ob)
	#pw.Dzf = pw.Dlin(z)	
	th = 54.5589364964
	#ngnorm = om*(pw.hf/pw.speedc)**2.*3.*1.686/pw.Dzf
	#print ngnorm
	fo = open('w'+nm+'lg.dat','w')
	for i in range(0,30):
		#w,w1,w2 = cl2om_fnl(th,nm)		
		#fo.write(str(th)+' '+str(w)+' '+str(2.*w1*ngnorm)+' '+str(w2*ngnorm**2.)+'\n')
		w = cl2om_camb(th,nm=nm)
		fo.write(str(th)+' '+str(w)+'\n')
		th = th*.826855
	fo.close()
	return True

def wxifile(nm,dz=0.05357,samp=1,min=.1,mult=1.2):
	#8.97164
	fo = open('wxi'+nm+'.dat','w')
	th = min
	while th < 10:
		w = calcwxi(nm,th,dz,samp)
		fo.write(str(th)+' '+str(w)+'\n')
		th = th*mult
		print th,w
	fo.close()
	return True

def mkwxifile(nm,dz=0.001,samp=1,max=8.97164,mult=0.794,min=.1,gam=-1.7):
	#8.97164
	fo = open('wxi'+nm+'.dat','w')
	th = max
	while th > min:
		w = calcwxi(nm,th,dz,samp,gam=gam)
		print th,w[0],w[1],w[2]
		fo.write(str(th)+' '+str(w[0])+' '+str(w[1])+' '+str(w[2])+'\n')
		th = th*mult
		
	fo.close()
	return True

def mkwxifile_thf(nm,dz=0.001,samp=1,thf='',min=.1,gam=-1.7):
	fo = open('wxi'+nm+'.dat','w')
	thf = load(thf+'.dat').transpose()#[0]
	for i in range(0,len(thf)):
		th = thf[i]
		w = calcwxi(nm,th,dz,samp,gam=gam)
		print th,w[0],w[1],w[2]
		fo.write(str(th)+' '+str(w[0])+' '+str(w[1])+' '+str(w[2])+'\n')
		
	fo.close()
	return True


def mkwxifilelinbin(nm,dz=0.001,samp=1,max=9,bs=.15,min=0,gam=-1.7,bao='',cosm='MICE',nzm='',md='three'):
	#8.97164
	fo = open('wxi'+nm+cosm+str(bs)+bao+md+'lb.dat','w')
	th = max-bs/2.
	while th > min:
		if md == 'one':
			w = calcwxi(nm,th,dz,samp,gam=gam,bao=bao,cosm=cosm,nzm=nzm)
			fo.write(str(th)+' '+str(w)+'\n')
		if md == 'three':	
		#print th,w[0],w[1],w[2]
			w = calcwxi3(nm,th,dz,samp,gam=gam,bao=bao,cosm=cosm,nzm=nzm)
			fo.write(str(th)+' '+str(w[0])+' '+str(w[1])+' '+str(w[2])+'\n')
		th = th-bs
		
	fo.close()
	return True


def mkwXxifile(nm,nm2,dz=0.001,samp=1,max=8.97164,mult=0.794):
	#8.97164
	fo = open('wxi'+nm+nm2+'.dat','w')
	th = max
	while th > .1:
		w = calcwXxi(nm,nm2,th,dz,samp)
		print th,w[0],w[1],w[2]
		fo.write(str(th)+' '+str(w[0])+' '+str(w[1])+' '+str(w[2])+'\n')
		th = th*mult
		
	fo.close()
	return True


def calcwxi3(nm,theta,dz=.001,samp=1,zm=.6,beta=.4,gam=-1.7,bao='',cosm='MICE',nzm='',dirz='/Users/ashleyross/DESY1/'):
	from Cosmo import distance
	from numpy import loadtxt
	#mf0 = open('/Users/ashleyross/DR12/xi0tkqpm0.413.03.57.015.0.dat').readlines()
	#mf2 = open('/Users/ashleyross/DR12/xi2tkqpm0.413.03.57.015.0.dat').readlines()
	#mf4 = open('/Users/ashleyross/DR12/xi4tkqpm0.413.03.57.015.0.dat').readlines()
	if cosm == 'MICE':
		xif = 'MICE_matterpower0.43.06.010.015.00'
		dir = ''
		dir = '/Users/ashleyross/LSSanalysis/BAOtemplates/'
	if cosm == 'Challenge':	
		xif = 'Challenge_matterpower0.406.010.015.00'
		dir = '/Users/ashleyross/LSSanalysis/BAOtemplates/'
	if bao == '':
		
		mf0 = open(dir+'xi0'+xif+'.dat').readlines()
		mf2 = open(dir+'xi2'+xif+'.dat').readlines()
		mf4 = open(dir+'xi4'+xif+'.dat').readlines()
	if bao == 'nobao':
		mf0 = open(dir+'xi0sm'+xif+'.dat').readlines()
		mf2 = open(dir+'xi2sm'+xif+'.dat').readlines()
		mf4 = open(dir+'xi4sm'+xif+'.dat').readlines()
	
	mf = []
	sump = 0
	sumpp = 0
	xi10 = float(mf0[10].split()[1])/(1.+2/3.*beta+beta**2./5.)
	for i in range(0,len(mf0)):
		xi0 = float(mf0[i].split()[1])/(1.+2/3.*beta+beta**2./5.)
		xi2 = float(mf2[i].split()[1])/(4./3.*beta+4./7.*beta**2.)
		xi4 = float(mf4[i].split()[1])/(beta**2.*8./35.)
		#if i > 10:
		#	sump += xi0*(10.+1.*i)**2.
		#	sumpp += xi0*(10.+1.*i)**4.
		#xi2 = xi0-3.*(10.+1.*i)**(-3.)*(sump+20.**3.*xi10/(3.+gam))		
		#xi4 = xi0+2.5*3.*(10.+1.*i)**(-3.)*(sump+20.**3.*xi10/(3.+gam))-3.5*5.*(10.+1.*i)**(-5.)*(sumpp+20.**5.*xi10/(5.+gam))
		mf.append((xi0,xi2,xi4))
	sumw = 0
	sumw2 = 0
	sumw4 = 0
	dl = []
	dl2 = []
	gl = []
	zl = []
	dir = ''
	if nzm == 'SM':
		dir = 'SMnz/'
		fz = loadtxt(dir+'MICE_WN_inter_V1.2.1_lightcone_corner000_photoz_mock0.zph'+nm+'_ztruedist').transpose()
	else:	
		fz = loadtxt(dirz+nm+'.dat').transpose()
	normz = sum(fz[1])*dz
	nz = fz[1]/normz
	#d = distance(.3,.7)
	if cosm == 'Challenge':
		d = distance(.31,.69)
	if cosm == 'MICE':
		d = distance(.25,.75)
	for i in range(0,len(fz[0]),samp):
		z = fz[0][i]
		#print z
		#if z > 0:
		zl.append(nz[i])		
		zd = d.dc(z)
		dl.append(zd)
		zd2 = d.dc(z)#+dz/2.)
		dl2.append(zd2)
		gl.append(d.D(z))
	#print sum(zl)*dz
	dang = d.dc(zm)*pi/180.*theta
	xid = wmod(dang,0,mf)
	#print dang,wmod(dang,0,mf)
	#frev = open('revcomp.dat','w')
	cthr = cos(theta*pi/180.)
	rmax = 250.**2.
	for i in range(1,len(zl),samp):
		sumi = 0
		sumi2 = 0
		sumi4 = 0
		if zl[i] > 0:
			for j in range(1,len(zl),samp):
				dzg = dz*gl[j]*zl[j]
				if i == j:
					zi = fz[0][i]
					sumk = 0
					sumk2 = 0
					sumk4 = 0
					reva = sqrt(d.dc(zi-dz/2.)**2.+dl[j]**2.-2.*dl[j]*d.dc(zi-dz/2.)*cthr)
					for k in range(0,10):
						zk = zi-.5*dz+.1*dz*k+.05*dz
						dk = d.dc(zk)
						revsq = dk**2.+dl[j]**2.-2.*dl[j]*dk*cthr
						if revsq > 1.e-4 and revsq < rmax:
							rev = sqrt(revsq)
							
							ra = dl[j]-dk
							mu = ra/rev
							xi = wmod3(rev,mu,mf,gam=gam)
							sumk += .1*dzg*xi[0]
							sumk2 += .1*dzg*xi[1]
							sumk4 += .1*dzg*xi[2]
					sumi += sumk
					sumi2 += sumk2
					sumi4 += sumk4

				else:	
					revsq = dl2[i]**2.+dl[j]**2.-2.*dl[j]*dl2[i]*cthr
					if revsq > 1.e-4:# and revsq < rmax:
						rev = sqrt(revsq)
						
						ra = dl[j]-dl[i]
						mu = ra/rev
						xi = wmod3(rev,mu,mf,gam=gam)
						sumi += dzg*xi[0]
						sumi2 += dzg*xi[1]
						sumi4 += dzg*xi[2]
			dzg = dz*gl[i]*zl[i]
			sumw += dzg*sumi
			sumw2 += dzg*sumi2
			sumw4 += dzg*sumi4
	#frev.close()
	return sumw,sumw2,sumw4

def calcwxi(nm,theta,dz=.001,samp=1,zm=.6,beta=.4,gam=-1.7,bao='',cosm='MICE',nzm=''):
	from Cosmo import distance
	from numpy import loadtxt
	#mf0 = open('/Users/ashleyross/DR12/xi0tkqpm0.413.03.57.015.0.dat').readlines()
	#mf2 = open('/Users/ashleyross/DR12/xi2tkqpm0.413.03.57.015.0.dat').readlines()
	#mf4 = open('/Users/ashleyross/DR12/xi4tkqpm0.413.03.57.015.0.dat').readlines()
	if cosm == 'MICE':
		xif = 'MICE_matterpower0.43.06.010.015.00'
		dir = ''
		dir = '/Users/ashleyross/LSSanalysis/BAOtemplates/'
	if cosm == 'BOSS':	
		xif = 'Challenge_matterpower0.563.04.07.015.00'
		dir = '/Users/ashleyross/LSSanalysis/BAOtemplates/'
	if bao == '':
		
		mf0 = open(dir+'xi0'+xif+'.dat').readlines()
		mf2 = open(dir+'xi2'+xif+'.dat').readlines()
		mf4 = open(dir+'xi4'+xif+'.dat').readlines()
	if bao == 'nobao':
		mf0 = open(dir+'xi0sm'+xif+'.dat').readlines()
		mf2 = open(dir+'xi2sm'+xif+'.dat').readlines()
		mf4 = open(dir+'xi4sm'+xif+'.dat').readlines()
	
	mf = []
	sump = 0
	sumpp = 0
	xi10 = float(mf0[10].split()[1])#/(1.+2/3.*beta+beta**2./5.)
	for i in range(0,len(mf0)):
		xi0 = float(mf0[i].split()[1])#/(1.+2/3.*beta+beta**2./5.)
		xi2 = float(mf2[i].split()[1])#/(4./3.*beta+4./7.*beta**2.)
		xi4 = float(mf4[i].split()[1])#/(beta**2.*8./35.)
		#if i > 10:
		#	sump += xi0*(10.+1.*i)**2.
		#	sumpp += xi0*(10.+1.*i)**4.
		#xi2 = xi0-3.*(10.+1.*i)**(-3.)*(sump+20.**3.*xi10/(3.+gam))		
		#xi4 = xi0+2.5*3.*(10.+1.*i)**(-3.)*(sump+20.**3.*xi10/(3.+gam))-3.5*5.*(10.+1.*i)**(-5.)*(sumpp+20.**5.*xi10/(5.+gam))
		mf.append((xi0,xi2,xi4))
	sumw = 0
	sumw2 = 0
	sumw4 = 0
	dl = []
	dl2 = []
	gl = []
	zl = []
	dir = ''
	if nzm == 'SM':
		dir = 'SMnz/'
		fz = loadtxt(dir+'MICE_WN_inter_V1.2.1_lightcone_corner000_photoz_mock0.zph'+nm+'_ztruedist').transpose()
	else:	
		fz = loadtxt(nm+'.dat').transpose()
	normz = sum(fz[1])*dz
	nz = fz[1]/normz
	#d = distance(.3,.7)
	if cosm == 'BOSS':
		d = distance(.31,.69)
	if cosm == 'MICE':
		d = distance(.25,.75)
	for i in range(0,len(fz[0]),samp):
		z = fz[0][i]
		#print z
		#if z > 0:
		zl.append(nz[i])		
		zd = d.dc(z)
		dl.append(zd)
		zd2 = d.dc(z)#+dz/2.)
		dl2.append(zd2)
		gl.append(d.D(z))
	#print sum(zl)*dz
	dang = d.dc(zm)*pi/180.*theta
	xid = wmod(dang,0,mf)
	#print dang,wmod(dang,0,mf)
	#frev = open('revcomp.dat','w')
	cthr = cos(theta*pi/180.)
	rmax = 250.**2.
	for i in range(1,len(zl),samp):
		sumi = 0
		if zl[i] > 0:
			for j in range(1,len(zl),samp):
				dzgj = dz*gl[j]*zl[j]
				if i == j:
					zi = fz[0][i]
					sumk = 0
					reva = sqrt(d.dc(zi-dz/2.)**2.+dl[j]**2.-2.*dl[j]*d.dc(zi-dz/2.)*cthr)
					for k in range(0,10):
						zk = zi-.5*dz+.1*dz*k+.05*dz
						dk = d.dc(zk)
						revsq = dk**2.+dl[j]**2.-2.*dl[j]*dk*cthr
						if revsq > 1.e-4 and revsq < rmax:
							rev = sqrt(revsq)
							ra = dl[j]-dk
							mu = ra/rev
							xi = wmod(rev,mu,mf,gam=gam)
							sumk += .1*dzgj*xi
					sumi += sumk

				else:	
					revsq = dl2[i]**2.+dl[j]**2.-2.*dl[j]*dl2[i]*cthr
					if revsq > 1.e-4 and revsq < rmax:
						rev = sqrt(revsq)
						ra = dl[j]-dl2[i]
						mu = ra/rev
						xi = wmod(rev,mu,mf,gam=gam)
						sumi += dzgj*xi
			dzg = dz*gl[i]*zl[i]
			sumw += dzg*sumi
	#frev.close()
	return sumw


def calcwXxi(nm,nm2,theta,dz=.001,samp=1,zm=.6,beta=.4,gam=-1.7):
	from Cosmo import distance
	#mf0 = open('/Users/ashleyross/DR12/xi0tkqpm0.413.03.57.015.0.dat').readlines()
	#mf2 = open('/Users/ashleyross/DR12/xi2tkqpm0.413.03.57.015.0.dat').readlines()
	#mf4 = open('/Users/ashleyross/DR12/xi4tkqpm0.413.03.57.015.0.dat').readlines()
	mf0 = open('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.04.85.dat').readlines()
	mf2 = open('/Users/ashleyross/DR12/xi2Challenge_matterpower0.43.04.85.dat').readlines()
	mf4 = open('/Users/ashleyross/DR12/xi4Challenge_matterpower0.43.04.85.dat').readlines()

	mf = []
	sump = 0
	sumpp = 0
	for i in range(0,len(mf0)):
		xi0 = float(mf0[i].split()[1])/(1.+2/3.*beta+beta**2./5.)
		if i > 10:
			sump += xi0*(10.+1.*i)**2.
			sumpp += xi0*(10.+1.*i)**4.
		#xi2
		mf.append((float(mf0[i].split()[1])/(1.+2/3.*beta+beta**2./5.),float(mf2[i].split()[1])/(4/3.*beta+4.*beta**2./7.),float(mf4[i].split()[1])/(8.*beta**2./35.)))
	sumw = 0
	sumw2 = 0
	sumw4 = 0
	dl = []
	dl2 = []
	gl = []
	zl = []
	zl2 = []
	fz = open(nm+'.dat').readlines()
	fz2 = open(nm2+'.dat').readlines()
	if len(fz) != len(fz2):
		return 'need to write something for differing n(z) formats'
	d = distance(.3,.7)
	for i in range(0,len(fz),samp):
		z = float(float(fz[i].split()[0]))
		#if z > 0:
		zl.append(float(fz[i].split()[1]))
		zl2.append(float(fz2[i].split()[1]))		
		zd = d.dc(z)
		dl.append(zd)
		zd2 = d.dc(z)#+dz/2.)
		dl2.append(zd2)
		gl.append(d.D(z))
	#print sum(zl)*dz
	dang = d.dc(zm)*pi/180.*theta
	xid = wmod(dang,0,mf)
	#print dang,wmod(dang,0,mf)
	#frev = open('revcomp.dat','w')
	cthr = cos(theta*pi/180.)
	rmax = 250.**2.
	for i in range(1,len(zl),samp):
		sumi = 0
		sumi2 = 0
		sumi4 = 0
		if zl[i] > 0:
			for j in range(1,len(zl),samp):
				dzg = dz*gl[j]*zl2[j]
				if i == j:
					zi = float(fz[i].split()[0])
					sumk = 0
					sumk2 = 0
					sumk4 = 0
					#sumka = 0
					reva = sqrt(d.dc(zi-dz/2.)**2.+dl[j]**2.-2.*dl[j]*d.dc(zi-dz/2.)*cthr)
					#print reva,d.dc(zi-dz/2.),dl[j]
					#revb = sqrt(d.dc(zi+dz/2.)**2.+dl[j]**2.-2.*dl[j]*d.dc(zi+dz/2.)*cthr)
					#print revb,d.dc(zi+dz/2.),dl[j]
					#revc = sqrt(dl2[i]**2.+dl[j]**2.-2.*dl[j]*dl2[i]*cthr)
					#sumka = dz*zl[j]*wmod(reva,0,mf)
					for k in range(0,10):
						zk = zi-.5*dz+.1*dz*k+.05*dz
						#print zk,zi,i,zl[i]
						dk = d.dc(zk)
						revsq = dk**2.+dl[j]**2.-2.*dl[j]*dk*cthr
						if revsq > 1.e-4 and revsq < rmax:
							rev = sqrt(revsq)
							#print rev
							
							ra = dl[j]-dk
							mu = ra/rev
							#print mu,dm,dk,ra,rev/2.
							#mu = 0
							xi = wmod3(rev,mu,mf)
							#if xi > 10:
							#	print rev,xi
							sumk += .1*dzg*xi[0]
							sumk2 += .1*dzg*xi[1]
							sumk4 += .1*dzg*xi[2]
					sumi += sumk
					sumi2 += sumk2
					sumi4 += sumk4
					#print sumk,sumka,reva,revb,revc		

				else:	
					revsq = dl2[i]**2.+dl[j]**2.-2.*dl[j]*dl2[i]*cthr
					#print revsq
					#revsq = dang**2.+(dl[i]-dl[j])**2.
					if revsq > 1.e-4 and revsq < rmax:
						rev = sqrt(revsq)
						#print rev
						
						ra = dl[j]-dl[i]
						mu = ra/rev
						#mu = 0
						xi = wmod3(rev,mu,mf)
						#if xi > 10:
						#	print rev,xi
						sumi += dz*zl2[j]*xi[0]
						sumi2 += dz*zl2[j]*xi[1]
						sumi4 += dz*zl2[j]*xi[2]
					#frev.write(str(i*dz+dz/2.)+' '+str(j*dz+dz/2.)+' '+str(rev)+' '+str(sqrt(revsq2))+'\n')
				#else:
				#	print zl[i],zl[j]	
			#print sumi,zl[i]
			dzg = dz*gl[i]*zl[i]
			sumw += dzg*sumi
			sumw2 += dzg*sumi2
			sumw4 += dzg*sumi4
	#frev.close()
	return sumw,sumw2,sumw4


def projxi(theta,dr=.01,zm=.6):
	from Cosmo import distance
	mf0 = open('/Users/ashleyross/DR12/xi0tkqpm0.413.03.57.015.0.dat').readlines()
	mf2 = open('/Users/ashleyross/DR12/xi2tkqpm0.413.03.57.015.0.dat').readlines()
	mf4 = open('/Users/ashleyross/DR12/xi4tkqpm0.413.03.57.015.0.dat').readlines()
	rmax = 300
	nr = int(rmax/dr)
	mf = []
	for i in range(0,len(mf0)):
		mf.append((float(mf0[i].split()[1]),float(mf2[i].split()[1]),float(mf4[i].split()[1])))
	sumw = 0
	dl = []
	zl = []
	d = distance(.3,.7)
	dang = d.dc(zm)*pi/180.*theta
	xid = wmod(dang,0,mf)
	print dang,wmod(dang,0,mf)
	#frev = open('revcomp.dat','w')
	cthr = cos(theta*pi/180.)
	for i in range(0,nr):
		r = dr/2.+i*dr
		rth = sqrt(dang**2.+r**2.)
		xi = wmod(rth,0,mf)
					#if xi > 10:
					#	print rev,xi
		sumw += dr*xi

	return sumw


def wmod(r,mu,mf,sp=1.,gam=-2.,beta=.4):
	p2 = legendre(2,mu)
	p4 = legendre(4,mu)
	xi10 = mf[10][0]
	if r < 21.:
		#xi10 = mf[0][0]+p2*mf[0][1]+p4*mf[0][2]
		xi0 = xi10*(r/20.)**(gam)
		xi2 = xi0-3.*xi0/(3.+gam)
		xi4 = xi0+2.5*3.*xi0/(3.+gam)-3.5*5.*xi0/(5.+gam)

		return xi0+p2*xi2+p4*xi4
	else:
		indd = int((r-10.)/sp)
		indu = indd + 1
		fac = (r-10.)/sp-indd
		if fac > 1.:
			print 'BAD FAC in wmod'
			return 'ERROR'
		if indu >= len(mf)-1:
			return 0
		xi0 =mf[indu][0]*fac+(1.-fac)*mf[indd][0] 
		xi2 = (mf[indu][1]*fac+(1.-fac)*mf[indd][1])
		xi4 = (mf[indu][2]*fac+(1.-fac)*mf[indd][2])

		return xi0+p2*xi2+p4*xi4

def wmodtest(r,beta=.4,gam=-1.7):
	mf0 = open('/Users/ashleyross/DR12/xi0Challenge_matterpower0.43.04.85.dat').readlines()

	mf = []
	sump = 0
	sumpp = 0
	xi10 = float(mf0[10].split()[1])/(1.+2/3.*beta+beta**2./5.)
	for i in range(0,len(mf0)):
		xi0 = float(mf0[i].split()[1])/(1.+2/3.*beta+beta**2./5.)
		if i > 10:
			sump += xi0*(10.+1.*i)**2.
			sumpp += xi0*(10.+1.*i)**4.
			
		xi2 = xi0-3.*(10.+1.*i)**(-3.)*(sump+20.**3.*xi10/(3.+gam))		
		if i == 11:
			print sump,xi2,20.**3.*xi10/(3.+gam),xi0
		xi4 = xi0+2.5*3.*(10.+1.*i)**(-3.)*(sump+20.**3.*xi10/(3.+gam))-3.5*5.*(10.+1.*i)**(-5.)*(sumpp+20.**5.*xi10/(5.+gam))
		mf.append((xi0,xi2,xi4))
	w = wmod3(r,0,mf)
	return w[0],w[1]/legendre(2,0),w[2]/legendre(4,0)

def wmod3(r,mu,mf,sp=1.,gam=-1.7):
	p2 = legendre(2,mu)
	p4 = legendre(4,mu)
	#p2 = 0
	#p4 = 0
	xi10 = mf[10][0]#,mf[10][1],mf[10][2]
	if r < 21.:
		#xi10 = mf[0][0]+p2*mf[0][1]+p4*mf[0][2]		
		xi0 = xi10*(r/20.)**(gam)
		#xi2 = -xi0/2./p2
		#xi4 = 0
		#xi2 = xi10[1]*(r/20.)**(gam)
		#xi4 = xi10[2]*(r/20.)**(gam)
		#if r < 10.:
		#	xi2 = mf[0][1]/(xi10[0]*(10/20.)**(gam))*(r/20.)**(gam)
		#	xi4 = mf[0][2]/(xi10[0]*(10/20.)**(gam))*(r/20.)**(gam)
		#else:
		#	indd = int((r-10.)/sp)
		#	xi2 = mf[indd][1]
		#	xi4 = mf[indd][2]
		xi2 = xi0-3.*xi0/(3.+gam)
		xi4 = xi0+2.5*3.*xi0/(3.+gam)-3.5*5.*xi0/(5.+gam)
		return xi0,p2*xi2,p4*xi4
		#return xi0,0,0
	else:
		indd = int((r-10.)/sp)
		indu = indd + 1
		fac = (r-10.)/sp-indd
		if fac > 1.:
			print 'BAD FAC in wmod'
			return 'ERROR'
		if indu >= len(mf)-1:
			return mf[-1][0],p2*mf[-1][1],p4*mf[-1][2]
		#return (mf[indu][0]+p2*mf[indu][1]+mf[indu][2]*p4)*fac+(1.-fac)*(mf[indu][0]+p2*mf[indu][1]+mf[indu][2]*p4)	
		xi0 =mf[indu][0]*fac+(1.-fac)*mf[indd][0] 
		xi2 = (mf[indu][1]*fac+(1.-fac)*mf[indd][1])
		xi4 = (mf[indu][2]*fac+(1.-fac)*mf[indd][2])
		#return xi0,p2*(mf[indu][1]*fac+(1.-fac)*mf[indd][1]),p4*(mf[indu][2]*fac+(1.-fac)*mf[indd][2])
		return xi0,p2*xi2,p4*xi4
		
def xip(mf,r):
	indmax = int((r-10)/1.)
	x = 0
	for i in range(10,indmax):
		rp = 20.+i*1.
		x += 1.*mf[i][0]*rp**2.
	return x

def xipp(mf,r):
	indmax = int((r-10)/1.)
	x = 0
	for i in range(10,indmax):
		rp = 20.+i*1.
		x += 1.*mf[i][0]*rp**4.
	return x

	
def mkw_lg(nm,z=.74,md='',om=.25,ob=.0459):
	fo = open('w'+nm+'melg.dat','w')
	th = 8.97164
	for i in range(0,20):
		w = cl2om_camb(th,nm=nm)
		fo.write(str(th)+' '+str(w)+'\n')
		th = th*0.7943
	fo.close()
	return True

def mkw_lgfnl(nm,z=.74,md='',om=.25,ob=.0459):
	from simulate import wp, simulate
	ol = 1.-om
	#s = simulate(omega=om,lamda=ol,h=.7,h0=1.,ombaryon=ob,sig8=.8)
	pw = wp(om=om,lam=ol,omb=ob)
	#pw.Dzf = pw.Dlin(z)	
	th = 54.5589364964
	#ngnorm = om*(pw.hf/pw.speedc)**2.*3.*1.686#/pw.Dzf
	ngnorm = om*(1./pw.speedc)**2.*3.*1.686
	print ngnorm
	#print pw.Dzf
	#ngnorm = ngnorm/pw.Dzf
	#print ngnorm
	fo = open('w'+nm+'mefnllgl1.dat','w')
	th = 54.5589364964
	for i in range(0,30):
		w,w1,w2 = cl2om_fnl(th,nm)
		#w = cl2om_tom(th)
		#fo.write(str(th)+' '+str(w)+'\n')
		fo.write(str(th)+' '+str(w)+' '+str(2.*w1*ngnorm)+' '+str(w2*ngnorm**2.)+'\n')
		th = th*.826855
	fo.close()
	return True


def mkw_cambsp(nm,sp=.01,maxl=2000):
	thi = 15.
	fo = open('w'+nm+str(sp)+'.dat','w')
	while thi > .2:
		w = cl2om_camb(thi,crx=7,nm=nm,maxl=maxl)
		fo.write(str(thi)+' '+str(w)+'\n')
		thi -= sp
		#if int(thi/1.-.01) == thi/1.:
		print thi
	fo.close()
	return True

	

def mkcov_camb(nm,ngal,bias=1.,thmax=12,fs=.035,bins=.15,nbin=40,wa='a'):
	th0 = thmax-bins/2.
	for i in range(0,nbin):
		thi = th0-bins*i
		cth1 = cos(thi*pi/180.)
		fo = open('Pl'+str(thi)+'.dat','w')
		for j in range(0,1000):
			p = legendre(j,cth1)
			fo.write(str(p)+'\n')
		fo.close()
	print 'made legendre files'
	fo = open('cov'+nm+wa+str(bins)+'.dat','w')
	fd = open('diag'+nm+wa+str(bins)+'.dat','w')
	for i in range(0,nbin):
		thi = th0-bins*i
		for j in range(0,nbin):			
			thj = th0-bins*j
			cov = cl2cov_camb(thi,thj,nm,ngal=ngal,fs=fs,bins=bins,bias=bias)
			fo.write(str(cov)+' ')
			if i == j:
				fd.write(str(thi)+' '+str(sqrt(cov))+'\n')
		fo.write('\n')
		#print i
	fo.close()
	return True

def mkcov_cambminmax(nm,ngal,bias=1.,thmin=0.075,fs=.035,bins=.15,nbin=40,wa='a'):
	th0 = thmin
# 	for i in range(0,nbin):
# 		thi = th0-bins*i
# 		cth1 = cos(thi*pi/180.)
# 		fo = open('Pl'+str(thi)+'.dat','w')
# 		for j in range(0,1000):
# 			p = legendre(j,cth1)
# 			fo.write(str(p)+'\n')
# 		fo.close()
	print 'made legendre files'
	fo = open('cov'+nm+wa+str(bins)+'.dat','w')
	fd = open('diag'+nm+wa+str(bins)+'.dat','w')
	for i in range(0,nbin):
		thi = th0+bins*i
		for j in range(0,nbin):			
			thj = th0+bins*j
			cov = cl2cov_camb(thi,thj,nm,ngal=ngal,fs=fs,bins=bins,bias=bias)
			fo.write(str(cov)+' ')
			if i == j:
				fd.write(str(thi)+' '+str(sqrt(cov))+'\n')
		fo.write('\n')
		#print i
	fo.close()
	return True


def mkcov_cambNzbin(nmlist,ngallist,biaslist,cfm,thmax=6,fs=.035,bins=.15,nbin=33,wa='a'):
	th0 = thmax-bins/2.
	for i in range(0,nbin):
		thi = th0-bins*i
		cth1 = cos(thi*pi/180.)
		fo = open('Pl'+str(thi)+'.dat','w')
		for j in range(0,1000):
			p = legendre(j,cth1)
			fo.write(str(p)+'\n')
		fo.close()
	print 'made legendre files'
	fo = open('covall'+wa+str(bins)+'.dat','w')
	for k in range(0,len(nmlist)):
		for i in range(0,nbin):
			thi = th0-bins*i
			for l in range(0,len(nmlist)):
				for j in range(0,nbin):			
					thj = th0-bins*j
					if k == l:
						cov = cl2cov_camb(thi,thj,nmlist[k],ngal=ngallist[k],fs=fs,bins=bins,bias = biaslist[k])
					else:
						cov = cl2cov_camb2bin(thi,thj,cfm[k][l],nm=nmlist[k],nm2=nmlist[l],fs=fs,bias=biaslist[k],bias2=biaslist[l])
					fo.write(str(cov)+' ')
			fo.write('\n')
		#print i
	fo.close()
	return True


def mkcov_camblg(nm,ngal,bias=1.,thmax=8.97164,nbin=20,fs=.0275):
	#16650764
	thi = thmax
	m = 0.7943
	for i in range(0,nbin):		
		cth1 = cos(thi*pi/180.)
		fo = open('Pl'+str(thi)+'.dat','w')
		for j in range(0,1000):
			p = legendre(j,cth1)
			fo.write(str(p)+'\n')
		thi = thi*m
		fo.close()
	print 'made legendre files'
	fo = open('cov'+nm+'lg.dat','w')
	fd = open('diag'+nm+'lg.dat','w')
	thi = thmax
	
	
	bins = thmax-thmax*m
	for i in range(0,nbin):
		thj = thmax
		for j in range(0,nbin):						
			cov = cl2cov_camb(thi,thj,nm,ngal=ngal,fs=fs,bins=bins,bias=bias)
			fo.write(str(cov)+' ')
			if i == j:
				fd.write(str(thi)+' '+str(sqrt(cov))+'\n')
			thj = thj*m
		oldi = thi
		thi = thi*m
		bins = oldi-thi
		fo.write('\n')
		print i
	fo.close()
	return True

def mkcov_lg(nm,ngal,bias,thmax=54.5589364964,nbin=20,fs=.125):
	thi = thmax
	m = 0.7943
	for i in range(0,nbin):		
		cth1 = cos(thi*pi/180.)
		fo = open('Pl'+str(thi)+'.dat','w')
		for j in range(0,1000):
			p = legendre(j,cth1)
			fo.write(str(p)+'\n')
		thi = thi*m
		fo.close()
	print 'made legendre files'
	fo = open('cov'+nm+'melg.dat','w')
	fd = open('diag'+nm+'melg.dat','w')
	thi = thmax
	
	
	bins = thmax-thmax*m
	for i in range(0,nbin):
		thj = thmax
		for j in range(0,nbin):						
			cov = cl2cov_d(thi,thj,nm,ngal=ngal,fs=fs,bins=bins,bias=bias)
			fo.write(str(cov)+' ')
			if i == j:
				fd.write(str(thi)+' '+str(sqrt(cov))+'\n')
			thj = thj*m
		oldi = thi
		thi = thi*m
		bins = oldi-thi
		fo.write('\n')
		print i
	fo.close()
	return True


def mkcov_cambrp(nm='BOSS0.450.5wall',ngal=863000.,nb=40):
	th0 = .163
	for i in range(0,nb):
		thi = th0+0.2452*i
		cth1 = cos(thi*pi/180.)		
		fo = open('Pl'+str(thi)+'.dat','w')
		for j in range(0,1000):
			p = legendre(j,cth1)
			fo.write(str(p)+'\n')
		fo.close()
	print 'made legendre files'
	fo = open('covBOSSrp.dat','w')
	fj = open('diagBOSSrp.dat','w')
	for i in range(0,nb):
		thi = th0+0.2452*i
		
		for j in range(0,nb):				
			thj = th0+0.2452*j
			cov = cl2cov_camb(thi,thj,nm,ngal=ngal/4.)/4.
			fo.write(str(cov)+' ')
			if i == j:
				fj.write(str(thj)+' '+str(sqrt(cov))+'\n')
		fo.write('\n')
		print i
	fo.close()
	fj.close()
	return True


def cl2cov_camb_Nsn(theta1,theta2,nm='',fs=.25,ngal=5000000.,minl=2,maxl=1000,bins=.2,bias=1.):
	f = open('/Users/ashleyross/CAMB_sources/counts_test_scalCls'+nm+'.dat').readlines()
	sum = 0
	innp = 0
	if theta1 == theta2:
		np = 2.*pi*theta1*(pi/180.)**2.*bins*.5*ngal*ngal/(fs*4.*pi)
		innp = 1./np	
	nbar = ngal/(4.*pi*fs)
	cth1 = cos(theta1*pi/180.)
	lth1 = open('Pl'+str(theta1)+'.dat').readlines()
	cth2 = cos(theta2*pi/180.)
	lth2 = open('Pl'+str(theta2)+'.dat').readlines()
	sumcl2 = 0
	sumcln = 0
	sumnb2 = 0
	for l in range(minl,maxl):
		if l < minl:
			c = float(f[0].split()[-2])*2.*pi/(1.*2*(2+1))*bias**2.
		else:
			c = float(f[l-minl].split()[-2])*2.*pi/(1.*l*(l+1.))*bias**2.
		pl1 = float(lth1[l])
		pl2 = float(lth2[l])
		errl = (2.*l+1.)*pl1*pl2*(c**2.)#+2.*c/nbar)#+1./nbar/nbar)
		sum += errl
	return 2./fs*sum/(4.*pi)**2.#+innp

def cl2cov_camb(theta1,theta2,nm='',fs=.25,ngal=5000000.,minl=2,maxl=1000,bins=.15,bias=1.):
	f = open('/Users/ashleyross/CAMB_sources/counts_test_scalCls'+nm+'.dat').readlines()
	sum = 0
	innp = 0
	nbar = ngal/(4.*pi*fs)
	Abin = 2.*pi*theta1*pi/180.*bins*pi/180.
	if theta1 == theta2:
		#np = 2.*pi*theta1*(pi/180.)**2.*bins*.5*ngal*ngal/(fs*4.*pi)
		np = 0.5*ngal*Abin*nbar
		
		innp = 1./np	
	
	cth1 = cos(theta1*pi/180.)
	lth1 = open('Pl'+str(theta1)+'.dat').readlines()
	cth2 = cos(theta2*pi/180.)
	lth2 = open('Pl'+str(theta2)+'.dat').readlines()
	sumcl2 = 0
	sumcln = 0
	sumnb2 = 0
	for l in range(minl,maxl):
		if l < minl:
			c = float(f[0].split()[-2])*2.*pi/(1.*2*(2+1))*bias**2.
		else:
			c = float(f[l-minl].split()[-2])*2.*pi/(1.*l*(l+1.))*bias**2.
		pl1 = float(lth1[l])
		pl2 = float(lth2[l])
		errl = (2.*l+1.)*pl1*pl2*(c**2.+2.*c/nbar+1./nbar/nbar)
		sum += errl
	return 2./fs*sum/(4.*pi)**2.+innp

def cl2cov_camb2bin(theta1,theta2,cf,nm='',nm2='',fs=.25,minl=2,maxl=1000,bins=.15,bias=1.,bias2=1.):
	f = open('/Users/ashleyross/CAMB_sources/counts_test_scalCls'+nm+'.dat').readlines()
	f2 = open('/Users/ashleyross/CAMB_sources/counts_test_scalCls'+nm2+'.dat').readlines()
	sum = 0
	
	cth1 = cos(theta1*pi/180.)
	lth1 = open('Pl'+str(theta1)+'.dat').readlines()
	cth2 = cos(theta2*pi/180.)
	lth2 = open('Pl'+str(theta2)+'.dat').readlines()
	sumcl2 = 0
	sumcln = 0
	sumnb2 = 0
	for l in range(minl,maxl):
		if l < minl:
			c = float(f[0].split()[-2])*2.*pi/(1.*2*(2+1))*bias**2.
			c2 = float(f2[0].split()[-2])*2.*pi/(1.*2*(2+1))*bias2**2.
		else:
			c = float(f[l-minl].split()[-2])*2.*pi/(1.*l*(l+1.))*bias**2.
			c2 = float(f2[l-minl].split()[-2])*2.*pi/(1.*l*(l+1.))*bias2**2.
		pl1 = float(lth1[l])
		pl2 = float(lth2[l])
		cx = cf*sqrt(c*c2)
		errl = (2.*l+1.)*pl1*pl2*(cx**2.)
		sum += errl
	return 2./fs*sum/(4.*pi)**2.




def cl2cov_d(theta1,theta2,nm,bias,fs=.25,ngal=5000000.,maxl=1000,bins=.2):
	f = open('clgnorm'+nm+'nodamp.dat').readlines()
	sum = 0
	innp = 0
	if theta1 == theta2:
		np = 2.*pi*theta1*(pi/180.)**2.*bins*.5*ngal*ngal/(fs*4.*pi)
		innp = 1./np	
	#print innp
	nbar = ngal/(4.*pi*fs)
	cth1 = cos(theta1*pi/180.)
	lth1 = open('Pl'+str(theta1)+'.dat').readlines()
	cth2 = cos(theta2*pi/180.)
	lth2 = open('Pl'+str(theta2)+'.dat').readlines()
	sumcl2 = 0
	sumcln = 0
	sumnb2 = 0
	for l in range(0,maxl):
		c = float(f[l].split()[1])*bias**2.
		pl1 = float(lth1[l])
		pl2 = float(lth2[l])
		#errl = (2.*l+1.)*pl1*pl2*(c+1./nbar)**2.
		errl = (2.*l+1.)*pl1*pl2*(c**2.+2.*c/nbar)
		sum += errl
	return 2./fs*sum/(4.*pi)**2.+innp


def cl2ompc(theta,zf='0.450.55',maxl=1000,md = 'rs'):
	om = 0
	omrs = 0
	rad = theta*pi/180.
	fpc = open('cl100pcDES'+zf+'.dat').readlines()
	fth = open('clTHzsig0.03'+zf+'nzDESnodamp.dat').readlines()
	mult = (float(fpc[97].split()[1])/float(fth[97].split()[1])+float(fpc[98].split()[1])/float(fth[98].split()[1])+float(fpc[99].split()[1])/float(fth[99].split()[1]))/3.
	try:
		frs = open('clrs'+zf+'.dat').readlines()
	except:
		print 'no z-dist file'
		md = 'n'
	oldpl = 1.
	cr = 0
	for l in range(0,maxl):
		#if l < 30:
		#	c = cl(l,zf)
		#else:
		#	c = cl_limb(l,zf)
		if l < 100:
			c = float(fpc[l].split()[1])
		else:
			c = float(fth[l].split()[1])*mult
		if md == 'n':
			crs = c
		else:
			if l < 30:
				crs = float(frs[l].split()[1])
			else:
				crs = c
		if l == 500:
			c1 = c
		if l == 999:
			c2 = c
		pl = legendre(l,cos(rad))
		if pl < 0 and oldpl >0:
			cr += 1
			print l
			if cr == 7:
				break
		om += c*pl*(2.*l+1.)/(4.*pi)
		omrs += crs*pl*(2.*l+1.)/(4.*pi)
		oldpl = pl
		#print l,om
	
	if l == 999:
		plin = (log(c2)-log(c1))/(log(999)-log(500))
		for l in range(1000,3000):
			pl = legendre(l,cos(rad))
			if pl < 0 and oldpl >0:
				cr += 1
				print l
				if cr == 7:
					break
			om += l**plin*c2/(999**plin)*pl*(2.*l+1.)/(4.*pi)
			omrs += l**plin*c2/(999**plin)*pl*(2.*l+1.)/(4.*pi)
			oldpl = pl
	if md == 'rs':		
		return omrs,om
	else:
		return om


def mkomfile(zf='THzsig0.030.450.55',thetmin=.5,thetap=.25,nbin=20):
	f = open('omcl'+zf+'.dat','w')
	theta = thetmin
	for i in range(0,nbin):
		om = cl2om(theta,zf+'nzDES')
		if len(om) > 1:
			f.write(str(theta)+' '+str(om[0])+' '+str(om[1])+'\n')
		else:
			f.write(str(theta)+' '+str(om)+'\n')
		theta += .25
	f.close()

def mkclfile(zf='TH0.40.60.05',maxl=1000,lphi=60,RS='RS',om=.285):
	if lphi != 0:
		f = open('cl'+zf+RS+'nodamp.dat','w')
		#f = open('cl'+zf+'.dat','w')
	else:
		f = open('clLO'+zf+'nodamp.dat','w')
	for l in range(0,maxl):
		if l < lphi:
			c = cl(l,zf,md=RS,om=om)
		else:
			c = cl_limb(l,zf,om=om)
		f.write(str(l)+' '+str(c)+'\n')
		#print l,c
	f.close()
	return True

def mkclfile_fnl(zf='TH0.40.60.05',maxl=1000,lphi=60):
	
	f = open('cl'+zf+'_fnl.dat','w')
	for l in range(0,maxl):
		if l < lphi:
			c = cl_fnl(l,zf)
		else:
			c = cl_limb_fnl(l,zf)
		f.write(str(l)+' '+str(c[0])+' '+str(c[1])+' '+str(c[2])+'\n')
	#	print l,c
	f.close()
	return True


def mkclfilecross(zf1,zf2='TH0.40.60.05',maxl=1000,lphi=30,mod='n'):
	if lphi != 0:
		f = open('cl'+zf1+zf2+'nodamp.dat','w')
	else:
		f = open('clLO'+zf1+zf2+'nodamp.dat','w')
	for l in range(0,maxl):
		if l < lphi:
			c = cl_cross(l,zf1,zf2,md=mod)
		else:
			c = cl_limbcross(l,zf1,zf2)
		f.write(str(l)+' '+str(c)+'\n')
		#print l,c
	f.close()




def mkclfilepc(mult,zf='zsig0.030.450.55',maxl=1000):
	f = open('clpc'+zf+'nodamp.dat','w')
	for l in range(0,maxl):
		if l < 30:
			c = cl_pc(l,zf)
		else:
			c = cl_limb(l,'TH'+zf+'nzDES')*mult
		f.write(str(l)+' '+str(c)+'\n')
		print l,c
	f.close()


def mkclfilers(zf='TH0.40.60.05',maxl=1000,bias=2.):
	f = open('clrs'+zf+str(bias)+'.dat','w')
	for l in range(0,maxl):
		#if l < maxl:
		c = cl(l,zf,md='RS',bias=bias)
		#else:
		#	c = cl_limb(l,zf)
		f.write(str(l)+' '+str(c)+'\n')
		#print l,c
	f.close()

def cl2cov(theta1,theta2,zf='TH0.40.6',fs=.125,ngal=5000000.,md='rs',maxl=1000,bins=.25,mult=1.,nb=1.,bias=2.,lmax=60,dir='/Users/ashleyr/BOSS/'):
	f = open(dir+'cl'+zf+'nodamp.dat').readlines()
	#f = open('cl'+zf+'.dat').readlines()
	if md == 'rs':
		frs = open(dir+'clrs'+zf+str(bias)+'.dat').readlines()
	sumrs = 0
	sum = 0
	innp = 0
	if theta1 == theta2:
		np = 2.*pi*theta1*(pi/180.)**2.*bins*.5*ngal*ngal/(fs*4.*pi)
		innp = 1./np	
	nbar = ngal/(4.*pi*fs)
	#nbar = nb
	#print nbar,sqrt(np)
	cth1 = cos(theta1*pi/180.)
	lth1 = []
	try:
		flth1 = open(dir+'lthfiles/lthfile'+str(theta1)+'.dat').readlines()
		for i in range(0,len(flth1)):
			lth1.append(float(flth1[i]))
	except:
		flth1 = open(dir+'lthfiles/lthfile'+str(theta1)+'.dat','w')
		for l in range(0,maxl):
			lth = legendre(l,cth1)
			lth1.append(lth)
			flth1.write(str(lth)+'\n')
		flth1.close()
	cth2 = cos(theta2*pi/180.)
	lth2 = []
	
	try:
		flth2 = open(dir+'lthfiles/lthfile'+str(theta2)+'.dat').readlines()
		for i in range(0,len(flth2)):
			lth2.append(float(flth2[i]))
	except:
		flth2 = open(dir+'lthfiles/lthfile'+str(theta2)+'.dat','w')
		for l in range(0,maxl):
			lth = legendre(l,cth2)
			lth2.append(lth)
			flth2.write(str(lth)+'\n')
	sumcl2 = 0
	sumcln = 0
	sumnb2 = 0
	for l in range(0,maxl):
		errl = (2.*l+1.)*lth1[l]*lth2[l]*(float(f[l].split()[1])*mult+1./nbar)**2.
		fac = (2.*l+1.)*lth1[l]*lth2[l]
		cl2 = (float(f[l].split()[1])*mult)**2.
		cln = 2.*(float(f[l].split()[1])*mult)/nbar
		nb2 = (1./nbar)**2.
		
		if l < lmax and md == 'rs':
			errlrs = (2.*l+1.)*lth1[l]*lth2[l]*(float(frs[l].split()[1])*mult+1./nbar)**2.
		else:
			errlrs = errl
		sum += errl
		sumrs += errlrs
		sumcl2 += fac*cl2
		sumcln += fac*cln
		sumnb2 += fac*nb2
		#print l, sum
	#if md == 'rs':
	#print 2./fs*sumcl2,2./fs*sumcln,2./fs*sumnb2
	return 2./fs*sumrs/(4.*pi)**2.+innp, 2./fs*sum/(4.*pi)**2.+innp

def cl2cov2bin(theta1,theta2,zf1,zf2,fs=.125,ngal=5000000.,md='n',maxl=1000,bins=.25,mult=1.,nb=1.,bias=2.,lmax=60,dir='/Users/ashleyr/BOSS/'):
	f = open(dir+'cl'+zf1+zf2+'nodamp.dat').readlines()
	#f = open('cl'+zf+'.dat').readlines()
	if md == 'rs':
		print 'STILL NEEDS TO BE DONE'
		frs = open(dir+'clrs'+zf+str(bias)+'.dat').readlines()
	sumrs = 0
	sum = 0
	#nbar = nb
	cth1 = cos(theta1*pi/180.)
	lth1 = []
	try:
		flth1 = open(dir+'lthfiles/lthfile'+str(theta1)+'.dat').readlines()
		for i in range(0,len(flth1)):
			lth1.append(float(flth1[i]))
	except:
		flth1 = open(dir+'lthfiles/lthfile'+str(theta1)+'.dat','w')
		for l in range(0,maxl):
			lth = legendre(l,cth1)
			lth1.append(lth)
			flth1.write(str(lth)+'\n')
		flth1.close()
	cth2 = cos(theta2*pi/180.)
	lth2 = []
	
	try:
		flth2 = open(dir+'lthfiles/lthfile'+str(theta2)+'.dat').readlines()
		for i in range(0,len(flth2)):
			lth2.append(float(flth2[i]))
	except:
		flth2 = open(dir+'lthfiles/lthfile'+str(theta2)+'.dat','w')
		for l in range(0,maxl):
			lth = legendre(l,cth2)
			lth2.append(lth)
			flth2.write(str(lth)+'\n')
	sumcl2 = 0
	sumcln = 0
	sumnb2 = 0
	for l in range(0,maxl):
		errl = (2.*l+1.)*lth1[l]*lth2[l]*(float(f[l].split()[1])*mult)**2.
		sum += errl
	return 2./fs*sum/(4.*pi)**2.


def cl2covcross(theta1,theta2,zf,zf2,fs=.125,ngal=5000000.,ngal2=5000000.,md='n',maxl=1000,bins=.25,mult=1.,nb=1.,bias=2.,bias2=2.,dir='/Users/ashleyr/BOSS/'):
	f = open(dir+'cl'+zf+'nodamp.dat').readlines()
	f2 = open(dir+'cl'+zf2+'nodamp.dat').readlines()
	#f = open('cl'+zf+'.dat').readlines()
	if md == 'rs':
		frs = open(dir+'clrs'+zf+str(bias)+'.dat').readlines()
		frs2 = open(dir+'clrs'+zf2+str(bias2)+'.dat').readlines()
	sumrs = 0
	sum = 0
	innp = 0
	if theta1 == theta2:
		np = 2.*pi*theta1*(pi/180.)**2.*bins*.5*ngal*ngal2/(fs*4.*pi)
		innp = 1./np	
	nbar = ngal/(4.*pi*fs)
	nbar2 = ngal2/(4.*pi*fs)
	#nbar = nb
	#print nbar,sqrt(np)
	cth1 = cos(theta1*pi/180.)
	lth1 = []
	try:
		flth1 = open(dir+'lthfiles/lthfile'+str(theta1)+'.dat').readlines()
		for i in range(0,len(flth1)):
			lth1.append(float(flth1[i]))
	except:
		flth1 = open(dir+'lthfiles/lthfile'+str(theta1)+'.dat','w')
		for l in range(0,maxl):
			lth = legendre(l,cth1)
			lth1.append(lth)
			flth1.write(str(lth)+'\n')
		flth1.close()
	cth2 = cos(theta2*pi/180.)
	lth2 = []
	
	try:
		flth2 = open(dir+'lthfiles/lthfile'+str(theta2)+'.dat').readlines()
		for i in range(0,len(flth2)):
			lth2.append(float(flth2[i]))
	except:
		flth2 = open(dir+'lthfiles/lthfile'+str(theta2)+'.dat','w')
		for l in range(0,maxl):
			lth = legendre(l,cth2)
			lth2.append(lth)
			flth2.write(str(lth)+'\n')
	sumcl2 = 0
	sumcln = 0
	sumnb2 = 0
	for l in range(0,maxl):
		errl = (2.*l+1.)*lth1[l]*lth2[l]*(float(f[l].split()[1])*mult+1./nbar)*(float(f2[l].split()[1])*mult+1./nbar2)
		fac = (2.*l+1.)*lth1[l]*lth2[l]
		cl2 = (float(f[l].split()[1])*mult)*(float(f2[l].split()[1])*mult)
		cln = (float(f[l].split()[1])*mult)/nbar2+(float(f2[l].split()[1])*mult)/nbar
		nb2 = (1./(nbar*nbar2))
		
		if l < 30 and md == 'rs':
			errlrs = (2.*l+1.)*lth1[l]*lth2[l]*(float(frs[l].split()[1])*mult+1./nbar)*(float(frs2[l].split()[1])*mult+1./nbar2)
		else:
			errlrs = errl
		sum += errl
		sumrs += errlrs
		sumcl2 += fac*cl2
		sumcln += fac*cln
		sumnb2 += fac*nb2
		if zf != zf2:
			innp = 0
		#print l, sum
	#if md == 'rs':
	#print sumcl2,sumcln,sumnb2,innp,2./fs*(sumcl2+sumcln+sumnb2)/(4.*pi)**2.+innp
	return 2./fs*sumrs/(4.*pi)**2.+innp, 2./fs*sum/(4.*pi)**2.+innp
	
def mkcovarmcl(file,zmed=.5,thetamin=.5,nbin=20,thetap=.25,ng=5000000.,m=1.,rsm='rs',fs=.125):
	d = distance()
	theta = thetamin
	fo = open('covarm'+str(thetamin)+'deg'+str(nbin)+'bin'+file+'mult'+str(m)+'cl.dat','w')
	fd = open('diagerr'+str(thetamin)+'deg'+str(nbin)+'bin'+file+'mult'+str(m)+'cl.dat','w')
	for i in range(0,nbin):
		print i
		theta1 = (thetamin + i*thetap)		
		for j in range(0,nbin):
			theta2 = (thetamin + j*thetap)
			c = cl2cov(theta1,theta2,zf=file,ngal=ng,md=rsm,bins=thetap,mult=m,fs=fs)[0]
			fo.write(str(c)+' ')
			if i ==j:
				fd.write(str(sqrt(c))+'\n')
			if j == nbin-1:
				fo.write('\n')
	fo.close()

def mkcovarmcl2bin(file1,file2,zmed=.5,thetamin=.5,nbin=20,thetap=.25,ng=5000000.,m=1.,rsm='n',fs=.125):
	d = distance()
	theta = thetamin
	fo = open('covarm2bin'+str(thetamin)+'deg'+str(nbin)+'bin'+file1+file2+'mult'+str(m)+'cl.dat','w')
	#fd = open('diagerr'+str(thetamin)+'deg'+str(nbin)+'bin'+file+'mult'+str(m)+'cl.dat','w')
	for i in range(0,2*nbin):
		print i
		if i < nbin:
			theta1 = (thetamin + i*thetap)
			im = 0
		else:
			theta1 = (thetamin + (i-nbin)*thetap)
			im = 1
		for j in range(0,2*nbin):
			if j < nbin:
				theta2 = (thetamin + j*thetap)
				jm = 0
			else:
				theta2 = (thetamin + (j-nbin)*thetap)
				jm = 1
			if im == jm:
				if im == 0:
					c = cl2cov(theta1,theta2,zf=file1,ngal=ng,md=rsm,bins=thetap,mult=m,fs=fs)[0]
				else:
					c = cl2cov(theta1,theta2,zf=file2,ngal=ng,md=rsm,bins=thetap,mult=m,fs=fs)[0]
			else:
				c = cl2cov2bin(theta1,theta2,file1,file2,ngal=ng,md=rsm,bins=thetap,mult=m,fs=fs)
			fo.write(str(c)+' ')
			if j == 2*nbin-1:
				fo.write('\n')
	fo.close()


def nzDESth(zmin,zmax):
	f = open('nzDES.dat')
	nzl = []
	zmind = int(1000*zmin)
	zmaxd = int(1000*zmax)
	for i in range(0,2000):
		nzl.append(0)
	c = 0
	for line in f:
		if c >= zmind and c < zmaxd:
			nzl[c] = float(line.split()[1])
		c += 1
	fo = open('nzDESTH'+str(zmin)+str(zmax)+'.dat','w')
	t = 0
	for i in range(0,len(nzl)):
		t += nzl[i]*.001
	div = t
	for i in range(0,len(nzl)):
		fo.write(str(i/1000.)+' '+str(nzl[i]/div)+'\n')
	fo.close()

def mkpcbinclsnophot(zmin,zmax,zw,zpcd=.475,zpcu=.525):
	nbins = int((zmax-zmin)/zw)
	z = zmin
	for i in range(0,nbins):
		zu = z+zw
		print z,zu
		try:
			open('nzDESTH'+str(z)+str(zu)+'.dat')
		except:
			nzDESth(z,zu)
		for i in range(0,30):
			try:
				open('phiL/phi'+str(i)+'DESTH'+str(z)+str(zu)+'.dat')
			except:
				mkphi(i,zf='DESTH'+str(z)+str(zu))
		try:
			open('phiL/phi50DESTH'+str(z)+str(zu))
		except:
			mkphi(50,zf='DESTH'+str(z)+str(zu))
		z = zu
	
	for i in range(0,nbins):
		for j in range(i,nbins):
			z1d = zmin+i*zw
			z1u = z1d+zw
			zm1 = (z1d+z1u)/2.
			z2d = zmin+j*zw
			z2u = z2d+zw
			zm2 = (z2d+z2u)/2.
			zmm = (zm1+zm2)/2.
			if zmm >= zpcd and zmm <= zpcu:
				if i != j:
					try:
						open('clLODESTH'+str(z1d)+str(z1u)+'DESTH'+str(z2d)+str(z2u)+'nodamp.dat')
					except:
						mkclfilecross('DESTH'+str(z1d)+str(z1u),zf2='DESTH'+str(z2d)+str(z2u),maxl=1000,lphi=0)
				else:
					try:
						open('clLODESTH'+str(z1d)+str(z1u)+'nodamp.dat')
					except:
						mkclfile('DESTH'+str(z1d)+str(z1u),maxl=1000,lphi=0)
					
	return True

def mkpcbincls(zmin,zmax,zw,zpcd=.475,zpcu=.525):
	nbins = int((zmax-zmin)/zw)
	z = zmin
	for i in range(0,nbins):
		zu = z+zw
		try:
			open('nzTHzsig0.03'+str(z)+str(zu)+'nzDES.dat')
		except:
			mknzdTH('nzDES.dat',z,zu)
		for i in range(0,30):
			try:
				open('phiL/phi'+str(i)+'THzsig0.03'+str(z)+str(zu)+'nzDES.dat')
			except:
				mkphi(i,zf='THzsig0.03'+str(z)+str(zu)+'nzDES')
		try:
			open('phiL/phi50THzsig0.03'+str(z)+str(zu)+'nzDES.dat')
		except:
			mkphi(50,zf='THzsig0.03'+str(z)+str(zu)+'nzDES')
		z = zu
	for i in range(0,(nbins+1)/2):
		z1d = zmin+i*zw
		z1u = z1d+zw
		z2d = zmin+(nbins-i-1.)*zw
		z2u = z2d+zw
		try:
			open('clTHzsig0.03'+str(z1d)+str(z1u)+'nzDESTHzsig0.03'+str(z2d)+str(z2u)+'nzDESnodamp.dat')
		except:
			mkclfilecross('THzsig0.03'+str(z1d)+str(z1u)+'nzDES',zf2='THzsig0.03'+str(z2d)+str(z2u)+'nzDES',maxl=1000,lphi=30)

def wtest(zmin,zmax,zw,zpcd=.475,zpcu=.525):
	w = .2
	nbins = int((zmax-zmin)/zw)
	sum = 0
	weight = 0
	nzf = open('nzDES.dat').readlines()
	norm1 = 0
	norm2 = 0
	bd = zpcd
	bu = zpcu
 	wt = 0
 	for i in range(int(zmin*1000),int(1000*zmax)):
 		norm1 += float(nzf[i].split()[1])*.001
	for i in range(0,nbins):
		for j in range(i,nbins):
			z1d = zmin+i*zw
			z1u = z1d+zw
			zm1 = (z1d+z1u)/2.
			z2d = zmin+j*zw
			z2u = z2d+zw
			zm2 = (z2d+z2u)/2.
			zmm = (zm1+zm2)/2.
			z1mind = int(zm1*1000)
			z2mind = int(zm2*1000)
			zf1 = 'DESTH'+str(z1d)+str(z1u)
			zf2 = 'DESTH'+str(z2d)+str(z2u)
			if zmm >= zpcd and zmm < zpcu:
				if i == j:
					zm = (float(nzf[z1mind].split()[1])/norm1)*(float(nzf[z2mind].split()[1])/norm1)
					print zmm,zm*zw
					wt += zm*w*zw*zw	

	return wt

def pccl(zmin,zmax,zw,zpcd=.475,zpcu=.525):
	nbins = int((zmax-zmin)/zw)
	sum = 0
	weight = 0
	nzf = open('nzDES.dat').readlines()
	norm1 = 0
	norm2 = 0
	bd = zpcd
	bu = zpcu
#  	for i in range(int(zmin*1000),int(1000*zmax)):
#  		norm1 += float(nzf[i].split()[1])*.001

	for i in range(int(zmin*1000),int(1000*bu)):
		norm1 += float(nzf[i].split()[1])*.001
	for i in range(int(1000*bd),int(zmax*1000)):
		norm2 += float(nzf[i].split()[1])*.001
# 	print bd,norm1,bu,norm2
	clpc = []
	for i in range(0,1000):
		clpc.append(0)
	for i in range(0,nbins):
		for j in range(i,nbins):
			z1d = zmin+i*zw
			z1u = z1d+zw
			zm1 = (z1d+z1u)/2.
			z2d = zmin+j*zw
			z2u = z2d+zw
			zm2 = (z2d+z2u)/2.
			zmm = (zm1+zm2)/2.
			z1mind = int(zm1*1000)
			z2mind = int(zm2*1000)
			zf1 = 'DESTH'+str(z1d)+str(z1u)
			zf2 = 'DESTH'+str(z2d)+str(z2u)
			if zmm >= zpcd and zmm <= zpcu:
				if i != j:
					zm = 2.*(float(nzf[z1mind].split()[1])/norm1)*(float(nzf[z2mind].split()[1])/norm2)
					f = open('clLO'+zf1+zf2+'nodamp.dat')
					for i in range(0,1000):
						cl = float(f.readline().split()[1])
						clpc[i] += zm*cl*zw*zw
				else:
					zm = (float(nzf[z1mind].split()[1])/norm1)*(float(nzf[z2mind].split()[1])/norm2)
					f = open('clLO'+zf1+'nodamp.dat')
					print zmm,zm*zw
					for i in range(0,1000):
						cl = float(f.readline().split()[1])
						clpc[i] += zm*cl*zw*zw
	#	weight += zm
	weight = 1.
	print weight
	for i in range(0,len(clpc)):
		clpc[i] = clpc[i]/weight
	fo = open('clpcDES'+str(zmin)+str(zmax)+str(zw)+'nodamp.dat','w')
	for i in range(0,1000):
		fo.write(str(i)+' '+str(clpc[i])+'\n')
	fo.close()


def pcbinerr(zmin,zmax,zw,theta=1.):
	nbins = int((zmax-zmin)/zw)
	sum = 0
	weight = 0
	nzf = open('nzDES.dat').readlines()
	norm = 0
	for i in range(0,len(nzf)):
		norm += float(nzf[i].split()[1])*.001
	sig2l = []
	for i in range(0,(nbins+1)/2):
		z1d = zmin+i*zw
		z1u = z1d+zw
		z2d = zmin+(nbins-i-1.)*zw
		z2u = z2d+zw
		zf1 = 'THzsig0.03'+str(z1d)+str(z1u)+'nzDES'
		zf2 ='THzsig0.03'+str(z2d)+str(z2u)+'nzDES'
		z1mind = int(1000*(z1d+z1u)/2.)
		z2mind = int(1000*(z2d+z2u)/2.)
		#print z1d,z1u,z2d,z2u,z1mind,z2mind,len(nzf)
		sig2 = cl2cov(theta,theta,zf=zf1+zf2,fs=.125,ngal=5000000.,md='n',maxl=1000,bins=.25)[0]
		sig2l.append(sig2)
		sum += sig2*(float(nzf[z1mind].split()[1])/norm)**2.*(float(nzf[z2mind].split()[1])/norm)**2.
		print sum
		weight += (float(nzf[z1mind].split()[1])/norm)*(float(nzf[z2mind].split()[1])/norm)
	print weight
	for i in range(0,(nbins+1)/2):
		z1d = zmin+i*zw
		z1u = z1d+zw		
		zud = zmin+(nbins-i-1.)*zw
		zuu = z2d+zw
		zumind = int(1000*(zud+zuu)/2.)
		for j in range(i+1,(nbins+1)/2):
			#print i,j
			z2d = zmin+j*zw
			z2u = z2d+zw
			z1mind = int(1000*(z1d+z1u)/2.)
			z2mind = int(1000*(z2d+z2u)/2.)
			clc = cl_cross(50,'THzsig0.03'+str(z1d)+str(z1u)+'nzDES',zf2='THzsig0.03'+str(z2d)+str(z2u)+'nzDES',zmed=.5,md='norm')
			cl1 = cl(50,'THzsig0.03'+str(z1d)+str(z1u)+'nzDES')
			cl2 = cl(50,'THzsig0.03'+str(z2d)+str(z2u)+'nzDES')
			cc = clc/sqrt(cl1*cl2)
			sig = cc*sqrt(sig2l[i]*sig2l[j])
			print i,j,sig,cc,clc,cl1,cl2
			sum += 2.*sig*(float(nzf[z1mind].split()[1])/norm)*(float(nzf[zumind].split()[1])/norm)**2.
			print sum
	return sum/weight
	

	
def phi(k,l,zf='TH0.40.60.05'):
	nzf = open('nz'+zf+'nzDES.dat').readlines()
	d = distance(.25,.75)
	dl = []
	rl = []
	sumz = 0
	for i in range(0,len(nzf)):
		z = i*.001
		dl.append(d.D(z))
		rl.append(d.dc(z))
		sumz += float(nzf[i].split()[1])
	norm = 1./sumz
	sum = 0
	for i in range(0,len(nzf)):
		sum += norm*float(nzf[i].split()[1])*dl[i]*bessj(l,k*abs(rl[i]))
		#print i,sum,dl[i],rl[i]
	return sum
	
def cl_delta(l,zf='TH0.40.60.05',zmed=.5,mult=.01):
	nzf = open('nz'+zf+'nzDES.dat').readlines()
	d = distance(.25,.75)
	dl = d.D(zmed)
	rl = d.dc(zmed)
	cfac = d.c
	sumz = 0
	for i in range(0,len(nzf)):
		z = i*.001
		sumz += float(nzf[i].split()[1])*.001
	norm = 1./sumz
	zfac = 0
	for i in range(0,len(nzf)):
		zfac += .001*(norm*float(nzf[i].split()[1]))**2.*d.Hz(i*.001)*d.D(i*.001)**2.
	zfac = zfac/d.c
	sum = 0
	s = simulate(.25,.75)
	sigv = 7.
	k = .001
	while k < 10.:
		#if k*sigv > 5.:
		#	break
		#dk = float(f[i+1].split()[0])-k
		dk = mult*k
		damp = exp(-1.*k**2.*sigv**2.)
		#damp = 1.
		pk = s.P(k,0)
		jkfac = k*k*sph_jn(l,rl*k)**2.
		sum += pk*damp*dk*jkfac
		k = k*(1.+mult)
		#print k, sum, jkfac, pk
	return 2/pi*sum*dl**2.#*zfac

def cl_delta_limb(l,zmed=.5):
	d = distance(.25,.75)
	dl = d.D(zmed)
	rl = d.dc(zmed)
	hz = d.Hz(zmed)
	s = simulate(.25,.75)
	return dl**2.*s.P((l+.5)/rl,0)*hz/rl**2.	

def pproj(kp,nzfile,zmed=.5):
	pkl = pkz(om=0.25,lam=0.75,sig8=0.8)
	s = simulate()
	f = open('Wsq'+nzfile+'.dat')
	ans = 0
	c = 0
	for line in f:
		ln = line.split()
		kz = float(ln[0])
		#if kz > 1.2:
		#	break
		#if kz > 2.85*10.**-5.:
		dk = kz*.01*.95
		kev = sqrt(kp**2.+kz**2.)
			#pk = pkl.DELsqlin(kev,zmed)/kev**3.*2.*pi**2
		pk = s.P(kev,zmed)
		wsq = float(ln[1])
		ans += wsq*dk/2./pi*pk#/-8.
		#print ans,kz,kev,wsq,dk,pk
		#c += 1.
	return ans

def mkpprojfile(nzf,zmed):
	fo = open('p'+nzf+'.dat','w')
	k = .001
	while k < 1000:
		p = pproj(k,'nz'+nzf,zmed)
		fo.write(str(k)+' '+str(p)+'\n')
		k = k*1.03
	fo.close()
	

def pprojcross(kp,nzfile1,nzfile2,zmed=.5):
	pkl = pkz(om=0.25,lam=0.75,sig8=0.8)
	s = simulate()
	f1 = open('W'+nzfile1+'.dat')
	f2 = open('W'+nzfile2+'.dat')
	ans = 0
	c = 0
	for line in f1:
		ln = line.split()
		ln2 = f2.readline().split()
		kz = float(ln[0])
		#if kz > 1.2:
		#	break
		#if kz > 2.85*10.**-5.:
		dk = kz*.01*.95
		kev = sqrt(kp**2.+kz**2.)
			#pk = pkl.DELsqlin(kev,zmed)/kev**3.*2.*pi**2
		pk = s.P(kev,zmed)
		w1 = float(ln[1])
		w2 = float(ln2[1])
		ans += w1*w2*dk/2./pi*pk#/-8.
		#print ans,kz,kev,wsq,dk,pk
		#c += 1.
	return ans

def mkpprojcrossfile(nzf,zmed):
	fo = open('ppc'+nzf+'.dat','w')
	k = .001
	while k < 1000:
		p = pprojcross(k,'nzppcd'+nzf,'nzppcu'+nzf,zmed)
		fo.write(str(k)+' '+str(p)+'\n')
		k = k*1.03
	fo.close()

		
def pproj2w2(r,pkfile='DR7.2.3_w'):
	#f = open('pkprojtemp.dat')
	f = open('p'+pkfile+'.dat')
	ans = 0
	for line in f:
		ln = line.split()
		k = float(ln[0])
		if k*r > 8.6535:
			break
		if k > .01:	
			dk = .03*k*.95
			ans += k*dk*bessj0(k*r)*float(ln[1])/(2.*pi)
	return ans

def erfint(t):
	#t = (x-x0)/(sqrt(2)sigma)
	return exp(-t**2.)

def erf(x):
	if x > 0:
		return 2./sqrt(pi)*rom(0,x,erfint)	
	else:
		return -1.*2./sqrt(pi)*rom(0,-1.*x,erfint)

def chisquare(bias,meas,model,thetamin=1.,thetamax=20.):
	fmeas = open('auto2pt'+meas+'p.75_wDR71.0Pix.dat').readlines()
	fmodel = open('omegawpmeas'+model+meas+'_w_mrss01.dat').readlines()
	chi = 0
	for i in range(0,len(fmeas)):
		lmeas = fmeas[i].split()
		theta = float(lmeas[0])
		if theta > thetamin and theta < thetamax:
			meas = float(lmeas[1])
			err = float(lmeas[2])
			lmod = fmodel[i].split()
			mod = float(lmod[1])*bias**2.
			chi += (mod-meas)**2./err**2.
			print theta,chi
	return chi

def com2rcom(file):
	m = open(file).readlines()
	fo = open(file.strip('.dat')+'r.dat','w')
	for i in range(0,len(m)):
		ln = m[i].split()
		for j in range(0,len(ln)):
			c = float(ln[j])
			di = float(ln[i])
			dj = float(m[j].split()[j])
			cr = c/sqrt(di*dj)
			fo.write(str(cr)+' ')
			if j == len(ln) -1:
				fo.write('\n')
	fo.close()

def chicalccons(cons,source1,source2,bins,min=.9,max=5.1):
	from numpy import zeros

	a = open('covarmTH'+bins+source1+'cohn.dat').readlines()
	b1 = open('om'+bins+'m0.25b0.045h0.7g0.557pcfzdTH'+source1+'.dat').readlines()
	b2 = open('om'+bins+'m0.25b0.045h0.7g0.557pcfzdpc'+source2+'.dat').readlines()
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(b1)):
		th = float(b1[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	print leng
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(a)):
		line = a[i].strip('\n').split()
		for j in range(0,len(a)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])/(float(b2[j].split()[1]))**2.
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					matrix[i-mini,j-mini] = float(line[j])/float(b2[j].split()[1])/float(b2[i].split()[1])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0

	for i in range(0,len(b2)):
		x1 = b1[i].strip('\n').split()
		y1 = b2[i].strip('\n').split()
		ra = float(x1[2])/float(y1[1])
		for j in range(0,len(b2)):	
			x2 = b1[j].strip('\n').split()
			y2 = b2[j].strip('\n').split()
			theta1 = float(x1[0])
			theta2 = float(x2[0])
			
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				rb = float(x2[2])/float(y2[1])
				moda = cons
				modb = cons
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				if i == j:
					print ChiSq,ra,rb,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
	return ChiSq

def chicalcratio(source,nzs,mult,min=.4,max=5.1,modf='temp',md='nm',gamma=0.557):
	from numpy import zeros

	mfile = open('covarratio'+nzs+'.dat').readlines()
	pcf = open('om'+source+'bias1.0m0.25b0.045h0.7g0.557pcfzdpczsig0.03'+nzs+'.dat').readlines()
	#THf = open('om'+source+'bias1.0m0.25b0.045h0.7g0.557pcfzdTHzsig0.03'+nzs+'.dat').readlines()
	THf = open('omfiles/om'+source+'bias1.0m0.25b0.045h0.7f0.849w1.0pcfzdTHzsig0.03'+nzs+'.dat').readlines()
	if md == 'mod':
		modTH = open('omfiles/om'+source+'bias1.0m0.25b0.045h0.7f'+str(gamma)+'w1.0pcfzdTHzsig0.03'+nzs+'.dat').readlines()
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(pcf)):
		th = float(pcf[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	print leng
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					matrix[i-mini,j-mini] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0

	for i in range(0,len(pcf)):
		x1 = THf[i].strip('\n').split()
		y1 = pcf[i].strip('\n').split()
		ra = float(x1[1])/float(y1[1])
		for j in range(0,len(pcf)):	
			x2 = THf[j].strip('\n').split()
			y2 = pcf[j].strip('\n').split()
			theta1 = float(x1[0])
			theta2 = float(x2[0])
			
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				rb = float(x2[1])/float(y2[1])
				moda = float(x1[2])/float(y1[2])
				modb = float(x2[2])/float(y2[2])
				if md == 'mod':
					moda = float(modTH[i].split()[1])/float(y1[1])
					modb = float(modTH[j].split()[1])/float(y2[1])
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				if i == j:
					print ChiSq,ra,moda,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
	return ChiSq


def fsigeig(binw,zs=.03,gamma=0.557):
	from numpy import zeros

	mfile = open('covarmDESfsig8'+str(zs)+str(binw)+'.dat').readlines()
	nbin = len(mfile)
	matrix = zeros((nbin,nbin),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			matrix[i,j] = float(line[j])
	return det(matrix),eig(matrix)

def fsiginvsum(binw,zs=.03):
	from numpy import zeros

	#mfile = open('covarmDESfsig8'+str(zs)+str(binw)+'.dat').readlines()
	mfile = open('covarmDESfsig8'+str(zs)+str(binw)+'bias.dat').readlines()
	nbin = len(mfile)
	matrix = zeros((nbin,nbin),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			matrix[i,j] = float(line[j])
	invC = inv(matrix)
	sum = 0
	for i in range(0,nbin):
		for j in range(0,nbin):
			sum += float(invC[i][j])
	return sum


def chicalcgamma(binw,zs=.03,gamma=0.557):
	from numpy import zeros

	mfile = open('covarmDESfsig8'+str(zs)+str(binw)+'.dat').readlines()
	fid = open('fsigfiles/fom0.25g0.557'+str(binw)+'.dat').readlines()
	mod = open('fsigfiles/fom0.25g'+str(gamma)+str(binw)+'.dat').readlines()
	nbin = len(fid)
	matrix = zeros((nbin,nbin),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
				if i == j:
					
					matrix[i,j] = float(line[j])
					if matrix[i,j] < 0:
						return 'something wrong, negative diagonal!'
				else:
					matrix[i,j] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	#print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0

	for i in range(0,len(mfile)):
		for j in range(0,len(mfile)):	
			ra = float(fid[i])
			rb = float(fid[j])
			moda = float(mod[i])
			modb = float(mod[j])
				#
			ChiSq += (ra-moda)*invCov[i,j]*(rb-modb)
			#if i == j:
			#	print ChiSq,ra,moda,invCov[i,j],matrix[i,j],mtest[i,j]
	return ChiSq

def mkffile8bin(gam):
	d = distance(.25,.75)
	fo = open('fom0.25g'+str(gam)+'8bin.dat','w')
	z = .5
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	z = .6
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	z = .705
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	z = .82
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	z = .94
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	z = 1.065
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	z = 1.195
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	z = 1.33
	f = d.omz(z)**gam
	fo.write(str(z)+' '+str(f)+'\n')
	fo.close()
	
def mkbins(binw,z0=.5,zs=0.03):
	
	d = distance(.25,.75)
	z = z0
	zmin = z-binw*(1.+z)/2.
	zmax = z+binw*(1.+z)/2.
	zmin = int(1000*zmin)/1000.
	zmax = int(1000*zmax)/1000.
	z = (zmax+zmin)/2.
	zbl = []
	zbl.append(zmin)
	zbl.append(zmax)
	zml = []
	zml.append(z)
	ngal = 15000000.*binw/.0666667
	fb = open('DESbins'+str(zs)+str(binw)+'.dat','w')
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat','w')
	fb.write(str(zmin)+'\n')
	fb.write(str(zmax)+'\n')
	fm.write(str(z)+'\n')
	while zmax < 1.4:
		print zmin,zmax
		#mknzdTH('nzDES.dat',zmin,zmax,zsig=zs)
		#print 'nz file done'
		for l in range(0,30):
			mkphi(l,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
		print 'phi made'
		mkclfile('THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
		print 'cl file made'
		thmin = 12./(d.dc(z)*pi/180.)
		thmin = int(100*thmin)/100.
		thplus = 6./(d.dc(z)*pi/180.)
		thplus = int(100*thplus)/100.
		#omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thmin,gam=.557,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
		#print 'w2 file made'
		mkcovarmcl('THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES',zmed=z,thetamin=thmin,nbin=20,thetap=thplus,ng=ngal,m=1.,rsm='nrs')
		print 'covar file made'
		zmin = zmax
		zmax = zmin+binw*(1.+zmin)
		z = (zmin+zmax)/2.
		zmax = zmin+binw*(1.+z)
		z = (zmin+zmax)/2.
		zmax = zmin+binw*(1.+z)
		z = (zmin+zmax)/2.
		zmin = int(1000*zmin)/1000.
		zmax = int(1000*zmax)/1000.
		zbl.append(zmax)		
		zml.append(z)
		fb.write(str(zmax)+'\n')
		fm.write(str(z)+'\n')
	print zmin,zmax	
	mknzdTH('nzDES.dat',zmin,zmax,zsig=zs)
	print 'nz file done'
	for l in range(0,30):
		mkphi(l,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
	print 'phi made'
	mkclfile('THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
	print 'cl file made'
	thmin = 12./(d.dc(z)*pi/180.)
	thmin = int(100*thmin)/100.
	thplus = 6./(d.dc(.5)*pi/180.)
	thplus = int(100*thplus)/100.
	omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thmin,gam=.557,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	print 'w2 file made'
	mkcovarmcl('THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES',zmed=z,thetamin=thmin,nbin=20,thetap=thplus,ng=ngal,m=1.,rsm='nrs')
	print 'covar file made'
	fm.close()
	fb.close()
	return zbl,zml
		
def findferrs(binw,zs=.03):
	d = distance(.25,.75)
	fb = open('DESbins'+str(zs)+str(binw)+'.dat').readlines()
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	fo = open('DESferr'+str(zs)+str(binw)+'a.dat','w')
	for i in range(0,len(fm)):
		z = float(fm[i])
		print z
		zmin = float(fb[i])
		zmax = float(fb[i+1])
		gamhi = findgamhi(zmin,zmax,z,zs=zs)
		gamlo = findgamlo(zmin,zmax,z,zs=zs)
		flo = d.omz(z)**gamhi
		fhi = d.omz(z)**gamlo
		err = abs(fhi-flo)/2.
		fo.write(str(err)+'\n')
	fo.close()

def mkferrssubfiles(binw,zs):
	d = distance(.25,.75)
	fb = open('DESbins'+str(zs)+str(binw)+'.dat').readlines()
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	ngal = 15000000.
	ngalth = 15000000.*binw/.0666
	for i in range(0,len(fm)):
		z = float(fm[i])
		print z
		zmin = float(fb[i])
		zmax = float(fb[i+1])
		#if i > 0:
		fc = d.omz(z)**.557
		fc = int(fc*1000)/1000.
		fileth = 'zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat'
		filepc = 'zdpczsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat'
		thmin = 12./(d.dc(z)*pi/180.)
		thmin = int(100*thmin)/100.
		thplus = 6./(d.dc(z)*pi/180.)
		thplus = int(100*thplus)/100.
		#if i > 0:
		fo = open('desfsub'+str(z)+'.sub','w')
		fo.write('universe = vanilla \n executable = dzz.tcsh \n'+'arguments = 0.557 '+str(thmin)+' '+str(thplus)+' '+filepc+' pc '+str(z)+'\n'+'output = dzz'+str(z)+'.out \n log = crosscall.log \n error = dzz'+str(z)+'.err \n queue \n\n')
		fo.write('universe = vanilla \n executable = dzz.tcsh \n'+'arguments = '+str(fc)+' '+str(thmin)+' '+str(thplus)+' '+filepc+' th '+str(z)+'\n'+'output = dzz'+str(z)+'.out \n log = crosscall.log \n error = dzz'+str(z)+'.err \n queue \n\n')
		f = .025
		while f < 1.5:
			fo.write('universe = vanilla \n executable = dzz.tcsh \n'+'arguments = '+str(f)+' '+str(thmin)+' '+str(thplus)+' '+filepc+' th '+str(z)+'\n'+'output = dzz'+str(z)+'.out \n log = crosscall.log \n error = dzz'+str(z)+'.err \n queue \n\n')

			if f > .5 and f < 1.1:
				f += .05
			else:
				f += .1
		fo.close()


def findferrsb(binw,zs=.03):
	d = distance(.25,.75)
	fb = open('DESbins'+str(zs)+str(binw)+'.dat').readlines()
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	fo = open('DESferr'+str(zs)+str(binw)+'bias.dat','w')
	ngal = 15000000.
	ngalth = 15000000.*binw/.0666
	start = 0
	if zs == .03:
		start = 1
	for i in range(start,len(fm)):
		z = float(fm[i])
		print z
		zmin = float(fb[i])
		zmax = float(fb[i+1])
		#if i > 0:
		mknzdpc('nzDES.dat',zmin,zmax,zerr=zs)
		fc = d.omz(z)**.557
		fc = int(fc*1000)/1000.
		fileth = 'zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat'
		filepc = 'zdpczsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat'
		thmin = 12./(d.dc(z)*pi/180.)
		thmin = int(100*thmin)/100.
		thplus = 6./(d.dc(z)*pi/180.)
		thplus = int(100*thplus)/100.
		#if i > 0:
		omegadzz(filepc,z,thmin=thmin,p=thplus,md='pc')
		omegadzz(fileth,z,thmin=thmin,p=thplus,gam=fc)
		fpc = open('omfiles/om'+str(thmin)+'deg20binbias1.0m0.25b0.045h0.7g0.557pcfzdpczsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat').readlines()
		fth = open('omfiles/om'+str(thmin)+'deg20binbias1.0m0.25b0.045h0.7f'+str(fc)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat').readlines()
		mult = float(fpc[0].split()[1])/float(fth[0].split()[2])
		mkcovarmcl('THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES',zmed=z,thetamin=thmin,nbin=20,thetap=thplus,ng=ngalth,rsm='nrs')
		#if i > 0:
		mkcovarmcl('THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES',zmed=z,thetamin=thmin,nbin=20,thetap=thplus,ng=ngal,m=mult,rsm='nrs')
		f = .025
		while f < 1.5:
			omegadzz(fileth,z,thmin=thmin,p=thplus,gam=f)
			if f > .5 and f < 1.1:
				f += .05
			else:
				f += .1
		omf = str(thmin)+'deg20bin'
		zmf = str(zmin)+str(zmax)+'nzDES'
		ferr = chicombggrid(omf,zmf,fc,min=2.*thmin,mt=mult,fmin=.025,fmax=1.5,zsig=zs)
		print 'ferr is '+str(ferr)+' at z = '+str(z)
		fo.write(str(z)+' '+str(ferr)+'\n')
	fo.close()

def mksig8vz(binw,zs='',sig8=.8):
	d = distance(.25,.75)
	f = open('DESbinsmedz'+zs+str(binw)+'.dat').readlines()
	fo = open('sig8vz'+str(binw)+'.dat','w')
	for i in range(0,len(f)):
		z = float(f[i])
		sig8z = sig8*d.D(z)
		fo.write(str(sig8z)+'\n')
	fo.close()

	

def mkfcovarm(binw,zs=.03):
	fb = open('DESbins'+str(zs)+str(binw)+'.dat').readlines()
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	cll = []
	for i in range(0,len(fm)):
		z = float(fm[i])
		zmin = float(fb[i])
		zmax = float(fb[i+1])
		try:
			c = cl(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
			cll.append(c)
		except:
			mkphi(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
			cll.append(cl(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES'))
	ferr = open('DESferr'+str(zs)+str(binw)+'.dat').readlines()
	fo = open('covarmDESf'+str(zs)+str(binw)+'.dat','w')
	for i in range(0,len(fm)):
		zmini = float(fb[i])
		zmaxi = float(fb[i+1])
		for j in range(0,len(fm)):
			zminj = float(fb[j])
			zmaxj = float(fb[j+1])
			if i == j:
				fo.write(str(float(ferr[i])**2.))
			else:
				clc = cl_cross(50,'THzsig'+str(zs)+str(zmini)+str(zmaxi)+'nzDES','THzsig'+str(zs)+str(zminj)+str(zmaxj)+'nzDES')
				rcov = clc**2./(cll[i]*cll[j])
				fcov = rcov*float(ferr[i])*float(ferr[j])
				fo.write(str(fcov))
			if j == len(fm)-1:
				fo.write('\n')
			else:
				fo.write(' ')
	fo.close()

def mkfsig8covarm(binw,zs=.03):
	if zs != 0.03:
		fb = open('DESbins'+str(zs)+str(binw)+'.dat').readlines()
		fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
		ferr = open('DESferr'+str(zs)+str(binw)+'.dat').readlines()
	else:
		fb = open('DESbins'+str(binw)+'.dat').readlines()
		fm = open('DESbinsmedz'+str(binw)+'.dat').readlines()
		ferr = open('DESferr'+str(binw)+'.dat').readlines()
	
	fsig8 = open('sig8vz'+str(binw)+'.dat').readlines()
	cll = []
	for i in range(0,len(fm)):
		z = float(fm[i])
		zmin = float(fb[i])
		zmax = float(fb[i+1])
		try:
			c = cl(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
			cll.append(c)
		except:
			mkphi(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
			cll.append(cl(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES'))
	fo = open('covarmDESfsig8'+str(zs)+str(binw)+'.dat','w')
	for i in range(0,len(fm)):
		zmini = float(fb[i])
		zmaxi = float(fb[i+1])
		for j in range(0,len(fm)):
			zminj = float(fb[j])
			zmaxj = float(fb[j+1])
			if i == j:
				fo.write(str((float(ferr[i])*float(fsig8[i]))**2.))
			else:
				clc = cl_cross(50,'THzsig'+str(zs)+str(zmini)+str(zmaxi)+'nzDES','THzsig'+str(zs)+str(zminj)+str(zmaxj)+'nzDES')
				rcov = clc**2./(cll[i]*cll[j])
				fcov = rcov*float(ferr[i])*float(fsig8[i])*float(ferr[j])*float(fsig8[j])
				fo.write(str(fcov))
			if j == len(fm)-1:
				fo.write('\n')
			else:
				fo.write(' ')
	fo.close()

def mkfsig8covarmbias(binw,zs=.03):
	fb = open('DESbins'+str(zs)+str(binw)+'.dat').readlines()
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	ferr = open('DESferr'+str(zs)+str(binw)+'bias.dat').readlines()
	
	fsig8 = open('sig8vz'+str(binw)+'.dat').readlines()
	cll = []
	for i in range(0,len(fm)):
		z = float(fm[i])
		zmin = float(fb[i])
		zmax = float(fb[i+1])
		try:
			c = cl(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
			cll.append(c)
		except:
			mkphi(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES')
			cll.append(cl(50,'THzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES'))
	fo = open('covarmDESfsig8'+str(zs)+str(binw)+'bias.dat','w')
	for i in range(0,len(fm)):
		zmini = float(fb[i])
		zmaxi = float(fb[i+1])
		zi = (zmaxi+zmini)/2.
		bi = 1. + zi-.5
		for j in range(0,len(fm)):
			zminj = float(fb[j])
			zmaxj = float(fb[j+1])
			zj = (zmaxj+zminj)/2.
			bj = 1. + zj-.5
			if i == j:
				fo.write(str((bi*float(ferr[i].split()[1])*float(fsig8[i]))**2.))
			else:
				clc = cl_cross(50,'THzsig'+str(zs)+str(zmini)+str(zmaxi)+'nzDES','THzsig'+str(zs)+str(zminj)+str(zmaxj)+'nzDES')
				rcov = clc**2./(cll[i]*cll[j])
				fcov = rcov*float(ferr[i].split()[1])*bi*float(fsig8[i])*bj*float(ferr[j].split()[1])*float(fsig8[j])
				fo.write(str(fcov))
			if j == len(fm)-1:
				fo.write('\n')
			else:
				fo.write(' ')
	fo.close()


def mkfsigzfile(gam,binw,zs=.03):
	d = distance(.25,.75,gam=gam)
	ld = d.mkD()
	if zs != 0.03:
		fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	else:
		fm = open('DESbinsmedz'+str(binw)+'.dat').readlines()
	fo = open('fsigfiles/fom0.25g'+str(gam)+str(binw)+'.dat','w')
	for i in range(0,len(fm)):
		z = float(fm[i])
		zind = int(1000*z)
		fsig = d.omz(z)**gam*.8*ld[zind]
		fo.write(str(fsig)+'\n')

def mkfsigfile(binw,zs,gam):
	d = distance(.25,.75,gam=gam)
	ld = d.mkD()	
	d = distance(.25,.75,gam=.557)
	ld0 = d.mkD()	
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	f0 = []
	fgl = []
	bl = []
	norm = 1.
	#norm = ld0[1088000]/ld[1088000]
	print norm
	fo = open('fsig'+str(binw)+str(zs)+str(gam)+'z0.dat','w')
	for i in range(0,len(fm)):
		z = float(fm[i])
		zind = int(1000*z)
		fsig = d.omz(z)**gam*.8*ld[zind]*norm
		fo.write(str(z)+' '+str(fsig)+'\n')
	fo.close()
	
	
def gamchi(binw,zs,gam,bm='n'):
	from numpy import zeros
	d = distance(.25,.75,gam=gam)
	ld = d.mkD()	
	d = distance(.25,.75,gam=.557)
	ld0 = d.mkD()	
	fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	f0 = []
	fgl = []
	bl = []
	norm = ld0[1088000]/ld[1088000]
	print norm
	for i in range(0,len(fm)):
		z = float(fm[i])
		zind = int(1000*z)
		fsig = d.omz(z)**gam*.8*ld[zind]*norm
		fgl.append(fsig)
		fsig = d.omz(z)**.557*.8*ld0[zind]
		f0.append(fsig)
		if bm == 'y':
			bl.append(1.+(z-.5))
		else:
			bl.append(1.)

	covf = open('covarmDESfsig8'+str(zs)+str(binw)+'bias.dat').readlines()
	m = zeros((len(covf),len(covf)),'f')
	for i in range(0,len(covf)):
		ln = covf[i].split()
		for j in range(0,len(ln)):
			m[i][j] = float(ln[j])*bl[i]*bl[j]
	invm = inv(m)
	chi2 = 0
	for i in range(0,len(f0)):
		for j in range(0,len(f0)):
			chi2 += (f0[i]-fgl[i])*invm[i][j]*(f0[j]-fgl[j])
	return chi2

	
def findgamhirange(binw,zs=.03,tol=.01):
	d = distance(.25,.75)
	if zs != .03:
		fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	else:
		fm = open('DESbinsmedz'+str(binw)+'.dat').readlines()
	gamhi = 1.2
	gamlo = .6
	gamg = .9
	mkfsigzfile(0.557,binw,zs)
	mkfsigzfile(gamhi,binw,zs)
	mkfsigzfile(gamlo,binw,zs)
	mkfsigzfile(gamg,binw,zs)
	chicalcgamma(binw,zs,gamma=0.557)
	chihi = chicalcgamma(binw,zs,gamhi)
	chilo = chicalcgamma(binw,zs,gamlo)
	chig = chicalcgamma(binw,zs,gamg)
	chihi = abs(chihi-1.)
	chilo = abs(chilo-1.)
	chig = abs(chig-1.)
	if chihi < chilo and chihi < chig:
		chimin = chihi
		gammin = gamhi
	if chilo < chihi and chilo < chig:
		chimin = chilo
		gammin = gamlo
	if chig < chilo and chig < chihi:
		chimin = chig
		gammin = gamg
	while chimin > tol:
		if chig == chimin:
			if chihi < chilo:
				gamlo = gamg
				gamg = (gamg+gamhi)/2.
				
			else:
				gamhi = gamg
				gamg = (gamg+gamlo)/2.
			newgam = gamg
		if chihi == chimin:
			gamlo = gamg
			gamg = gamhi
			gamhi = gamhi + (gamg-gamlo)
			newgam = gamhi
		if chilo == chimin:
			gamhi = gamg
			gamg = gamlo
			gamlo = gamlo - (gamhi-gamg)
			if gamlo < .557:
				gamlo = .557
				print 'gamlo = .557, possible problem!!!'
			newgam = gamlo
		mkfsigzfile(newgam,binw,zs)
		chihi = chicalcgamma(binw,zs,gamhi)
		chilo = chicalcgamma(binw,zs,gamlo)
		chig = chicalcgamma(binw,zs,gamg)
		chihi = abs(chihi-1.)
		chilo = abs(chilo-1.)
		chig = abs(chig-1.)
		if chihi < chilo and chihi < chig:
			chimin = chihi
			gammin = gamhi
		if chilo < chihi and chilo < chig:
			chimin = chilo
			gammin = gamlo
		if chig < chilo and chig < chihi:
			chimin = chig
			gammin = gamg
	return gammin
		
def findgamlorange(binw,zs=.03,tol=.01):
	d = distance(.25,.75)
	if zs != .03:
		fm = open('DESbinsmedz'+str(zs)+str(binw)+'.dat').readlines()
	else:
		fm = open('DESbinsmedz'+str(binw)+'.dat').readlines()
	gamhi = .5
	gamlo = -.3
	gamg = 0
	mkfsigzfile(0.557,binw,zs)
	mkfsigzfile(gamhi,binw,zs)
	mkfsigzfile(gamlo,binw,zs)
	mkfsigzfile(gamg,binw,zs)
	chihi = chicalcgamma(binw,zs,gamhi)
	chilo = chicalcgamma(binw,zs,gamlo)
	chig = chicalcgamma(binw,zs,gamg)
	chihi = abs(chihi-1.)
	chilo = abs(chilo-1.)
	chig = abs(chig-1.)
	if chihi < chilo and chihi < chig:
		chimin = chihi
		gammin = gamhi
	if chilo < chihi and chilo < chig:
		chimin = chilo
		gammin = gamlo
	if chig < chilo and chig < chihi:
		chimin = chig
		gammin = gamg
	while chimin > tol:
		print gamlo, gamg,gamhi
		if chig == chimin:
			if chihi < chilo:
				gamlo = gamg
				gamg = (gamg+gamhi)/2.
				
			else:
				gamhi = gamg
				gamg = (gamg+gamlo)/2.
			newgam = gamg
		if chihi == chimin:
			gamlo = gamg
			gamg = gamhi
			gamhi = gamhi + (gamg-gamlo)
			if gamhi > .557:
				gamhi = .557
				print 'gamhi = .557, possible problem!!!'
			newgam = gamhi
		if chilo == chimin:
			gamhi = gamg
			gamg = gamlo
			gamlo = gamlo - (gamhi-gamg)
			newgam = gamlo
		mkfsigzfile(newgam,binw,zs)
		chihi = chicalcgamma(binw,zs,gamhi)
		chilo = chicalcgamma(binw,zs,gamlo)
		chig = chicalcgamma(binw,zs,gamg)
		chihi = abs(chihi-1.)
		chilo = abs(chilo-1.)
		chig = abs(chig-1.)
		if chihi < chilo and chihi < chig:
			chimin = chihi
			gammin = gamhi
		if chilo < chihi and chilo < chig:
			chimin = chilo
			gammin = gamlo
		if chig < chilo and chig < chihi:
			chimin = chig
			gammin = gamg
	return gammin
	
	
def findgamhi(zmin,zmax,z,zs=0.03,tol=.03):
	d = distance(.25,.75)
	thm = 12./(d.dc(z)*pi/180.)
	thm = int(100*thm)/100.
	thplus = 6./(d.dc(.5)*pi/180.)
	thplus = int(100*thplus)/100.
	gamhi = 3.
	gamlo = 1.5
	gamg = 2.25
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g0.557pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=0.557,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(gamhi)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=gamhi,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(gamlo)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=gamlo,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(gamg)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=gamg,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	chihi = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamhi,h=.7,bias=1.0,mm='mod',zsig=zs)
	chilo = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamlo,h=.7,bias=1.0,mm='mod',zsig=zs)
	chig = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamg,h=.7,bias=1.0,mm='mod',zsig=zs)
	chihi = abs(chihi-1.)
	chilo = abs(chilo-1.)
	chig = abs(chig-1.)
	if chihi < chilo and chihi < chig:
		chimin = chihi
		gammin = gamhi
	if chilo < chihi and chilo < chig:
		chimin = chilo
		gammin = gamlo
	if chig < chilo and chig < chihi:
		chimin = chig
		gammin = gamg
	while chimin > tol:
		if chig == chimin:
			if chihi < chilo:
				gamlo = gamg
				gamg = (gamg+gamhi)/2.
				
			else:
				gamhi = gamg
				gamg = (gamg+gamlo)/2.
			newgam = gamg
		if chihi == chimin:
			gamlo = gamg
			gamg = gamhi
			gamhi = gamhi + (gamg-gamlo)
			newgam = gamhi
		if chilo == chimin:
			gamhi = gamg
			gamg = gamlo
			gamlo = gamlo - (gamhi-gamg)
			if gamlo < .557:
				gamlo = .557
				print 'gamlo = .557, possible problem!!!'
			newgam = gamlo
		try:
			open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(newgam)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
		except:
			omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=newgam,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
		chihi = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamhi,h=.7,bias=1.0,mm='mod',zsig=zs)
		chilo = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamlo,h=.7,bias=1.0,mm='mod',zsig=zs)
		chig = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamg,h=.7,bias=1.0,mm='mod',zsig=zs)
		chihi = abs(chihi-1.)
		chilo = abs(chilo-1.)
		chig = abs(chig-1.)
		if chihi < chilo and chihi < chig:
			chimin = chihi
			gammin = gamhi
			if chilo < chig:
				chi2 = chilo
				gam2 = gamlo
			else:
				chi2 = chig
				gam2 = gamg
		if chilo < chihi and chilo < chig:
			chimin = chilo
			gammin = gamlo
			if chihi < chig:
				chi2 = chihi
				gam2 = gamhi
			else:
				chi2 = chig
				gam2 = gamg
		if chig < chilo and chig < chihi:
			chimin = chig
			gammin = gamg
			if chilo < chihi:
				chi2 = chilo
				gam2 = gamlo
			else:
				chi2 = chihi
				gam2 = gamhi
	#gamr = (gammin/chimin+gam2/chi2)/(1./chim+1./chi2)
	return gammin
		
def findgamlo(zmin,zmax,z,zs=0.03,tol=.03):
	d = distance(.25,.75)
	thmin = 12./(d.dc(z)*pi/180.)
	thmin = int(100*thmin)/100.
	thm = thmin
	thplus = 6./(d.dc(.5)*pi/180.)
	thplus = int(100*thplus)/100.
	gamhi = 0
	gamlo = -.5
	gamg = -.25
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g0.557pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=0.557,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(gamhi)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=gamhi,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(gamlo)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=gamlo,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	try:
		open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(gamg)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
	except:
		omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=gamg,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
	chihi = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamhi,h=.7,bias=1.0,mm='mod',zsig=zs)
	chilo = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamlo,h=.7,bias=1.0,mm='mod',zsig=zs)
	chig = chicalcTH(str(thm)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thm+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamg,h=.7,bias=1.0,mm='mod',zsig=zs)
	chihi = abs(chihi-1.)
	chilo = abs(chilo-1.)
	chig = abs(chig-1.)
	if chihi < chilo and chihi < chig:
		chimin = chihi
		gammin = gamhi
	if chilo < chihi and chilo < chig:
		chimin = chilo
		gammin = gamlo
	if chig < chilo and chig < chihi:
		chimin = chig
		gammin = gamg
	while chimin > tol:
		if chig == chimin:
			if chihi < chilo:
				gamlo = gamg
				gamg = (gamg+gamhi)/2.
			else:
				gamhi = gamg
				gamg = (gamg+gamlo)/2.
			newgam = gamg
		if chihi == chimin:
			gamlo = gamg
			gamg = gamhi
			gamhi = gamhi + (gamg-gamlo)
			if gamhi > .557:
				gamhi = .557
				print 'gamhi = .557, possible problem!!!'
			newgam = gamhi
		if chilo == chimin:
			gamhi = gamg
			gamg = gamlo
			gamlo = gamlo - (gamhi-gamg)
			newgam = gamlo
		try:
			open('omfiles/om'+str(thm)+'deg20binbias1.0m0.25b0.045h0.7g'+str(newgam)+'pcfzdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat')
		except:
			omegadzz('zdTHzsig'+str(zs)+str(zmin)+str(zmax)+'nzDES.dat',z,of='om',nbin=20,thmin=thm,gam=newgam,m=.25,b=.045,hub=.7,bias=1.0,p=thplus)
		chihi = chicalcTH(str(thmin)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thmin+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamhi,h=.7,bias=1.0,mm='mod',zsig=zs)
		chilo = chicalcTH(str(thmin)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thmin+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamlo,h=.7,bias=1.0,mm='mod',zsig=zs)
		chig = chicalcTH(str(thmin)+'deg20bin',str(zmin)+str(zmax)+'nzDES',min=thmin+.01,max=5.1,md='f',b=0.045,m=0.25,g=gamg,h=.7,bias=1.0,mm='mod',zsig=zs)
		chihi = abs(chihi-1.)
		chilo = abs(chilo-1.)
		chig = abs(chig-1.)
		if chihi < chilo and chihi < chig:
			chimin = chihi
			gammin = gamhi
			if chilo < chig:
				chi2 = chilo
				gam2 = gamlo
			else:
				chi2 = chig
				gam2 = gamg
		if chilo < chihi and chilo < chig:
			chimin = chilo
			gammin = gamlo
			if chihi < chig:
				chi2 = chihi
				gam2 = gamhi
			else:
				chi2 = chig
				gam2 = gamg
		if chig < chilo and chig < chihi:
			chimin = chig
			gammin = gamg
			if chilo < chihi:
				chi2 = chilo
				gam2 = gamlo
			else:
				chi2 = chihi
				gam2 = gamhi
	#gamr = (gammin/chimin+gam2/chi2)/(1./chim+1./chi2)
	return gammin
			
	
def chicalcratiocosmo(source,nzs,min=.4,max=5.1,gamma1=0.557,g2=0.557,bias=1.,m=0.25,h=0.7,b=.045,mult=0.116279069767,zsig=.03):
	from numpy import zeros

	mfile = open('covarratio'+nzs+str(mult)+'.dat').readlines()
	pcf = open('omfiles/om'+source+'bias1.0'+'m'+str(m)+'b'+str(b)+'h'+str(h)+'f0.7pcfzdpczsig'+str(zsig)+nzs+'.dat').readlines()
	THf = open('omfiles/om'+source+'bias1.0'+'m'+str(m)+'b'+str(b)+'h'+str(h)+'f'+str(gamma1)+'pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	modTH = open('omfiles/om'+source+'bias'+str(bias)+'m0.25b0.045h0.7f'+str(g2)+'pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	modpc = open('omfiles/om'+source+'bias'+str(bias)+'m0.25b0.045h0.7f0.7pcfzdpczsig'+str(zsig)+nzs+'.dat').readlines()
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(pcf)):
		th = float(pcf[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	print leng,mini,maxi
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					matrix[i-mini,j-mini] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0

	for i in range(0,len(pcf)):
		#print i
		x1 = THf[i].strip('\n').split()
		y1 = pcf[i].strip('\n').split()
		ra = float(x1[1])/float(y1[1])
		for j in range(0,len(pcf)):	
			x2 = THf[j].strip('\n').split()
			y2 = pcf[j].strip('\n').split()
			theta1 = float(x1[0])
			theta2 = float(x2[0])
			
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				rb = float(x2[1])/float(y2[1])
				moda = float(modTH[i].split()[1])/float(modpc[i].split()[1])
				modb = float(modTH[j].split()[1])/float(modpc[j].split()[1])
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				if i == j:
					print ChiSq,ra,moda,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
	return ChiSq

def testcosmoratio(bs,g0,b=0.054,m=0.25,g1=0.557,h=.7):
	omegadzz('zdTHzsig0.031.01.13nzDES.dat',1.065,nbin=20,thmin=.3,p=.15,bias=bs,gam=g0)
	omegadzz('zdpczsig0.031.01.13nzDES.dat',1.065,nbin=20,thmin=.3,p=.15,bias=bs,gam=g0)
	chi = chicalcratiocosmo('0.3deg20bin','zsig0.031.01.13nzDES',min=.29,max=5.1,b=b,m=m,gamma1=g1,g2=g0,h=h,bias=bs)
	return chi

def chol(zf1='0.4850.515',zf2='0.5150.545',th=0.5):
	from numpy import zeros
	from numpy.linalg import *
	from random import gauss
	#mt = bias**2.
	#mt = 1.
	mfile = open('covarm2bin'+str(th)+'deg20binTHzsig0.03'+zf1+'nzDESTHzsig0.03'+zf2+'nzDESmult1.0cl.dat').readlines()

	leng= len(mfile)
	print leng
	matrix = zeros((leng,leng),'f')
	for i in range(0,leng):
		line = mfile[i].strip('\n').split()
		for j in range(0,leng):
			matrix[i,j] = float(line[j])#*.001
	a = cholesky(matrix)
	#print a
			
	return a

def pchol(zf1='0.4850.515',zf2='0.5150.545',n=100,fmin1=0.702,fmin2=.712,thm1=0.5,thm2=0.48,mth=.9):
	#n = 100
	fo = open('dfbins'+zf1+zf2+'.dat','w')
	from random import gauss
	a = chol(zf1,zf2,thm1)
	fw1 = open('omfiles/om'+str(thm1)+'deg20binbias1.0m0.25b0.045h0.7f'+str(fmin1)+'pcfzdTHzsig0.03'+zf1+'nzDES.dat').readlines()
	#fsig1 = open('diagerr0.5deg20binTHzsig0.03'+zf+'nzDESmult1.0cl.dat').readlines()
	if zf2 == '0.5250.575':
		fw2 = open('omfiles/om'+str(thm1)+'deg20binbias1.0m0.25b0.045h0.7f'+str(fmin2)+'w1.0pcfzdTHzsig0.03'+zf2+'nzDES.dat').readlines()
	else:
		fw2 = open('omfiles/om'+str(thm2)+'deg20binbias1.0m0.25b0.045h0.7f'+str(fmin2)+'pcfzdTHzsig0.03'+zf2+'nzDES.dat').readlines()
	#fsig2 = open('diagerr0.5deg20binTHzsig0.030.5250.575nzDESmult1.0cl.dat').readlines()
	
	for k in range(0,n):
		print k
		dx = []
		for i in range(0,40):
			if i < 20:
				dx.append(gauss(0,1.))
			else:
				dx.append(gauss(0,1.))
		pw1 = []
		pw2 = []
		for i in range(0,40):
			dxc = 0
			for j in range(0,40):
				dxc += a[i][j]*dx[j]
			if i < 20:
				pw1.append(dxc)
			else:
				pw2.append(dxc)
		chimin1 = 1000
		chimin2 = 1000
		chiu1 = 1000
		chiu2 = 1000
		fv = .025
		fv2 = .025
		fm1 = 0
		fm2 = 0
		while fv < 1.45:
			chi1 = chicalcTH(str(thm1)+'deg20bin',zf1+'nzDES',min=mth,fmin=fmin1,g=fv,pw=pw1)
			if chi1 < chimin1:
				chiu1 = chimin1
				fo1 = fm1
				chimin1 = chi1
				fm1 = fv
			else:
				if chi1 < chiu1:
					chiu1 = chi1
					fo1 = fv
			chi2 = chicalcTH(str(thm2)+'deg20bin',zf2+'nzDES',min=mth,fmin=fmin2,g=fv2,pw=pw2)
			if chi2 < chimin2:
				chiu2 = chimin2
				fo2 = fm2
				chimin2 = chi2
				fm2 = fv2
			else:
				if chi2 < chiu2:
					chiu2 = chi2
					fo2 = fv
			if fv > .45 and fv < 1.1:
				fv += .05
				fv2 += .05
			else:
				fv += .1
				fv2 += .1
		fmm1 = (fm1/chimin1+fo1/chiu1)/(1./chimin1+1./chiu1)
		fmm2 = (fm2/chimin2+fo2/chiu2)/(1./chimin2+1./chiu2)
		fo.write(str(fmin1-fmm1)+' '+str(fmin2-fmm2)+'\n')
	fo.close()
	f = open('dfbins'+zf1+zf2+'.dat').readlines()
	sum1 = 0
	sum2 = 0
	sum3 = 0
	for i in range(0,len(f)):
		ln = f[i].split()
		d1 = float(ln[0])
		d2 = float(ln[1])
		sum1 += d1*d1
		sum2 += d2*d2
		sum3 += d2*d1
	sig1 = sqrt(sum1/float(n))
	sig2 = sqrt(sum2/float(n))
	cov = sum3/float(n)
	return sig1,sig2,cov, cov/(sig1*sig2)

def chicalcTH(source,nzs,fmin=0.557,min=.4,max=5.1,md='f',b=0.045,m=0.25,g=0.557,h=.7,bias=1.0,mm='mod',zsig=0.03,mt=1.,pw='n'):
	from numpy import zeros
	#mt = bias**2.
	#mt = 1.
	mfile = open('covarm'+source+'THzsig'+str(zsig)+nzs+'mult'+str(bias**2.)+'cl.dat').readlines()
	#mfile = open('covarm'+source+'THzsig0.03'+nzs+'cl.dat').readlines()
	#THf = open('omfiles/om'+source+'bias'+str(bias)+'m0.25b0.045h0.7g0.557pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	#THf = open('omfiles/om'+source+'bias'+str(bias)+'m0.25b0.045h0.7f'+str(fmin)+'pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	#THf = open('omfiles/om'+source+'bias'+str(bias)+'m'+str(m)+'b'+str(b)+'h'+str(h)+'f'+str(fmin)+'pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	try:
		THf = open('omfiles/om'+source+'bias'+str(bias)+'m'+str(m)+'b'+str(b)+'h'+str(h)+'f'+str(fmin)+'pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	except:
		THf = open('omfiles/om'+source+'bias'+str(bias)+'m'+str(m)+'b'+str(b)+'h'+str(h)+'f'+str(fmin)+'w1.0pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	#THf = open('om'+source+'bias1.0'+'m0.25b0.045h0.7g0.557pcfzdTHzerr0.9'+nzs+'.dat').readlines()
	if mm == 'mod':
		try:
			modf = open('omfiles/om'+source+'bias'+str(bias)+'m0.25b0.045h0.7f'+str(g)+'pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
		except:
			modf = open('omfiles/om'+source+'bias'+str(bias)+'m0.25b0.045h0.7f'+str(g)+'w1.0pcfzdTHzsig'+str(zsig)+nzs+'.dat').readlines()
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(THf)):
		th = float(THf[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	#print leng
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					if md != 'd':
						matrix[i-mini,j-mini] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	#print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0
	oldchi = 0
	for i in range(0,len(THf)):
		x1 = THf[i].strip('\n').split()
		ra = float(x1[1])
		if pw != 'n':
			ra = ra+pw[i]
		for j in range(0,len(THf)):	
			x2 = THf[j].strip('\n').split()
			theta1 = float(x1[0])
			theta2 = float(x2[0])
			
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				rb = float(x2[1])
				if pw != 'n':
					rb = rb+pw[j]
				moda = float(x1[2])
				modb = float(x2[2])
				if mm == 'mod':
					moda = float(modf[i].split()[1])*mt
					modb = float(modf[j].split()[1])*mt
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				if i == j:
					chidif = ChiSq - oldchi
					#print chidif,ra-moda,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
					oldchi = ChiSq
	return ChiSq

def chicombggrid(source,nzs,fmc,min=.5,max=5.1,mt=0.116279069767,zsig=.03,bmin=.8,bmax=1.2,bstep=.01,fmin=.05,fmax=1.6,fstep=.05,bmode='y'):
	b = bmin	
	fo = open(nzs+'fbgrid.dat','w')
	bthmin = min/2.
	while b < bmax:
		f = fmin
		if bmode == 'y':
			chib = chicalcpcbias(source,nzs,min=bthmin,max=max,bias=b,zsig=zsig,mt=mt)
		else:
			if b == 1.:
				chib = .1
			else:
				chib = 100.
		while f < fmax:
			feff = f/b
			if feff < fmin or feff > fmax:
				fo.write('0 ')
			else:
				ftest = fmin
				oldftest = ftest
				while ftest < feff:
					oldftest = ftest
					if ftest > .5 and ftest < 1.1:
						ftest += fstep
					else:
						ftest += 2.*fstep
				if abs(feff-oldftest) < abs(feff-ftest) or ftest > fmax:
					ftest = oldftest
				chig = chicalcTH(source,nzs,fmc,min=min,max=max,g=ftest,mm='mod',zsig=zsig)
				chi = chib+chig
				fo.write(str(exp(-1.*chi/2.))+' ')
			if f > .5 and f < 1.1:
				f += fstep
			else:
				f += 2.*fstep
			
		fo.write('\n')
		b += bstep
		#print b
	findd = -1000
	findu = 2000
	findd = int(0.5/(fstep*2))
	findu = int(findd+.6/fstep)+1
	fo.close()
	f = open(nzs+'fbgrid.dat').readlines()
	sum = 0
	for i in range(0,len(f)):
		ln = f[i].split()
		for j in range(0,len(ln)):
			if j > findd and j < findu:
				sum += float(ln[j])
			else:
				sum += float(ln[j])*2.
	print sum
	find = int((fmc-fmin+.001)/fstep)
	find = int(findd+(fmc-0.5+.001)/fstep)
	
	per = 0
	ran = 0
	while per < .68:
		oldper = per
		per = 0
		ran += 1		
		for i in range(0,len(f)):
			ln = f[i].split()
			for j in range(find-ran,find+1+ran):
				per += float(ln[j])/sum
		print per,ran,find
	if ran <= findu - findd:
		f1 = ran*fstep
	else:
		#f1 =fstep*(2.*ran-(findu-find))
		f1 = (findu-find)*fstep+(ran-(findu-find))*2.*fstep
	if ran-1 <= findu - findd:
		f2 = (ran-1)*fstep
	else:
		#f1 =fstep*(2.*ran-(findu-find))
		f2 = (findu-find)*fstep+(ran-1.-(findu-find))*2.*fstep

	err = (f1/abs(per-.68)+f2/abs(oldper-.68))/(1./abs(per-.68)+1./abs(oldper-.68))
	print per, oldper
	return err	

def chicombggridsimp(source,nzs,fmc,min=.5,max=5.1,mt=0.116279069767,zsig=.03,bmin=.8,bmax=1.2,bstep=.01,fmin=.05,fmax=1.6,fstep=.05,bmode='y'):
	b = bmin	
	fo = open(nzs+'fbgridsimp.dat','w')
	bthmin = min/2.
	while b < bmax:
		f = fmin
		if bmode == 'y':
			chib = chicalcpcbias(source,nzs,min=bthmin,max=max,bias=b,zsig=zsig,mt=mt)
		else:
			if b == 1.:
				chib = .1
			else:
				chib = 100.
		#print chib
		while f < fmax:
			#if bmode == 'y':
			#	feff = f/b
			#else:
			feff = f
			if feff < fmin or feff > fmax:
				fo.write('0 ')
			else:
				ftest = fmin
				oldftest = ftest
				while ftest < feff:
					oldftest = ftest
					ftest += fstep
				if abs(feff-oldftest) < abs(feff-ftest) or ftest > fmax:
					ftest = oldftest
				chig = chicalcTH(source,nzs,fmc,min=min,max=max,g=ftest,mm='mod',zsig=zsig)
				if bmode == 'y':
					chi = chib+chig
				else:
					chi = chig
				fo.write(str(exp(-1.*chi/2.))+' ')
			f += fstep
			
		fo.write('\n')
		b += bstep
		#print b
	findd = -1000
	findu = 2000
	findd = int(0.5/(fstep*2))
	findu = int(findd+.6/fstep)+1
	fo.close()
	f = open(nzs+'fbgridsimp.dat').readlines()
	sum = 0
	for i in range(0,len(f)):
		ln = f[i].split()
		for j in range(0,len(ln)):
			#if j > findd and j < findu:
			sum += float(ln[j])
			#else:
			#	sum += float(ln[j])*2.
	print sum
	find = int((fmc-fmin+.001)/fstep)
	#find = int(findd+(fmc-0.5+.001)/fstep)
	
	per = 0
	ran = 0
	while per < .68:
		oldper = per
		per = 0
		ran += 1		
		for i in range(0,len(f)):
			ln = f[i].split()
			for j in range(find-ran,find+1+ran):
				per += float(ln[j])/sum
		print per,ran,find
	#if ran <= findu - findd:
	f1 = ran*fstep
	#else:
	#	#f1 =fstep*(2.*ran-(findu-find))
	#	f1 = (findu-find)*fstep+(ran-(findu-find))*2.*fstep
	#if ran-1 <= findu - findd:
	f2 = (ran-1)*fstep
	#else:
	#	#f1 =fstep*(2.*ran-(findu-find))
	#	f2 = (findu-find)*fstep+(ran-1.-(findu-find))*2.*fstep

	err = (f1/abs(per-.68)+f2/abs(oldper-.68))/(1./abs(per-.68)+1./abs(oldper-.68))
	print per, oldper
	return err	

	
def chicalcpcbias(source,nzs,min=.4,max=5.1,md='f',b=0.045,m=0.25,g=0.557,h=.7,bias=1.0,zsig=0.03,mt=0.116279069767):
	from numpy import zeros
	#mfile = open('covarm'+source+'THzsig'+str(zsig)+nzs+'mult'+str(mt)+'cl.dat').readlines()
	mfile = open('pccov'+nzs+str(zsig)+'th.dat').readlines()
	THf = open('omfiles/om'+source+'bias1.0m0.25b0.045h0.7g0.557pcfzdpczsig'+str(zsig)+nzs+'.dat').readlines()
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(THf)):
		th = float(THf[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	#print leng
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					if md != 'd':
						matrix[i-mini,j-mini] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	#print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0
	oldchi = 0
	for i in range(0,len(THf)):
		x1 = THf[i].strip('\n').split()
		ra = float(x1[1])
		for j in range(0,len(THf)):	
			x2 = THf[j].strip('\n').split()
			theta1 = float(x1[0])
			theta2 = float(x2[0])
			
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				rb = float(x2[1])
				moda = float(x1[1])*bias**2.
				modb = float(x2[1])*bias**2.
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				if i == j:
					chidif = ChiSq - oldchi
					#print chidif,ra-moda,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
					oldchi = ChiSq
	return ChiSq


def chiminfindmock(ind,source='0.2deg40bin',nzs='MICE0.05',min=.8,max=5.01,gammax=2.11,gammin=-.6):
	if nzs == 'MICE0.05':
		mockf = '1345.83_width-89.12_i-'
	if nzs == 'MICE0.15':
		mockf = '1338.77_width-267.42_i-'
	if nzs == 'zsig0.040.46250.5375nzDES':
		mockf = '1345.83_width-89.12_i-'
	if nzs == 'zsig0.020.46250.5375nzDES':
		mockf = '1345.83_width-89.12_i-'
	if nzs == 'zsig0.040.47750.5225nzDES':
		mockf = '1346.40_width-53.47_i-'
	if nzs == 'zsig0.020.47750.5225nzDES':
		mockf = '1346.40_width-53.47_i-'
	if nzs == 'zsig0.030.4850.515nzDES':
		mf = 'Outputphoto2/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(ind)+'_zm-0.50_dz-0.03'
	#mf = 'RedshiftSpace-Photoz0p03-dNdz-wtheta/Realizations/CorrFile_zspace_photozspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-'+mockf+str(ind)+'-0.20'
	mf = 'RedshiftSpace-Photoz0p03-dNdz-wtheta/Realizations/CorrFile_zspace_photozspace_0p03_dNdz_Lbox-7680.00_Nmocks-125_rmean-'+mockf+str(ind)+'-0.20'
	chimin = 1000.
	f = open(mf).readlines()
	ind = int(min/.2)
	#fm = open('omeganl'+nzs+'g0.55701.dat').readlines()
	#mt = float(fm[ind].split()[1])/float(f[ind].split()[0])
	gam = gammin
	while gam < gammax:
		if abs(gam) < .000001:
			gam = 0.0
		chi = chidiagmice(source,nzs,g=gam,THfile=mf,min=min,max=max)
		#print chi
		if chi < chimin:
			chimin = chi
			gammin = gam
		if gam < -.6:
			gam += .1
		else:
			gam += .1
# 	gam = 2.5
# 	chi = chidiagmice(source,nzs,g=gam,THfile=mf,min=min,max=max)
# 		#print chi
# 	if chi < chimin:
# 		chimin = chi
# 		gammin = gam
# 	gam = 3.5
# 	chi = chidiagmice(source,nzs,g=gam,THfile=mf,min=min,max=max)
# 		#print chi
# 	if chi < chimin:
# 		chimin = chi
# 		gammin = gam
# 	
# 	chi = chidiagmice(source,nzs,g=gam,THfile=mf,min=min,max=max,md='nf')
# 	if chi < chimin:
# 		chimin = chi
# 		gammin = 999
# 	
# 	if gammin == 999:
# 		fmin = 0
# 	else:
# 		from Cosmo import *
# 		d = distance(.25,.75)
# 		fmin = d.omz(.5)**gammin
	fmin = gammin
	print gammin, chimin, fmin
	return gammin, chimin, fmin

def chiminmultmock(nzs,g,min=.15,max=5.1,TH='n',chid=.001):
	mo = 1.
	mu = 1.2
	md = .8
	chi = chicalcmock('0.2deg40bin',nzs,min=min,max=max,md='f',g=g,THfile=TH,mult=mo)
	chi1 = chicalcmock('0.2deg40bin',nzs,min=min,max=max,md='f',g=g,THfile=TH,mult=mu)
	chi2 = chicalcmock('0.2deg40bin',nzs,min=min,max=max,md='f',g=g,THfile=TH,mult=md)
	chif = abs(chi/chi1-1.)
	while chif > chid:
		if chi1 < chi and chi1 < chi2:
			chimin = chi1
			md = mo
			mo = mu
			mu = mu + (mu-md)
		else:
			if chi2 < chi and chi2 < chi1:
				chimin = chi1
				mu = mo
				mo = md
				md = md + (md-mu)
			else:
				if chi1 < chi2:
					md = mo
					mo = mo + (mu-mo)/2.
				else:
					mu = mo
					mo = mo + (md-mo)/2.
		chi = chicalcmock('0.2deg40bin',nzs,min=min,max=max,md='f',g=g,THfile=TH,mult=mo)
		chi1 = chicalcmock('0.2deg40bin',nzs,min=min,max=max,md='f',g=g,THfile=TH,mult=mu)
		chi2 = chicalcmock('0.2deg40bin',nzs,min=min,max=max,md='f',g=g,THfile=TH,mult=md)
		chif = abs(chi/chi1-1.)
		print md,mo,mu,chi
	return mo
		
def chidiagmice(source,nzs,THfile='n',min=.4,max=5.1,g=0.557,h=.7,bias=1.0,mult=1.,md='g'):
	from numpy import zeros
	#fm = open('omfiles/om0.2deg40binbias1.0m0.25b0.044h0.7g0.557pcfzdTH'+nzs+'.dat').readlines()
	fm = open('omfiles/om0.2deg40binbias1.0m0.25b0.044h0.7f0.7pcfzdTH'+nzs+'.dat').readlines()
	#modf = open('omfiles/om0.2deg40binbias1.0m0.25b0.044h0.7g'+str(g)+'pcfzdTH'+nzs+'.dat').readlines()
	modf = open('omfiles/om0.2deg40binbias1.0m0.25b0.044h0.7f'+str(g)+'pcfzdTH'+nzs+'.dat').readlines()
	modff = open('omfiles/om0.2deg40binbias1.0m0.25b0.044h0.7f0.7pcfzdTH'+nzs+'.dat').readlines()
	THf = open(THfile).readlines()
	#errf = open('errTH'+nzs+'.dat').readlines()
	ind = int(min/.2)
	#print ind,len(modf),len(THf)
	avef = open('w2MICE'+nzs.strip('nzDES')+'.dat').readlines()
	mt = float(modff[ind].split()[1])/float(avef[ind].split()[1])
	#mt = float(modf[ind].split()[1])/float(THf[ind].split()[0])
	ChiSq = 0
	errl = []
	for i in range(0,len(THf)):
		errl.append((float(THf[i].split()[1])*mt)**2.)
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(THf)):
		th = float(fm[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	#print leng
	#leng = len(b2)

# 	matrix = zeros((leng,leng),'f')
# 	for i in range(0,leng):
# 		e1 = errl[i+mini]
# 		for j in range(0,leng):
# 			e2 = errl[j+mini]
# 
# 			matrix[i,j] = .95**abs(i-j)*sqrt(e1*e2)
# 					#if matrix[i,j] < 0:# or abs(i-j) > 7:
# 					#	matrix[i,j] = 0.1
# 	#print 'made matrix'		
# 	#print matrix
# 	invCov = inv(matrix)
# 
# 	ChiSq = 0
# 	oldchi = 0
# 	for i in range(0,len(THf)):
# 		x1 = THf[i].strip('\n').split()
# 		ra = float(x1[0])*mt
# 		for j in range(0,len(THf)):	
# 			x2 = THf[j].strip('\n').split()
# 			theta1 = float(modf[i].split()[0])
# 			theta2 = float(modf[j].split()[0])
# 			
# 			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
# 				
# 				rb = float(x2[0])*mt
# 
# 				moda = float(modf[i].split()[1])
# 				modb = float(modf[j].split()[1])
# 				#
# 				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				#if i == j:
				#	print moda,modb,ra,rb,invCov[i-mini,j-mini]
	#print ChiSq



	for i in range(0,len(THf)):
		theta1 = float(fm[i].split()[0])
		ra = float(THf[i].split()[0])*mt
		err = float(THf[i].split()[1])*mt
		#err = float(errf[i].split()[1])*mt
		if theta1 > min and theta1 < max:
				
			moda = float(modf[i].split()[1])
				#
			if md == 'nf':
				moda = float(modf[i].split()[2])
			#print(ra,moda,err)
			ChiSq += (ra-moda)**2./err**2.
	return ChiSq



def chicalcmock(source,nzs,min=.4,max=5.1,md='f',g=0.557,h=.7,bias=1.0,THfile='n',mult=1.,fm='n'):
	from numpy import zeros
	mt = bias**2.
	mfile = open('covarm'+source+nzs+'mult'+str(mt)+'cl.dat').readlines()
	#mfile = open('covarm'+nzs+'mock.dat').readlines()
	fm = open('omeganl'+nzs+'g0.55701.dat').readlines()
	modf = open('omeganl'+nzs+'g'+str(g)+'01.dat').readlines()
	ind = int(min/.2)
	if THfile == 'n':
		THf = open('omave'+nzs+'.dat').readlines()
		mt = float(modf[ind].split()[1])/float(THf[ind].split()[1])
	else:
		reffile = open('omave'+nzs+'.dat').readlines()
		THf = open(THfile).readlines()
		mt = float(fm[ind].split()[1])/float(THf[ind].split()[0])
	if fm != 'n':
		THf=fm
	#mult = mt
	#if mm != 'mt':
	#	mult = 1.
	#print mult
	
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(THf)):
		if THfile == 'n':
			th = float(THf[i].split()[0])
		else:
			th = float(reffile[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	#print leng
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					if md != 'd':
						matrix[i-mini,j-mini] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	#print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0

	for i in range(0,len(THf)):
		if THfile == 'n':
			x1 = THf[i].strip('\n').split()
			ra = float(x1[1])*mult
		else:
			theta1 = float(reffile[i].split()[0])
			ra = float(THf[i].split()[0])*mult
		for j in range(0,len(THf)):	
			if THfile == 'n':
				x2 = THf[j].strip('\n').split()
				theta1 = float(x1[0])
				theta2 = float(x2[0])
				rb = float(x2[1])*mult	
			else:
				theta2 = float(reffile[j].split()[0])
				rb = float(THf[j].split()[0])*mult
				
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				moda = float(modf[i].split()[1])
				modb = float(modf[j].split()[1])
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				#if i == j:
				#	print theta1,mult,ChiSq,ra-moda,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
	return ChiSq



def testcosmoTH(bs,g0,b=0.054,m=0.25,g1=0.557,h=.7):
	omegadzz('zdTHzsig0.031.01.13nzDES.dat',1.065,nbin=20,thmin=.3,p=.15,bias=bs,gam=g0)
	chi = chicalcTHcosm('0.3deg20bin','zsig0.031.01.13nzDES',min=.29,max=5.1,b=b,m=m,g0=g0,g1=g1,h=h,bias=bs)
	return chi


def chicalcTHcosm(source,nzs,min=.4,max=5.1,md='f',b=0.045,m=0.25,g0=0.557,g1=0.557,h=.7,bias=1.0):
	from numpy import zeros

	mfile = open('covarm'+source+'TH'+nzs+'cl.dat').readlines()
	THf = open('om'+source+'bias'+str(bias)+'m0.25b0.045h0.7g'+str(g0)+'pcfzdTH'+nzs+'.dat').readlines()
	
	modf = open('om'+source+'bias1.0m'+str(m)+'b'+str(b)+'h'+str(h)+'g'+str(g1)+'pcfzdTH'+nzs+'.dat').readlines()
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(THf)):
		th = float(THf[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	print leng
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					if md != 'd':
						matrix[i-mini,j-mini] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0

	for i in range(0,len(THf)):
		x1 = THf[i].strip('\n').split()
		ra = float(x1[1])
		for j in range(0,len(THf)):	
			x2 = THf[j].strip('\n').split()
			theta1 = float(x1[0])
			theta2 = float(x2[0])
			
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				rb = float(x2[1])
				moda = float(modf[i].split()[1])
				modb = float(modf[j].split()[1])
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				if i == j:
					print ChiSq,ra,moda,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
	return ChiSq



def chicalcf(source,f,nzs,min=.4,max=5.1):
	from numpy import zeros

	mfile = open('covarmTH'+source+nzs+'cohn.dat').readlines()
	THf = open('om'+source+'h0.7g0.557pcfzdTH'+nzs+'.dat').readlines()
	ff = open('om'+source+'h0.7g'+str(f)+'pcfzdTH'+nzs+'.dat').readlines()
	leng = 0
	mini = 1000
	maxi = 0
	for i in range(0,len(THf)):
		th = float(THf[i].split()[0])
		#print th,min,m
		if th > min and th < max:
			
			if i < mini:
				mini = i
			if i > maxi:
				maxi = i
			leng += 1
	print leng
	#leng = len(b2)

	matrix = zeros((leng,leng),'f')
	for i in range(0,len(mfile)):
		line = mfile[i].strip('\n').split()
		for j in range(0,len(mfile)):
			if i >= mini and i <= maxi and j>=mini and j <= maxi:
				if i == j:
					
					matrix[i-mini,j-mini] = float(line[j])
					if matrix[i-mini,j-mini] < 0:
						return 'something wrong, negative diagonal!'
				else:
					matrix[i-mini,j-mini] = float(line[j])#*.001
					#if matrix[i,j] < 0:# or abs(i-j) > 7:
					#	matrix[i,j] = 0.1
	print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	#print invCov
	ChiSq = 0

	for i in range(0,len(THf)):
		x1 = THf[i].strip('\n').split()
		f1 = ff[i].strip('\n').split()
		ra = float(x1[1])
		for j in range(0,len(THf)):	
			x2 = THf[j].strip('\n').split()
			f2 = ff[j].strip('\n').split()
			theta1 = float(x1[0])
			theta2 = float(x2[0])
			
			if theta1 > min and theta1 < max and theta2 > min and theta2 < max:
				
				rb = float(x2[1])
				moda = float(f1[1])
				modb = float(f2[1])
				#
				ChiSq += (ra-moda)*invCov[i-mini,j-mini]*(rb-modb)
				if i == j:
					print ChiSq,ra,rb,invCov[i-mini,j-mini],matrix[i-mini,j-mini],mtest[i-mini,j-mini]
	return ChiSq



def nzerr(file,zerrt,zc=3):
	f = open(file)
	c = 0
	zerrsum = 0
	for line in f:
		zerr = float(line.split()[3])
		if zerr < zerrt:
			c += 1.
			zerrsum += zerr
	if c > 0:
		zerrave = zerrsum/c
	else:
		zerrave = 0
	return c,zerrave

def mknvzerr(file,of='temp'):
	d = distance(h=.7)
	solidang = 6000.
	vol = solidang/(360.*360./pi)*4.*pi/3.*(d.dc(.3)**3.-d.dc(.2)**3.)
	zt = .1
	fo = open('nvzerr'+of+'.dat','w')
	while zt > 0:
		n,zea = nzerr(file,zt)
		n = n/vol
		fo.write(str(zt)+' '+str(n)+' '+str(zea)+'\n')
		zt -= .001
	fo.close()

def pperCV(k,n,z1,z2):
	
	pkl = pkz(om=0.3,lam=0.7,sig8=0.8,h=.7)
	pk = pkl.DELsqlin(k,(z1+z2)/2.)/k**3.*2.*pi**2
	print pk
	d = distance(h=1.)
	#v = d.covol(z1,z2)
	v = 4.*pi/3.*(d.dc(z2)**3.-d.dc(z1)**3.)
	print v
	dk = k*.1
	vk = 4.*pi*k**2.*dk#/(2.*pi)**3.
	nbar = n/v
	print nbar,pk,sqrt(v),sqrt(vk)
	return (2.*pi)**1.5/sqrt(vk)/sqrt(v)*(1.+1./(nbar*pk))
	
def pperCA(k,zmed,nzfile,n=5000000.,f=.125):
	
	#pkl = pkz(om=0.3,lam=0.7,sig8=0.8,h=.7)
	#pk = pkl.DELsqlin(k,(z1+z2)/2.)/k**3.*2.*pi**2
	pk = pproj(k,nzfile,zmed)
	#print pk
	d = distance(h=1.)
	#v = d.covol(z1,z2)
	A = 4.*pi*(d.dc(zmed)**2.)*f
	#print A
	dk = k*.1
	Ak = 2.*pi*k**dk#/(2.*pi)**3.
	nbar = n/A
	#print nbar,pk,sqrt(v),sqrt(vk)
	return (2.*pi)/sqrt(Ak)/sqrt(A)*(1.+1./(nbar*pk))*pk

def sigmaw2sum(r,dr,zmed=.5,n=5000000.,f=.125,pkfile='projcrosspc'):
	d = distance(h=1.)
	pkfile = open('p'+pkfile+'.dat').readlines()
	sig2 = 0
	#k = .001

	#while k < 1000.:
	m = 1.-float(pkfile[1].split()[0])/float(pkfile[0].split()[0])
	for i in range(0,len(pkfile)):
		ln = pkfile[i].split()
		k = float(ln[0])
		pk = float(ln[1])
		A = 4.*pi*(d.dc(zmed)**2.)*f
		dk = k*.95*m
		dka = 2.*pi/dr
		Ak = 2.*pi*k*dka/(2.*pi)**2.
		nbar = n/A
		pksig = (2.*pi)/sqrt(Ak)/sqrt(A)*(1.+1./(nbar*pk))*pk		
		sig2 += pksig**2.*(1./(2.*pi**2.)*k*dk*pk*bessj0(k*r))**2.
		
	return sig2

def covarw2sum(r,r2,dr,dr2,zmed=.5,n=5000000.,f=.125,pkfile='projcrosspc'):
	d = distance(h=1.)
	pkfile = open('p'+pkfile+'.dat').readlines()
	covarp = 0
	covarc = 0
	covarn = 0
	m = 1.-float(pkfile[1].split()[0])/float(pkfile[0].split()[0])
	A = 4.*pi*(d.dc(zmed)**2.)*f
	nbar = n/A
	print A,nbar,m
	for i in range(0,len(pkfile)):
		ln = pkfile[i].split()
		k = float(ln[0])
		pk = float(ln[1])
		dk = abs(k*.95*m)
		dka = 2.*pi/dr
		dka2 = 2.*pi/dr2
		#Ak = 2.*pi*k*dka/(2.*pi)**2.
		#Ak2 = 2.*pi*k*dka2/(2.*pi)**2.
		Ak = 2.*pi*k*dk
		Ak2 = Ak
		
		pksig = (2.*pi)/sqrt(Ak)/sqrt(A)#*(1.+1./(nbar*pk))*pk
		pksig2 = (2.*pi)/sqrt(Ak2)/sqrt(A)#*(1.+1./(nbar*pk))*pk
		sig2r = pksig*(1./(2.*pi**2.)*k*dk*bessj0(k*r))#(pksig*(1./(2.*pi**2.)*k*dk*pk*bessj0(k*r))+pksig*(1./(2.*pi**2.)*k*dk*pk*bessj0(k*ru))+pksig*(1./(2.*pi**2.)*k*dk*pk*bessj0(k*rd)))/3.
		sig2r2 = pksig2*(1./(2.*pi**2.)*k*dk*bessj0(k*r2))#(pksig*(1./(2.*pi**2.)*k*dk*pk*bessj0(k*r2))+pksig*(1./(2.*pi**2.)*k*dk*pk*bessj0(k*r2u))+pksig*(1./(2.*pi**2.)*k*dk*pk*bessj0(k*r2d)))/3.
		frac = sig2r/sig2r2
		frac = 1.
		#if frac > 1.:
		#	frac = 1./frac
		covarp += frac*sig2r*sig2r2*pk**2.
		covarc += 2.*sig2r*sig2r2*pk/nbar
		covarn += sig2r*sig2r2*(1./nbar)**2.
	covar = covarp + covarc + covarn
	print covarp,covarc,covarn
	return covar

def jdt0(r,k,dr):
	return 1./(r*k*dr)*((r+dr/2.)*bessj1(k*(r+dr/2.))-(r-dr/2.)*bessj1(k*(r-dr/2.)))

def covarw2cohn(r,r2,dr,dr2,zmed=.5,n=5000000.,f=.125,pkfile='kprojpcDES.4.6'):
	d = distance(h=1.)
	pkfile = open('p'+pkfile+'.dat').readlines()
	covar = 0
	covarp = 0
	covarc = 0
	
	sa = 4.*pi*f
	A = (d.dc(zmed)**2.)*sa
	innp = 0.
	theta = r/d.dc(zmed)
	dth = dr/d.dc(zmed)
	if r == r2:
		#np = 2.*pi*theta*dth*.5*n*n/(4.*pi)
		np = 2.*pi*r*dr*.5*n*n/A
		innp = 1./np
	
	#m = 1./A/pi/pi/2.
	m = 1./pi/A
	nbar = n/A
	
	mk = abs(1.-float(pkfile[1].split()[0])/float(pkfile[0].split()[0]))
	#print innp,nbar,A,mk
	for i in range(0,len(pkfile)):
		ln = pkfile[i].split()
		k = float(ln[0])
		pk = float(ln[1])#*pi/2.
		dk = k*.95*mk		
		#covar += k*dk*(pk**2.+2./nbar*pk)*jdt0(r,k,dr)*jdt0(r2,k,dr2)
		covarp += k*dk*(pk**2.)*jdt0(r,k,dr)*jdt0(r2,k,dr2)
		covarc += k*dk*(2./nbar*pk)*jdt0(r,k,dr)*jdt0(r2,k,dr2)
		#print k, dk, pk, jdt0(r,k,dr)*jdt0(r2,k,dr2), covar
	covar = covarp*m+covarc*m+innp
	print covarp*m,covarc*m,innp
	return covar




def mkcovarm(file,zmed=.5,thetamin=.5,nbin=10.,thetap=.5,ng=5000000.):
	d = distance()
	theta = thetamin
	fo = open('covarm'+file+'.dat','w')
	rs = d.dc(zmed)
	#dr = rs*pi/180.*thetap
	dr = rs*pi/180.*.5
	for i in range(0,nbin):
		theta1 = thetamin + i*thetap
		r1 = rs*pi/180.*theta1
		
		for j in range(0,nbin):
			theta2 = thetamin + j*thetap
			r2 = rs*pi/180.*theta2
			
			c = covarw2sum(r1,r2,dr,dr,zmed=zmed,n=ng,pkfile=file)
			fo.write(str(c)+' ')
			if j == nbin-1:
				fo.write('\n')
	fo.close()

def mkcovarmc(file,zmed=.5,thetamin=.5,nbin=10.,thetap=.5,ng=5000000.):
	d = distance()
	theta = thetamin
	fo = open('covarm'+file+'cohn.dat','w')
	rs = d.dc(zmed)
	#dr = rs*pi/180.*thetap
	dr = rs*pi/180.*.5
	for i in range(0,nbin):
		theta1 = thetamin + i*thetap
		r1 = rs*pi/180.*theta1
		
		for j in range(0,nbin):
			theta2 = thetamin + j*thetap
			r2 = rs*pi/180.*theta2
			
			c = covarw2cohn(r1,r2,dr,dr,zmed=zmed,n=ng,pkfile=file)
			fo.write(str(c)+' ')
			if j == nbin-1:
				fo.write('\n')
	fo.close()

def mkcovarmratcl(file,zf,mult=1.,cv=1./7.,m=.25,b=.045,h=.7,g=.849,mcon='n',con=.2,thmin=.3):
	#d1 = open('om'+str(thmin)+file+'bias1.0m0.25b0.045h0.7f'+str(g)+'pcfzdTHzsig0.03'+zf+'.dat').readlines()
	d1 = open('om0.3deg20binbias1.0m0.25b0.045h0.7f0.849pcfzdTHzsig0.031.01.13nzDES.dat').readlines()
	#if mcon != 'y':
	#	d2 = open('om'+str(thmin)+file+'bias1.0m0.25b0.045h0.7f0.7pcfzdpczsig0.03'+zf+'.dat').readlines()
	#else:
	#	d2 = open('om'+str(thmin)+'g'+str(g)+'pcfzdpczsig0.03cons'+str(con)+zf+'.dat').readlines()
	d2 = open('om0.3deg20binbias1.0m0.25b0.045h0.7f0.849pcfzdpczsig0.031.01.13nzDES.dat').readlines()
	#mult = (float(d2[0].split()[1])/float(d1[0].split()[1]))**2.
	#print mult
	if len(d1) != len(d2):
		return 'ERROR, mismatched files!'
	dl1 = []
	dl2 = []
	for i in range(0,len(d1)):
		dl1.append(float(d1[i].split()[1]))
		dl2.append(float(d2[i].split()[1]))
	cf1 = open('covarm'+str(thmin)+file+'THzsig0.03'+zf+'mult1.0'+'cl.dat').readlines()
	#cf2 = open('covarm'+str(thmin)+file+'THzsig0.03'+zf+'mult'+str(mult)+'cl.dat').readlines()
	cf2 = open('pccov1.01.13th.dat').readlines()
	if len(cf1) != len(dl1):
		return 'ERROR, mismatched covariance/data files'
	#fo = open('covarratio'+zf+str(mult)+'.dat','w')
	#fdiag = open('ratiodiagerrs'+zf+str(mult)+'.dat','w')
	fo = open('covarratio'+zf+'.dat','w')
	fdiag = open('ratiodiagerrs'+zf+'.dat','w')

	for i in range(0,len(dl1)):
		x1 = dl1[i]
		y1 = dl2[i]
		lnc1 = cf1[i].split()
		lnc2 = cf2[i].split()
		for j in range(0,len(dl1)):
			x2 = dl1[j]
			y2 = dl2[j]
			cx = float(lnc1[j])
			cy = float(lnc2[j])
			cxy = sqrt(abs(cx*cy))*cv
			if cx*cy < 0:
				cxy = -1.*cxy
			ctot = cx/y1/y2+cy*x1*x2/y1**2./y2**2.-1.*cxy*x1/y1**2./y2-1.*cxy*x2/y1/y2**2.
			if i == j:
				print cx/y1/y2,cy*x1*x2/y1**2./y2**2.,-1.*cxy*x1/y1**2./y2,-1.*cxy*x2/y1/y2**2.,cx,cy,cxy,ctot
				fdiag.write(str(sqrt(ctot))+'\n')
			fo.write(str(ctot)+' ')
			if j == len(dl1) -1:
				fo.write('\n')
	fo.close()



def mkcovarmr(f1,f2,cm1,cm2,cv=.125):
	d1 = open(f1).readlines()
	d2 = open(f2).readlines()
	if len(d1) != len(d2):
		return 'ERROR, mismatched files!'
	dl1 = []
	dl2 = []
	for i in range(0,len(d1)):
		dl1.append(float(d1[i].split()[1]))
		dl2.append(float(d2[i].split()[1]))
	cf1 = open(cm1).readlines()
	cf2 = open(cm2).readlines()
	if len(cf1) != len(cf2):
		return 'ERROR, mismatched covariance files'
	if len(cf1) != len(dl1):
		return 'ERROR, mismatched covariance/data files'
	fo = open('covarratiotemp.dat','w')
	for i in range(0,len(dl1)):
		x1 = dl1[i]
		y1 = dl2[i]
		lnc1 = cf1[i].split()
		lnc2 = cf2[i].split()
		for j in range(0,len(dl1)):
			x2 = dl1[j]
			y2 = dl2[j]
			cx = float(lnc1[j])
			cy = float(lnc2[j])
			cxy = sqrt(abs(cx*cy))*cv
			if cx*cy < 0:
				cxy = -1.*cxy
			ctot = cx/y1/y2+cy*x1*x2/y1**2./y2**2.-2.*cxy/y1/y2
			if i == j:
				print x1,x2,y1,y2,cx,cy,cxy,ctot
			fo.write(str(ctot)+' ')
			if j == len(dl1) -1:
				fo.write('\n')
		
		
def pkerr(sigz,n,k=.01,zl=.2,zh=.3,zmed=.25,solidang=6000.):
	d = distance(h=.7)
	sigd = d.cHz(zmed)*sigz
	pkl = pkz(om=0.3,lam=0.7,sig8=0.8,h=.7)
	pk = pkl.DELsqlin(k,zmed)/k**3.*2.*pi**2.
	vol = solidang/(360.*360./pi)*4.*pi/3.*(d.dc(zh)**3.-d.dc(zl)**3.)
	print vol
	#sigfz = exp(k**2.*sigd**2.)
	#sigsq = 1./pi*((2.*sigd*k)/erf(k*sigd))**2.*(2.*pi)**3./vol*(1.+1./(n*pk**2.))
	sigfz = 1./sqrt(pi)*((2.*sigd*k)/erf(k*sigd))
	sig = sqrt((2.*pi)**3./vol)*(1.+1./(n*pk))*sigfz
	return sig
	
def pkerrfile(file):
	f = open('nvzerr'+file+'.dat')
	fo = open('pkerrzerr'+file+'.dat','w')
	for line in f:
		ln = line.split()
		n = float(ln[1])
		zerr = float(ln[2])
		perr = pkerr(zerr,n)
		fo.write(str(zerr)+' '+str(perr)+' '+ln[0]+'\n')
	fo.close()

def mksratiofile(file,thetamin,thetamax):
	gam = .4
	fo = open('ratio3.73'+file,'w')
	while gam < .8:
		f1 = open('om'+str(gam)+'pcfzdTH'+file).readlines()
		for i in range(0,len(f1)):
			theta = float(f1[i].split()[0])
			if theta > thetamin and theta < thetamax:
				ind = i
				break
		omTH = float(f1[ind].split()[1])
		f2 = open('om'+str(gam)+'pcfzdpc'+file).readlines()
		ompc = float(f2[ind].split()[1])
		ratio = omTH/ompc
		fo.write(str(.25**gam)+' '+str(ratio)+'\n')
		gam += .025
	fo.close()

def mkfakecovarm(nbin):
	nbin = int(nbin)
	fo = open('fakecovar'+str(nbin)+'.dat','w')
	for i in range(0,nbin):
		for j in range(0,nbin):
			od = abs(j-i)
			frac = float(od)/float(nbin)
			covar = 1.-frac
			fo.write(str(covar)+' ')
			if j == nbin-1:
				fo.write('\n')
	fo.close()

def chicalcfake(nbin,md='f'):
	from numpy import zeros

	mfile = open('covarfake'+str(nbin)+'.dat').readlines()
	leng = len(mfile)

	matrix = zeros((leng,leng),'f')
	for i in range(0,leng):
		line = mfile[i].strip('\n').split()
		for j in range(0,leng):
				if i == j:
					
					matrix[i,j] = float(line[j])
					if matrix[i,j] < 0:
						return 'something wrong, negative diagonal!'
				else:
					if md != 'd':
						matrix[i,j] = float(line[j])#*.001
	print 'made matrix'		
	invCov = inv(matrix)
	mtest = inv(invCov)
	ChiSq = 0

	for i in range(0,leng):
		moda = 1.
		ra = 2.
		for j in range(0,leng):	
			modb = 1.				
			rb = 2.
			ChiSq += (ra-moda)*invCov[i,j]*(rb-modb)
			if i == j:
				print ChiSq,ra,moda,invCov[i,j],matrix[i,j],mtest[i,j]
	return ChiSq

def covarTHpc(thfile,pcfile):
	nzpcu = open('nzppcu'+pcfile+'nzDES.dat').readlines()
	nzpcd = open('nzppcd'+pcfile+'nzDES.dat').readlines()
	nzTH = open('nzTH'+thfile+'nzDES.dat').readlines()
	nzpc = []
	pcnormu = 0
	pcnormd = 0
	for i in range(0,len(nzpcu)):
		#nzpc.append(float(nzpcu[i].split()[1])*float(nzpcd[i].split()[1]))
		pcnormu += float(nzpcu[i].split()[1])*.001
		pcnormd += float(nzpcu[i].split()[1])*.001
	thnorm = 0
	for i in range(0,len(nzTH)):
		thnorm += float(nzTH[i].split()[1])*.001
	pcTHint = 0
	pcint = 0
	thint = 0
	for i in range(0,len(nzTH)):
		pcTHint += sqrt(float(nzpcu[i].split()[1])*.001/pcnormu*float(nzpcd[i].split()[1])*.001/pcnormd)*float(nzTH[i].split()[1])*.001/thnorm
		pcint += float(nzpcu[i].split()[1])*.001/pcnormu*float(nzpcd[i].split()[1])*.001/pcnormd
		thint += (float(nzTH[i].split()[1])*.001/thnorm)**2.
	print pcTHint,pcint,thint
	return pcTHint**2./pcint/thint

def davecalc(file):
	f = open(file).readlines()
	sum = 0
	norm = 0
	for i in range(0,len(f)):
		line = f[i].split()
		sum += .01*float(line[0])*float(line[1])
		norm += .01*float(line[1])
	
	return sum/norm

def mknzpc(zmin,zmax,nzf='nzDES.dat',err=.03):
	f = open(nzf).readlines()
	from random import gauss
	pcl = []
	zindmin = int(1000*zmin)
	zindmax = int(1000*zmax)
	for i in range(0,len(f)):
		pcl.append(0)
	for i in range(0,len(f)):
		indmin = 2*zindmin-i
		indmax = 2*zindmax-i
		if indmin < 0:
			indmin = 0
		if indmax < 0:
			indmax = 0
		if indmax > len(f):
			indmax = len(f)
		n = 0
		nz = float(f[i].split()[1])
		z = .001*i
		zerr = (1.+z)*err
		for j in range(indmin,indmax):
			n += sqrt(nz*float(f[j].split()[1]))
		n = int(n/1000)
		for j in range(0,n+1):
			zg = gauss(z,zerr)
			indg = int(z*1000)
			if indg >= 0 and indg < len(f):
				pcl[indg] += 1
		print i,n
	fo = open('nzpc'+str(zmin)+str(zmax)+'DES.dat','w')
	for i in range(0,len(pcl)):
		fo.write(str(i*.001)+' '+str(pcl[i])+'\n')
			
			
def zbinnerb(dNdZfile="nzDR51821Aa.dat",zmax=.9,zmin=0):
	
	file = open(dNdZfile+'.dat', 'r')
	lines = file.readlines()
	#print(len(lines))
	Nz = []
	x = []
	xmax = 0
	for i in range(1,len(lines)):
		
		parsed = float(lines[i].rstrip('\n').split()[0])
		#print(parsed)
		if parsed < zmax:
			x.append(parsed)
			if parsed > .01 and parsed > zmin:
				Nz.append(float(lines[i].rstrip('\n').split()[1]))
			else:
				Nz.append(0)
		if parsed > xmax:
			xmax = parsed
	step = x[1]-x[0]
	#print(dNdZfile,step)
	file.close()
	#print(Nz)
	
	x.reverse()
	x.append(0.)
	x.reverse()
	Nz.reverse()
	Nz.append(0.)
	Nz.reverse()

	b = array(Nz)
	
	norm = b/(sum(b)*step)
	#print(len(norm))
	f = open(dNdZfile+'normtest','w')
	for i in range(0,len(norm)):
		f.write(str(i*step)+' '+str(norm[i])+'\n')
	return norm, zmin, zmax, step
	
def mdif():
	f = open('mdifdiag0.40.6.dat').readlines()
	fd = open('nloopAnna0.40.6.dat').readlines()
	dl = []
	nl = []
	for i in range(0,20):
		dl.append(0)
		nl.append(0)
	for i in range(0,len(f)):
		d = float(f[i])
		zd = float(fd[i].split()[1])-float(fd[i].split()[0])
		zdind = int(zd*100)
		if zdind >= 0 and zdind < 20:
			nl[zdind] += 1.
			dl[zdind] += d
		else:
			print zd,zdind
	fo = open('mdifavediag0.40.6.dat','w')
	for i in range(0,20):
		fo.write(str(i/100.)+' '+str(dl[i]/nl[i])+'\n')
	fo.close()

def thdif():
	fd = open('nloopAnna0.40.6.dat').readlines()
	dl = []
	nl = []
	ml = []
	pl = []
	for i in range(0,20):
		dl.append(0)
		nl.append(0)
		ml.append(0)
		pl.append(0)
	ml[19] = .194
	ml[18] = .243
	ml[17] = .296
	ml[16] = .347
	ml[15] = .397
	ml[14] = .452
	ml[13] = .502
	ml[12] = .55
	ml[11] = .595
	ml[10] = .65
	ml[9] = .702
	ml[8] = .75
	ml[7] = .8
	ml[6] = .85
	ml[5] = .88
	ml[4] = .91
	ml[3] = .94
	ml[2] = .97
	ml[1] = .99
	ml[0] = .999
	for i in range(0,223):
		f = open('waveAnnaslice'+str(i+1)+'0.03.dat').readlines()
		dm = float(f[9].split()[1])
		 
		zd = float(fd[i].split()[1])-float(fd[i].split()[0])
		zdind = int(zd*100)
		dth = 0.002625*ml[zdind]
		d = dth - dm
		if zdind >= 0 and zdind < 20:
			nl[zdind] += 1.
			dl[zdind] += d
			pl[zdind] += abs(d)/dth
		else:
			print zd,zdind
	fo = open('thdifavediag0.40.6.dat','w')
	for i in range(0,20):
		fo.write(str(i/100.)+' '+str(dl[i]/nl[i])+' '+str(pl[i]/nl[i])+'\n')
	fo.close()

def errvzdif(theta=2.,ngal=500000.):
	fo = open('varvzdif'+str(theta)+'.dat','w')
	#zf0 = 'TH0.4960.4060.045'
	z1 = .496
	z1a = .496
	while z1 < .6 and z1a >= .4:
		z2 = z1+.006
		z2a = z1a +.006
		zf0 = 'TH'+str(z1a)+str(z2a)+'0.045'
		zf2 = 'TH'+str(z1)+str(z2)+'0.045'
		ng1 = ngal*z1/.5
		ng2 = ngal*z1a/.5
		var = cl2covcross(theta,theta,zf0,zf2,ngal=ng1,ngal2=ng2)
		zd = z1-z1a
		fo.write(str(zd)+' '+str(var[0])+'\n')
		z1 += .006
		z1a -= .006
	fo.close()
	

def mknloop4TH(zmin,zmax,step=.006):
	fo = open('nloop4th'+str(zmin)+str(zmax)+str(step)+'.dat','w')
	nbin = int((zmax-zmin)/step)
	for i in range(0,nbin):
		for j in range(i,nbin):
			z1 = zmin+step/2.+i*step
			z2 = zmin+step/2.+j*step
			fo.write(str(z1)+' '+str(z2)+'\n')
	fo.close()

def mknloop4pc(pcmin,pcmax,zmin=0,zmax=1.4,step=.006):
	fo = open('nloop4pc'+str(pcmin)+str(pcmax)+str(step)+'.dat','w')
	nbin = int((zmax-zmin)/step)
	for i in range(0,nbin):
		for j in range(i,nbin):
			z1 = zmin+step/2.+i*step
			z2 = zmin+step/2.+j*step
			zave = (z1+z2)/2.
			if zave > pcmin and zave < pcmax:
				fo.write(str(z1)+' '+str(z2)+'\n')
	fo.close()

def diagpcerr(zmin,zmax,thmin=.5,thmax=5.1,md='mock',step=.01,nbin=20):
	if md == 'mock':
		fo = open('pcdiag'+str(zmin)+str(zmax)+'mock.dat','w')
		astep = .2
		thmin = .2
		thmax = 6.1
		ng = 900000.
		rmin = .2
		rmax = .8
	if md == 'th':
		fo = open('pcdiag'+str(zmin)+str(zmax)+'th.dat','w')
		astep = thmin/2.
		ng = 1200000.
		rmin = 0
		rmax = 1.4
	ang = thmin
	#while ang < thmax:		
	for i in range(0,20):
		sum = 0
		cov = mkcov4commock(ang,ang,zmin=zmin,zmax=zmax,nsl=ng,step=step,zm=(zmin+zmax)/2.,zrmin=rmin,zrmax=rmax,zsig=.03)

		print ang,cov
		fo.write(str(ang)+' '+str(cov)+'\n')
		ang += astep
	fo.close()
	return True

def mkpccovm(zmin,zmax,thmin=.5,thmax=5.1,nbin=20,astep=.25,md='th',step=.01,zsig=.03):
	from numpy import zeros
	if md == 'mock':
		fo = open('pccov'+str(zmin)+str(zmax)+'mock.dat','w')
		astep = .2
		thmin = .2
		thmax = 6.1
		ng = 900000.
		rmin = .2
		rmax = .8
	if md == 'th':
		fo = open('pccov'+str(zmin)+str(zmax)+'nzDES'+str(zsig)+'th.dat','w')
		#astep = thmin/2.
		ng = 1200000.
		rmin = 0
		rmax = 1.4
	ang = thmin
	#nbin = int(thmax/astep)
	m = zeros((nbin,nbin),'f')
	for i in range(0,nbin):
		for j in range(i,nbin):		
			ang = i*astep+thmin
			ang2 = j*astep+thmin
			cov = mkcov4commock(ang,ang2,zmin=zmin,zmax=zmax,nsl=ng,step=step,zm=(zmin+zmax)/2.,zrmin=rmin,zrmax=rmax,zsig=.03)
			m[i][j] = cov
		print ang,m[i][i]
		

	for i in range(0,nbin):
		for j in range(0,nbin):
			cov = m[i][j]
			if cov == 0.0:
				cov = m[j][i]
			fo.write(str(cov)+' ')
		fo.write('\n')
			
	fo.close()
	return True

		

def mkcov4commock(ang,ang2,zmin=.4,zmax=.6,md='pc',nsl=500000.,step=.006,zm=.48,zrmin=.2,zrmax=.8,zsig=.03):
	from nz import mknzdTH
	from numpy import zeros
	from time import time
	t0 = time()
	#lf = open('nloop.3.z.0.50.dz.0.03.zreal.0.40-0.60.dat').readlines()
	mknloop4pc(zmin,zmax,zmin=zrmin,zmax=zrmax,step=step)
	if md == 'pc':
		lf = open('nloop4pc'+str(zmin)+str(zmax)+str(step)+'.dat').readlines()
	else:
		lf = open('nloop.3.z.0.50.dz.0.03.zreal.0.40-0.60.dat').readlines()
	try:
		clf = open('clTHzsig'+str(zsig)+str(zm)+str(zm+step)+'nzDESnodamp.dat')
	except:
		mknzdTH('nzDES.dat',zm,zm+step,zsig=zsig)
		for i in range(0,30):
			mkphi(i,'THzsig'+str(zsig)+str(zm)+str(zm+step)+'nzDES')
		mkclfile('THzsig'+str(zsig)+str(zm)+str(zm+step)+'nzDES')
		clf = open('clTHzsig'+str(zsig)+str(zm)+str(zm+step)+'nzDESnodamp.dat')
	diagerrl = []
	diagnsl = []
	zerr = zsig*(1.+zm)+.03
	gl = []
	ml = []
	norm = 0
	leng = len(lf)
	mat = zeros((leng,leng),'f')
	for i in range(0,1500):
		g = 1./sqrt(2*pi)/zerr*exp(-1.*(i/1000.-zm)**2./(2.*zerr**2.))
		norm += g
		gl.append(g)
	m0 = 0
	for i in range(0,1500):
		m0 += gl[i]**2.*norm
	zdl = []
	for i in range(0,100):
		zdl.append(0)
	mm = 0
	try:
		mf = open('m'+str(zmin)+str(zmax)+str(zrmin)+str(zrmax)+'.dat').readlines()
	except:
		mm = 1
		mf = open('m'+str(zmin)+str(zmax)+str(zrmin)+str(zrmax)+'.dat','w')
	#print mm
	ng = nsl

	covm1 = cl2cov(ang,ang2,zf='THzsig0.03'+str(zm)+str(zm+step)+'nzDES',fs=.125,ngal=ng,md='n',maxl=1000,bins=.25,mult=1.)[0]
#	print covm1
	covns = cl2cov(ang,ang2,zf='THzsig0.03'+str(zm)+str(zm+step)+'nzDES',fs=.125,ngal=10000000000.,md='n',maxl=1000,bins=.25,mult=1.)[0]/2.
#	print covns
	for k in range(0,len(lf)):
		lni = lf[k].split()
		z1 = float(lni[0])
		z2 = float(lni[1])
		#print z1,z2
		m = 0
		inds1 = int((zm-z1)*1000)
		inds2 = int((zm-z2)*1000)

		if mm == 1:
			for i in range(0,1500):
				ind1 = i+inds1
				ind2 = i+inds2
				if ind1 >= 0 and ind1 < 1500 and ind2 >= 0 and ind2 < 1500:
					m += gl[ind1]*gl[ind2]*norm
			m = m/m0
			ml.append(m)
			mf.write(str(m)+'\n')
			#print z1,z2,m
			if m > 1.:
				print m
		else:
			m = float(mf[k])

		zdind = abs((inds2-inds1)/10)
		if zdind < 100:
			zdl[zdind] = m

		cov = (covm1+cl2cov(ang,ang2,zf='THzsig0.03'+str(zm)+str(zm+step)+'nzDES',fs=.125,ngal=10000000000.,md='n',maxl=1000,bins=.25,mult=m)[0])/2.
		mat[k][k] = cov
		diagerrl.append(cov)
#	print zdl
	t1 = time()
	#print t1-t0
	for k in range(0,len(lf)):
	#for k in range(0,1):
		lni = lf[k].split()
 		#print k
		z1i = float(lni[0])
		z2i = float(lni[1])
		for j in range(0,len(lf)):			
			lnj = lf[j].split()
			z1j = float(lnj[0])
			z2j = float(lnj[1])
				
			if j != k:
						
				m1 = 0
				m2 = 0
				inds1i = int((zm-z1i)*1000)
				inds2i = int((zm-z2i)*1000)
				inds1j = int((zm-z1j)*1000)
				inds2j = int((zm-z2j)*1000)
				zdind1 = abs((inds1j-inds1i)/10)
				zdind2 = abs((inds2j-inds2i)/10)
				if zdind1 < 100 and zdind2 < 100:
					m1 = zdl[zdind1]
					m2 = zdl[zdind2]
					mt = m1*m2
				else:
					mt = 0
				#cov1 = mt*sqrt(diagnsl[k]*diagnsl[j])
				cov1 = mt*covns
				mat[k][j] = cov1

	sum = 0
	#ft = open('temp.dat','w')
	for i in range(0,len(mat)):
		for j in range(0,len(mat)):
			sum += mat[i][j]
		#	ft.write(str(mat[i][j])+' ')
		#ft.write('\n')
	#print time()-t1
	return sum/float(len(mat))**2.

def mkcovpcAnna(bins,nmock=60,res='n'):
	from numpy import zeros
	thl = []
	pcl = []
	for i in range(0,30):
		pcl.append(0)
	for i in range(1,nmock+1):
		if res == 'n':
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		else:
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			pcl[j] += float(line[2])/float(nmock)
	covm = zeros((30,30),'f')
	for i in range(1,nmock+1):
		if res == 'n':
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		else:
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			linej = f[j].split()
			for k in range(0,30):
				linek = f[k].split()
				covm[k][j] += (pcl[j]-float(linej[2]))*(pcl[k]-float(linek[2]))/float(nmock)
	if res == 'n':
		fo = open('covpcAnna'+str(bins)+'.dat','w')
	else:
		fo = open('covpcAnna0.4-0.6'+str(bins)+'.dat','w')
	for i in range(0,30):
		for j in range(0,30):
			fo.write(str(covm[i][j])+' ')
		fo.write('\n')
	fo.close()

def mdifpc():
	fm = open('covpcAnna0.03.dat').readlines()
	ft = open('pccov0.4750.525mock.dat').readlines()
	fo = open('mdifpc.dat','w')
	for i in range(0,30):
		ang1 = .2+i*.2
		for j in range(0,30):
			ang2 = .2+j*.2
			d = float(fm[i].split()[j])-float(ft[i].split()[j])
			fd = d/float(ft[i].split()[j])
			fo.write(str(ang1)+' '+str(ang2)+' '+str(fd)+' '+str(d)+'\n')
	fo.close()

def mkaveAnna(bins,nmock=60,res='n'):
	thl = []
	pcl = []
	for i in range(0,30):
		thl.append(0)
		pcl.append(0)
	for i in range(1,nmock+1):
		if res == 'n':
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		else:
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			thl[j] += float(line[1])/float(nmock)
			pcl[j] += float(line[2])/float(nmock)
	stdpc = []
	stdth = []
	cov = []
	for i in range(0,30):
		stdpc.append(0)
		stdth.append(0)
		cov.append(0)
	for i in range(1,nmock+1):
		if res == 'n':
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		else:
			f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			stdth[j] += (thl[j]-float(line[1]))**2./float(nmock)
			stdpc[j] += (pcl[j]-float(line[2]))**2./float(nmock)
			cov[j] += (thl[j]-float(line[1]))*(pcl[j]-float(line[2]))/float(nmock)
	if res == 'n':
		fo = open('waveAnna'+str(bins)+'125.dat','w')
	else:
		fo = open('waveAnna0.4-0.6'+str(bins)+'.dat','w')
	for i in range(0,30):
		fo.write(f[i].split()[0]+' '+str(thl[i])+' '+str(stdth[i])+' '+str(pcl[i])+' '+str(stdpc[i])+' '+str(cov[i]/sqrt(stdpc[i]*stdth[i]))+'\n')
	fo.close()

def mkaveAnnaN(bins,nmock=82):
	thl = []
	pcl = []
	for i in range(0,30):
		thl.append(0)
		pcl.append(0)
	for i in range(1,nmock+1):
		if i != 42:
			f = open('Output3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-0.03').readlines()
			for j in range(0,30):
				line = f[j].split()
				thl[j] += float(line[1])/float(nmock-1)
				pcl[j] += float(line[2])/float(nmock-1)
	stdpc = []
	stdth = []
	cov = []
	for i in range(0,30):
		stdpc.append(0)
		stdth.append(0)
		cov.append(0)
	for i in range(1,nmock+1):
		if i != 42:
			f = open('Output3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-0.03').readlines()
			for j in range(0,30):
				line = f[j].split()
				stdth[j] += (thl[j]-float(line[1]))**2./float(nmock-1)
				stdpc[j] += (pcl[j]-float(line[2]))**2./float(nmock-1)
				cov[j] += (thl[j]-float(line[1]))*(pcl[j]-float(line[2]))/float(nmock-1)
	fo = open('waveAnnaN'+str(bins)+'.dat','w')
	for i in range(0,30):
		fo.write(f[i].split()[0]+' '+str(thl[i])+' '+str(stdth[i])+' '+str(pcl[i])+' '+str(stdpc[i])+' '+str(cov[i]/sqrt(stdpc[i]*stdth[i]))+'\n')
	fo.close()


def mkaveAnnaslices(bins,nmock=60,col=1):
	thl = []
	for i in range(0,30):
		thl.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			thl[j] += float(line[col])/float(nmock)
	stdth = []
	for i in range(0,30):
		stdth.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			stdth[j] += (thl[j]-float(line[col]))**2./float(nmock)
	fo = open('waveAnnaslice'+str(col)+str(bins)+'.dat','w')
	for i in range(0,30):
		fo.write(f[i].split()[0]+' '+str(thl[i])+' '+str(stdth[i])+'\n')
	fo.close()

def mkaveAnnaslices1deg(bins,nmock=60,row=1):
	f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-1_zm-0.50_dz-'+str(bins))
	nc = len(f.readline().split())
	th= []
	for i in range(1,nc):
		th.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		line = f[row].split()
		for j in range(1,nc):			
			th[j-1] += float(line[j])/float(nmock)
	stdth = []
	for i in range(1,nc):
		stdth.append(0)
	from numpy import zeros
	leng = nc-1
	covarm = zeros((leng,leng),'f')
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		line = f[row].split()
		for j in range(1,nc):
			for k in range(1,nc):
			#stdth[j-1] += (th[j-1]-float(line[j]))**2./float(nmock)
				covarm[j-1,k-1] += (th[j-1]-float(line[j]))*(th[k-1]-float(line[k]))/float(nmock)
	tht = 0
	stdtht = 0
	fo = open('mockslicescovarm.dat','w')
	for i in range(0,len(th)):
		tht += th[i]/float(nc)
		for j in range(0,len(th)):
			stdtht += covarm[i][j]
			if j == leng-1:
				fo.write(str(covarm[i][j])+'\n')
			else:
				fo.write(str(covarm[i][j])+' ')
	fo.close()
	return tht,stdtht/float(nc)/float(nc)

def mkaveAnnasliceszr1deg(zmin,zmax,bins='0.03',nmock=60,row=1):
	f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-1_zm-0.50_dz-'+str(bins))
	th= []
	lf = open('nloop.3.z.0.50.dz.0.03.zreal.0.40-0.60.dat').readlines()
	nc = 0
	for i in range(0,len(lf)):
		ln = lf[i].split()
		z1 = float(ln[0])
		z2 = float(ln[1])
		if z1 > zmin and z2 < zmax and z1 < zmax and z2 > zmin:
			nc += 1
	print nc
	for i in range(0,nc):
		th.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		line = f[row].split()
		cind = 0
		for j in range(1,len(line)):
			zind = j-1
			lnz = lf[zind].split()
			z1 = float(lnz[0])
			z2 = float(lnz[1])
			if z1 > zmin and z2 < zmax and z1 < zmax and z2 > zmin:
				th[cind] += float(line[j])/float(nmock)
				cind += 1
	stdth = []
	for i in range(0,nc):
		stdth.append(0)
	from numpy import zeros
	leng = nc
	covarm = zeros((leng,leng),'f')
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		line = f[row].split()
		cjind = 0
		for j in range(1,len(line)):
			zindj = j-1
			lnzj = lf[zindj].split()
			z1j = float(lnzj[0])
			z2j = float(lnzj[1])
			if z1j > zmin and z2j < zmax and z1j < zmax and z2j > zmin:
				ckind = 0
				for k in range(1,len(line)):
					zindk = k-1
					lnzk = lf[zindk].split()
					z1k = float(lnzk[0])
					z2k = float(lnzk[1])
					if z1k > zmin and z2k < zmax and z1k < zmax and z2k > zmin:
	
	#stdth[j-1] += (th[j-1]-float(line[j]))**2./float(nmock)
						covarm[cjind,ckind] += (th[cjind]-float(line[j]))*(th[ckind]-float(line[k]))/float(nmock)
						ckind += 1
				cjind += 1
	tht = 0
	stdtht = 0
	fo = open('mockslicescovarm'+str(zmin)+str(zmax)+'.dat','w')
	for i in range(0,len(th)):
		tht += th[i]/float(nc)
		for j in range(0,len(th)):
			stdtht += covarm[i][j]
			if j == leng-1:
				fo.write(str(covarm[i][j])+'\n')
			else:
				fo.write(str(covarm[i][j])+' ')
	fo.close()
	return tht,stdtht/float(nc)/float(nc)

	
def mkaveAnna2slicescov(bins,nmock=60,col1=1,col2=2):
	thl1 = []
	thl2 = []
	for i in range(0,30):
		thl1.append(0)
		thl2.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			thl1[j] += float(line[col1])/float(nmock)
			thl2[j] += float(line[col2])/float(nmock)
	stdth1 = []
	stdth2 = []
	cov = []
	for i in range(0,30):
		stdth1.append(0)
		stdth2.append(0)
		cov.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/SliceCorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			stdth1[j] += (thl1[j]-float(line[col1]))**2./float(nmock)
			stdth2[j] += (thl2[j]-float(line[col2]))**2./float(nmock)
			cov[j] += (thl1[j]-float(line[col1]))*(thl2[j]-float(line[col2]))/float(nmock)
	fo = open('waveAnnaslicecov'+str(col1)+str(col2)+str(bins)+'.dat','w')
	for i in range(0,30):
		fo.write(f[i].split()[0]+' '+str(thl1[i])+' '+str(thl2[i])+' '+str(stdth1[i])+' '+str(stdth2[i])+' '+str(cov[i]/sqrt(stdth1[i]*stdth2[i]))+'\n')
	fo.close()


def mkaveAnnacon(bins,nmock=60,con='0.4-0.6'):
	pcl = []
	for i in range(0,30):
		pcl.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			pcl[j] += float(line[1])/float(nmock)
	stdpc = []
	for i in range(0,30):
		stdpc.append(0)
	for i in range(1,nmock+1):
		f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz0.4-0.6_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		for j in range(0,30):
			line = f[j].split()
			stdpc[j] += (pcl[j]-float(line[1]))**2./float(nmock)
	fo = open('waveAnna0.4-0.6'+str(bins)+'.dat','w')
	for i in range(0,30):
		fo.write(f[i].split()[0]+' '+str(pcl[i])+' '+str(stdpc[i])+'\n')
	fo.close()


def mkhistAnna(bins,nmock=60,ind=10):
	thave = 0
	pcave = 0
	for i in range(1,nmock+1):
		f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		line = f[ind].split()
		thave += float(line[1])/float(nmock)
		pcave += float(line[2])/float(nmock)
	fhth = open('thhist.dat','w')
	fhpc = open('pchist.dat','w')
	thl = []
	pcl =  []
	nbin = 20
	for i in range(0,nbin):
		thl.append(0)
		pcl.append(0)
	thm = .001
	pcm = .0006
	
	for i in range(1,nmock+1):
		f = open('Outputphoto3/CorrFile_photoz_zspace_dNdz_Lbox-7680.00_Nmocks-125_rmean-1284.96_width-750.10_sigz-0.045_imock-'+str(i)+'_zm-0.50_dz-'+str(bins)).readlines()
		line = f[ind].split()
		dth = thave - float(line[1])
		thind = int((dth+thm)*float(nbin)/(2.*thm))
		if thind >= 0:
			if thind < nbin:
				thl[thind] += 1
			else:
				thl[nbin-1] += 1
		else:
			thl[0] += 1
		dpc = pcave - float(line[2])
		pcind = int((dpc+pcm)*float(nbin)/(2.*pcm))
		if pcind >= 0:
			if pcind < nbin:
				pcl[pcind] += 1
			else:
				pcl[nbin-1] += 1
		else:
			pcl[0] += 1
		#print dth,dpc
	for i in range(0,nbin):
		fhth.write(str(-1.*thm+i*2.*thm/float(nbin))+' '+str(thl[i])+'\n')
		fhpc.write(str(-1.*pcm+i*2.*pcm/float(nbin))+' '+str(pcl[i])+'\n')
	return True