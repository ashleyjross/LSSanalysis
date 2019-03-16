Hdir = '/Users/ashleyross/Dropbox/eBOSS/'
import numpy as np
from numpy import zeros,dot
from numpy.linalg import pinv
from math import *

def P2(mu):
	return .5*(3.*mu**2.-1.)
	
def P4(mu):
	return .125*(35.*mu**4.-30.*mu**2.+3.)

def findPolya(H,ci,d):
	ht = H.transpose()
	onei = pinv(dot(dot(H,ci),ht))
	comb = dot(dot(onei,H),ci)
	return dot(comb,d)

class baofit3D_ellFull_1cov:
	def __init__(self,dv,ic,mod,rl):
		self.xim = dv
		xinl = []
		xisl = []
		xinli = []
		xisli = []
		self.rl = rl
		s = 0
		m2 = 1.
		#self.B0fac = (1.+2/3.*B0+.2*B0**2.)
		#self.B2fac = (4/3.*B0+4/7.*B0**2.)
		#self.B4fac = 8/35.*B0**2.

		self.nbin = len(self.rl)
		
		print(self.nbin)
		print(self.xim)		
		print(self.rl)
		self.invt = ic
		if self.nbin != len(self.invt):
			print('vector matrix mismatch!')
			return 'vector matrix mismatch!'

		self.ximodmin = 10.

		self.x0 = [] #empty lists to be filled for model templates
		self.x2 = []
		self.x4 = []
		self.x0sm = []
		self.x2sm = []
		self.x4sm = []
		
		mf0 = open('BAOtemplates/xi0'+mod).readlines()
		mf2 = open('BAOtemplates/xi2'+mod).readlines()
		mf4 = open('BAOtemplates/xi4'+mod).readlines()
		mf0sm = open('BAOtemplates/xi0sm'+mod).readlines()
		mf2sm = open('BAOtemplates/xi2sm'+mod).readlines()
		mf4sm = open('BAOtemplates/xi4sm'+mod).readlines()
		for i in range(0,len(mf0)):
			ln0 = mf0[i].split()
			self.x0.append(2.1*(float(ln0[1])))
			self.x0sm.append(2.1*float(mf0sm[i].split()[1]))
			ln2 = mf2[i].split()
			m2 = 2.1*(float(ln2[1]))
			self.x2.append(m2)
			self.x2sm.append(2.1*float(mf2sm[i].split()[1]))
			ln4 = mf4[i].split()
			m4 = 2.1*(float(ln4[1]))
			self.x4.append(m4)
			self.x4sm.append(2.1*float(mf4sm[i].split()[1]))
		self.at = 1.
		self.ar = 1.
		self.b0 = 1.
		self.b2 = 1.
		self.b4 = 1.
		r = 20.
		self.H = zeros((6,self.nbin))
		for i in range(0,self.nbin):
			if i < self.nbin/2:
				self.H[0][i] = 1.
				self.H[1][i] = 1./self.rl[i]
				self.H[2][i] = 1./self.rl[i]**2.
			if i >= self.nbin/2:
				self.H[3][i] = 1.
				self.H[4][i] = 1./self.rl[i]
				self.H[5][i] = 1./self.rl[i]**2.
		
	def wmod(self,r,sp=1.):
		self.sp = sp
		sum = 0
		sum2 = 0
		nmu = 100
		dmu = 1./float(nmu)
		for i in range(0,nmu):
			mu = i*dmu+dmu/2.
			al = sqrt(mu**2.*self.ar**2.+(1.-mu**2.)*self.at**2.)
			mup = mu*self.ar/al
			rp = r*al
			#ximu = self.b0*self.lininterp(self.x0,rp)+self.b2*P2(mup)*self.lininterp(self.x2,rp)+self.b4*P4(mup)*self.lininterp(self.x4,rp)
			ximu = self.lininterp(self.x0,rp)+P2(mup)*self.lininterp(self.x2,rp)+P4(mup)*self.lininterp(self.x4,rp)
			sum += ximu
			#sum2 += P2(mu)*ximu
			sum2 += mu**2.*ximu
		return dmu*sum,1.5*dmu*sum2#3.*dmu*sum2

	def wmodsm(self,r,sp=1.):
		self.sp = sp
		sum = 0
		sum2 = 0
		nmu = 100
		dmu = 1./float(nmu)
		for i in range(0,nmu):
			mu = i*dmu+dmu/2.
			al = sqrt(mu**2.*self.ar**2.+(1.-mu**2.)*self.at**2.)
			mup = mu*self.ar/al
			rp = r*al
			ximu = self.lininterp(self.x0sm,rp)+P2(mup)*self.lininterp(self.x2sm,rp)+P4(mup)*self.lininterp(self.x4sm,rp)
			sum += ximu
			#sum2 += P2(mu)*ximu
			sum2 += mu**2.*ximu
		return dmu*sum,1.5*dmu*sum2#3.*dmu*sum2

	def wmodW(self,r,sp=1.):
		self.sp = sp
		sum = 0
		sum2 = 0
		nmu = 100
		dmu = 1./float(nmu)
		nspl = int(nmu*self.Wsp)
		for i in range(0,nspl):
			mu = i*dmu+dmu/2.
			al = sqrt(mu**2.*self.ar**2.+(1.-mu**2.)*self.at**2.)
			mup = mu*self.ar/al
			rp = r*al
			#ximu = self.b0*self.lininterp(self.x0,rp)+self.b2*P2(mup)*self.lininterp(self.x2,rp)+self.b4*P4(mup)*self.lininterp(self.x4,rp)
			ximu = self.lininterp(self.x0,rp)+P2(mup)*self.lininterp(self.x2,rp)+P4(mup)*self.lininterp(self.x4,rp)
			sum += ximu
		for i in range(nspl,nmu):
			mu = i*dmu+dmu/2.
			al = sqrt(mu**2.*self.ar**2.+(1.-mu**2.)*self.at**2.)
			mup = mu*self.ar/al
			rp = r*al
			#ximu = self.b0*self.lininterp(self.x0,rp)+self.b2*P2(mup)*self.lininterp(self.x2,rp)+self.b4*P4(mup)*self.lininterp(self.x4,rp)
			ximu = self.lininterp(self.x0,rp)+P2(mup)*self.lininterp(self.x2,rp)+P4(mup)*self.lininterp(self.x4,rp)
			sum2 += ximu
		n0 = 1./float(self.Wsp)
		n2 = 1./(1.-self.Wsp)
		return dmu*sum*n0,dmu*sum2*n2

			
	def lininterp(self,f,r):
		indd = int((r-self.ximodmin)/self.sp)
		indu = indd + 1
		fac = (r-self.ximodmin)/self.sp-indd
		if fac > 1.:
			print('BAD FAC in wmod')
			return 'ERROR'
		if indu >= len(f)-1:
			return 0
		return f[indu]*fac+(1.-fac)*f[indd]	

	def mkxi(self):
		self.xia = []
		for i in range(0,len(self.rl)//2):
			xi0,xi2 = self.wmod(self.rl[i])
			self.xia.append(xi0)
			#self.xi2a.append(xi2)
		for i in range(len(self.rl)//2,self.nbin):
			xi0,xi2 = self.wmod(self.rl[i])
			self.xia.append(xi2)
			#self.xi2a.append(xi2)
		return True	

	def mkxism(self):
		self.xiasm = []
		for i in range(0,len(self.rl)//2):
			xi0,xi2 = self.wmodsm(self.rl[i])
			self.xiasm.append(xi0)
			#self.xi2a.append(xi2)
		for i in range(len(self.rl)//2,self.nbin):
			xi0,xi2 = self.wmodsm(self.rl[i])
			self.xiasm.append(xi2)
			#self.xi2a.append(xi2)
		return True	

	def mkxiW(self):
		self.xi0a = []
		self.xi2a = []
		for i in range(0,len(self.rl)):
			xi0,xi2 = self.wmodW(self.rl[i])
			self.xi0a.append(xi0)
			self.xi2a.append(xi2)
		return True	
		
	def chi_templ_alphfXX(self,list,wo='n',fw='',v='n'):
		from time import time
		t = time()
		BB = list[0]
		if BB < 0:
			return 1000
		A0 = list[1]
		A1 = list[2]
		A2 = list[3]
		Beta = list[4]
		if Beta < 0:
			return 1000
		A02 = list[5]
		A12 = list[6]
		A22 = list[7]
		#self.b0 = BB*(1.+2/3.*Beta+.2*Beta**2.)/self.B0fac
		#self.b2 = BB*(4/3.*Beta+4/7.*Beta**2.)/self.B2fac
		#self.b4 = 8/35.*BB*Beta**2./self.B4fac
		modl = []
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin//2):
			r = self.rl[i]
			
			mod0 = BB*self.xia[i]+A0+A1/r+A2/r**2.
			modl.append(mod0)
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod0)+'\n')
		
		for i in range(self.nbin//2,self.nbin):
			r = self.rl[i]
			mod2 = 5.*(Beta*self.xia[i]-BB*0.5*self.xia[i-self.nbin//2])+A02+A12/r+A22/r**2.
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod2)+'\n')			
			modl.append(mod2)
		if wo == 'y':
			fo.close()		
			

		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.xim[i]-modl[i])
		chit = dot(dot(dl,self.invt),dl)
		if v == 'y':
			print(dl,chit)
		BBfac = (log(BB/self.BB)/self.Bp)**2.
		#Btfac = ((Beta-self.B0)/self.Bt)**2.
		Btfac = (log(Beta/self.B0)/self.Bt)**2.
		return chit+BBfac+Btfac

	def chi_templ_alphfXX_an(self,list,wo='n',fw='',v='n'):
		from time import time
		t = time()
		BB = list[0]
		if BB < 0:
			return 1000
		#A0 = list[1]
		#A1 = list[2]
		#A2 = list[3]
		Beta = list[1]
		if Beta < 0:
			return 1000
		#A02 = list[5]
		#A12 = list[6]
		#A22 = list[7]
		#self.b0 = BB*(1.+2/3.*Beta+.2*Beta**2.)/self.B0fac
		#self.b2 = BB*(4/3.*Beta+4/7.*Beta**2.)/self.B2fac
		#self.b4 = 8/35.*BB*Beta**2./self.B4fac
		nrbin = self.nbin//2
		modl = []
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		pv = []
		for i in range(0,self.nbin//2):
			pv.append(self.xim[i]-BB*self.xia[i])
		for i in range(self.nbin//2,self.nbin):
			pv.append(self.xim[i]-(5.*(Beta*self.xia[i]-BB*0.5*self.xia[i-self.nbin//2])))
		 
		Al = findPolya(self.H,self.invt,pv)
		A0,A1,A2,A02,A12,A22 = Al[0],Al[1],Al[2],Al[3],Al[4],Al[5]
		for i in range(0,self.nbin//2):
			r = self.rl[i]
			
			mod0 = BB*self.xia[i]+A0+A1/r+A2/r**2.
			modl.append(mod0)
			if wo == 'y':
				mod0sm = BB*self.xiasm[i]+A0+A1/r+A2/r**2.
				fo.write(str(self.rl[i])+' '+str(mod0)+' '+str(mod0sm)+'\n')
		
		for i in range(self.nbin//2,self.nbin):
			r = self.rl[i]
			mod2 = 5.*(Beta*self.xia[i]-BB*0.5*self.xia[i-nrbin])+A02+A12/r+A22/r**2.
			if wo == 'y':
				mod2sm = 5.*(Beta*self.xiasm[i]-BB*0.5*self.xiasm[i-nrbin])+A02+A12/r+A22/r**2.
				fo.write(str(self.rl[i])+' '+str(mod2)+' '+str(mod2sm)+'\n')			
			modl.append(mod2)
		if wo == 'y':
			fo.close()		
			

		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.xim[i]-modl[i])	
		chit = dot(dot(dl,self.invt),dl)
		if v == 'y':
			print(dl,chit)
		BBfac = (log(BB/self.BB)/self.Bp)**2.
		#Btfac = ((Beta-self.B0)/self.Bt)**2.
		Btfac = (log(Beta/self.B0)/self.Bt)**2.
		return chit+BBfac+Btfac


	def chi_templ_alphfXXnA(self,list,wo='n'):
		from time import time
		t = time()
		BB = list[0]
		if BB < 0:
			return 1000
		#A0 = list[1]
		#A1 = list[2]
		#A2 = list[3]
		Beta = list[1]
		if Beta < 0:
			return 1000
		#A02 = list[5]
		#A12 = list[6]
		#A22 = list[7]
		#self.b0 = BB*(1.+2/3.*Beta+.2*Beta**2.)/self.B0fac
		#self.b2 = BB*(4/3.*Beta+4/7.*Beta**2.)/self.B2fac
		#self.b4 = 8/35.*BB*Beta**2./self.B4fac
		modl = []
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin/2):
			r = self.rl[i]
			
			mod0 = BB*self.xia[i]#+A0+A1/r+A2/r**2.
			modl.append(mod0)
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod0)+'\n')
		
		for i in range(self.nbin/2,self.nbin):
			r = self.rl[i]
			mod2 = 5.*(Beta*self.xia[i]-BB*0.5*self.xia[i-self.nbin/2])#+A02+A12/r+A22/r**2.
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod2)+'\n')			
			modl.append(mod2)
		if wo == 'y':
			fo.close()		
			

		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.xim[i]-modl[i])
		chit = dot(dot(dl,self.invt),dl)
		BBfac = (log(BB/self.BB)/self.Bp)**2.
		#Btfac = ((Beta-self.B0)/self.Bt)**2.
		Btfac = (log(Beta/self.B0)/self.Bt)**2.
		return chit+BBfac+Btfac

	def chi_templ_alphfXXnA0(self,list,wo='n'):
		from time import time
		t = time()
		BB = list[0]
		if BB < 0:
			return 1000
		#A0 = list[1]
		#A1 = list[2]
		#A2 = list[3]
		Beta = list[1]
		if Beta < 0:
			return 1000
		A02 = list[2]
		A12 = list[3]
		A22 = list[4]
		#self.b0 = BB*(1.+2/3.*Beta+.2*Beta**2.)/self.B0fac
		#self.b2 = BB*(4/3.*Beta+4/7.*Beta**2.)/self.B2fac
		#self.b4 = 8/35.*BB*Beta**2./self.B4fac
		modl = []
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin/2):
			r = self.rl[i]
			
			mod0 = BB*self.xia[i]#+A0+A1/r+A2/r**2.
			modl.append(mod0)
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod0)+'\n')
		
		for i in range(self.nbin/2,self.nbin):
			r = self.rl[i]
			mod2 = 5.*(Beta*self.xia[i]-BB*0.5*self.xia[i-self.nbin/2])+A02+A12/r+A22/r**2.
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod2)+'\n')			
			modl.append(mod2)
		if wo == 'y':
			fo.close()		
			

		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.xim[i]-modl[i])
		chit = dot(dot(dl,self.invt),dl)
		BBfac = (log(BB/self.BB)/self.Bp)**2.
		#Btfac = ((Beta-self.B0)/self.Bt)**2.
		Btfac = (log(Beta/self.B0)/self.Bt)**2.
		return chit+BBfac+Btfac

	def chi_templ_alphfXXnA2(self,list,wo='n'):
		from time import time
		t = time()
		BB = list[0]
		if BB < 0:
			return 1000
		A0 = list[1]
		A1 = list[2]
		A2 = list[3]
		Beta = list[4]
		if Beta < 0:
			return 1000
		#A02 = list[5]
		#A12 = list[6]
		#A22 = list[7]
		#self.b0 = BB*(1.+2/3.*Beta+.2*Beta**2.)/self.B0fac
		#self.b2 = BB*(4/3.*Beta+4/7.*Beta**2.)/self.B2fac
		#self.b4 = 8/35.*BB*Beta**2./self.B4fac
		modl = []
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin/2):
			r = self.rl[i]
			
			mod0 = BB*self.xia[i]+A0+A1/r+A2/r**2.
			modl.append(mod0)
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod0)+'\n')
		
		for i in range(self.nbin/2,self.nbin):
			r = self.rl[i]
			mod2 = 5.*(Beta*self.xia[i]-BB*0.5*self.xia[i-self.nbin/2])#+A02+A12/r+A22/r**2.
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod2)+'\n')			
			modl.append(mod2)
		if wo == 'y':
			fo.close()		
			

		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.xim[i]-modl[i])
		chit = dot(dot(dl,self.invt),dl)
		BBfac = (log(BB/self.BB)/self.Bp)**2.
		#Btfac = ((Beta-self.B0)/self.Bt)**2.
		Btfac = (log(Beta/self.B0)/self.Bt)**2.
		return chit+BBfac+Btfac


	def chi_templ_alphfXXW(self,list,wo='n'):
		from time import time
		t = time()
		BB = list[0]
		if BB < 0:
			return 1000
		A0 = list[1]
		A1 = list[2]
		A2 = list[3]
		Beta = list[4]
		if Beta < 0:
			return 1000
		A02 = list[5]
		A12 = list[6]
		A22 = list[7]
		#self.b0 = BB*(1.+2/3.*Beta+.2*Beta**2.)/self.B0fac
		#self.b2 = BB*(4/3.*Beta+4/7.*Beta**2.)/self.B2fac
		#self.b4 = 8/35.*BB*Beta**2./self.B4fac
		mod0l = []
		mod2l = []
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin):
			r = self.rl[i]
			
			mod0 = BB*self.xi0a[i]+A0+A1/r+A2/r**2.
			
			mod2 = Beta*self.xi2a[i]+A02+A12/r+A22/r**2.	

			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod0)+' '+str(mod2)+'\n')
			mod0l.append(mod0)
			mod2l.append(mod2)
		if wo == 'y':
			fo.close()		
			

		ChiSq = 0
		chid = 0
		Chioff = 0
		dl0 = []
		dl2 = []
		for i in range(0,self.nbin):
			dl0.append(self.xim[i][0]-mod0l[i])
			dl2.append(self.xim[i][1]-mod2l[i])
		chit0 = dot(dot(dl0,self.invt0),dl0)
		chit2 = dot(dot(dl2,self.invt2),dl2)
		BBfac = (log(BB/self.BB)/self.Bp)**2.
		Btfac = (log(Beta/self.Bta)/self.Bt)**2.
		return chit0+chit2+BBfac+Btfac



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


def Xism_arat_1C_an(dv,icov,rl,mod,dvb,icovb,rlb,B0=1.,spat=.003,spar=.006,mina=.8,maxa=1.2,nobao='n',Bp=.4,Bt=.4,meth='Powell',bs=8,fout=''):
	from time import time
	import numpy	
	#from optimize import fmin
	from random import gauss, random
	from scipy.optimize import minimize 
	print('try meth = "Nelder-Mead" if does not work or answer is weird')
	bb = baofit3D_ellFull_1cov(dvb,icovb,mod,rlb) #initialize for bias prior
	b = baofit3D_ellFull_1cov(dv,icov,mod,rl) #initialize for fitting
	b.B0 = B0
	b.Bt = Bt	
	
	bb.Bp = 100.
	bb.BB = 1.
	bb.B0 = B0
	bb.Bt = 100.
	bb.mkxi()
	b.bins = bs
	t = time()	
	
	B = .1
	chiBmin = 1000
	while B < 4.:
		chiB = bb.chi_templ_alphfXX((B,0,0,0,1,0,0,0))
		if chiB < chiBmin:
			chiBmin = chiB			
			BB = B
			print(BB,chiBmin)
		B += .01	
 	#bb.chi_templ_alphfXX((BB,0,0,0,1.,0,0,0),v='y')	
	print(BB,chiBmin)
	b.BB = BB
	b.B0 = BB		
	b.Bp= Bp
	b.Bt = Bt
	fo = open(Hdir+'2Dbaofits/arat'+fout+'1covchi.dat','w')
	fg = open(Hdir+'2Dbaofits/arat'+fout+'1covchigrid.dat','w')
	chim = 1000
	nar = int((maxa-mina)/spar)
	nat = int((maxa-mina)/spat)
	grid = numpy.zeros((nar,nat))
	pt = 0
	A0 = 0
	A1 = 0
	A2 = 0
	A02 = 0
	A12 = 0
	A22 = 0
	fac = (998.-float(len(dv)))/999.

	for i in range(0,nar):		
		b.ar = mina+spar*i+spar/2.
		print(b.ar)
		for j in range(0,nat):
			b.at = mina+spat*j+spat/2.
			b.mkxi()
			inl = (B,B0)
			(B,B0) = minimize(b.chi_templ_alphfXX_an,inl,method=meth,options={'disp': False}).x
			chi = b.chi_templ_alphfXX_an((B,B0))*fac
			grid[i][j] = chi
			fo.write(str(b.ar)+' '+str(b.at)+' '+str(chi)+'\n')
			fg.write(str(chi)+' ')
			if chi < chim:
				print(b.ar,b.at,chi)	
				chim = chi
				alrm = b.ar
				altm = b.at
				Bm = B
				Betam = B0
		fg.write('\n')
	b.ar = alrm
	b.at = altm	
	b.mkxi()
	b.mkxism()
	chi = b.chi_templ_alphfXX_an((Bm,Betam),wo='y',fw=fout) #writes out best-fit model
	print(alrm,altm,chim)#,alphlk,likm
	alph = (alrm*altm**2.)**(1/3.)
	b.ar = alph
	b.at = alph	
	b.mkxi()
	b.mkxism()
	chi = b.chi_templ_alphfXX_an((Bm,Betam),wo='y',fw=fout+'ep0')
	#print chi
	fo.close()
	fg.close()
	ans = sigreg_2dme(Hdir+'2Dbaofits/arat'+fout+'1covchi',spar=spar,spat=spat)
	return ans


if __name__ == '__main__':
	import sys
	#This is setup to run the data in the Ross_2016_COMBINEDDR12 folder
	min = 50.
	max = 150. #the minimum and maximum scales to be used in the fit
	maxb = 80. #the maximum scale to be used to set the bias prior
	dir = '/Users/ashleyross/DR12/Ross_2016_COMBINEDDR12/' #change to wherever the data is
	ft = 'Ross_2016_COMBINEDDR12_'
	zb = 'zbin3_' #change number to change zbin
	binc = 0 #change number to change bin center
	bs = 5. #the bin size
	bc = 'post_recon_bincent'+str(binc)+'.dat' 
	c = np.loadtxt(dir+ft+zb+'covariance_monoquad_'+bc)
	d0 = np.loadtxt(dir+ft+zb+'correlation_function_monopole_'+bc).transpose()[1]
	d2 = np.loadtxt(dir+ft+zb+'correlation_function_quadrupole_'+bc).transpose()[1]
	if len(c) != len(d0)*2:
		print('MISMATCHED data and cov matrix!')
	dv = [] #empty list to become data vector
	dvb = [] #empty list to become data vector for setting bias prior
	rl = [] #empty list to become list of r values to evaluate model at	
	rlb  = [] #empty list to become list of r values to evaluate model at for bias prior
	mini = 0
	for i in range(0,len(d0)):
		r = i*bs+bs/2.+binc
		if r > min and r < max:
			dv.append(d0[i])
			rbc = .75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.) #correct for pairs should have slightly larger average pair distance than the bin center
			rl.append(rbc) 
			if mini == 0:
				mini = i #minimum index number to be used for covariance matrix index assignment
			if r < maxb:
				dvb.append(d0[i])
				rlb.append(rbc)
	for i in range(0,len(d2)):
		r = i*bs+bs/2.+binc
		if r > min and r < max:
			dv.append(d2[i])
			rbc = .75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.)
			rl.append(rbc)
			if r < maxb:
				dvb.append(d2[i])
				rlb.append(rbc)

	dv = np.array(dv)
	print(len(dv))
	covm = zeros((len(dv),len(dv))) #will become covariance matrix to be used with data vector
	#need to cut it to correct size
	for i in range(0,len(c)):
		if i < len(d0):
			ri = i*bs+bs/2.+binc
			indi = i-mini
		else:
			ri = (i-len(d0))*bs+bs/2.+binc
			indi = len(dv)/2+i-mini-len(d0)	
		for j in range(0,len(c)):		
			if j < len(d0):
				rj = j*bs+bs/2.+binc
				indj = j-mini
			else:
				rj = (j-len(d0))*bs+bs/2.+binc
				indj = len(dv)/2+j-mini-len(d0)
			if ri > min and ri < max and rj > min and rj < max:
				#print ri,rj,i,j,indi,indj
				covm[indi][indj] = c[i][j]
	invc = pinv(covm) #the inverse covariance matrix to pass to the code
	covmb = zeros((len(dvb),len(dvb)))
	for i in range(0,len(dvb)):
		if i < len(dvb)/2:
			indi = i
		else:
			indi = i-len(dvb)/2+len(covm)/2
		for j in range(0,len(dvb)):
			if j < len(dvb)/2:
				indj = j
			else:
				indj = j-len(dvb)/2+len(covm)/2
			covmb[i][j] = covm[indi][indj]
	invcb = pinv(covmb)
	mod = 'Challenge_matterpower0.44.02.54.015.01.0.dat' #BAO template used		
	fout = ft+zb+bc
	spa = .001
	mina = .8
	maxa = 1.2
	Xism_arat_1C_an(dv,invc,rl,mod,dvb,invcb,rlb)
	#cl = doxi_isolike(d,ct,mod,r,spa=spa,mina=mina,maxa=maxa)
