import numpy as np
import numpy.linalg as linalg
from math import *

class baofit_iso:
	def __init__(self,xid,covd,modl,rl,rmin=50,rmax=150,sp=1.,cov2 = ''):
		#xid and covd are the data vector and the covariance matrix and should be matched in length
		#modl is the BAO template, assumed to have spacing sp between 10 and 300 mpc/h
		#rl is the list of r values matching xid and covd
		self.xim = [] #these lists will be filled based on the r limits
		self.rl = []
		self.sp = sp
		bs = rl[1]-rl[0] #assumes linear bins
		s = 0
		nxib = len(xid)
		rsw = 0
		Bnode = 50.
		print 'total available xi(s) bins '+str(nxib)
		for i in range(0,nxib):
			r = rl[i]
			if r > rmin and r < rmax:				
				if s == 0:
					mini = i
					s = 1
				rbc = .75*((r+bs/2.)**4.-(r-bs/2.)**4.)/((r+bs/2.)**3.-(r-bs/2.)**3.)
				#print rbc
				#this is the properly weighted bin center assume spherical symmetry
				self.rl.append(rbc)
				self.xim.append(xid[i])
		self.nbin = len(self.rl)
		print 'using '+ str(self.nbin)+' xi(s) bins'
		#mt = zeros((self.nbin,self.nbin)) #this will be the trimmed covariance matrix
		#for i in range(mini,mini+self.nbin):
		#	for j in range(mini,mini+self.nbin):
		#		mt[i-mini][j-mini] = covd[i][j]
		mt = covd[mini:mini+self.nbin,mini:mini+self.nbin]
		self.invt = linalg.pinv(mt)
		if cov2 != '':
			mt2 = cov2[mini:mini+self.nbin,mini:mini+self.nbin]
			self.invt += linalg.pinv(mt2)
		self.ximodmin = 10. #minimum of template
		self.modl = modl
						
	def wmod(self,r):
		sp = self.sp
		indd = int((r-self.ximodmin)/sp)
		indu = indd + 1
		fac = (r-self.ximodmin)/sp-indd
		if fac > 1.:
			print 'BAD FAC in wmod'
			return 'ERROR, BAD FAC in wmod'
		if indu >= len(self.modl)-1:
			return -5.47608128044e-05,5.7422824622e-06
		a = self.modl[indu]*fac+(1.-fac)*self.modl[indd]
		return a
		
		

	def chi_templ_alphfXX(self,list,Bnode=50.,wo='n',fw=''):
		from time import time
		t = time()
		B = list#[0]
		if B < 0:
			return 1000
		alph = self.alph
		#A0 = list[1]
		#A1 = list[2]
		#A2 = list[3]
		pv = []
		for i in range(0,self.nbin):
			r = self.rl[i]*alph
			wm = self.wmod(r)
			pv.append(self.xim[i]-B*wm)
		 
		Al = findPolya(self.H,self.invt,pv)
		A0,A1,A2 = Al[0],Al[1],Al[2]
		self.A0 = A0
		self.A1 = A1
		self.A2 = A2
		modl = np.zeros((self.nbin))
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin):
			r = self.rl[i]
			ply = A0+A1/r+A2/r**2.
			r = self.rl[i]*alph
			wm = self.wmod(r)
			mod = B*wm+ply
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod[0])+'\n')
				#print mod,self.xim[i]
			modl[i] = mod
		if wo == 'y':
			fo.close()		
			
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = np.zeros((self.nbin))
		
		for i in range(0,self.nbin):
			dl[i] = self.xim[i]-modl[i]
		chit = np.dot(np.dot(dl,self.invt),dl)
		BBfac = (log(B/self.BB)/self.Bp)**2. #bias prior
		return chit+BBfac

	def chi_templ_alphfXXn(self,list,Bnode=50.,wo='n',fw=''):
		from time import time
		t = time()
		B = list[0]
		if B < 0:
			return 1000
		alph = self.alph
		modl = np.zeros((self.nbin))
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin):
			r = self.rl[i]*alph
			wm = self.wmod(r)
			mod = B*wm
			modl[i] = mod
		if wo == 'y':
			fo.close()		
		
		dl = np.zeros((self.nbin))	
		for i in range(0,self.nbin):
			dl[i] = self.xim[i]-modl[i]
		chit = np.dot(np.dot(dl,self.invt),dl)	
		return chit

class baofit_isoN:
	def __init__(self,N,xid,covd,modl,rl,sp=1.):
		#should take in N data vectors with corresponding models
		#xid is list of data vectors; should already be cut to correct scales, given in rl 
		#covd is full covariance matrix
		#modl is list of BAO templates
		#rl is the list of r values matching xid and covd
		self.xim = [] #these lists will be filled based on the r limits
		self.rl = []
		self.sp = sp
		bs = rl[1]-rl[0] #assumes linear bins
		s = 0
		nxib = len(xid)
		rsw = 0
		Bnode = 50.
		self.rl = rl
		self.xim = []
		for j in range(0,N):
			for i in range(0,len(xid[0])):
				self.xim.append(xid[j][i])
		print self.xim,len(self.xim)		
		
		self.nbin = len(self.rl)
		print 'using '+ str(self.nbin)+' xi(s) bins'
		#mt = zeros((self.nbin,self.nbin)) #this will be the trimmed covariance matrix
		#for i in range(mini,mini+self.nbin):
		#	for j in range(mini,mini+self.nbin):
		#		mt[i-mini][j-mini] = covd[i][j]
		mt = covd#[self.nbin,mini:mini+self.nbin]
		self.invt = linalg.pinv(mt)
		self.ximodmin = 10. #minimum of template
		self.modl = modl
		self.N = N
		print self.wmod(100.)
						
	def wmod(self,r):
		sp = self.sp
		indd = int((r-self.ximodmin)/sp)
		indu = indd + 1
		fac = (r-self.ximodmin)/sp-indd
		if fac > 1.:
			print 'BAD FAC in wmod'
			return 'ERROR, BAD FAC in wmod'
		if indu >= len(self.modl[0])-1:
			return -5.47608128044e-05,5.7422824622e-06
		a = []
		for i in range(0,self.N): 
			ans = self.modl[i][indu]*fac+(1.-fac)*self.modl[i][indd]
			a.append(ans)
		return a
		
		

	def chi_templ_alphfXX(self,list,Bnode=50.,wo='n',fw=''):
		from time import time
		t = time()
		B = list#[0]
		if B[0] < 0 or B[1] < 0:
			return 1000
		alph = self.alph
		#A0 = list[1]
		#A1 = list[2]
		#A2 = list[3]
		pv = []
		for j in range(0,self.N):
			for i in range(0,self.nbin):
				r = self.rl[i]*alph
				wm = self.wmod(r)
				#print B[j],wm[j],self.xim[j]
				pv.append(self.xim[j*self.nbin+i]-B[j]*wm[j])
		 
		Al = findPolya(self.H,self.invt,pv)
		A0,A1,A2 = [],[],[]
		for i in range(0,self.N):
			A0.append(Al[0+3*i])
			A1.append(Al[1+3*i])
			A2.append(Al[2+3*i])
		#A0,A1,A2 = Al[0],Al[1],Al[2]
		self.A0 = A0
		self.A1 = A1
		self.A2 = A2
		modl = np.zeros((self.N*self.nbin))
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for j in range(0,self.N):
			for i in range(0,self.nbin):
				r = self.rl[i]				
				ply = A0[j]+A1[j]/r+A2[j]/r**2.
				r = self.rl[i]*alph
				wm = self.wmod(r)[j]
				mod = B[j]*wm+ply
				if wo == 'y':
					#fo.write(str(self.rl[i])+' '+str(mod[0])+'\n')
					#print mod,self.xim[i]
					pass
				modl[i+self.nbin*j] = mod
		if wo == 'y':
			fo.close()		
			
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = np.zeros((self.N*self.nbin))
		
		for i in range(0,self.N*self.nbin):
			dl[i] = self.xim[i]-modl[i]
		chit = np.dot(np.dot(dl,self.invt),dl)
		BBfac = 0
		for i in range(0,self.N):
			BBfac += (log(B[i]/self.BB)/self.Bp)**2. #bias prior
		return chit+BBfac

	def chi_templ_alphfXXna(self,list,Bnode=50.,wo='n',fw=''):
		from time import time
		t = time()
		B = (list[0],list[1])
		if B[0] < 0 or B[1] < 0:
			return 1000
		alph = self.alph
		A0 = (list[2],list[3])
		A1 = (list[4],list[5])
		A2 = (list[6],list[7])
		pv = []
		#A0,A1,A2 = Al[0],Al[1],Al[2]
		self.A0 = A0
		self.A1 = A1
		self.A2 = A2
		modl = np.zeros((self.N*self.nbin))
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for j in range(0,self.N):
			for i in range(0,self.nbin):
				r = self.rl[i]				
				ply = A0[j]+A1[j]/r+A2[j]/r**2.
				r = self.rl[i]*alph
				wm = self.wmod(r)[j]
				mod = B[j]*wm+ply
				if wo == 'y':
					#fo.write(str(self.rl[i])+' '+str(mod[0])+'\n')
					#print mod,self.xim[i]
					pass
				modl[i+self.nbin*j] = mod
		if wo == 'y':
			fo.close()		
			
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = np.zeros((self.N*self.nbin))
		
		for i in range(0,self.N*self.nbin):
			dl[i] = self.xim[i]-modl[i]
		chit = np.dot(np.dot(dl,self.invt),dl)
		BBfac = 0
		for i in range(0,self.N):
			BBfac += (log(B[i]/self.BB)/self.Bp)**2. #bias prior
		return chit+BBfac


	def chi_templ_alphfXXn(self,list,Bnode=50.,wo='n',fw=''):
		from time import time
		t = time()
		B = list#[0]
		if B < 0:
			return 1000
		alph = self.alph
		modl = np.zeros((self.N*self.nbin))
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for j in range(0,self.N):
			for i in range(0,self.nbin):
				r = self.rl[i]*alph
				wm = self.wmod(r)[j]
				mod = B[j]*wm
				modl[i+j*self.nbin] = mod
		if wo == 'y':
			fo.close()		
		
		dl = np.zeros((self.nbin))	
		for i in range(0,self.nbin):
			dl[i] = self.xim[i]-modl[i]
		chit = np.dot(np.dot(dl,self.invt),dl)	
		return chit


def doxi_isolike(xid,covd,modl,rl,rmin=50,rmax=150,sp=1.,Bp=.4,rminb=50.,rmaxb=80.,spa=.001,mina=.8,maxa=1.2,chi2fac=1.,v='n',wo='',cov2=''):
	#chi2fac should be hartlap factor
	from time import time
	from optimize import fmin
	b = baofit_iso(xid,covd,modl,rl,rmin=rmin,rmax=rmax,sp=sp,cov2=cov2)
	b.Bp = Bp
	b.H = np.zeros((3,b.nbin))
	for i in range(0,b.nbin):
		b.H[0][i] = 1.
		b.H[1][i] = 1./b.rl[i]
		b.H[2][i] = 1./b.rl[i]**2.

	bb = baofit_iso(xid,covd,modl,rl,rmin=rmin,rmax=rmaxb,sp=sp,cov2=cov2)
	#bb is to set bias prior
	bb.H = np.zeros((3,bb.nbin))
	for i in range(0,bb.nbin):
		bb.H[0][i] = 1.
		bb.H[1][i] = 1./bb.rl[i]
		bb.H[2][i] = 1./bb.rl[i]**2.

	alphl= []
	chil = []
	likl = []
	chim = 1000
	na = int((maxa-mina)/spa)
	likm = -1
	pt = 0
	A0 = 0
	A1 = 1.
	A2 = 0
	b.alph = 1.
	bb.alph = b.alph
	#b.alph = 0
	B = .1
	chiBmin = 1000
	while B < 2.:
		bl = [B]
		chiB = bb.chi_templ_alphfXXn(bl)*chi2fac
		if chiB < chiBmin:
			chiBmin = chiB
			BB = B
		B += .01	
	print 'best-fit bias factor is '+str(BB)+' '+str(chiBmin)
	b.BB = BB		
	#b.BB = 1. #switch to this to make bias prior centered on input rather than fit value
	B = BB
	for i in range(0,na):		
		b.alph = mina+spa*i+spa/2.
		#inl = np.array([B,A0,A1,A2])
		#inl = np.array(B)
		inl = B
		#(B,A0,A1,A2) = fmin(b.chi_templ_alphfXX,inl,disp=False)
		B = fmin(b.chi_templ_alphfXX,inl,disp=False)
		#B = fmin(b.chi_templ_alphfXXn,inl,disp=False)
		#chi = b.chi_templ_alphfXX((B,A0,A1,A2))*chi2fac
		chi = b.chi_templ_alphfXX((B))*chi2fac
		if v == 'y':
			print b.alph,chi,B[0],b.A0[0],b.A1[0],b.A2[0] #single values getting output as arrays, silly, but works so not worrying about it
		alphl.append(b.alph)
		chil.append(chi)
		if chi < chim:
			chim = chi
			alphm = b.alph
			Bm = B
			A0m = b.A0
			A1m = b.A1
			A2m = b.A2
	print alphm,chim,Bm,A0m,A1m,A2m
	fo = open('BAOisobestfit'+wo+'.dat','w')
	b.alph = alphm	
	b.chi_templ_alphfXX((Bm),wo='y',fw=wo)
	fo.write(str(alphm)+' '+str(chim)+' '+str(Bm[0])+' '+str(A0m[0])+' '+str(A1m[0])+' '+str(A2m[0])+'\n')
	fo.close()
	return chil

def doxi_isolikeN(N,xid,covd,modl,rl,rmin=50,rmax=150,sp=1.,Bp=.4,rminb=50.,rmaxb=80.,spa=.001,mina=.8,maxa=1.2,chi2fac=1.,v='n',wo=''):
	#chi2fac should be hartlap factor
	#function is to 
	from time import time
	from optimize import fmin
	b = baofit_isoN(N,xid,covd,modl,rl,sp=sp)
	b.Bp = Bp
	b.H = np.zeros((3*N,N*b.nbin))
	print b.nbin
	for j in range(0,N):
		for i in range(0,N*b.nbin):
			if i > j*b.nbin and i < (j+1)*b.nbin:
				#print i+j*nbin
				b.H[0+3*j][i] = 1.
				b.H[1+3*j][i] = 1./b.rl[i-j*b.nbin]
				b.H[2+3*j][i] = 1./b.rl[i-j*b.nbin]**2.

	bb = baofit_iso(xid[0],covd,modl[0],rl,rmin=rminb,rmax=rmaxb,sp=sp)
	#bb is to set bias prior
	bb.H = np.zeros((3,bb.nbin))
	#for j in range(0,N):
	for i in range(0,bb.nbin):
		bb.H[0][i] = 1.
		bb.H[1][i] = 1./bb.rl[i]
		bb.H[2][i] = 1./bb.rl[i]**2.

	alphl= []
	chil = []
	likl = []
	chim = 1000
	na = int((maxa-mina)/spa)
	likm = -1
	pt = 0
	A0 = 0
	A1 = 1.
	A2 = 0
	b.alph = 1.
	bb.alph = b.alph
	#b.alph = 0
	B = .1
	chiBmin = 1000
	while B < 2.:
		bl = [B]
		#for i in range(0,N):
		#	bl.append(B)
		chiB = bb.chi_templ_alphfXXn(bl)*chi2fac
		if chiB < chiBmin:
			chiBmin = chiB
			BB = B
		B += .01	
	print 'best-fit bias factor is '+str(BB)+' '+str(chiBmin)
	b.BB = BB		
	#b.BB = 1. #switch to this to make bias prior centered on input rather than fit value
	B = [BB,BB]
	for i in range(0,na):		
		b.alph = mina+spa*i+spa/2.
		#inl = np.array([B,A0,A1,A2])
		#inl = np.array(B)
		inl = B
		#(B,A0,A1,A2) = fmin(b.chi_templ_alphfXX,inl,disp=False)
		B = fmin(b.chi_templ_alphfXX,inl,disp=False)
		#B = fmin(b.chi_templ_alphfXXn,inl,disp=False)
		#chi = b.chi_templ_alphfXX((B,A0,A1,A2))*chi2fac
		chi = b.chi_templ_alphfXX((B))*chi2fac
		if v == 'y':
			print b.alph,chi,B[0],b.A0[0],b.A1[0],b.A2[0] #single values getting output as arrays, silly, but works so not worrying about it
		alphl.append(b.alph)
		chil.append(chi)
		if chi < chim:
			chim = chi
			alphm = b.alph
			Bm = B
			A0m = b.A0
			A1m = b.A1
			A2m = b.A2
	print alphm,chim,Bm,A0m,A1m,A2m
	fo = open('BAOisobestfit'+wo+'.dat','w')
	b.alph = alphm	
	b.chi_templ_alphfXX((Bm),wo='y',fw=wo)
	fo.write(str(alphm)+' '+str(chim)+' '+str(Bm[0])+' '+str(A0m[0])+' '+str(A1m[0])+' '+str(A2m[0])+'\n')
	fo.close()
	return chil

def doxi_isolikeNna(N,xid,covd,modl,rl,rmin=50,rmax=150,sp=1.,Bp=.4,rminb=50.,rmaxb=80.,spa=.001,mina=.8,maxa=1.2,chi2fac=1.,v='n',wo=''):
	#chi2fac should be hartlap factor
	from time import time
	from optimize import fmin
	b = baofit_isoN(N,xid,covd,modl,rl,sp=sp)
	b.Bp = Bp

	bb = baofit_iso(xid[0],covd,modl[0],rl,rmin=rminb,rmax=rmaxb,sp=sp)
	#bb is to set bias prior
	bb.H = np.zeros((3,bb.nbin))
	#for j in range(0,N):
	for i in range(0,bb.nbin):
		bb.H[0][i] = 1.
		bb.H[1][i] = 1./bb.rl[i]
		bb.H[2][i] = 1./bb.rl[i]**2.

	alphl= []
	chil = []
	likl = []
	chim = 1000
	na = int((maxa-mina)/spa)
	likm = -1
	pt = 0
	A0 = 0
	A1 = 1.
	A2 = 0
	b.alph = 1.
	bb.alph = b.alph
	#b.alph = 0
	B = .1
	chiBmin = 1000
	while B < 2.:
		bl = [B]
		#for i in range(0,N):
		#	bl.append(B)
		chiB = bb.chi_templ_alphfXXn(bl)*chi2fac
		if chiB < chiBmin:
			chiBmin = chiB
			BB = B
		B += .01	
	print 'best-fit bias factor is '+str(BB)+' '+str(chiBmin)
	b.BB = BB		
	#b.BB = 1. #switch to this to make bias prior centered on input rather than fit value
	B = [BB,BB,0,0,0,0,0,0]
	for i in range(0,na):		
		b.alph = mina+spa*i+spa/2.
		#inl = np.array([B,A0,A1,A2])
		#inl = np.array(B)
		inl = B
		#(B,A0,A1,A2) = fmin(b.chi_templ_alphfXX,inl,disp=False)
		B = fmin(b.chi_templ_alphfXXna,inl,disp=False)
		#B = fmin(b.chi_templ_alphfXXn,inl,disp=False)
		#chi = b.chi_templ_alphfXX((B,A0,A1,A2))*chi2fac
		chi = b.chi_templ_alphfXXna((B))*chi2fac
		if v == 'y':
			print b.alph,chi,B[0],b.A0[0],b.A1[0],b.A2[0] #single values getting output as arrays, silly, but works so not worrying about it
		alphl.append(b.alph)
		chil.append(chi)
		if chi < chim:
			chim = chi
			alphm = b.alph
			Bm = B
			A0m = b.A0
			A1m = b.A1
			A2m = b.A2
	print alphm,chim,Bm
	#fo = open('BAOisobestfit'+wo+'.dat','w')
	b.alph = alphm	
	#b.chi_templ_alphfXX((Bm),wo='y',fw=wo)
	#fo.write(str(alphm)+' '+str(chim)+' '+str(Bm[0])+' '+str(A0m[0])+' '+str(A1m[0])+' '+str(A2m[0])+'\n')
	#fo.close()
	return chil


def findPolya(H,ci,d):
	#analytically solves for best-fit polynomial nuisance terms, given BAO and bias parameters
	ht = H.transpose()
	onei = linalg.pinv(np.dot(np.dot(H,ci),ht))
	comb = np.dot(np.dot(onei,H),ci)
	return np.dot(comb,d)

def sigreg_c12(chilist,mina=.8,maxa=1.2,astep=.001,fac=1.):
	#returns maximum likelihod, +/-1 chi2, +/-4 chi2, chimin
	dir = ''
	chil = []
	chim = 1000
	#na = int(1/astep)
	na = len(chilist)
	fl = []
	for i in range(0,na):
		a = mina + astep/2.+astep*i
		chiv = chilist[i]
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


if __name__ == '__main__':
	import sys
	#This is setup to run the data in the exampledata folder, which is BOSS DR11 data
	c1 = np.loadtxt('exampledata/cov0CMASSreconNScomb1DR118st2.dat')
	c2 = np.loadtxt('exampledata/cov0CMASSreconNScomb2DR118st2.dat')
	ct = (c1+c2)/2. #two independent sets of mocks are averaged for the DR11 covariance matrix
	d = np.loadtxt('exampledata/xi0ACbossv1NScombreconbs8st2.dat').transpose()[1]
	r = np.loadtxt('exampledata/xi0ACbossv1NScombreconbs8st2.dat').transpose()[0]
	mod = np.loadtxt('exampledata/xi0camb_Nacc0n3.00.051.9_And.dat').transpose()[1]
	spa = .001
	mina = .8
	maxa = 1.2
	cl = doxi_isolike(d,ct,mod,r,spa=spa,mina=mina,maxa=maxa)
	al = [] #list to be filled with alpha values
	for i in range(0,len(cl)):
		a = .8+spa/2.+spa*i
		al.append(a)
	#below assumes you have matplotlib to plot things, if not, save the above info to a file or something
	from matplotlib import pyplot as plt
	plt.plot(al,cl-min(cl),'k-')
	plt.show()