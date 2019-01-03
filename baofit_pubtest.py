import numpy as np
import numpy.linalg as linalg
from math import *

class baofit_iso:
	def __init__(self,xid,covd,modl,modsmoothl,rl,rmin=50,rmax=150,sp=1.,cov2 = '',facc=1.):
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
		self.modsmoothl = modsmoothl
						
	def wmod(self,r):
		sp = self.sp
		indd = int((r-self.ximodmin)/sp)
		indu = indd + 1
		fac = (r-self.ximodmin)/sp-indd
		if fac > 1.:
			print 'BAD FAC in wmod'
			return 'ERROR, BAD FAC in wmod'
		if indu >= len(self.modl)-1:
			#return -5.47608128044e-05,5.7422824622e-06
			return self.modl[-1]
		a = self.modl[indu]*fac+(1.-fac)*self.modl[indd]
		return a

	def wmodsmooth(self,r):
		sp = self.sp
		indd = int((r-self.ximodmin)/sp)
		indu = indd + 1
		fac = (r-self.ximodmin)/sp-indd
		if fac > 1.:
			print 'BAD FAC in wmod'
			return 'ERROR, BAD FAC in wmod'
		if indu >= len(self.modsmoothl)-1:
			#return -5.47608128044e-05,5.7422824622e-06
			return self.modsmoothl[-1]
		a = self.modsmoothl[indu]*fac+(1.-fac)*self.modsmoothl[indd]
		return a
		
		

	def chi_templ_alphfXX(self,list,Bnode=50.,wo='n',fw=''):
		from time import time
		t = time()
		B = list[0]
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
			#if len(self.xim[i]-B*wm) > 1:
			#print self.xim[i]-B*wm,self.xim[i],B,wm
			pv.append(self.xim[i]-B*wm)
		 
		Al = findPolya(self.H,self.invt,pv)
		#A0,A1,A2 = Al[0],Al[1],Al[2]
		#self.A0 = A0
		#self.A1 = A1
		#self.A2 = A2
		modl = np.zeros((self.nbin))
		if wo == 'y':
			fo = open('ximod'+fw+'.dat','w')
		for i in range(0,self.nbin):
			r = self.rl[i]
			#ply = A0+A1/r+A2/r**2.
			ply = 0
			for j in range(0,self.np):
				ply += Al[j]/r**j
			r = self.rl[i]*alph
			wm = self.wmod(r)
			wsm = self.wmodsmooth(r)
			mod = B*wm+ply
			modsm = B*wsm+ply
			if wo == 'y':
				fo.write(str(self.rl[i])+' '+str(mod)+' '+str(modsm)+' '+str(self.xim[i])+'\n')
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
		#A0,A1,A2 = [],[],[]
		#for i in range(0,self.N):
		#	for j in range(0,self.np):
		#		A0.append(Al[0+3*i])
		#	A1.append(Al[1+3*i])
		#	A2.append(Al[2+3*i])
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
				ply = 0
				for k in range(0,self.np):
					ply += Al[k+self.np*j]/r**k				
				#ply = A0[j]+A1[j]/r+A2[j]/r**2.
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

class baofitPk:
	def __init__(self,pkl,kl,cov,min=.02,max=.3,damp=6.,pktemp='Challenge_matterpower'):
		from EH import simulate
		self.Pm = []
		self.kl = []
		self.damp = damp
		s = 0
		for i in range(0,len(kl)):
			if kl[i] > min and kl[i] < max:
				if s == 0:
					mini = i
					s = 1
				self.kl.append(kl[i])
				self.Pm.append(pkl[i])
		self.nbin = len(self.kl)
		print 'using '+ str(self.nbin)+' P(k) bins'
		mt = cov[mini:mini+self.nbin,mini:mini+self.nbin]
		self.invt = linalg.pinv(mt)
		#setup p wiggle
		ptemp = np.loadtxt('powerspectra/'+pktemp+'.dat').transpose()
		self.ptkl = ptemp[0]		
		if pktemp=='Challenge_matterpower' or pktemp == 'TSPT_out':
			om = 0.31
			lam = 0.69
			h = .676
			nindex = .963
			ombhh = .022
		s = simulate(omega=om,lamda=lam,h=h,nindex=nindex,ombhh=ombhh) #smooth
		norm = float(ptemp[1][0])/s.Psmooth(ptemp[0][0],0)
		self.pwigl = []
		self.psl = []
		self.ks = .001
		kt = self.ks
		fo = open('psmoothtemp.dat','w')
		while kt < .6:
			ps = s.Psmooth(kt,0)*norm
			fo.write(str(kt)+' '+str(ps)+'\n')
			pw = np.interp(kt,self.ptkl,ptemp[1])
			self.pwigl.append(pw/ps)
			kt += self.ks
		for k in self.kl:
			self.psl.append(s.Psmooth(k,0)*norm)

		#for i in range(0,len(self.ptkl)):
		#	k = self.ptkl[i]
		#	ps = s.Psmooth(k,0)*norm
		#	pw = ptemp[1][i]/(ps)
		#	self.psl.append(ps)
		#	self.pwigl.append(pw)
	
	def Pwig(self,k):
		indd = int((k-self.ks)/self.ks)
		indu = indd+1
		fac = (k-(indd*self.ks+self.ks))/self.ks
		if fac > 1:
			print 'factor greater than 1'
			return 'ERROR'
		#return ((1.-fac)*self.Pkcl[indd][0]+fac*self.Pkcl[indu][0])/((1.-fac)*self.Pkcl[indd][1]+fac*self.Pkcl[indu][1])	
		return ((1.-fac)*self.pwigl[indd]+fac*self.pwigl[indu])
		
	def Pconv(self,B,alph,A0,A1,A2,snl=8.,k0=0,k1=0,w='n',fl=''):
		from time import time
		#t = time()
		pthl = zeros((self.nkin))
		#if w == 'y':
		#	fo = open('pkrat.dat','w')
		for ktb in range(0,self.nkin):
			k = self.kin[ktb]
			#psm = self.s.Psmooth(k,0)+A0/k+A1/k**2.+A2/k**3.
			psm = self.Psml[ktb]+A0/k+A1/k**2.+A2/k**3.+k0+k1*k
			#pc = self.Pkc(k/alph)/(self.Psml[ktb])
			if alph != 0:
				pc = self.Pwig(k/alph)#/(self.Psml[ktb])
				pm = (pc-1.)*exp(-(k/alph)**2.*snl**2./2.)+1.
				mod = B*pm*psm
			else:
				mod = B*psm
			#mod = 1.
			#if w == 'y':
			#	fo.write(str(k)+' '+str(pc/psm)+'\n')
			pthl[ktb] = mod
		#if w == 'y':
		#	fo.close()	
		#print time()-t
		pkconv0 = 0
		kmaxmax = 2.
		#for kb in range(0,len(winf)-3):
		#	pkconv0 += float(winf[kb+3].split()[1])*pthl[kb]
		winfac = pkconv0#/float(self.winf[2].split()[0]) no zero bin subtraction currently!!!
		#pconv = []
		pconv = zeros((self.nbin))
		for kb in range(0,self.nbin):
			pc = 0
			#sub = winfac*float(winf[2].split()[kb-1])
			vb = self.wmat[kb]
			pc = dot(vb,pthl)		
			#for ktb in range(0,self.nkin):
			#	pth =pthl[ktb]
			#	pc += self.wmat[ktb][kb]*pth
			#if self.lg == 'lg':
			#pconv.append(pc)
			pconv[kb] = pc
			#else:
			#	pconv.append(pc)	
		#print time()-t			
		if w == 'y':
			fo = open('pconv'+fl+'.dat','w')
			for i in range(0,len(pconv)):
				fo.write(str(pconv[i])+'\n')
			fo.close()			

		return pconv

	def Pconv_alph(self,alph,snl=8.,w='n'):
		from time import time
		t = time()
		pthl = []
		if w == 'y':
			fo = open('pkrat.dat','w')
		for ktb in range(0,self.nkin):
			k = self.kin[ktb]
			#pc = self.Pkc(k/alph)/(self.Psml[ktb])
			pc = self.Pwig(k/alph)/(self.Psml[ktb])
			pm = (pc-1)*exp(-k**2.*snl**2./2.)+1
			mod = pm
			#mod = 1.
			if w == 'y':
				fo.write(str(k)+' '+str(pm)+'\n')
			pthl.append(mod)
		if w == 'y':
			fo.close()	
		#print time()-t
		pkconv0 = 0
		kmaxmax = 2.
		#for kb in range(0,len(winf)-3):
		#	pkconv0 += float(winf[kb+3].split()[1])*pthl[kb]
		winfac = pkconv0#/float(self.winf[2].split()[0]) no zero bin subtraction currently!!!
		pconv = []
		for kb in range(0,self.nbin):
			pc = 0
			#sub = winfac*float(winf[2].split()[kb-1])		
			for ktb in range(0,self.nkin):
				pth =pthl[ktb]
				pc += self.wmat[ktb][kb]*pth		
			pconv.append(pc)
		#print time()-t			
		return pconv


	def Pconv_sm(self,B,A0,A1,A2,w='n',fl=''):
		from time import time
		t = time()
		pthl = []
		for ktb in range(0,self.nkin):
			k = self.kin[ktb]
			psm = self.Psml[ktb]+A0/k+A1/k**2.+A2/k**3.
			mod = B*psm
			pthl.append(mod)
		#print time()-t
		pkconv0 = 0
		kmaxmax = 2.
		#for kb in range(0,len(winf)-3):
		#	pkconv0 += float(winf[kb+3].split()[1])*pthl[kb]
		winfac = pkconv0#/float(self.winf[2].split()[0]) no zero bin subtraction currently!!!
		pconv = []
		for kb in range(0,self.nbin):
			pc = 0
			#sub = winfac*float(winf[2].split()[kb-1])		
			for ktb in range(0,self.nkin):
				pth =pthl[ktb]
				pc += self.wmat[ktb][kb]*pth
		
			pconv.append(pc)
		#print time()-t	
		if w == 'y':
			fo = open('pconvsm'+fl+'.dat','w')
			for i in range(0,len(pconv)):
				fo.write(str(pconv[i])+'\n')
			fo.close()			
		return pconv

	def chi_temp_sm(self,B,A0,A1,A2):
		from time import time
		t = time()
		modl = self.Pconv_sm(B,A0,A1,A2)
		#print time()-t	
		for i in range(0,self.nbin):
			mod = modl[i]
			if mod > 0:
				modl[i] = (log(mod))
			else:
				return 1000.
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.Pm[i]-modl[i])
		#print dl
		#print time()-t
		for i in range(0,self.nbin):
			chid += (dl[i])**2*self.invt[i,i]
			for j in range(i+1,self.nbin):	
				Chioff += dl[i]*self.invt[i,j]*dl[j]
			#print chid+2.*Chioff,dl[i],sqrt(1./self.invt[i,i])
			#for j in range(0,self.nbin):
			#	ChiSq += dl[i]*self.invt[i,j]*dl[j]		
		#print time()-t
		return 2.*Chioff+chid		
		
		
	def chi_temp(self,B,alph,A0,A1,A2,snl=8.):
		from time import time
		t = time()
		modl = self.Pconv(B,alph,A0,A1,A2,snl)
		#print time()-t	
		for i in range(0,self.nbin):
			mod = modl[i]
			if mod > 0:
				modl[i] = (log(mod))
			else:
				return 1000.
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.Pm[i]-modl[i])
		#print dl
		#print time()-t
		for i in range(0,self.nbin):
			chid += (dl[i])**2*self.invt[i,i]
			for j in range(i+1,self.nbin):	
				Chioff += dl[i]*self.invt[i,j]*dl[j]
			#print chid+2.*Chioff,dl[i],sqrt(1./self.invt[i,i])
			#for j in range(0,self.nbin):
			#	ChiSq += dl[i]*self.invt[i,j]*dl[j]		
		#print time()-t
		return 2.*Chioff+chid		
		#return ChiSq

	def chi_templ_alphf(self,list,snl=4.3,k0m='y',k1m='y'):
		from time import time
		t = time()

		B = list[0]
		A0 = list[1]
		A1 = list[2]
		A2 = list[3]
		if k0m == 'y':
			k0 = list[4]
		else:
			k0 = 0
		if k1m == 'y':
			k1 = list[5]
		else:
			k1 = 0
		modl = self.Pconv(B,self.alph,A0,A1,A2,snl,k0,k1)
		#print time()-t	
		for i in range(0,self.nbin):
			mod = modl[i]
			if mod > 0:
				if self.lg == 0:
					modl[i] = log(mod)
				if self.lg == 10:
					modl[i] = log(mod,10)	
			else:
				return 1000.
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.Pm[i]-modl[i])
		#print dl
		#print time()-t
		for i in range(0,self.nbin):
			chid += (dl[i])**2*self.invt[i,i]
			for j in range(i+1,self.nbin):	
				Chioff += dl[i]*self.invt[i,j]*dl[j]
			#print chid+2.*Chioff,dl[i],sqrt(1./self.invt[i,i])
			#for j in range(0,self.nbin):
			#	ChiSq += dl[i]*self.invt[i,j]*dl[j]		
		#print time()-t
		return 2.*Chioff+chid		

	def chi_temp_noconv(self,list):
		B = list[0]
		A0 = list[1]
		A1 = list[2]
		A2 = list[3]
		try:
			A3 = list[4]
			A4 = list[5]
		except:
			pass
		#snl = list[4]
		snl = self.snl
		modl = []

		for i in range(0,self.nbin):
			k = self.kl[i]
			ka = k/self.alph
			#pnw = np.interp(k,self.ptkl,self.psl)
			#pw = np.interp(ka,self.ptkl,self.pwigl)
			pw = self.Pwig(ka)
			psm = B*self.psl[i]+A0*k+A1+A2/k
			try:
				psm += A3/k**2.+A4/k**3.
			except:
				pass
			mod = psm*(1.+(pw-1.)*exp(-.5*snl**2.*k**2.))
			modl.append(mod)
		dl = np.zeros((self.nbin))
		for i in range(0,self.nbin):
			dl[i] = (self.Pm[i]-modl[i])
		chit = np.dot(np.dot(dl,self.invt),dl)
		return chit
			
	def chi_templ_alphf_snl(self,list,v='n'):
		from time import time
		t = time()
		B = list[0]
		A0 = list[1]
		A1 = list[2]
		A2 = list[3]
		snl = list[4]
		if snl < 0 or snl > 12:
			return 1000.
		k0 = 0
		k1 = 0
		if self.k0m == 'y':
			k0 = list[5]
		if self.k1m == 'y':
			k1 = list[6]	
		modl = self.Pconv(B,self.alph,A0,A1,A2,snl,k0,k1)
		#print time()-t	
		for i in range(0,self.nbin):
			mod = modl[i]
			if self.lg == 'lg':
				if mod > 0:					
					modl[i] = log(mod)
					if self.lg == 10:
						modl[i] = log(mod,10)	
				else:
					return 1000.
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = zeros((self.nbin))
		for i in range(0,self.nbin):
			dl[i] = (self.Pm[i]-modl[i])
		#print dl
		#print time()-t
		#print self.pm,modl
		chit = dot(dot(dl,self.invt),dl)
		#for i in range(0,self.nbin):
		#	chid += (dl[i])**2*self.invt[i,i]
		#	if v == 'y':
		#		print dl[i],sqrt(self.cov[i,i])
		#	for j in range(i+1,self.nbin):	
		#		Chioff += dl[i]*self.invt[i,j]*dl[j]
			#print chid+2.*Chioff,dl[i],sqrt(1./self.invt[i,i])
			#for j in range(0,self.nbin):
			#	ChiSq += dl[i]*self.invt[i,j]*dl[j]		
		#print time()-t
		#sigsnl = 2.
		snlfac =0
		snlfac = ((self.damp-snl)/self.sigsnl)**2.
		#sigsnl = 1.
		#snlfac = -2.*log(exp(-.5*((4.3-snl)/sigsnl)**2.))
		#return 2.*Chioff+chid+snlfac
		return chit+snlfac		


	def chi_temp2(self,B,A0,A1,A2):
		from time import time
		t = time()
		#print time()-t	
		modl = []
		for i in range(0,self.nbin):
			k = self.kl[i]
			mod = B*self.Palph[i]*(self.Psm0l[i]+A0/k+A1/k**2.+A2/k**3.)
			if mod > 0:
				if self.lg == 0:
					modl.append(log(mod))
				else:
					modl.append(log(mod,10.))
			else:
				return 100000.
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.Pm[i]-modl[i])
		#print dl
		#print time()-t
		for i in range(0,self.nbin):
			chid += (dl[i])**2*self.invt[i,i]
			for j in range(i+1,self.nbin):	
				Chioff += dl[i]*self.invt[i,j]*dl[j]
			#print chid+2.*Chioff,dl[i],sqrt(1./self.invt[i,i])
			#for j in range(0,self.nbin):
			#	ChiSq += dl[i]*self.invt[i,j]*dl[j]		
		#print time()-t
		return 2.*Chioff+chid		

	def chi_temp2l(self,list):
		from time import time
		t = time()
		B = list[0]
		A0 = list[1]
		A1 = list[2]
		A2 = list[3]
		#print time()-t	
		modl = []
		for i in range(0,self.nbin):
			k = self.kl[i]
			mod = B*self.Palph[i]*(self.Psm0l[i]+A0/k+A1/k**2.+A2/k**3.)
			if mod > 0:
				if self.lg == 0:
					modl.append(log(mod))
				else:
					modl.append(log(mod,10.))
			else:
				return 100000.
		ChiSq = 0
		chid = 0
		Chioff = 0
		dl = []
		for i in range(0,self.nbin):
			dl.append(self.Pm[i]-modl[i])
		#print dl
		#print time()-t
		for i in range(0,self.nbin):
			chid += (dl[i])**2*self.invt[i,i]
			for j in range(i+1,self.nbin):	
				Chioff += dl[i]*self.invt[i,j]*dl[j]
			#print chid+2.*Chioff,dl[i],sqrt(1./self.invt[i,i])
			#for j in range(0,self.nbin):
			#	ChiSq += dl[i]*self.invt[i,j]*dl[j]		
		#print time()-t
		return 2.*Chioff+chid		


	def chi_smooth(self,B,A0,A1,A2):
		ChiSq = 0
		for i in range(0,self.nkin):
			k = self.kin[i]
			psm = B*(self.Psml[i]+A0/k+A1/k**2.+A2/k**3.)
			pc = self.Pc0l[i]
			ChiSq += (pc-psm)**2./(pc*.1)**2.
		return ChiSq		

	def Pksm_TH(self,Berr=.001,A0err=.01,A1err=.0001,A2err=.00001,N=100):
		from random import gauss, random
		A0 = random()*A0err
		A1 = random()*A1err
		B = 1. +random()*Berr
		A2 = random()*A2err
		chim = 1000
		bs = .2
		nbin = int((9.)/bs)
		dof = float(nbin-5)
		nchitot = 0
		nfA0 = 0
		nfA1 = 0
		nfA2 = 0
		nfB = 0
		nfal = 0
		alph = 1.
		for i in range(0,N):
			s = 0
			sb = 0
			sg = 0
			ssig =0
			sC= 0
			sth = 0
			chi = self.chi_smooth(B,A0,A1,A2)
			#print chi
			if chi < chim:
				chim = chi
				print i,chim
				bestA0 = A0
				bestB = B
				bestA1 = A1
				bestA2 = A2
			chid = chi-chim
			nchitot += 1.
			lc = exp(-.5*chi)
			while s == 0:
				Bn = exp(gauss(log(B),Berr))			
				#while Bn < 0:
				#	Bn = gauss(B,Berr)
				#	nchitot += 1.
				ran = random()
				chin =  self.chi_smooth(B,A0,A1,A2)
				nchitot += 1.
				#print chi,chin
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					s = 1
					B = Bn
				else:
					nfB += 1.

			while sg == 0:
				A0n = gauss(A0,A0err)
				ran = random()
				chin =  self.chi_smooth(B,A0,A1,A2)
				nchitot += 1.
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					sg = 1
					A0 = A0n
				else:
					nfA0 += 1.
			while sC == 0:
				A1n = gauss(A1,A1err)
				ran = random()
				chin =  self.chi_smooth(B,A0,A1,A2)
				nchitot += 1.
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					sC = 1
					A1 = A1n
				else:
					nfA1 += 1.
			while sth == 0:
				A2n = gauss(A2,A2err)
				while abs(A2n) > 40:
					A2n = gauss(A2,A2err)
				ran = random()
				chin =  self.chi_smooth(B,A0,A1,A2)
				nchitot += 1.
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					sth = 1
					A2 = A2n
				else:
					nfA2 += 1.
		print bestB,bestA0,bestA1,bestA2,chim
		fo = open('smoothfit.dat','w')
		fo.write(str(bestB)+' '+str(bestA0)+' '+str(bestA1)+' '+str(bestA2)+' '+str(chim)+'\n')
		fo.close()
		return bestB,bestA0,bestA1,bestA2	
	
	def smfind(self,N,B,A0,A1,A2,Berr=.1,A0err=10,A1err=1.,A2err=.1):
		from random import gauss, random
		chim = 10000

		for i in range(0,N):
			s = 0
			sb = 0
			sg = 0
			ssig =0
			sC= 0
			sth = 0
			chi = self.chi_temp2(B,A0,A1,A2)
			if chi < chim:
				chim = chi
				print i,chim,B,A0,A1,A2
				bestA0 = A0
				bestB = B
				bestA1 = A1
				bestA2 = A2
			B = bestB
			A0 = bestA0
			A1 = bestA1
			A2 = bestA2	
			chid = chi-chim
			lc = exp(-.5*chi)
			while s == 0:
				Bn = exp(gauss(log(B),Berr))			
				ran = random()
				chin =  self.chi_temp2(Bn,A0,A1,A2)
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					s = 1
					B = Bn

			while sg == 0:
				A0n = gauss(A0,A0err)
				ran = random()
				chin =  self.chi_temp2(B,A0n,A1,A2)
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					sg = 1
					A0 = A0n
			while sC == 0:
				A1n = gauss(A1,A1err)
				ran = random()
				chin =  self.chi_temp2(B,A0,A1n,A2)
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					sC = 1
					A1 = A1n
			while sth == 0:
				A2n = gauss(A2,A2err)
				#while abs(A2n) > 20:
				#	A2n = gauss(A2,A2err)
				ran = random()
				chin =  self.chi_temp2(B,A0,A1,A2n)
				ln = exp(-.5*(chin-chi))
				if ln > ran:
					sth = 1
					A2 = A2n
		return chim,bestB,bestA0,bestA1,bestA2

def doPk_isolike_noconv(pkd,kl,cov,snl=6.,kmin=0.02,kmax=.3,npar=3,sp=1.,spa=.001,mina=.8,maxa=1.2,chi2fac=1.,Nmock=1000,v='n',wo=''):
	#chi2fac should be hartlap factor
	from time import time
	from optimize import fmin
	print np
	
	b = baofitPk(pkd,kl,cov,min=kmin,max=kmax)
	b.snl = snl
	#b.Bp = Bp
	#b.np = npar
	#b.H = np.zeros((npar,b.nbin))
	chi2fac = (Nmock-b.nbin-2.)/(Nmock-1.)
	print b.nbin,chi2fac
	#for i in range(0,b.nbin):
	#	for j in range(0,npar):
	#		b.H[j][i] = 1./b.rl[i]**j

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
	A3 = 0
	A4 = 0
	b.alph = 1.
	B = 1.
	for i in range(0,na):		
		b.alph = mina+spa*i+spa/2.
		inl = np.array([B,A0,A1,A2,A3,A4])
		#inl = np.array(B)
		#inl = B
		#(B,A0,A1,A2) = fmin(b.chi_templ_alphfXX,inl,disp=False)
		#if npar > 0:
		(B,A0,A1,A2,A3,A4) = fmin(b.chi_temp_noconv,inl,disp=False)
		chi = b.chi_temp_noconv((B,A0,A1,A2,A3,A4))*chi2fac
		
		if v == 'y':
			print b.alph,chi,B#[0],b.A0[0],b.A1[0],b.A2[0] #single values getting output as arrays, silly, but works so not worrying about it
		alphl.append(b.alph)
		chil.append(chi)
		if chi < chim:
			chim = chi
			alphm = b.alph
			Bm = B
			#A0m = b.A0
			#A1m = b.A1
			#A2m = b.A2
	#print alphm,chim,Bm,A0m,A1m,A2m
	#fo = open('BAOisobestfit'+wo+'.dat','w')
	b.alph = alphm
	#if npar > 0:	
	#	b.chi_templ_alphfXX((Bm),wo='y',fw=wo)
	#fo.write(str(alphm)+' '+str(chim)+' '+str(Bm[0])+' '+str(A0m[0])+' '+str(A1m[0])+' '+str(A2m[0])+'\n')
	#fo.close()
	return chil


def doxi_isolike(xid,covd,modl,modsmoothl,rl,rmin=50,rmax=150,npar=3,sp=1.,Bp=.4,rminb=50.,rmaxb=50.,spa=.001,mina=.8,maxa=1.2,chi2fac=1.,Nmock=1000,v='n',wo='',cov2=''):
	#chi2fac should be hartlap factor
	from time import time
	from optimize import fmin
	print np
	b = baofit_iso(xid,covd,modl,modsmoothl,rl,rmin=rmin,rmax=rmax,sp=sp,cov2=cov2)
	b.Bp = Bp
	b.np = npar
	b.H = np.zeros((npar,b.nbin))
	chi2fac = (Nmock-b.nbin-2.)/(Nmock-1.)
	print b.nbin,chi2fac
	for i in range(0,b.nbin):
		for j in range(0,npar):
			b.H[j][i] = 1./b.rl[i]**j
	if rmin == rmaxb:
		rmaxb += (b.rl[1]-b.rl[0])*1.1 #increase rmaxb by one bin size if set poorly
	bb = baofit_iso(xid,covd,modl,modsmoothl,rl,rmin=rmin,rmax=rmaxb,sp=sp,cov2=cov2)
	#bb is to set bias prior
	bb.np = npar
	bb.H = np.zeros((npar,bb.nbin))
	for i in range(0,bb.nbin):
		for j in range(0,npar):
			bb.H[0][i] = 1./bb.rl[i]**j

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
	Bmax = 10.
	while B < Bmax:
		bl = [B]
		chiB = bb.chi_templ_alphfXXn(bl)*chi2fac
		if chiB < chiBmin:
			chiBmin = chiB
			BB = B
		B += .01	
	print 'best-fit bias factor is '+str(BB)+' '+str(chiBmin)
	#print BB,Bmax,Bmax-.01
	if BB >= Bmax-.011:
		print 'WARNING, best-fit bias is at max tested value'
	#else:
	#	print BB,Bmax,Bmax-.01	
	b.BB = BB		
	#b.BB = 1. #switch to this to make bias prior centered on input rather than fit value
	B = BB
	for i in range(0,na):		
		b.alph = mina+spa*i+spa/2.
		#inl = np.array([B,A0,A1,A2])
		#inl = np.array(B)
		inl = B
		#(B,A0,A1,A2) = fmin(b.chi_templ_alphfXX,inl,disp=False)
		if npar > 0:
			B = fmin(b.chi_templ_alphfXX,inl,disp=False)
			chi = b.chi_templ_alphfXX((B))*chi2fac
		else:
			B = fmin(b.chi_templ_alphfXXn,inl,disp=False)
			chi = b.chi_templ_alphfXXn((B))*chi2fac
		#B = fmin(b.chi_templ_alphfXXn,inl,disp=False)
		#chi = b.chi_templ_alphfXX((B,A0,A1,A2))*chi2fac
		
		if v == 'y':
			print b.alph,chi,B[0],b.A0[0],b.A1[0],b.A2[0] #single values getting output as arrays, silly, but works so not worrying about it
		alphl.append(b.alph)
		chil.append(chi)
		if chi < chim:
			chim = chi
			alphm = b.alph
			Bm = B
			#A0m = b.A0
			#A1m = b.A1
			#A2m = b.A2
	#print alphm,chim,Bm,A0m,A1m,A2m
	#fo = open('BAOisobestfit'+wo+'.dat','w')
	b.alph = alphm
	if npar > 0:	
		b.chi_templ_alphfXX((Bm),wo='y',fw=wo)
	#fo.write(str(alphm)+' '+str(chim)+' '+str(Bm[0])+' '+str(A0m[0])+' '+str(A1m[0])+' '+str(A2m[0])+'\n')
	#fo.close()
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
	#http://pdg.lbl.gov/2016/reviews/rpp2016-rev-statistics.pdf eq. 39.22 and surrounding
	ht = H.transpose()
	onei = linalg.pinv(np.dot(np.dot(H,ci),ht))
	comb = np.dot(np.dot(onei,H),ci)
	#print comb,d
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