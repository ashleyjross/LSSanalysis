from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from math import *

def plotfnl(fnl,b):
	fac = 10.e-8
	pk = np.loadtxt('powerspectra/Challenge_matterpower.dat').transpose()
	pl = []
	pln = []
	plm = []
	kl = []
	for i in range(0,len(pk[0])):
		k = pk[0][i]
		if k > 0.0005 and k < 0.1:
			bk = (b+fnl*(b-1.)*fac/k**2.)**2.
			pl.append(bk*pk[1][i])
			pln.append(b*b*pk[1][i])
			plm.append(pk[1][i])
			kl.append(k)
	kl = np.array(kl)
	pl = np.array(pl)		
	plm = np.array(plm)
	#plt.loglog(k,pl)
	#plt.xlim(1.e-4,1.e-1)
	plt.plot(kl*1000.,(pl-pln)/plm)
	plt.xlim(.5,1.e2)
	plt.xscale('log')
	#plt.ylim(1.e-6,3e-5)
	plt.xlabel(r'$k$ ($h$Gpc$^{-1}$)')
	plt.ylabel(r'$[P(k)-P(k,f_{NL}=0)]/P_{\rm matter}(k)$')

	#plt.yscale('log')
	plt.show()
	return True

def plotfnl_b(fnl,kmin=.0008,kmax=.1):
	fac = 10.e-8
	pk = np.loadtxt('powerspectra/Challenge_matterpower.dat').transpose()
	bl = [.2,1.,1.5,2.]
	#b = .2
	#while b < 4:
	for b in bl:
		pl = []
		pln = []
		kl = []
		plm = []
		for i in range(0,len(pk[0])):
			k = pk[0][i]
			if k > kmin and k < kmax:
				bk = (b+fnl*(b-1.)*fac/k**2.)**2.
				pl.append(bk*pk[1][i])
				pln.append(b*b*pk[1][i])
				plm.append(pk[1][i])
				kl.append(k)
		kl = np.array(kl)
		pl = np.array(pl)	
		plm = np.array(plm)
		if b == 1:
			plt.plot(kl*1000.,1.e-2*(kl/0.001)*(pl-pln),'k--')
			#plt.plot(kl*1000.,(pl-pln),'k--')	
		else:
			plt.plot(kl*1000.,1.e-2*(kl/0.001)*(pl-pln))
			#plt.plot(kl*1000.,(pl-pln))
		#b += .2
	#plt.xlim(.5,1.e2)
	plt.xscale('log')
	#plt.ylim(1.e-6,3e-5)
	plt.xlabel(r'$k$ ($h$Gpc$^{-1}$)')
	#plt.ylabel(r'$[P(k)-P(k,f_{NL}=0)]/P_{\rm matter}(k)$')
	plt.ylabel(r'$10^{-2}k[P(k)-P(k,f_{NL}=0)]$ ($h^{-2}$Gpc$^2$)')
	plt.legend(labels=[r'$b=0.2$',r'$b=1$',r'$b=1.5$',r'$b=2$'])
	plt.title(r'$f_{NL,local}=10$')
	#plt.yscale('log')
	plt.show()
	return True


def plotfnlx(fnl,b1,b2):
	fac = 5.e-8
	pk = np.loadtxt('powerspectra/Challenge_matterpower.dat').transpose()
	pl = []
	for i in range(0,len(pk[0])):
		k = pk[0][i]
		bk = (b1+fnl*(b1-1.)*fac/k**2.)*(b2+fnl*(b2-1.)*fac/k**2.)
		pl.append(bk*pk[1][i])
	plt.loglog(pk[0],pl)
	plt.xlim(1.e-4,1.e-1)
	#plt.yscale('log')
	plt.show()
	return True

def plotfnlpsys(fnl,b1,b2):
	fac = 5.e-8
	pk = np.loadtxt('powerspectra/Challenge_matterpower.dat').transpose()
	pl = []
	for i in range(0,len(pk[0])):
		k = pk[0][i]
		b1t = (b1+fnl*(b1-1.)*fac/k**2.)
		bk = (b2+fnl*(b2-1.)*fac/k**2.)
		pl.append(bk*pk[1][i])
	plt.loglog(pk[0],pl)
	plt.xlim(1.e-4,1.e-1)
	#plt.yscale('log')
	plt.show()
	return True
	
def plotGPC():
	pp = PdfPages('PkMpc.pdf')

	pk = np.loadtxt('powerspectra/Challenge_matterpower.dat').transpose()
	plt.loglog(pk[0],pk[1])
	plt.xlim(1.e-3,1.e-1)
	plt.ylim(1.e3,3e4)
	plt.xlabel(r'$k$ ($h$Mpc$^{-1}$)')
	plt.ylabel(r'$P(k)$ ($h^{-3}$Mpc$^3$)')
	pp.savefig()
	pp.close()
	pp = PdfPages('PkGpc.pdf')
	plt.clf()
	plt.loglog(pk[0]*1000.,pk[1]*1.e-9)
	plt.xlim(1.,1.e2)
	plt.ylim(1.e-6,3e-5)
	plt.xlabel(r'$k$ ($h$Gpc$^{-1}$)')
	plt.ylabel(r'$P(k)$ ($h^{-3}$Gpc$^3$)')
	pp.savefig()
	pp.close()
	return True	
		
def plotvolz():
	pp = PdfPages('volz.pdf')
	from Cosmo import distance
	d = distance(.3,.7)
	vl = []
	zl = []
	dz = .01
	z = .01
	while z < 10:
		v = pi*(d.dc(z)/1000.)**3.
		vl.append(v)
		zl.append(z)
		z += dz	
	#plt.loglog(zl,vl)
	plt.plot(zl,vl)
	plt.ylim(1.,1000)
	plt.yscale('log')
	plt.xlabel('Maximum Redshift')
	plt.ylabel(r'Volume for 3/4 of sky (Gpc$h^{-1}$)$^3$')
	pp.savefig()
	pp.close()
	return True	
	
def pksnr(k,v):
	kvol = 4.*pi*k**3.
	return kvol*v/(2.*pi)**3.	
