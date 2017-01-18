from math import *
from numpy import loadtxt as load
import numpy as np
diro = '/Users/ashleyross/DESY1/'

def baoerr(zmin,zmax,sigz,area,num,bias,recon_fac=1.,sig8=0.8,dampz='y',keff=.15,cosm='Planck'):
	#based on Seo & Eisenstein 2007
	from Cosmo import distance
	from numpy import ones
	wo = False
	if num < 1:
		fo = open('Fmufiles/FdaHvsmu_z'+str(zmin)+str(zmax)+'_zerr'+str(sigz)+'_10e3n'+str(10**3.*num)+'_b'+str(bias)+dampz+'.dat','w')
		wo = True
	
	Pbao_list =  [ 9.034, 14.52, 12.63, 9.481, 7.409, 6.397, 5.688, 4.804, 3.841, 3.108,
    2.707, 2.503, 2.300, 2.014, 1.707, 1.473, 1.338, 1.259, 1.174, 1.061,
    0.9409, 0.8435, 0.7792, 0.7351, 0.6915, 0.6398, 0.5851, 0.5376, 0.5018, 0.4741,
    0.4484, 0.4210, 0.3929, 0.3671, 0.3456, 0.3276, 0.3112, 0.2950, 0.2788, 0.2635,
    0.2499, 0.2379, 0.2270, 0.2165, 0.2062, 0.1965, 0.1876, 0.1794, 0.1718, 0.1646]
	if cosm == 'Planck':
		BAO_POWER = 0.18961E+04   # /* The power spectrum at k=0.2h Mpc^-1 for sigma8=0.8 and Planck cosmo */
		BAO_SILK = 7.50#based on footnote in Seo & Eisenstein and 2015 Plik from Table 1 of Planck 2015 cosmo paper
		BAO_AMP = 0.39 #approximate, check this
		om = 0.3
	if cosm == 'WMAP3': #as given in Seo & Eisenstein 2007
		BAO_POWER = 0.17713E+04 #from the MICE_matterpower.dat file
		BAO_SILK = 8.38
		BAO_AMP = 0.5817
		om = 0.25
	d = distance(om,1.-om)
	mustep = 0.01
	KSTEP = .01
	fsky = area/(360*360./pi)
	dtot = 0
	neff = 0
	kl = []
	keffl = []
	for i in range(0,len(Pbao_list)):
		k=0.5*KSTEP+KSTEP*i
		kl.append(k)
		keffl.append(0)

	z = (zmax+zmin)/2.
	print z
	z1 = zmin
	z2 = zmax
	sigzdampl = BAOdampsigz(z,sigz)
	dr = d.dc(z2)-d.dc(z1)
	volume= 4./3.*pi*fsky*(d.dc(z2)**3.-d.dc(z1)**3.)
	if num > 1:
		num = num/volume
		print num
	if dampz == 'n':
		sigzdampl = ones((len(Pbao_list)))
	Dg = d.D(z)
	f = d.omz(z)**.557
	#print f
	beta = f/bias
	#Sig0 = 12.4*sig8/0.9*Dg*.758*recon_fac
	Sig0 = 10.4*sig8*Dg*recon_fac
	Sigma_perp = Sig0
	Sigma_par = Sig0*(1.+f)
	Sigma_perp2 = Sigma_perp*Sigma_perp
	Sigma_par2 = Sigma_par*Sigma_par
#	print Sigma_perp,Sigma_par

	Sigma_z = d.cHz(z)*sigz
	Sigma_zb = Sigma_z/d.dc(z)*105. #percentage distance error multiplied by BAO scale
	#print Sigma_zb
	Sigma_z2 = Sigma_z*Sigma_z
	print Sigma_z2
	sigma8 = bias*Dg
	#print sigma8
	
	power = sigma8*sigma8*BAO_POWER
	#print power,sigma8**2.
	nP = num*power
	#print nP
	Silk_list  = []

	for i in range(0,len(Pbao_list)):
		k=0.5*KSTEP+KSTEP*i
		Silk_list.append(exp(-2.0*pow(k*BAO_SILK,1.40))*k*k*sigzdampl[i]**2.)
	mu = .5*mustep
	sumt = 0
	sumW1 = 0
	sumW2 = 0
	monosum = 0
	Fdd = Fdh = Fhh = 0.0
	while mu<1:
		mu2 = mu*mu
		redshift_distort = (1.+beta*mu2)*(1.+beta*mu2)
		tmp = 1.0/(nP*redshift_distort)
		Sigma2_tot = Sigma_perp2*(1.-mu2)+Sigma_par2*mu2#+Sigma_zb*Sigma_zb/2.
		sum = 0
		for i in range(0,len(Pbao_list)):
			k=0.5*KSTEP+KSTEP*i
			try:
				tmpz = Pbao_list[i]+tmp*exp(k*k*Sigma_z2*mu2)
				#print redshift_distort/exp(k*k*Sigma_z2*mu2)
				Fmu = Silk_list[i]*exp(-k*k*Sigma2_tot)/tmpz/tmpz
				sum += Fmu
				keffl[i] += Fmu
			except:
				pass
				#print k,mu
		neff += exp(-1.*keff**2.*Sigma_z2*mu2)*num*mustep
		if wo:
			fo.write(str(mu)+' '+str(sum)+'\n')
		Fdd += sum*(1.-mu2)*(1.-mu2)
		Fdh += sum*(1.-mu2)*mu2
		Fhh += sum*mu2*mu2
		sumt += sum
		monosum += 1./sum
		if mu < 0.5:
			sumW1 += sum
		if mu > 0.5:
			sumW2 += sum
		mu += mustep
	r = Fdh/sqrt(Fhh*Fdd)
	Fdd *= BAO_AMP*BAO_AMP*KSTEP*mustep*volume 
	Fhh *= BAO_AMP*BAO_AMP*KSTEP*mustep*volume
	sumt *= BAO_AMP*BAO_AMP*KSTEP*mustep*volume
	Drms = 1.0/sqrt(Fdd*(1.0-(r)*(r)))
	Hrms = 1.0/sqrt(Fhh*(1.0-(r)*(r)))
	Rrms = (Drms)*sqrt((1-(r)*(r))/(1+(Drms)/(Hrms)*(2*(r)+(Drms)/(Hrms))))
	print Drms,Hrms,Rrms,r,z1,z2,volume,z
	dtot = sumt
	keff = 0
	wkeff = 0
	for i in range(0,len(keffl)):
		keff += kl[i]*keffl[i]
		wkeff += keffl[i]
	print neff,keff/wkeff
	#print 'total BAO error '+str(sqrt(1.0/dtot))
	return sqrt(1.0/dtot)
	#return Drms,Hrms,Rrms,r,1./sqrt(sumt),sqrt(monosum),1./sqrt(sumW1),1./sqrt(sumW2)

def baoerr_full(zmin,zmax,sigz,area,num,bias,recon_fac=1.,sig8=0.8,vis='n',dampz='y',keff=.15,cosm='Challenge',kmin=.02,kmax=.3):
	#based on Seo & Eisenstein 2007, but using input linear power spectrum as the signal
	#inputs are stored in the powerspectra folder; Challenge is the fiducial BOSS DR12 power spectrum and MICE is the MICE one
	from Cosmo import distance
	from numpy import ones
	from EH import simulate
	#fo = open('Fmufiles/FdaHvsmu_z'+str(zmin)+str(zmax)+'_zerr'+str(sigz)+'_10e3n'+str(10**3.*num)+'_b'+str(bias)+dampz+'.dat','w')
	
	pf = load('powerspectra/'+cosm+'_matterpower.dat').transpose()
	k_list = pf[0]
	P_list = pf[1]
	if cosm == 'Challenge':
		om = 0.31
		lam = 0.69
		h = .676
		nindex = .963
		ombhh = .022
	if cosm == 'MICE': 
		om = 0.25
		lam = .75
		h = .7
		ombhh = .044*0.7*.7	
		nindex = .949
	if cosm == 'WMAP3': 
		om = 0.24
		lam = .76
		h = .73
		ombhh = .0223	
		nindex = .949
	s = simulate(omega=om,lamda=lam,h=h,nindex=nindex,ombhh=ombhh)
	BAO_list = []
	for i in range(0,len(k_list)):
		k = k_list[i]
		dpk = P_list[i]-s.Psmooth(k,0)
		BAO_list.append(dpk)
	if vis == 'y':
		from matplotlib import pyplot as plt
		plt.plot(k_list,BAO_list)
		plt.xlim(kmin,kmax)
		plt.show()	
	BAO_list = np.array(BAO_list)	
	d = distance(om,1.-om)
	mustep = 0.01
	fsky = area/(360*360./pi)
	dtot = 0
	neff = 0
	keffl = np.zeros((len(k_list)))
	z = (zmax+zmin)/2.
	print z
	z1 = zmin
	z2 = zmax
	#sigzdampl = BAOdampsigz(z,sigz)
	volume= 4./3.*pi*fsky*(d.dc(z2)**3.-d.dc(z1)**3.)
	#if dampz == 'n':
	#	sigzdampl = ones((len(Pbao_list)))
	Dg = d.D(z)
	f = d.omz(z)**.557
	beta = f/bias

	#Sig0 = 12.4*sig8/0.9*Dg*.758*recon_fac
	Sig0 = 10.4*sig8*Dg*recon_fac
	Sigma_perp = Sig0
	Sigma_par = Sig0*(1.+f)
	print Sigma_perp,Sigma_par
	Sigma_perp2 = Sigma_perp*Sigma_perp
	Sigma_par2 = Sigma_par*Sigma_par

	Sigma_z = d.cHz(z)*sigz
	Sigma_z2 = Sigma_z*Sigma_z
	print Sigma_z2

	Pkamp = bias*Dg
	power = Pkamp*Pkamp
	BAO_list = BAO_list*power
	P_list = P_list*power
	nP_list = num**P_list
	
	mu = .5*mustep
	sumt = 0
	sumW1 = 0
	sumW2 = 0
	monosum = 0
	Fdd = Fdh = Fhh = 0.0


	k0 = k_list[0]
	k1 = k_list[1]
	ldk = log(k1)-log(k0)

	while mu<1:
		mu2 = mu*mu
		redshift_distort = (1.+beta*mu2)*(1.+beta*mu2)
		#redshift_distort = 1.
		Sigma2_tot = Sigma_perp2*(1.-mu2)+Sigma_par2*mu2#+Sigma_zb*Sigma_zb/2.
		sum = 0
		for i in range(0,len(BAO_list)):		
			k=k_list[i]
			if k > kmin and k < kmax:
				Rmu = redshift_distort*exp(-1.*k*k*Sigma_z2*mu2)
				Rmup = redshift_distort*exp(-1.*k_list[i+1]*k_list[i+1]*Sigma_z2*mu2)
				dk = ldk*k
				kvol = dk*k*k/(4.*pi**2.)
				sig = BAO_list[i]*Rmu*exp(-.5*k*k*Sigma2_tot)
				dsigdlk = (sig-BAO_list[i+1]*Rmup*exp(-.5*k_list[i+1]*k_list[i+1]*Sigma2_tot))/ldk 
				noiseV = (P_list[i]*Rmu+1./num)#/kvol
				#Fmu = (sig/noiseV)**2.*kvol
				Fmu = (dsigdlk/noiseV)**2.*kvol
				sum += Fmu
				keffl[i] += Fmu 
		neff += exp(-1.*keff**2.*Sigma_z2*mu2)*num*mustep
#		fo.write(str(mu)+' '+str(sum)+'\n')
		Fdd += sum*(1.-mu2)*(1.-mu2)
		Fdh += sum*(1.-mu2)*mu2
		Fhh += sum*mu2*mu2
		sumt += sum
		monosum += 1./sum
		if mu < 0.5:
			sumW1 += sum
		if mu > 0.5:
			sumW2 += sum
		mu += mustep
	
	r = Fdh/sqrt(Fhh*Fdd)
	Fdd *= mustep*volume
	Fhh *= mustep*volume
	sumt *= mustep*volume
	print sqrt(sumt)
	Drms = 1.0/sqrt(Fdd*(1.0-(r)*(r)))
	Hrms = 1.0/sqrt(Fhh*(1.0-(r)*(r)))
	Rrms = (Drms)*sqrt((1-(r)*(r))/(1+(Drms)/(Hrms)*(2*(r)+(Drms)/(Hrms))))
	print Drms,Hrms,Rrms,r,z1,z2,volume,z
	dtot = sumt
	keff = 0
	wkeff = 0
	for i in range(0,len(keffl)):
		k = k_list[i]
		if k > kmin and k < kmax:
			keff += k_list[i]*keffl[i]
			wkeff += keffl[i]
	print neff,keff/wkeff
	#print 'total BAO error '+str(sqrt(1.0/dtot))
	return sqrt(1.0/dtot)


def errvsig(num,sigmin=.01,sigmax=.06,sigstep=.001,zmin=.75,zmax=.85,bias=1.8,dampz='y'):
	fo = open('Fmufiles/errvsig'+str(zmin)+str(zmax)+'_10e3n'+str(10**3.*num)+'_b'+str(bias)+dampz+'.dat','w')
	nsig = int((sigmax+sigstep*.01-sigmin)/sigstep)
	for i in range(0,nsig):
		sigz = sigmin+i*sigstep+sigstep/2.
		err = baoerr(zmin,zmax,.1,sigz,1400.,num,bias,dampz=dampz)
		fo.write(str(sigz)+' '+str(1./err**2.)+' '+str(err)+'\n')
	fo.close()
	return True


def baoerr_comp2(z,sigz,sigz2,num,num2,bias,bias2,sig8=0.8,recon_fac=1.):
	from Cosmo import distance
	fo = open('FdaHvsmu_z'+str(z)+'_zerr'+str(sigz)+str(sigz2)+'_10e3n'+str(10**3.*num)+str(10**3.*num2)+'_b'+str(bias)+str(bias2)+'.dat','w')
	d = distance(.3,.7)
	Dg = d.D(z)
	f = d.omz(z)**.557
	beta = f/bias
	beta2 = f/bias2
	Pbao_list =  [ 9.034, 14.52, 12.63, 9.481, 7.409, 6.397, 5.688, 4.804, 3.841, 3.108,
    2.707, 2.503, 2.300, 2.014, 1.707, 1.473, 1.338, 1.259, 1.174, 1.061,
    0.9409, 0.8435, 0.7792, 0.7351, 0.6915, 0.6398, 0.5851, 0.5376, 0.5018, 0.4741,
    0.4484, 0.4210, 0.3929, 0.3671, 0.3456, 0.3276, 0.3112, 0.2950, 0.2788, 0.2635,
    0.2499, 0.2379, 0.2270, 0.2165, 0.2062, 0.1965, 0.1876, 0.1794, 0.1718, 0.1646]
	BAO_POWER = 2875.0   # /* The power spectrum at k=0.2h Mpc^-1 for sigma8=1 */
	BAO_SILK = 7.76
	BAO_AMP = 0.04024
	mustep = 0.01
	Sig0 = 12.4*sig8/0.9*Dg*.758*recon_fac
	Sigma_perp = Sig0
	Sigma_par = Sig0*(1.+f)
	print Sigma_perp,Sigma_par
	volume=1.0
	Sigma_perp2 = Sigma_perp*Sigma_perp
	Sigma_par2 = Sigma_par*Sigma_par
	#Sigma_z = (d.dc(z+0.01)-d.dc(z))/.01*sigz
	Sigma_z = d.cHz(z)*sigz
	Sigma_zn = d.cHz(z)*sigz2
	Sigma_z2 = Sigma_z*Sigma_z
	Sigma_z22 = Sigma_zn*Sigma_zn
	Sigma_zb = Sigma_z/d.dc(z)*105.
	Sigma_zb2 = Sigma_zn/d.dc(z)*105.
	print Sigma_z2
	sigma8 = sig8*bias*Dg
	sigma82 = sig8*bias2*Dg
	KSTEP = .01
	nP = num*sigma8*sigma8*BAO_POWER
	nP2 = num2*sigma82*sigma82*BAO_POWER
	Silk_list  = []
	
	for i in range(0,len(Pbao_list)):
		k=0.5*KSTEP+KSTEP*i
		Silk_list.append(exp(-2.0*pow(k*BAO_SILK,1.40))*k*k)
	mu = .5*mustep
	sumt = 0
	sumt2 = 0
	sumtt = 0

	while mu<1:
		mu2 = mu*mu
		redshift_distort = (1+beta*mu2)*(1+beta*mu2)
		redshift_distort2 = (1+beta2*mu2)*(1+beta2*mu2)
		tmp = 1.0/(nP*redshift_distort)
		tmp2 = 1.0/(nP2*redshift_distort2)
		Sigma2_tot = Sigma_perp2*(1-mu2)+Sigma_par2*mu2+Sigma_zb*Sigma_zb/2.
		Sigma2_tot2 = Sigma_perp2*(1-mu2)+Sigma_par2*mu2+Sigma_zb2*Sigma_zb2/2.
		sum = 0
		sum2 = 0
		sumkt = 0
		for i in range(0,len(Pbao_list)):
			k=0.5*KSTEP+KSTEP*i
			try:
				tmpz = Pbao_list[i]+tmp*exp(k*k*Sigma_z2*mu2)
				Fmu = Silk_list[i]*exp(-k*k*Sigma2_tot)/tmpz/tmpz
				sum += Fmu
			except:
				Fmu = 0
			try:
				tmpz2 = Pbao_list[i]+tmp2*exp(k*k*Sigma_z22*mu2)
				Fmu2 = Silk_list[i]*exp(-k*k*Sigma2_tot2)/tmpz2/tmpz2
				sum2 += Fmu2
			except:
				Fmu2 = 0
				#print k,mu
			try:
				C12 = 1./(1.+tmp*exp(k*k*Sigma_z2*mu2))/(1.+tmp2*exp(k*k*Sigma_z22*mu2))
			except:
				C12 = 0
			sumkt += 1./(1.-C12**2.)*(Fmu+Fmu2-2.*sqrt(Fmu*Fmu2)*C12)		
		fo.write(str(mu)+' '+str(sum)+' '+str(sum2)+' '+str(sumkt)+'\n')
		sumt += sum
		sumt2 += sum2
		sumtt += sumkt
		mu += mustep
	sumt *= BAO_AMP*BAO_AMP/8.0/pi/pi*1.0e9*KSTEP*mustep*volume
	sumt2 *= BAO_AMP*BAO_AMP/8.0/pi/pi*1.0e9*KSTEP*mustep*volume
	sumtt *= BAO_AMP*BAO_AMP/8.0/pi/pi*1.0e9*KSTEP*mustep*volume
	return 1./sqrt(sumt),1./sqrt(sumt2),1./sqrt(sumtt),1./sqrt(sumt+sumt2)

def BAOdampsigz(z,sigz,nthbin=50,dz=.01,zmax=3.,vis='n'):
	from Cosmo import distance
	from numpy import ones
	d = distance(.3,.7)
	disz = d.dc(z)
	cl  =[]
	cl0 = []
	thl  = []
	dth = 8*2.*pi/float(nthbin)
	for i in range(0,nthbin):
		cl.append(0)
		cl0.append(cos(dth/2.+dth*i))
		thl.append(dth/2.+dth*i)
	nzbin = int((zmax)/dz)
	sum = 0
	
	for i in range(0,nthbin):
		th = dth/2.+dth*i
		c = 0
		for j in range(0,nzbin):
			zj = dz/2.+dz*j
			c += dz/sqrt(2.*pi*sigz**2.)*exp(-1.*(z-zj)**2./(2.*sigz**2.))*cos(th*d.dc(zj)/disz)
		cl[i] = c
		sum += c*c
	sum = sum*dth/pi/3.
	if vis == 'y':
		import matplotlib.pyplot as plt
		plt.plot(thl,cl0,'k-',thl,cl,'r-')
		plt.show()
	outl = []	
	for i in range(0,nthbin):
		outl.append(cl[i]/cl0[i])
	if sigz < 0.01 and sum < .1:
		outl = ones((nthbin))
	return outl
	

def plot_Fvmu():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('Fvmusigz_simBOSS.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$\mu$',size=16)
	plt.ylabel(r'$F(\mu)$',size=16)
	d0 = load('FdaHvsmu_z0.6_zerr0.001_10e3n0.2_b2.0.dat').transpose()
	d1 = load('FdaHvsmu_z0.6_zerr0.01_10e3n0.2_b2.0.dat').transpose()
	d2 = load('FdaHvsmu_z0.6_zerr0.02_10e3n0.2_b2.0.dat').transpose()
	d3 = load('FdaHvsmu_z0.6_zerr0.03_10e3n0.2_b2.0.dat').transpose()
	d5 = load('FdaHvsmu_z0.6_zerr0.05_10e3n0.2_b2.0.dat').transpose()
	d10 = load('FdaHvsmu_z0.6_zerr0.1_10e3n0.2_b2.0.dat').transpose()
	plt.yscale('log')
	plt.plot(d0[0],d0[1]/d0[1][0],'k-')
	plt.plot(d0[0],1./d1[1][0]*d1[1],'r-')
	plt.plot(d0[0],1./d2[1][0]*d2[1],'-',color='orange')
	plt.plot(d0[0],1./d3[1][0]*d3[1],'-',color='yellow')  
	plt.plot(d0[0],1./d5[1][0]*d5[1],'g-')  
	plt.plot(d0[0],1./d10[1][0]*d10[1],'b-')
	plt.show()
	return True

def mudistobszerr(r=100.,z=0.8,sigz=0.029,bias=1.5,rmin=10.,rmax=300,muww='0',a='',v='y',gam=-1.7,file='MICE_matterpower',dir='',mun=0,beta=0.4,sfog=3.0,amc=0.0,sigt=6.,sigr=10.,mult=1.,sigs=15.,mumin=0,mumax=1):
	from random import gauss
	from numpy import zeros
	spf = 1.
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
	nz = 40000
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
	fmuw = load('Fmufiles/FdaHvsmu_z0.750.85_zerr0.03_10e3n1.0_b1.8n.dat').transpose()
	#muwt = sum(fmuw[1])
	muwt = 0
	#for i in range(0, len(fmuw[1])):
	#	if i > int(mumin*100) and i < int(mumax*100):
	mmu = zeros((100*(mumax-mumin),100*(mumax-mumin)))
	for i in range(int(100*mumin),int(100*mumax)):	
		if muww == 'muw':
			muwt += fmuw[1][i]
		else:
			muwt += .01
				
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	fo = open('/Users/ashleyross/DESY1/muobsmutrue'+muww+file+muw+str(sigz)+'.dat','w')
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
			rtp = rt*rzl[i]
			rp = sqrt(rrp**2.+rtp**2.)
			mup = abs(rrp/rp)
			muind = int(mup*100)
			mmu[m][muind] += wl[i]
	for i in range(0,len(mmu)):
		for j in range(0,len(mmu[i])):
			fo.write(str(mmu[i][j])+' ')
		fo.write('\n')		
	fo.close()
	fo = open('Fmufiles/FmuobsdaHvsmu_z0.750.85_zerr0.03_10e3n1.0_b1.8n.dat','w')
	for i in range(0,len(mmu)):
		mu = .005+.01*i
		smu = 0
		for j in range(0,len(mmu[i])):
			smu += fmuw[1][j]*mmu[i][j]
		fo.write(str(mu)+' '+str(smu)+'\n')
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


def plot_Fvmusig3():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('/Users/ashleyross/Dropbox/zerrBAOpaper/Fvmusigz3prmld.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$\mu$',size=16)
	plt.ylabel(r'$F(\mu)$',size=16)
	#d0 = load('FdaHvsmu_z0.60.7_zerr0.03_10e3n1.5_b1.8n.dat').transpose()
	#dc = load('FdaHvsmu_z0.60.7_zerr0.001_10e3n1.5_b1.8n.dat').transpose()
	d0 = load('Fmufiles/FdaHvsmu_z0.750.85_zerr0.03_10e3n1.0_b1.8n.dat').transpose()
	dc = load('Fmufiles/FdaHvsmu_z0.750.85_zerr0.001_10e3n1.0_b1.8n.dat').transpose()	
	ds = load('Fmufiles/FmuobsdaHvsmu_z0.750.85_zerr0.03_10e3n1.0_b1.8n.dat').transpose()
	dr = load('Fmufiles/FdaHvsmu_z0.750.85_zerr0.02_10e3n0.4_b2.0n.dat').transpose()
	dl = load('Fmufiles/FdaHvsmu_z0.750.85_zerr0.001_10e3n0.01_b1.8n.dat').transpose()	
	#plt.yscale('log')
	plt.plot(d0[0],d0[1]/d0[1][0],'k-',linewidth=4)
	plt.text(.012,.28,r'$\sigma_{zf}=0.03$',color='k',size=14)
	plt.text(.012,.33,r'$n=10^{-3}$',color='k',size=14)
	plt.plot(dr[0],dr[1]/d0[1][0],'r-',linewidth=4)
	plt.text(.2,.33,r'$\sigma_{zf}=0.02$',color='r',size=14)
	plt.text(.2,.38,r'$n=4\times10^{-4}$',color='r',size=14)
	plt.plot(dl[0],dl[1]/dl[1][-1]*1.15,'b--',linewidth=4)
	plt.text(.83,.95,r'$\sigma_{zf}=0.001$',color='b',size=14)
	plt.text(.85,1.,r'$n=10^{-5}$',color='b',size=14)
	plt.plot(d0[0],dc[1]/d0[1][0],'k--',linewidth=4)
	plt.text(.25,1.05,r'$\sigma_{zf}=0.001$',color='k',size=14)
	plt.text(.25,1.1,r'$n=10^{-3}$',color='k',size=14)
	plt.plot(d0[0],ds[1]*sum(d0[1]/d0[1][0])/sum(ds[1]),'k:',linewidth=4)
	plt.text(.55,.2,r'Apparent $\mu$, $n=10^{-3}$',size=14)
	plt.text(.73,.15,r'$\sigma_{zf}=0.03$',color='k',size=14)
	plt.ylabel('Relative BAO Information', fontsize=18)
	plt.xlabel(r'$\mu$ (Cosine of angle to line-of-sight)',fontsize=16)
	pp.savefig()
	pp.close()
	return True

def xisigmuplot():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('/Users/ashleyross/DESY1/xisigmuscale.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s_{\perp}$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$10^3(\xi_{\rm BAO} - \xi_{\rm no BAO})$',size=16)
	d0 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumax0.20.406.010.00.029sp1.0.dat').transpose()
	#dt = load('xizconvmuwMICE_matterpower0.406.010.00.029sp5.0.dat').transpose()
	drp = load('/Users/ashleyross/DESY1/xirpzconvmuwMICE_matterpower0.406.010.00.029sp5.0.dat').transpose()
	d1 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.2mumax0.40.406.010.00.029sp1.0.dat').transpose()
	d2 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.4mumax0.60.406.010.00.029sp1.0.dat').transpose()
	d3 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.6mumax0.80.406.010.00.029sp1.0.dat').transpose()
	d4 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.80.406.010.00.029sp1.0.dat').transpose()
 	#plt.plot(d0[0]*sin(acos(.1)),(d0[1]-d0[2])*1000,'k-',linewidth=3)
 	#plt.plot(d0[0]*sin(acos(.3)),(d1[1]-d1[2])*1000,'k--',linewidth=3)
 	#plt.plot(d0[0]*sin(acos(.5)),(d2[1]-d2[2])*1000,'k:',linewidth=3)
 	#plt.plot(d0[0]*sin(acos(.7)),(d3[1]-d3[2])*1000,':',color='r',linewidth=3)
 	#plt.plot(d0[0]*sin(acos(.9)),(d4[1]-d4[2])*1000,':',color='b',linewidth=3)
	plt.plot(d0[0]*sqrt(1.-.1**2.),(d0[1]-d0[2])*1000,'k-',linewidth=3)
	xl = [104,104]
	yl = [-.5,.5]
	plt.plot(xl,yl,'k:')
	xl = [130,140]
	yl = [.36,.36]
	plt.plot(xl,yl,'k-',linewidth=3)
	plt.text(150,.35,r'$\mu < 0.2$',color='k')
	#plt.plot(dt[0],(dt[1]-dt[2])*1000,'y-',linewidth=3)
	#plt.plot(drp[0],(drp[1]-drp[2])*1000,'-',color='purple',linewidth=3)
	plt.plot(d0[0]*sqrt(1.-.3**2.),(d1[1]-d1[2])*1000,'k--',linewidth=3)
	yl = [.33,.33]
	plt.plot(xl,yl,'k--',linewidth=3)
	plt.text(150,.32,r'$0.2 < \mu < 0.4$',color='k')
	plt.plot(d0[0]*sqrt(1.-.5**2.),(d2[1]-d2[2])*1000,'k:',linewidth=3)
	yl = [.29,.29]
	plt.plot(xl,yl,'k:',linewidth=3)
	plt.text(150,.28,r'$0.4 < \mu < 0.6$',color='k')
	plt.plot(d0[0]*sqrt(1.-.7**2.),(d3[1]-d3[2])*1000,'-',color='r',linewidth=3)
	yl = [.26,.26]
	plt.plot(xl,yl,'r-',linewidth=3)
	plt.text(150,.25,r'$0.6 < \mu < 0.8$',color='k')
	plt.plot(d0[0]*sqrt(1.-.9**2.),(d4[1]-d4[2])*1000,'-',color='b',linewidth=3)
	yl = [.23,.23]
	plt.plot(xl,yl,'b-',linewidth=3)
	plt.text(150,.22,r'$0.8 < \mu $',color='k')

	plt.ylim(-.2,.4)
	plt.xlim(30,200)
	pp.savefig()
	pp.close()
	return True

def xisigmuplotzerr(zerr):
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('/Users/ashleyross/DESY1/xisigmu'+str(zerr)+'.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s$ ($h^{-1}$ Mpc)',size=16)	
	plt.ylabel(r'$10^3(\xi_{\rm BAO} - \xi_{\rm no BAO})$',size=16)
	if zerr == 0.029:
		d0 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumax0.20.406.010.00.029sp1.0.dat').transpose()
		#dt = load('xizconvmuwMICE_matterpower0.406.010.00.029sp5.0.dat').transpose()
		drp = load('/Users/ashleyross/DESY1/xirpzconvmuwMICE_matterpower0.406.010.00.029sp5.0.dat').transpose()
		d1 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.2mumax0.40.406.010.00.029sp1.0.dat').transpose()
		d2 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.4mumax0.60.406.010.00.029sp1.0.dat').transpose()
		d3 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.6mumax0.80.406.010.00.029sp1.0.dat').transpose()
		d4 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.80.406.010.00.029sp1.0.dat').transpose()
	else:	
		d0 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumax0.20.406.010.0combzsiglsp1.0.dat').transpose()
		#dt = load('xizconvmuwMICE_matterpower0.406.010.00.029sp5.0.dat').transpose()
		d1 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.2mumax0.40.406.010.0combzsiglsp1.0.dat').transpose()
		d2 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.4mumax0.60.406.010.0combzsiglsp1.0.dat').transpose()
		d3 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.6mumax0.80.406.010.0combzsiglsp1.0.dat').transpose()
		d4 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.80.406.010.0combzsiglsp1.0.dat').transpose()
 	plt.plot(d0[0],(d0[1]-d0[2])*1000,'k-',linewidth=3)
 	plt.plot(d0[0],(d1[1]-d1[2])*1000,'k--',linewidth=3)
 	plt.plot(d0[0],(d2[1]-d2[2])*1000,'k:',linewidth=3)
 	plt.plot(d0[0],(d3[1]-d3[2])*1000,'-',color='r',linewidth=3)
 	plt.plot(d0[0],(d4[1]-d4[2])*1000,'-',color='b',linewidth=3)
	plt.ylim(-.2,.4)
	if zerr == 0.02:
		plt.ylim(-.2,.5)
	plt.xlim(30,300)

	xl = [190,200]
	yl = [.8,.8]
	plt.plot(xl,yl,'k-',linewidth=3)
	plt.text(210,.79,r'$\mu < 0.2$',color='k')
	#plt.plot(dt[0],(dt[1]-dt[2])*1000,'y-',linewidth=3)
	#plt.plot(drp[0],(drp[1]-drp[2])*1000,'-',color='purple',linewidth=3)
	yl = [.74,.74]
	plt.plot(xl,yl,'k--',linewidth=3)
	plt.text(210,.73,r'$0.2 < \mu < 0.4$',color='k')
	
	yl = [.68,.68]
	plt.plot(xl,yl,'k:',linewidth=3)
	plt.text(210,.67,r'$0.4 < \mu < 0.6$',color='k')
	
	yl = [.62,.62]
	plt.plot(xl,yl,'r-',linewidth=3)
	plt.text(210,.61,r'$0.6 < \mu < 0.8$',color='k')
	
	yl = [.56,.56]
	plt.plot(xl,yl,'b-',linewidth=3)
	plt.text(210,.55,r'$0.8 < \mu $',color='k')
	plt.title(r'$\sigma_z/(1+z)$ = '+str(zerr))
	pp.savefig()
	pp.close()
	d0 = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumax0.20.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	#dt = load('xizconvmuwMICE_matterpower0.406.010.00.029sp5.0.dat').transpose()
	d1 = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.2mumax0.40.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d2 = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.4mumax0.60.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d3 = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.6mumax0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d4 = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()

	pp = PdfPages('/Users/ashleyross/DESY1/xisigmuscale'+str(zerr)+'.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s_{\perp}$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$10^3(\xi_{\rm BAO} - \xi_{\rm no BAO})$',size=16)

# 	plt.plot(d0[0]*sqrt(1.-.1**2.),(d0[1]-d0[2])*1000,'k-',linewidth=3)
# 	plt.plot(d0[0]*sqrt(1.-.3**2.),(d1[1]-d1[2])*1000,'k--',linewidth=3)
# 	plt.plot(d0[0]*sqrt(1.-.5**2.),(d2[1]-d2[2])*1000,'k:',linewidth=3)
# 	plt.plot(d0[0]*sqrt(1.-.7**2.),(d3[1]-d3[2])*1000,'-',color='r',linewidth=3)
# 	plt.plot(d0[0]*sqrt(1.-.9**2.),(d4[1]-d4[2])*1000,'-',color='b',linewidth=3)
	plt.plot(d0[0],(d0[1]-d0[2])*1000,'k-',linewidth=3)
	plt.plot(d0[0],(d1[1]-d1[2])*1000,'k--',linewidth=3)
	plt.plot(d0[0],(d2[1]-d2[2])*1000,'k:',linewidth=3)
	plt.plot(d0[0],(d3[1]-d3[2])*1000,'-',color='r',linewidth=3)
	plt.plot(d0[0],(d4[1]-d4[2])*1000,'-',color='b',linewidth=3)

	xl = [104,104]
	yl = [-20,20]
	plt.plot(xl,yl,'k:')

	plt.ylim(-.2,.4)
	if zerr == 0.02:
		plt.ylim(-.2,.5)
	plt.xlim(30,200)
	pp.savefig()
	pp.close()
	return True

def xisigmuplot6pan():
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	from pylab import *
	from numpy import loadtxt as load
	import numpy as np
	pp = PdfPages('/Users/ashleyross/DESY1/xisigmu6pan.pdf')
	
	fig = plt.figure(figsize=(8.5,8))
	zerr = 0.01
	#fig.tick_params(axis='both', which='major', labelsize=12)
	mpl.rcParams['xtick.labelsize'] = 10
	mpl.rcParams['ytick.labelsize'] = 9
	d0 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumax0.20.406.010.0combzsiglsp1.0.dat').transpose()
	d1 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.2mumax0.40.406.010.0combzsiglsp1.0.dat').transpose()
	d2 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.4mumax0.60.406.010.0combzsiglsp1.0.dat').transpose()
	d3 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.6mumax0.80.406.010.0combzsiglsp1.0.dat').transpose()
	d4 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.80.406.010.0combzsiglsp1.0.dat').transpose()
	d0p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumax0.20.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d1p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.2mumax0.40.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d2p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.4mumax0.60.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d3p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.6mumax0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d4p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()

	xl = [104,104]
	yl = [-20,20]

	ax = fig.add_subplot(2,3,1)
	cnr = [145,289,.35,.9]
	cxl = [cnr[0],cnr[1]]
	cyl = [cnr[2],cnr[2]]
	ax.plot(cxl,cyl,'k-')
	cyl = [cnr[3],cnr[3]]
	ax.plot(cxl,cyl,'k-')
	cxl = [cnr[0],cnr[0]]
	cyl = [cnr[2],cnr[3]]
	ax.plot(cxl,cyl,'k-')
	cxl = [cnr[1],cnr[1]]
	ax.plot(cxl,cyl,'k-')
	lxl = [150,169]
	lxll = [145,175]
	ax.text(173,.8,r'$\mu<0.2$',size=10)
	lyl = [.83,.83]
	ax.plot(lxl,lyl,'k-',linewidth=3)
	ax.text(173,.7,r'$0.2<\mu<0.4$',size=10)
	lyl = [.73,.73]
	ax.plot(lxll,lyl,'k--',linewidth=3)
	ax.text(173,.6,r'$0.4<\mu<0.6$',size=10)
	lyl = [.63,.63]
	ax.plot(lxl,lyl,'k:',linewidth=3)
	ax.text(173,.5,r'$0.6<\mu<0.8$',size=10)
	lyl = [.53,.53]
	ax.plot(lxl,lyl,'r-',linewidth=3)
	ax.text(173,.4,r'$\mu>0.8$',size=10)
	lyl = [.43,.43]
	ax.plot(lxl,lyl,'b-',linewidth=3)

	
	ax.set_xlim(40,300)
	ax.set_ylim(-.5,1.)
	ax.minorticks_on()
	ax.plot(xl,yl,'k:')
 	ax.plot(d0[0],(d0[1]-d0[2])*1000,'k-',linewidth=3)
 	ax.plot(d0[0],(d1[1]-d1[2])*1000,'k--',linewidth=3)
 	ax.plot(d0[0],(d2[1]-d2[2])*1000,'k:',linewidth=3)
 	ax.plot(d0[0],(d3[1]-d3[2])*1000,'-',color='r',linewidth=3)
 	ax.plot(d0[0],(d4[1]-d4[2])*1000,'-',color='b',linewidth=3)
	
	ax.set_xlabel(r'$s$ ($h^{-1}$ Mpc)',size=12)
	ax.set_title(r'$\sigma_z /(1+z) = 0.01$')

	ax2 = fig.add_subplot(2,3,4)
	ax2.minorticks_on()
# 	for tick in ax2.xaxis.get_major_ticks():
	ax2.set_ylabel(r'                                                                        $10^3(\xi_{\rm BAO} - \xi_{\rm no BAO})$',size=12)
	ax2.set_xlim(30,200)
	ax2.set_ylim(-.5,1.)
	ax2.plot(xl,yl,'k:')
 	ax2.plot(d0p[0],(d0p[1]-d0p[2])*1000,'k-',linewidth=3)
 	ax2.plot(d0p[0],(d1p[1]-d1p[2])*1000,'k--',linewidth=3)
 	ax2.plot(d0p[0],(d2p[1]-d2p[2])*1000,'k:',linewidth=3)
 	ax2.plot(d0p[0],(d3p[1]-d3p[2])*1000,'-',color='r',linewidth=3)
 	ax2.plot(d0p[0],(d4p[1]-d4p[2])*1000,'-',color='b',linewidth=3)	
	ax2.set_xlabel(r'$s_{\perp}$ ($h^{-1}$ Mpc)',size=12)
	start, end = ax2.get_xlim()
	ax2.xaxis.set_ticks(np.arange(start, end, 40))

	zerr = 0.02
	d0 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumax0.20.406.010.0combzsiglsp1.0.dat').transpose()
	d1 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.2mumax0.40.406.010.0combzsiglsp1.0.dat').transpose()
	d2 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.4mumax0.60.406.010.0combzsiglsp1.0.dat').transpose()
	d3 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.6mumax0.80.406.010.0combzsiglsp1.0.dat').transpose()
	d4 = load('/Users/ashleyross/DESY1/xizconvcsigz'+str(zerr)+'MICE_matterpowermumin0.80.406.010.0combzsiglsp1.0.dat').transpose()
	d0p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumax0.20.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d1p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.2mumax0.40.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d2p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.4mumax0.60.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d3p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.6mumax0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d4p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()

	ax3 = fig.add_subplot(2,3,2)
	ax3.minorticks_on()
	ax3.set_xlim(40,300)
	ax3.set_ylim(-.2,.5)
	ax3.plot(xl,yl,'k:')
 	ax3.plot(d0[0],(d0[1]-d0[2])*1000,'k-',linewidth=3)
 	ax3.plot(d0[0],(d1[1]-d1[2])*1000,'k--',linewidth=3)
 	ax3.plot(d0[0],(d2[1]-d2[2])*1000,'k:',linewidth=3)
 	ax3.plot(d0[0],(d3[1]-d3[2])*1000,'-',color='r',linewidth=3)
 	ax3.plot(d0[0],(d4[1]-d4[2])*1000,'-',color='b',linewidth=3)	
	ax3.set_xlabel(r'$s$ ($h^{-1}$ Mpc)',size=12)
	ax3.set_title(r'$\sigma_z /(1+z) = 0.02$')

	ax4 = fig.add_subplot(2,3,5)
	ax4.minorticks_on()
	ax4.set_xlim(30,200)
	ax4.set_ylim(-.2,.5)
	ax4.plot(xl,yl,'k:')
 	ax4.plot(d0p[0],(d0p[1]-d0p[2])*1000,'k-',linewidth=3)
 	ax4.plot(d0p[0],(d1p[1]-d1p[2])*1000,'k--',linewidth=3)
 	ax4.plot(d0p[0],(d2p[1]-d2p[2])*1000,'k:',linewidth=3)
 	ax4.plot(d0p[0],(d3p[1]-d3p[2])*1000,'-',color='r',linewidth=3)
 	ax4.plot(d0p[0],(d4p[1]-d4p[2])*1000,'-',color='b',linewidth=3)	
	ax4.set_xlabel(r'$s_{\perp}$ ($h^{-1}$ Mpc)',size=12)
	start, end = ax4.get_xlim()
	ax4.xaxis.set_ticks(np.arange(start, end, 40))


	zerr = 0.029
	d0 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumax0.20.406.010.00.029sp1.0.dat').transpose()
	d1 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.2mumax0.40.406.010.00.029sp1.0.dat').transpose()
	d2 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.4mumax0.60.406.010.00.029sp1.0.dat').transpose()
	d3 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.6mumax0.80.406.010.00.029sp1.0.dat').transpose()
	d4 = load('/Users/ashleyross/DESY1/xizconv0MICE_matterpowermumin0.80.406.010.00.029sp1.0.dat').transpose()
	d0p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumax0.20.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d1p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.2mumax0.40.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d2p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.4mumax0.60.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d3p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.6mumax0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()
	d4p = load('/Users/ashleyross/DESY1/xizconvcrpMICE_matterpowermumin0.80.406.010.0combzsigl'+str(zerr)+'sp1.0.dat').transpose()


	ax5 = fig.add_subplot(2,3,3)
	ax5.minorticks_on()
	ax5.set_xlim(40,300)
	ax5.set_ylim(-.15,.35)
	ax5.plot(xl,yl,'k:')
 	ax5.plot(d0[0],(d0[1]-d0[2])*1000,'k-',linewidth=3)
 	ax5.plot(d0[0],(d1[1]-d1[2])*1000,'k--',linewidth=3)
 	ax5.plot(d0[0],(d2[1]-d2[2])*1000,'k:',linewidth=3)
 	ax5.plot(d0[0],(d3[1]-d3[2])*1000,'-',color='r',linewidth=3)
 	ax5.plot(d0[0],(d4[1]-d4[2])*1000,'-',color='b',linewidth=3)	
	ax5.set_xlabel(r'$s$ ($h^{-1}$ Mpc)',size=12)
	ax5.set_title(r'$\sigma_z /(1+z) = 0.029$')

	ax6 = fig.add_subplot(2,3,6)
	ax6.minorticks_on()
	ax6.set_xlim(30,200)
	ax6.set_ylim(-.15,.35)
	ax6.plot(xl,yl,'k:')
 	ax6.plot(d0p[0],(d0p[1]-d0p[2])*1000,'k-',linewidth=3)
 	ax6.plot(d0p[0],(d1p[1]-d1p[2])*1000,'k--',linewidth=3)
 	ax6.plot(d0p[0],(d2p[1]-d2p[2])*1000,'k:',linewidth=3)
 	ax6.plot(d0p[0],(d3p[1]-d3p[2])*1000,'-',color='r',linewidth=3)
 	ax6.plot(d0p[0],(d4p[1]-d4p[2])*1000,'-',color='b',linewidth=3)	
	ax6.set_xlabel(r'$s_{\perp}$ ($h^{-1}$ Mpc)',size=12)
	start, end = ax6.get_xlim()
	ax6.xaxis.set_ticks(np.arange(start, end, 40))
	fig.subplots_adjust(wspace=.25)
	pp.savefig()
	pp.close()
	return True


def xipzcompspeczplot():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xipzcompspecz.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$10^3(\xi_{\rm BAO} - \xi_{\rm no BAO})$',size=16)
	d0 = load('xizconv0MICE_matterpowermumax0.20.406.010.00.029sp1.0.dat').transpose()
	ds = load('xi0MICE_matterpower0.406.010.015.00.dat').transpose()
	dsm = load('xi0smMICE_matterpower0.406.010.015.00.dat').transpose()
	#dt = load('xizconvmuwMICE_matterpower0.406.010.00.029sp5.0.dat').transpose()
	plt.plot(d0[0],(d0[1]-d0[2])*1000,'k-',linewidth=3)
	plt.plot(ds[0],(ds[1]-dsm[1])*200+4./ds[0],'r-',linewidth=3)
	xl = [150,160]
	yl = [.36,.36]
	plt.plot(xl,yl,'k-',linewidth=3)
	plt.text(170,.35,r'$\mu < 0.2$',color='k')
	plt.text(170,.3,r'1/5 spec z',color='r')
	plt.ylim(-.2,.4)
	plt.xlim(30,200)
	pp.savefig()
	pp.close()
	return True


def xisigmumockplot():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xisigmumock.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$\xi_{\rm BAO} - \xi_{\rm no BAO}$',size=16)
	d0 = load('xizconvmuwMICE_matterpowermumax0.10.43.06.010.00.029sp5.0.dat').transpose()
	d1 = load('xizconvmuwMICE_matterpowermumin0.1mumax0.20.43.06.010.00.029sp5.0.dat').transpose()
	d2 = load('xizconvmuwMICE_matterpowermumin0.2mumax0.30.43.06.010.00.029sp5.0.dat').transpose()
	d3 = load('xizconvmuwMICE_matterpowermumin0.3mumax0.40.43.06.010.00.029sp5.0.dat').transpose()
	d4 = load('xizconvmuwMICE_matterpowermumin0.4mumax0.50.43.06.010.00.029sp5.0.dat').transpose()
	dm0 = load('xiavemuwmumax0.1SM.dat').transpose()
	dm1 = load('xiavemuwmumin0.1mumax0.2SM.dat').transpose()
	dm2 = load('xiavemuwmumin0.2mumax0.3SM.dat').transpose()
	dm3 = load('xiavemuwmumin0.3mumax0.4SM.dat').transpose()
	dm4 = load('xiavemuwmumin0.4mumax0.5SM.dat').transpose()
	plt.errorbar(d0[0][:40],dm0[1]-d0[2][:40]*dm0[1][10]/d0[1][10],dm0[2]/25.,fmt='k-',linewidth=3)
	plt.errorbar(d0[0][:40],dm1[1]-d1[2][:40]*dm0[1][10]/d0[1][10],dm1[2]/25.,fmt='k--',linewidth=3)
	plt.errorbar(d0[0][:40],dm2[1]-d2[2][:40]*dm0[1][10]/d0[1][10],dm2[2]/25.,fmt='k:',linewidth=3)
	plt.errorbar(d0[0][:40],dm3[1]-d3[2][:40]*dm0[1][10]/d0[1][10],dm3[2]/25.,fmt=':',color='r',linewidth=3)
	plt.errorbar(d0[0][:40],dm4[1]-d4[2][:40]*dm0[1][10]/d0[1][10],dm4[2]/25.,fmt=':',color='b',linewidth=3)
# 	plt.plot(d0[0][:40],dm0[1]-d0[2][:40],'k-',linewidth=3)
# 	plt.plot(d0[0][:40],dm1[1]-d1[2][:40],'k--',linewidth=3)
# 	plt.plot(d0[0][:40],dm2[1]-d2[2][:40],'k:',linewidth=3)
# 	plt.plot(d0[0][:40],dm3[1]-d3[2][:40],':',color='r',linewidth=3)
# 	plt.plot(d0[0][:40],dm4[1]-d4[2][:40],':',color='b',linewidth=3)
	plt.ylim(-.001,.001)
	pp.savefig()
	pp.close()
	return True


def xisigmumockplotscale():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xisigmumockscale.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$\xi_{\rm BAO} - \xi_{\rm no BAO}$',size=16)
	d0 = load('xizconvmuwMICE_matterpowermumax0.10.43.06.010.00.029sp5.0.dat').transpose()
	d1 = load('xizconvmuwMICE_matterpowermumin0.1mumax0.20.43.06.010.00.029sp5.0.dat').transpose()
	d2 = load('xizconvmuwMICE_matterpowermumin0.2mumax0.30.43.06.010.00.029sp5.0.dat').transpose()
	d3 = load('xizconvmuwMICE_matterpowermumin0.3mumax0.40.43.06.010.00.029sp5.0.dat').transpose()
	d4 = load('xizconvmuwMICE_matterpowermumin0.4mumax0.50.43.06.010.00.029sp5.0.dat').transpose()
	dm0 = load('xiavemuwmumax0.1SM.dat').transpose()
	dm1 = load('xiavemuwmumin0.1mumax0.2SM.dat').transpose()
	dm2 = load('xiavemuwmumin0.2mumax0.3SM.dat').transpose()
	dm3 = load('xiavemuwmumin0.3mumax0.4SM.dat').transpose()
	dm4 = load('xiavemuwmumin0.4mumax0.5SM.dat').transpose()
	fac = (1.1,1.3,1.5,1.7,1.9)
	plt.errorbar(d0[0][:40]/fac[0],dm0[1]-d0[2][:40]*dm0[1][10]/d0[1][10],dm0[2]/25.,fmt='k-',linewidth=3)
	plt.errorbar(d0[0][:40]/fac[1],dm1[1]-d1[2][:40]*dm0[1][10]/d0[1][10],dm1[2]/25.,fmt='k--',linewidth=3)
	plt.errorbar(d0[0][:40]/fac[2],dm2[1]-d2[2][:40]*dm0[1][10]/d0[1][10],dm2[2]/25.,fmt='k:',linewidth=3)
	plt.errorbar(d0[0][:40]/fac[3],dm3[1]-d3[2][:40]*dm0[1][10]/d0[1][10],dm3[2]/25.,fmt=':',color='r',linewidth=3)
	plt.errorbar(d0[0][:40]/fac[4],dm4[1]-d4[2][:40]*dm0[1][10]/d0[1][10],dm4[2]/25.,fmt=':',color='b',linewidth=3)
# 	plt.plot(d0[0][:40],dm0[1]-d0[2][:40],'k-',linewidth=3)
# 	plt.plot(d0[0][:40],dm1[1]-d1[2][:40],'k--',linewidth=3)
# 	plt.plot(d0[0][:40],dm2[1]-d2[2][:40],'k:',linewidth=3)
# 	plt.plot(d0[0][:40],dm3[1]-d3[2][:40],':',color='r',linewidth=3)
# 	plt.plot(d0[0][:40],dm4[1]-d4[2][:40],':',color='b',linewidth=3)
	plt.ylim(-.001,.001)
	pp.savefig()
	pp.close()
	return True


def xisigmumockcompthplot():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xisigmumockcompthnohimu.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$s^2\xi$',size=16)
	d0 = load('xizconv0MICE_matterpowermumax0.20.406.010.00.029sp1.0.dat').transpose()
	d1 = load('xizconv0MICE_matterpowermumin0.2mumax0.40.406.010.00.029sp1.0.dat').transpose()
	d2 = load('xizconv0MICE_matterpowermumin0.4mumax0.60.406.010.00.029sp1.0.dat').transpose()
	d3 = load('xizconv0MICE_matterpowermumin0.6mumax0.80.406.010.00.029sp1.0.dat').transpose()
	d4 = load('xizconv0MICE_matterpowermumin0.80.406.010.00.029sp1.0.dat').transpose()
	dm0 = load('xiave0mumax0.2SM.dat').transpose()
	dm1 = load('xiave0mumin0.2mumax0.4SM.dat').transpose()
	dm2 = load('xiave0mumin0.4mumax0.6SM.dat').transpose()
	dm3 = load('xiave0mumin0.6mumax0.8SM.dat').transpose()
	dm4 = load('xiave0mumin0.8SM.dat').transpose()
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm1[1],dm0[0][:40]**2.*dm1[2]/25.,fmt='rd',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm2[1],dm0[0][:40]**2.*dm2[2]/25.,fmt='b^',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0]**2.*dm3[1],dm0[0][:40]**2.*dm3[2]/25.,fmt='gs',linewidth=3)
	#plt.errorbar(dm0[0][:40],dm0[0]**2.*dm4[1],dm0[0][:40]**2.*dm4[2]/25.,fmt='y<',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm0[1],dm0[0][:40]**2.*dm0[2]/25.,fmt='ko',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d1[1]*1.4,'r-',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d2[1]*1.4,'b-',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d3[1]*1.4,'g-',linewidth=3)
 	#plt.plot(d0[0],d0[0]**2.*d4[1]*1.4,'y-',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d0[1]*1.4,'k-',linewidth=3)
 	plt.text(40,5,r'$\mu < 0.2$',color='k')
 	plt.text(40,3,r'$0.2 < \mu < 0.4$',color='r')
 	plt.text(40,1,r'$0.4 < \mu < 0.6$',color='b')
 	plt.text(40,-1,r'$0.6 < \mu < 0.8$',color='g')
	plt.ylim(-5,25)
	plt.xlim(30,200)
	pp.savefig()
	pp.close()
	return True

def xirpsigmumockcompthplot():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages(diro+'xirpsigmumockcompthnohimusig55.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s_{\perp}$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$s^2\xi_{\perp}$',size=16)
# 	d0 = load(diro+'xizconvcrpMICE_matterpowermumax0.20.406.010.0combzsiglsp1.0.dat').transpose()
# 	d1 = load(diro+'xizconvcrpMICE_matterpowermumin0.2mumax0.40.406.010.0combzsiglsp1.0.dat').transpose()
# 	d2 = load(diro+'xizconvcrpMICE_matterpowermumin0.4mumax0.60.406.010.0combzsiglsp1.0.dat').transpose()
# 	d3 = load(diro+'xizconvcrpMICE_matterpowermumin0.6mumax0.80.406.010.0combzsiglsp1.0.dat').transpose()
# 	d4 = load(diro+'xizconvcrpMICE_matterpowermumin0.80.406.010.0combzsiglsp1.0.dat').transpose()
	d0 = load(diro+'xizconvcrpMICE_matterpowermumax0.20.405.05.0combzsigl0.029sp1.0.dat').transpose()
	d1 = load(diro+'xizconvcrpMICE_matterpowermumin0.2mumax0.40.405.05.0combzsigl0.029sp1.0.dat').transpose()
	d2 = load(diro+'xizconvcrpMICE_matterpowermumin0.4mumax0.60.405.05.0combzsigl0.029sp1.0.dat').transpose()
	d3 = load(diro+'xizconvcrpMICE_matterpowermumin0.6mumax0.80.405.05.0combzsigl0.029sp1.0.dat').transpose()
	d4 = load(diro+'xizconvcrpMICE_matterpowermumin0.80.405.05.0combzsigl0.029sp1.0.dat').transpose()
	dm0 = load(diro+'xiaverpsqmumax0.2.dat').transpose()
	dm1 = load(diro+'xiaverpsqmumin0.2mumax0.4.dat').transpose()
	dm2 = load(diro+'xiaverpsqmumin0.4mumax0.6.dat').transpose()
	dm3 = load(diro+'xiaverpsqmumin0.6mumax0.8.dat').transpose()
	dm4 = load(diro+'xiaverpsqmumin0.8.dat').transpose()
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm1[1]+1.,dm0[0][:40]**2.*dm1[2]/sqrt(504.),fmt='rd',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm2[1],dm0[0][:40]**2.*dm2[2]/sqrt(504.),fmt='b^',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0]**2.*dm3[1],dm0[0][:40]**2.*dm3[2]/sqrt(504.),fmt='gs',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0]**2.*dm4[1],dm0[0][:40]**2.*dm4[2]/sqrt(504.),fmt='y<',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm0[1]+2.,dm0[0][:40]**2.*dm0[2]/sqrt(504.),fmt='ko',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d1[1]*1.48+1.,'r-',linewidth=2)
 	plt.plot(d0[0],d0[0]**2.*d2[1]*1.48,'b-',linewidth=2)
 	plt.plot(d0[0],d0[0]**2.*d3[1]*1.48,'g-',linewidth=2)
 	plt.plot(d0[0],d0[0]**2.*d4[1]*1.48,'y-',linewidth=2)
 	plt.plot(d0[0],d0[0]**2.*d0[1]*1.48+2.,'k-',linewidth=2)
 	plt.text(40,2.5,r'$\mu < 0.2$',color='k',size=16)
 	plt.text(40,1,r'$0.2 < \mu < 0.4$',color='r',size=16)
 	plt.text(40,-.5,r'$0.4 < \mu < 0.6$',color='b',size=16)
 	plt.text(40,-2,r'$0.6 < \mu < 0.8$',color='g',size=16)
 	plt.text(40,-3.5,r'$0.8 < \mu$',color='y',size=16)
	plt.ylim(-5,17)
	plt.xlim(30,200)
	pp.savefig()
	pp.close()
	return True


def xisigmulampcompthplot():
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	pp = PdfPages('xisigmumocklampcompth.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$s^2\xi$',size=16)
	d0 = load('xizconvcMICE_matterpowermumax0.20.4010.010.0combzsiglsp1.0.dat').transpose()
	d1 = load('xizconvcMICE_matterpowermumin0.2mumax0.40.4010.010.0combzsiglsp1.0.dat').transpose()
	d2 = load('xizconvcMICE_matterpowermumin0.4mumax0.60.4010.010.0combzsiglsp1.0.dat').transpose()
	d3 = load('xizconvcMICE_matterpowermumin0.6mumax0.80.4010.010.0combzsiglsp1.0.dat').transpose()
	d4 = load('xizconv0MICE_matterpowermumin0.80.406.010.00.029sp1.0.dat').transpose()
	dm0 = load('xiSMlampavemu00.25.dat').transpose()
	dm1 = load('xiSMlampavemu0.20.45.dat').transpose()
	dm2 = load('xiSMlampavemu0.40.65.dat').transpose()
	dm3 = load('xiSMlampavemu0.60.85.dat').transpose()
	#dm4 = load('xiave0mumin0.8SM.dat').transpose()
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm1[1],dm0[0][:40]**2.*dm1[2]/25.,fmt='rd',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm2[1],dm0[0][:40]**2.*dm2[2]/25.,fmt='b^',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0]**2.*dm3[1],dm0[0][:40]**2.*dm3[2]/25.,fmt='gs',linewidth=3)
	#plt.errorbar(dm0[0][:40],dm0[0]**2.*dm4[1],dm0[0][:40]**2.*dm4[2]/25.,fmt='y<',linewidth=3)
	plt.errorbar(dm0[0][:40],dm0[0][:40]**2.*dm0[1],dm0[0][:40]**2.*dm0[2]/25.,fmt='ko',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d1[1]*1.85,'r-',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d2[1]*1.85,'b-',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d3[1]*1.85,'g-',linewidth=3)
 	#plt.plot(d0[0],d0[0]**2.*d4[1]*1.4,'y-',linewidth=3)
 	plt.plot(d0[0],d0[0]**2.*d0[1]*1.85,'k-',linewidth=3)
 	plt.text(40,5,r'$\mu < 0.2$',color='k')
 	plt.text(40,3,r'$0.2 < \mu < 0.4$',color='r')
 	plt.text(40,1,r'$0.4 < \mu < 0.6$',color='b')
 	plt.text(40,-1,r'$0.6 < \mu < 0.8$',color='g')
	plt.ylim(-5,25)
	plt.xlim(30,200)
	pp.savefig()
	pp.close()
	return True


def xisigmumockcompthrsdplot(mumin,mumax,b=1.5,sfog=0,alph=1.):
	import matplotlib.pyplot as plt
	from matplotlib import rc
	from matplotlib.backends.backend_pdf import PdfPages
	#pp = PdfPages('xisigmumockcompth.pdf')
	plt.clf()
	plt.minorticks_on()
	plt.xlabel(r'$s$ ($h^{-1}$ Mpc)',size=16)
	plt.ylabel(r'$s^2\xi$',size=16)
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1:
		muw += 'mumax'+str(mumax)
	
	dt = load('xizconv0MICE_matterpower'+muw+'0.4'+str(sfog)+'6.010.00.029sp5.0.dat').transpose()
	#dnr = load('xizconv0MICE_matterpower'+muw+'0.4'+str(sfog)+'6.010.00.029norsdsp5.0.dat').transpose()
	dm = load('xiave0'+muw+'SM.dat').transpose()
	plt.errorbar(dm[0][:40],dm[0][:40]**2.*dm[1],dm[0][:40]**2.*dm[2]/25.,fmt='ko',linewidth=3)
 	plt.plot(dt[0]*alph,dt[0]**2.*dt[1]*b,'k-',linewidth=3)
 	#plt.plot(dnr[0],dnr[0]**2.*dnr[1]*b,'k--',linewidth=3)
	plt.show()
	#plt.ylim(-.001,.001)
#	pp.savefig()
#	pp.close()
	return True

def rbaomu():
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.axes as ax
	import numpy as np
	from Cosmo import *
	from numpy import loadtxt as load

	pp = PdfPages('rbaomu.pdf')
	d = load('/Users/ashleyross/DESY1/xizconvrbao504MICE_matterpower0.406.010.0combzsiglsp1.0.dat').transpose()
	plt.plot(d[0],d[1],'k-',linewidth=5)
	r0 = d[1][0]
	rpp = r0/(1.-d[0]**2.)**.5
	plt.plot(d[0],rpp,'k:',linewidth=2)
	plt.xlabel(r'$\mu$',size=16)
	plt.ylabel(r'$r_{\rm BAO}$ ($h^{-1}$Mpc)',size=16)
	plt.xlim(0,0.8)
	plt.ylim(100,200)
	pp.savefig()
	pp.close()
	return True

def BAOerrplot(wo='test',BOSS=True,MGS=False,wz=False,sdss=False,df6=False,des=True,eboss=False,desi=False):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.axes as ax
	import numpy as np
	from Cosmo import *
	from numpy import loadtxt as load

	pp = PdfPages('BAOerr'+wo+'.pdf')
	plt.clf()
	fig=plt.figure()
	ax0=fig.add_subplot(1,1,1)
	ax.Axes.axis(ax0,'auto')
	plt.tick_params(axis='both', which='major', labelsize=15)
	#x = pe[0]
	pe = np.loadtxt('/Users/ashleyross/DR7VAC/DVordPlanck.txt').transpose()
	x = pe[0]
	y=np.ones(len(x))
	yu = np.zeros(len(x))
	yd = np.zeros(len(x))

	x=np.arange(0.02,2.,0.01)
	y=np.ones(len(x))
	if df6:
		plt.text(0.11,.95,'6dFGS',fontsize=18)
		x6 = [0.1]
		e6 = [0.045]
		plt.plot(x6,e6,'s',markersize=7,color='green')
	if BOSS:
		plt.text(0.1,.021,'BOSS measured',fontsize=18)
		xl = [0.32,0.57]
		el = [0.020,0.01]
		plt.plot(xl,el,'-o',markersize=7,color='k')
	#plt.text(0.26,.9,'LOWZ',fontsize=18)
	#plt.text(0.51,.92,'BOSS',fontsize=18)
	#plt.text(0.5,.9,'CMASS',fontsize=18)
	if MGS:
		plt.text(0.16,1.05,'SDSS MGS',fontsize=18,color='r')
		xl = [0.15]
		el = [0.037]
		plt.plot(xl,el,'D',markeredgecolor='r',markersize=7,color='r')
	if wz:
		plt.text(0.6,.047,'WiggleZ measured',fontsize=18,color='.5')
		xwl = [.44,.6,.73] 
		ewl = [0.048,0.045,0.034]
		plt.plot(xwl,ewl,'-s',markersize=6,color='.5')

	if sdss:		
		xsl = [0.275,0.35]
		ysl = [0.972,0.969]
		esl = [0.027,0.018]
		plt.plot(xsl,esl,'s',markersize=6,elinewidth=1.,mfc='w',color='k')
	#plt.text(0.95,1.07,r'$\bar{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_avg)+'\pm'+'{:.3f}'.format(delta_alpha_std)+'$')
	#plt.text(0.95,1.04,r'$\tilde{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_med)+'^{'+'{:+.3f}'.format(delta_alpha_upp)+'}_{'+'{:+.3f}'.format(delta_alpha_low)+'}$')
	if eboss:
		plt.text(1.,.022,'eBOSS projected',fontsize=18,color='r')
		xl = [0.75,0.8,1.5]
		yl = [.99,.99,.99]
		el = [0.013,0.02,0.017]
		plt.plot(xl,el,'-D',markeredgecolor='r',markersize=7,color='r')
	if des:
		plt.text(.85,.0155,'DES projected',fontsize=18,color='b')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		#xld = [0.875,0.925]
		#yld = [1.0,1.0]
		#eld = [0.019,0.019]
		#yu = [1.019,1.019]
		#yd = [0.981,0.981]
		xl = [0.9]
		el = [0.019]
		#plt.errorbar(xld,yld,eld,fmt='-s',markeredgecolor='b',markersize=7,elinewidth=1.75,color='b')
		#plt.fill_between(xld,yd,yu,color='b')
		plt.plot(xl,el,'-^',markeredgecolor='b',markersize=7,color='b')
	if desi:
		plt.text(.9,.003,'DESI projected',fontsize=18,color='purple')		
		d = load('/Users/ashleyross/BAOpredictionscode/BAO_BigBOSS.dat').transpose()
		xl = d[0]
		yl = np.ones((len(d[0])))
		el = d[-1]
		plt.plot(xl,el,'-*',markeredgecolor='purple',markersize=7,color='purple')
	plt.xlim ( 0.0, 2. )
	plt.ylim ( 0.001, 0.055 )
	plt.xlabel ('Redshift', fontsize=18)
	#plt.ylabel (r'$(D_{\rm V}/r_{\rm d})/(D_{\rm V}/r_{\rm d})_{\rm Planck}$', fontsize=18)
	plt.ylabel ('Distance Uncertainty', fontsize=18)
	pp.savefig()
	pp.close()
	return True

if __name__ == '__main__':
	xisigmuplot6pan()