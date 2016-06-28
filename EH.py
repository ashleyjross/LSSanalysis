from Cosmo import distance as flatdistance
from math import *
from scipy.integrate import romberg

class simulate:
	def __init__(self,omega=0.3,lamda=0.7,h=0.7,h0=1.,ombhh=0.0224,CMBtemp=2.725,sig8=.8,c=2997,hmode=0,sbao=8.,nindex=.95):
		#self.cosmogend = gendistance(omega,lamda,h)
		ombaryon = ombhh/h/h
		self.Dlin = flatdistance(omega,lamda,h).D
		self.om_m0 = float(omega)
		self.om_v0 = float(lamda)
		self.hub = float(h)
		self.speedc = float(c) # speed of light in km/centi-second
		self.h0 = float(h0) # dimensionless hubble constant 
		self.codist = flatdistance(omega=self.om_m0,lamda=self.om_v0,h=self.h0).dc
		self.HubPar = flatdistance(omega=self.om_m0,lamda=self.om_v0,h=self.h0).Hz
		self.step = .001
		self.sbao = sbao
	# This is the baryon fraction in total matter
		self.omb = 1.*ombaryon/omega
		#print(self.omb)
	# This is the (cold) dark matter fraction in total matter
		self.omc = 1.-self.omb
	# Some frequently used quantitites
		self.omhh = 1.*omega*h*h
		self.ombhh = 1.*ombaryon*h*h
	# True CMB temperature relative to 2.7K
		self.qCMB = CMBtemp/2.7
	# This is the particle horizon scale at the equality epoch per Mpc
		self.keq = 0.0746*self.omhh*(self.qCMB**-2.)
	# This is the equality epoch
		zeq = 25000*self.omhh*(self.qCMB**-4.)
	# This is the drag epoch
		bee1 = 0.313*(self.omhh**-0.419)*(1.+(0.607*(self.omhh**0.674)))
		bee2 = 0.238*(self.omhh**0.223)
		zd = 1291.*(self.omhh**0.251)/(1.+(0.659*(self.omhh**0.828)))
		zd *= 1.+(bee1*(self.ombhh**bee2))
	# This is the sound horizon at the drag epoch (Mpc)
		Req = self.R(zeq)
		Rd = self.R(zd)        
		self.s = (2./(3.*self.keq))*sqrt(6./Req)*log((sqrt(1.+Rd)+sqrt(Rd+Req))/(1.+sqrt(Req)))
		self.sapprox = 44.5*log(9.83/self.omhh)/sqrt(1.+(10.*self.ombhh**0.75))
	# This is the silk damping scale per Mpc
		self.ksilk =1.6*(self.ombhh**0.52)*(self.omhh**0.73)*(1.+((10.4*self.omhh)**-0.95))
	# This is the CDM supression factor
		a1 = ((46.9*self.omhh)**0.670)*(1.+((32.1*self.omhh)**-0.532))
		a2 = ((12.0*self.omhh)**0.424)*(1.+((45.*self.omhh)**-0.582))
		self.ac = (a1**(-self.omb))*(a2**(-1.*(self.omb)**3.))
	# This is the CDM log shift
		b1 = 0.944/(1.+((458*self.omhh)**-0.708))
		b2 = (0.395*self.omhh)**(-0.0266)
		self.bc = 1./(1.+(b1*((self.omc**b2)-1.)))
	# This is the baryon supression factor
		y = (1.+zeq)/(1.+zd)
		self.ab = 2.07*self.keq*self.s*((1.+Rd)**-0.75)*HyperG(y)
	# This is the baryon envelope shift
		self.bb = 0.5+self.omb+((3.-(2.*self.omb))*sqrt(((17.2*self.omhh)**2.)+1.))
	# This is the node shift parameter
		self.bnode = 8.41*(self.omhh)**0.435
	# This is the suppression of the oft-used gamma term
		self.agam = 1.-(0.328*log(431.*self.omhh)*self.omb)+(0.38*log(22.3*self.omhh)*(self.omb**2.))
	# And that should be everything you need nominally for BAOs
	
	# Harrison-Peebles-Zeldovich case 'cos it's close enough these days
	# Now find normalization based on sigma_8 in a particular cosmology
		self.nindex = nindex 
		self.P0 = 1.
		self.P0smooth = 1.
	# units of h-1Mpc
		correc = (sig8/self.siglog())**2.
		#correc = sig8
		correcsmooth = (sig8/self.siglogsmooth())**2.
		keff = self.h0/8.
		#print self.siglog()
		#print self.P(keff,0)*keff**3./2./pi/pi
		self.P0 = correc
		#print correc
		self.P0smooth = correcsmooth
		
		#print self.P(keff,0)*keff**3./2./pi/pi
		if hmode == 1:
			from Halonoscim import HODint
			self.hdint = HODint(Mmin=10.,no = .0005,M1 = 14.3, al =1.5,sigMd=0,A=.32222 ,p=.3,q=.75,om_v=self.om_v0,om_m=self.om_m0,z=.5,sig8=sig8,h=self.hub)



########################
# P(k), sig_8 stuff....
########################
	def Wm(self,k):
		#the window function for a spherical tophat
		kr = k*self.rscale
		return 3./kr**3.*(sin(kr)-kr*cos(kr))
		#return 3.*(sin(kr)-kr*cos(kr)) #since the k^3 comes out in the next step
	
	def sigintfunc(self,lrk):
		#linear power spectrum convolved with the window function
		rk = 10**lrk
		wsq = self.Wm(rk)**2.
		pk = self.P(rk,0)
		return wsq*pk*rk**3.

	def siglog(self,r=8.,z=0.,acc=1e-10,lims=[-10,10]):
		"""rms mass fluctuations (integrated in logspace) for sig8 etc.
		
		NOTES: formalism is all real-space distance in Mpc/h, all k in
		h/Mpc
		
		v1.0 Adam D. Myers Jan 2007
		"""
		
		self.zee = z
		self.rscale = r
		#return sqrt((log(10.)/2./pi/pi)*(rom(-10,-5,self.intefunclog)+rom(-5,0,self.intefunclog)+rom(0,5,self.intefunclog)+rom(5,10,self.intefunclog)))
		return sqrt((log(10.)/2./pi/pi)*romberg(self.intefunclog,-10,10,tol=acc))
		#return sqrt((log(10.)/2./pi/pi)*rom(-10,10,self.intefunclog,eps=1e-10)) #changed to match Halonoscim.py

	def intefunclog(self,logk):
		"""Function to integrate to get rms mass fluctuations (logspace)
	
		NOTES: Because the k's here are all expressed in h/Mpc the
		value of P0 that arises using this integration is in
		(Mpc**3)/(h**3.) and therefore values of P(k) derived from this
		integration are in (Mpc**3)/(h**3.), i.e. h^-3Mpc^3 and need
		divided by h^3 to get the 'absolute' P(k) in Mpc^3.
	
		v1.0 Adam D. Myers Jan 2007
		"""
	
		k = 10.**logk
		kr = k*self.rscale
	
		return k*k*k*self.P(k,self.zee)*(w(kr)**2.)

	def siglogsmooth(self,r=8.,z=0.,acc=1e-10,lims=[-10,10]):
		"""rms mass fluctuations (integrated in logspace) for sig8 etc.
	
		NOTES: formalism is all real-space distance in Mpc/h, all k in
		h/Mpc
		
		v1.0 Adam D. Myers Jan 2007
		"""
	
		self.zee = z
		self.rscale = r
		
		#return sqrt((log(10.)/2./pi/pi)*rom(-10,10,self.intefunclogsmooth,eps=1.e-10))
		return sqrt((log(10.)/2./pi/pi)*romberg(self.intefunclogsmooth,-10,10,tol=acc))
	
	def intefunclogsmooth(self,logk):
		"""Function to integrate to get rms mass fluctuations (logspace)
	
		NOTES: Because the k's here are all expressed in h/Mpc the
		value of P0 that arises using this integration is in
		(Mpc**3)/(h**3.) and therefore values of P(k) derived from this
		integration are in (Mpc**3)/(h**3.), i.e. h^-3Mpc^3 and need
		divided by h^3 to get the 'absolute' P(k) in Mpc^3.
	
		v1.0 Adam D. Myers Jan 2007
		"""
	
		k = 10.**logk
		kr = k*self.rscale
	
		return k*k*k*self.Psmooth(k,self.zee)*(w(kr)**2.)
	
	
	
	def P(self,k,z):
		"""Full baryonic Linear power spectrum
	
		NOTES: Pass k in h/Mpc
	
		NOTES: Technically the k**n would be a problem depending on
		whether k has a factor of h in it or not, which, as we express
		it, it does. It gets absorbed, though....i.e, it all comes out
		in the normalization of PO for sigma8...so I wouldn't worry
		about it....but hey, try it if you like.
	
		v1.0 Adam D. Myers Jan 2007
		"""
	
		return self.P0*(self.T(k)**2)*(k**self.nindex)*self.Dlin(z)*self.Dlin(z)
	
	def Psmooth(self,k,z):
		"""Baryonic Linear power spectrum SHAPE 
	
		NOTES: Pass k in h/Mpc
	
		NOTES: Just to be clear for any student readers, there are
		three P(k) that are typically useful. (i) The full, baryonic
		P(k)..call function P(k). The P(k) in the absence of
		baryons..call P(k) with ombaryon=1e-16....now, removing baryons
		removes baryonic features but also augments power, because
		baryons suppress P(k) relative to CDM. Finally, there is the
		SHAPE of the full P(k), removing the oscillations that the
		baryons cause without removing the power suppression that the
		baryons cause. This is what I'm calling Psmooth.
		
		v1.0 Adam D. Myers Jan 2007
		"""
	
		return self.P0smooth*(self.T0(k)**2)*(k**self.nindex)*self.Dlin(z)*self.Dlin(z)

		
###Everything from here down codifies the Eisenstein & Hu Transfer
###functions. At some point I'll give this it's own module.
###ADM Jan 2007        

	def T(self,k):
		"""Full transfer function from Eisenstein & Hu 1998
	
		NOTES: k is always expressed in h/Mpc in the transfer functions
		as written in this code and should be passed in units of h/Mpc
	
		v1.0 Adam D. Myers Jan 2007"""
	
		k *= self.hub # correct from h/Mpc to Mpc
	   
		self.q = k/(13.41*self.keq)
	
		return (self.omb*self.Tb(k))+(self.omc*self.Tc(k))
	
	def Tb(self,k):
		"""Baryon component of full transfer function from E&H98
	
		v1.0 Adam D. Myers Jan 2007"""
	
		ks = self.s*k
		stilda = self.s/((1.+(self.bnode/ks)**3.)**(1/3.))
		
		return ((self.Ttilda(k,1.,1.)/(1.+(ks/5.2)**2.))+(e**(-1.*((k/self.ksilk)**1.4))*self.ab/(1.+((self.bb/ks)**3.))))*(sin(k*stilda)/(k*stilda))
	
	def Tc(self,k):
		"""CDM component of full transfer function from E&H98
	
		v1.0 Adam D. Myers Jan 2007"""
	
		f = 1./(1.+((k*self.s/5.4)**4.))
	
		return (f*self.Ttilda(k,1.,self.bc))+((1.-f)*self.Ttilda(k,self.ac,self.bc))
	
	def Ttilda(self,k,alph,bet):
		"""pressureless form of transfer function from E&H98
	
		v1.0 Adam D. Myers Jan 2007"""
		#print k,self.q
		fac = log(e+(1.8*bet*self.q))
		C = (14.2/alph)+(386./(1.+(69.9*(self.q**1.08))))
	
		return fac/(fac+(C*self.q*self.q))
	
	def T0(self,k):
		"""Zero-baryon transfer function shape from E&H98
	
		v1.0 Adam D. Myers Jan 2007"""
		
		k *= self.hub # correct from h/Mpc to Mpc
		gamefftimeshub = self.omhh*(self.agam+((1.-self.agam)/(1.+(0.43*k*self.sapprox)**4.)))
		q = k*(self.qCMB**2.)/gamefftimeshub
	
		C = 14.2+(731./(1.+(62.5*q)))
		L = log((2.*e)+1.8*q)
		return L/(L+(C*q*q))
	
	
	def R(self,z):
		"""Ratio between baryon and photon momentum from E&H98
	
		v1.0 Adam D. Myers Jan 2007"""
	
		return 31500.*self.ombhh*(self.qCMB**-4.)/z

def HyperG(y):
    """Hypergeometric function"""

    return y*((-6.*sqrt(1.+y))+((2.+(3.*y))*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.))))

def w(kr):
    """FT of a spherical top hat"""
    return 3.*(sin(kr)-(kr*cos(kr)))/kr/kr/kr
