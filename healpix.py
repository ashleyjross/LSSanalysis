from math import *
from legendre import legendre
twopi = 2.*pi
piover2 = .5*pi
ns_max = 8192
ee = exp(1.)


def pix2ang_ring(nside,ipix):

#c=======================================================================
#c     gives theta and phi corresponding to pixel ipix (RING) 
#c     for a parameter nside
#c=======================================================================


#		pi = 3.1415926535897932384626434
  
	ns_max=8192
	npix = 12*nside*nside


	ipix1 = ipix + 1
	nl2 = 2*nside
	nl4 = 4*nside
	ncap = 2*nside*(nside-1)#// ! points in each polar cap, =0 for nside =1
	fact1 = 1.5*nside
	fact2 = 3.0*nside*nside
	
	if ipix1 <= ncap :# {  //! North Polar cap -------------
		hip   = ipix1/2.
		fihip = floor(hip)
		iring = int(floor( sqrt( hip - sqrt(fihip) ) ) + 1)#;// ! counted from North pole
		iphi  = ipix1 - 2*iring*(iring - 1)
		
		theta = acos( 1. - iring*iring / fact2 )
		phi   = (1.*iphi - 0.5) * pi/(2.*iring)

	else:
		if ipix1 <= nl2*(5*nside+1):# {//then ! Equatorial region ------

			ip    = ipix1 - ncap - 1
			iring = int(floor( ip / nl4 ) + nside)#;// ! counted from North pole
			iphi  = int(fmod(ip,nl4) + 1)
			
			fodd  = 0.5 * (1 + fmod(float(iring+nside),2))#//  ! 1 if iring+nside is odd, 1/2 otherwise
			theta = acos( (nl2 - iring) / fact1 )
			phi   = (1.*iphi - fodd) * pi /(2.*nside)
		else:# {//! South Polar cap -----------------------------------

			ip    = npix - ipix1 + 1
			hip   = ip/2.
#/* bug corrige floor instead of 1.* */
			fihip = floor(hip)
			iring = int(floor( sqrt( hip - sqrt(fihip) ) ) + 1)#;//     ! counted from South pole
			iphi  = int((4.*iring + 1 - (ip - 2.*iring*(iring-1))))
			
			theta = acos( -1. + iring*iring / fact2 )
			phi   = (1.*iphi - 0.5) * pi/(2.*iring)
	
	return theta,phi

def ang2pix_ring(nside,theta,phi):
#    c=======================================================================
#   c     gives the pixel number ipix (RING) 
#    c     corresponding to angles theta and phi
#    c=======================================================================
  
	z0=2.0/3.0
	ns_max=8192
  
	z = cos(theta)
	za = fabs(z)
	if phi >= twopi:
		phi = phi - twopi
	if phi < 0.:
		phi = phi + twopi
	tt = phi / piover2#;//  ! in [0,4)
  
	nl2 = 2*nside
	nl4 = 4*nside
	ncap  = nl2*(nside-1)#// ! number of pixels in the north polar cap
	npix  = 12*nside*nside
  
	if za <= z0:# {
    
		jp = int(floor(nside*(0.5 + tt - z*0.75)))#; /*index of ascending edge line*/
		jm = int(floor(nside*(0.5 + tt + z*0.75)))#; /*index of descending edge line*/
    
		ir = nside + 1 + jp - jm#;// ! in {1,2n+1} (ring number counted from z=2/3)
		kshift = 0
		if fmod(ir,2)==0.:
			kshift = 1#;// ! kshift=1 if ir even, 0 otherwise
		ip = int(floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1)#;// ! in {1,4n}
		if ip>nl4:
			ip = ip - nl4
    
		ipix1 = ncap + nl4*(ir-1) + ip

	else:
    
		tp = tt - floor(tt)#;//      !MOD(tt,1.d0)
		tmp = sqrt( 3.*(1. - za) )
		
		jp = int(floor( nside * tp * tmp ))#;// ! increasing edge line index
		jm = int(floor( nside * (1. - tp) * tmp ))#;// ! decreasing edge line index
		
		ir = jp + jm + 1#;//        ! ring number counted from the closest pole
		ip = int(floor( tt * ir ) + 1)#;// ! in {1,4*ir}
		if ip>4*ir:
			ip = ip - 4*ir
		
		ipix1 = 2*ir*(ir-1) + ip
		if z<=0.:
			ipix1 = npix - 2*ir*(ir+1) + ip

	return ipix1 - 1

def SDSSDR8runmap(Nside=256):
	f = open('radecrun.csv')
	f.readline()
	h = healpix()
	n = Nside*Nside*12
	a = []
	for i in range(0,n):
		a.append([-999])
	for line in f:
		ln = line.split(',')
		ra,dec = float(ln[0]),float(ln[1])
		th,phi = radec2thphi(ra,dec)
		pix = int(h.ang2pix_nest(Nside,th,phi))
		run = float(ln[2])
		if a[pix][0] == -999:
			a[pix][0] = run
		else:
			if a[pix][0] != run:
				a[pix].append(run)
	f = open('radecrun19.119.5.csv')
	f.readline()
	for line in f:
		ln = line.split(',')
		ra,dec = float(ln[0]),float(ln[1])
		th,phi = radec2thphi(ra,dec)
		pix = int(h.ang2pix_nest(Nside,th,phi))
		run = float(ln[2])
		if a[pix][0] == -999:
			a[pix][0] = run
		else:
			if a[pix][0] != run:
				a[pix].append(run)

	fo = open('pixrun'+str(Nside)+'.dat','w')
	for i in range(0,len(a)):
		for j in range(0,len(a[i])):
			fo.write(str(a[i][j])+' ')
		fo.write('\n')
	
	fo.close()

def SDSSDR8datemap(Nside=256):
	f = open('radecmjdrun.csv')
	f.readline()
	h = healpix()
	n = Nside*Nside*12
	a = []
	for i in range(0,n):
		a.append([-999])
	for line in f:
		ln = line.split(',')
		ra,dec = float(ln[0]),float(ln[1])
		th,phi = radec2thphi(ra,dec)
		pix = int(h.ang2pix_nest(Nside,th,phi))
		date = float(ln[2])
		if a[pix][0] == -999:
			a[pix][0] = date
		else:
			if a[pix][0] != date:
				a[pix].append(date)

	fo = open('pixdate'+str(Nside)+'.dat','w')
	for i in range(0,len(a)):
		for j in range(0,len(a[i])):
			fo.write(str(a[i][j])+' ')
		fo.write('\n')
	
	fo.close()

def DR7seemap160(Nside=512):
	from stripeDR7 import pix2etalama
	h = healpix()
	n = Nside*Nside*12
	nl = []
	seel = []
	for i in range(0,n):
		nl.append(0)
		seel.append(0)
	f = open('/Users/ashleyr/masks/DR7seepix160com.dat')
	ps = 0
	for line in f:
		see = float(line)
		if see > 0 and see < 10.:
			eta,lam = pix2etalama(160,ps) 
			th,phi = le2thetaphi(lam,eta)
			p = int(h.ang2pix_nest(Nside,th,phi))
			nl[p] += 1.
			seel[p] += see
		ps += 1	
	fo = open('DR7seemap_healp'+str(Nside)+'.dat','w')
	for i in range(0,len(nl)):
		if nl[i] == 0:
			fo.write('0\n')
		else:
			fo.write(str(seel[i]/nl[i])+'\n')
	fo.close()
	return True


def DR7seemapp(Nside=256):
	h = healpix()
	n = Nside*Nside*12
	nl = []
	seel = []
	for i in range(0,n):
		nl.append(0)
		seel.append(0)
	f = open('/Users/ashleyr/masks/DR7pixelsRS.dat')
	for line in f:
		ln = line.split()
		see = float(ln[-2])
		if see > 0 and see < 10.: 
			lam = float(ln[0])
			eta = float(ln[1])
			th,phi = le2thetaphi(lam,eta)
			p = int(h.ang2pix_nest(Nside,th,phi))
			nl[p] += 1.
			seel[p] += see
	fo = open('DR7seemap_healp'+str(Nside)+'.dat','w')
	for i in range(0,len(nl)):
		if nl[i] == 0:
			fo.write('-999\n')
		else:
			fo.write(str(seel[i]/nl[i])+'\n')
	fo.close()
	return True

def DR7skymap(Nside=256):
	h = healpix()
	n = Nside*Nside*12
	nl = []
	seel = []
	for i in range(0,n):
		nl.append(0)
		seel.append(0)
	f = open('/Users/ashleyr/QSO/DR7gal1919.5skyr.csv')
	f.readline()
	for line in f:
		ln = line.split(',')
		sky = float(ln[2])
		ra = float(ln[0])
		dec = float(ln[1])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(Nside,th,phi))
		nl[p] += 1.
		seel[p] += sky
	fo = open('DR7skymap_healp'+str(Nside)+'.dat','w')
	for i in range(0,len(nl)):
		if nl[i] == 0:
			fo.write('-999\n')
		else:
			fo.write(str(seel[i]/nl[i])+'\n')
	fo.close()
	return True


		
	

def mkoffmap(Nside=256):
	f = open('runoff_grri.dat').readlines()
	fp = open('pixrun'+str(Nside)+'.dat')
	fo = open('offmap'+str(Nside)+'.dat','w')
	for line in fp:
		runl = line.split()
		n = 0
		rioff = 0
		groff = 0
		for i in range(0,len(runl)):
			run = float(runl[i])
			if run > 0 and run < 7943:								
				for j in range(0,len(f)):
					lnf = f[j].split()
					runo = float(lnf[0])
					if runo == run:
						rioff += float(lnf[1])
						groff += float(lnf[2])
						n += 1.
						break
		if n > 0:
			fo.write(str(rioff/float(len(runl)))+' '+str(groff/float(len(runl)))+'\n')
		else:
			fo.write('-99 -99 \n')

	fo.close()
	return True


def hpixrd(nside=1024):
	h = healpix()
	fo = open('hpixrd'+str(nside)+'.dat','w')
	np = int(nside**2.*12)
	for i in range(0,np):
		th,phi = h.pix2ang_nest(nside,i)
		ra,dec = thphi2radec(th,phi)
		fo.write(str(ra)+' '+str(dec)+'\n')
	fo.close()
	return True

def stars2healpix(nside=1024):
	h = healpix()
	f = open('allstars17.519.9.dat')
	fp = open('allstarp.dat','w')
	pixl =[]
	np = int(nside**2.*12)
	for i in range(0,np):
		pixl.append(0)
	for line in f:
		ln = line.split()
		ra,dec = float(ln[0]),float(ln[1])
		th,phi = radec2thphi(ra,dec)
		px = int(h.ang2pix_nest(nside,th,phi))
		pixl[px] += 1.
	f.close()
	print 'file done'
	for i in range(0,len(pixl)):
		fp.write(str(pixl[i])+'\n')
	fp.close()
	c = 0
	fp = open('healpix_star17.519.9.dat','w')
	for i in range(0,len(pixl)):
		c += 1
		if c == nside:
			fp.write(str(pixl[i])+'\n')
			c = 0
		else:
			fp.write(str(pixl[i])+' ')
	fp.close()
	return True
		
	

def healpixm2radecm(nside,tol=.9):
	f = open('shirleyhpixNorth1024.dat').readlines()
	fo = open('shirleyhpixNorthm.dat','w')
	frd = open('shirleyhpixNorthrd.dat','w')
#	h = healpix()
	for i in range(0,len(f)):
		row = f[i].split()
		for j in range(0,len(row)):
			#pix = (j/nside)*(12*nside)+fmod(j,nside)+i*nside
			pix = 1024*i+j
			if float(row[j]) < tol:
				fo.write(str(pix)+'\n')
#			else:
#				th,phi = h.pix2ang_nest(nside, pix)
#				ra,dec = thphi2radec(th,phi)
				frd.write(str(ra)+' '+str(dec)+'\n')
	fo.close()

def pixl2odensNS(file,res=256,t=.74,seec=100,skyc=100,redc=100,stm=10000,md='n',ast=4.5):
	f = open(file+'hpix'+str(res)+'.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	aveS = 0
	ntS = 0
	aveN = 0
	ntN = 0
	ng = 0
	np = 0
	arc2 = 169956.96315526127
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		stln = stf.readline().split()
		st = float(stln[2])/p
		eta = float(stln[1])
		if md == 'sub':
			pm = st*pi*(ast)**2./arc2*p
			p = p-pm
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			ng += float(ln[2])
			np += float(ln[3])
			if eta > 70.:
				aveS += float(ln[2])
				ntS += p
			else:
				aveN += float(ln[2])
				ntN += p
	aveS = aveS/ntS
	aveN = aveN/ntN
	aver = ng/np
	print np,ng,aveN,aveS,aver
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	f = open(file+'hpix'+str(res)+'.dat')
	if md == 'sub':
		file += 'sub'
	foS = open(file+ro+so+sko+sto+'hSodens64.dat','w')
	foN = open(file+ro+so+sko+sto+'hNodens64.dat','w')

	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		stln = stf.readline().split()
		eta = float(stln[1])
		ln = line.split()
		p = float(ln[3])
		st = float(stln[2])/p
		if md == 'sub':
			pm = st*pi*(ast)**2./arc2*p
			p = p-pm
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			#od = float(ln[2])/ave - p
			if eta > 70.:
				od = float(ln[2])/p/aveS - 1.
				foS.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			else:
				od = float(ln[2])/p/aveN - 1.
				foN.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	foN.close()
	foS.close()
	print sump
	return True


def pixl2odens_simp(file,res=1024,t=.9,md='g'):
	f = open(file+'hpix'+str(res)+'.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		if p > t:
			ng += float(ln[2])
			np += p
			ave += float(ln[2])/p
			nt += 1.
	
	aver = ng/nt
	ave = ng/np
	print np,ng,ave,aver
	f = open(file+'hpix'+str(res)+'.dat')
	if res != '1024':
		fo = open(file+str(res)+'hodens64.dat','w')
	else:
		fo = open(file+'hodens64.dat','w')

	sump = 0
	sumn = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])

		if p > t:
			if md == 'g':
				od = float(ln[2])/p/ave - 1.
			if md == 'sys':
				od = float(ln[2])/aver-1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	fo.close()
	print sump
	return True

def pixl2abs_simp(file,res=1024,t=.9,md='g'):
	f = open(file+'hpix'+str(res)+'.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		if p > t:
			ng += float(ln[2])
			np += p
			ave += float(ln[2])/p
			nt += 1.
	
	aver = ng/nt
	ave = ng/np
	print np,ng,ave,aver
	f = open(file+'hpix'+str(res)+'.dat')
	if res != '1024':
		fo = open(file+str(res)+'hodens64.dat','w')
	else:
		fo = open(file+'hodens64.dat','w')

	sump = 0
	sumn = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])

		if p > t:
			if md == 'sys':
				od = float(ln[2])-aver
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	fo.close()
	print sump
	return True


def healRDfullsky(res):
	np = 12*res*res
	h = healpix()
	fo = open('healRDfullsky'+str(res)+'.dat','w')
	for i in range(0,np):
		th,phi = h.pix2ang_nest(res,i)
		ra,dec = thphi2radec(th,phi)
		fo.write(str(ra)+' '+str(dec)+'\n')
	fo.close()
	return True


def pixl2odens(file,res=256,t=.74,seec=100,skymin=0,skyc=100,redc=100,stm=10000,schc=100,bmin=-100,bmax=100,abmin=0,dir='pc2pt',md='n',ast=4.5,wts='n'):
	f = open(file+'hpix'+str(res)+'.dat')
	seef = open(dir+'/healsee256.dat')
	skyf = open(dir+'/healsky256.dat')
	redf = open(dir+'/healred256.dat')
	#stf = open(dir+'/star_infohpix256.dat')
	stf = open(dir+'/allstars17.519.9n_mHealhpix256.dat')
	schf = open(dir+'/healschug256.dat')
	bf = open(dir+'/lbhpixl256.dat')
	arc2 = 169956.96315526127
	ave = 0
	nt = 0
	ng = 0
	np = 0
	if wts == 'y':
		wtsf = open(dir+'/wl'+file.strip('/pc2pt')+'.dat').readlines()
	i = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		sch = float(schf.readline())
		b = float(bf.readline().split()[1])
		st = float(stf.readline().split()[2])/p
		if md == 'sub':
			pm = st*pi*(ast)**2./arc2*p
			p = p-pm
		if p > t:# and 
			if see < seec and sky > skymin and sky < skyc and red < redc and st < stm and sch < schc and b > bmin and b < bmax and abs(b) > abmin:
				n = float(ln[2])
				if wts == 'y':
					ng += n*float(wtsf[i])
				else:
					ng += n
				np += p
				#ave += ng/p
				nt += 1.
		i += 1
			#else:
			#	print see,sky,red,st,sch,b
	#ave = ave/nt
	ave = ng/np
	print np,ng,ave
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if skymin == 0:
		skmo =''
	else:
		skmo = 'skym'+str(skymin)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
	if schc == 100:
		scho = ''
	else:
		scho = 'sch'+str(schc)
	if bmin == -100:
		bo = ''
	else:
		bo = 'bm'+str(bmin)
	if bmax == 100:
		bx = ''
	else:
		bx = 'bx'+str(bmax)

	if abmin == 0:
		abo = ''
	else:
		abo = 'abm'+str(abmin)
	f = open(file+'hpix'+str(res)+'.dat')
	seef = open(dir+'/healsee256.dat')
	skyf = open(dir+'/healsky256.dat')
	redf = open(dir+'/healred256.dat')
	#stf = open(dir+'/star_infohpix256.dat')
	stf = open(dir+'/allstars17.519.9n_mHealhpix256.dat')
	schf = open(dir+'/healschug256.dat')
	#bf = open(dir+'/LRGcmassV4sgpixlb.dat')
	bf = open(dir+'/lbhpixl256.dat')
	if md == 'sub':
		file += 'sub'
	if wts == 'y':
		file += 'wts'
	fo = open(file+ro+so+skmo+sko+sto+scho+bo+bx+abo+'hodens64.dat','w')

	sump = 0
	sumn = 0
	i = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		sch = float(schf.readline())
		b = float(bf.readline().split()[1])
		if md == 'sub':
			pm = st*pi*(ast)**2./arc2*p
			p = p-pm

		if p > t and see < seec and sky > skymin and sky < skyc and red < redc and st < stm and sch < schc and b > bmin and b < bmax and abs(b) > abmin:
			#od = float(ln[2])/ave - p
			if wts == 'y':
				od = float(ln[2])*float(wtsf[i])/p/ave - 1.
			else:
				od = float(ln[2])/p/ave - 1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
		i += 1
	fo.close()
	print sump
	return True

def pixl2odensstar(file,res=256,t=.74,seec=100,skyc=100,redc=100,stm=100):
	f = open(file+'hpix'+str(res)+'.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/star_infohpix256.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			ng += float(ln[2])
			np += float(ln[3])
			ave += float(ln[2])/float(ln[3])
			nt += 1.
	
	ave = ave/nt
	aver = ng/np
	print np,ng,ave,aver
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if stm == 100:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	fo = open(file+ro+so+sko+sto+'hstarodens64.dat','w')
	f = open(file+'hpix'+str(res)+'.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/star_infohpix256.dat')
	avef = open(file+ro+so+sko+sto+'nvsysH.dat').readlines()
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		ln = line.split()
		p = float(ln[3])
		st = float(stf.readline().split()[2])/p
		stdiv = stm/10.
		starind = int(st/stdiv)
		if starind > 9:
			starind = 9
		
		
		mave = float(avef[starind].split()[26])
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			#od = float(ln[2])/ave - p
			od = float(ln[2])/p/(ave*mave) - 1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+'\n')
			sump += p
	fo.close()
	print sump,
	return True

def pixl2odens_ext(res=256,t=.2,seec=100,skyc=100,redc=100,stm=10000):
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			ng += red*p
			nt += 1.*p
	ave = ng/nt
	#print np,ng,ave,aver
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	fo = open('pc2pt/ext'+ro+so+sko+sto+'hodens64.dat','w')
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		ln = line.split()
		p = float(ln[3])
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			#od = float(ln[2])/ave - p
			od = red/ave - 1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	fo.close()
	print sump,
 	from pcsourcezw import pcsourcezw
 	a = pcsourcezw.createSourcesLE_zw('pc2pt/exthodens64.dat')
 	for i in range(0,20):
 		a = pcsourcezw.createSourcesLE_jackzw('pc2pt/exthodens64.dat',i,jackt=20)
	return True
	
def pixl2odens_sky(res=256,t=.2,seec=100,skymin=0,skyc=100,redc=100,stm=10000):
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm:
			ng += sky*p
			nt += 1.*p
	ave = ng/nt
	#print np,ng,ave,aver
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if skymin == 0:
		skmo = ''
	else:
		skmo = 'skym'+str(skymin)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	fo = open('pc2pt/sky'+ro+so+skmo+sko+sto+'hodens64.dat','w')
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		ln = line.split()
		p = float(ln[3])
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm:
			#od = float(ln[2])/ave - p
			od = sky/ave - 1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	fo.close()
	print sump,
 	from pcsourcezw import pcsourcezw
 	a = pcsourcezw.createSourcesLE_zw('pc2pt/skyhodens64.dat')
 	for i in range(0,20):
 		a = pcsourcezw.createSourcesLE_jackzw('pc2pt/skyhodens64.dat',i,jackt=20)
	return True

def pixl2odens_skyS(res=256,t=.2,seec=100,skymin=0,skyc=100,redc=100,stm=10000):
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		eta = float(ln[1])
		if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm and eta > 90.:
			ng += sky*p
			nt += 1.*p
	ave = ng/nt
	#print np,ng,ave,aver
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if skymin == 0:
		skmo = ''
	else:
		skmo = 'skym'+str(skymin)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	fo = open('pc2pt/sky'+ro+so+skmo+sko+sto+'hSodens64.dat','w')
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		ln = line.split()
		p = float(ln[3])
		st = float(stf.readline().split()[2])/p
		eta = float(ln[1])
		if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm and eta > 90:
			#od = float(ln[2])/ave - p
			od = sky/ave - 1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	fo.close()
	print sump,
 	from pcsourcezw import pcsourcezw
 	a = pcsourcezw.createSourcesLE_zw('pc2pt/skyhSodens64.dat')
 	for i in range(0,20):
 		a = pcsourcezw.createSourcesLE_jackzw('pc2pt/skyhSodens64.dat',i,jackt=20)
	return True


def pixl2odens_air(res=256,t=.2,seec=100,skymin=0,skyc=100,redc=100,stm=10000):
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	airf = open('pc2pt/healair256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		air = float(airf.readline())
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm:
			ng += air*p
			nt += 1.*p
	ave = ng/nt
	#print np,ng,ave,aver
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if skymin == 0:
		skmo = ''
	else:
		skmo = 'skym'+str(skymin)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	fo = open('pc2pt/air'+ro+so+skmo+sko+sto+'hodens64.dat','w')
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	airf = open('pc2pt/healair256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		air = float(airf.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		ln = line.split()
		p = float(ln[3])
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm:
			#od = float(ln[2])/ave - p
			od = air/ave - 1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	fo.close()
	print sump,
 	from pcsourcezw import pcsourcezw
 	a = pcsourcezw.createSourcesLE_zw('pc2pt/airhodens64.dat')
 	for i in range(0,20):
 		a = pcsourcezw.createSourcesLE_jackzw('pc2pt/airhodens64.dat',i,jackt=20)
	return True


def pixl2odens_sch(res=256,t=.2,seec=100,skymin=0,skyc=100,redc=100,stm=10000,schC='ri'):
	f = open('pc2pt/LRGcasBwredab1.6_0.60.65sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	schf = open('pc2pt/healsch'+schC+'256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		sch = float(schf.readline())
		if sch != -1:
			if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm and sch > 0:
				ng += sch*p
				nt += 1.*p
	ave = ng/nt
	print np,ng,ave
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if skymin == 0:
		skmo = ''
	else:
		skmo = 'skym'+str(skymin)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	fo = open('pc2pt/sch'+schC+ro+so+skmo+sko+sto+'hodens64.dat','w')
	f = open('pc2pt/LRGcasBwredab1.6_0.60.65sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	schf = open('pc2pt/healsch'+schC+'256.dat')
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		sch = float(schf.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		ln = line.split()
		p = float(ln[3])
		st = float(stf.readline().split()[2])/p
		if sch != -1:
			if p > t and see < seec and sky < skyc and sky > skymin and red < redc and st < stm and sch > 0:
				#od = float(ln[2])/ave - p
				od = sch/ave - 1.
				fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
				sump += p
	fo.close()
	print sump
	from pcsourcezw import pcsourcezw
	a = pcsourcezw.createSourcesLE_zw('pc2pt/sch'+schC+'hodens64.dat')
	for i in range(0,20):
		a = pcsourcezw.createSourcesLE_jackzw('pc2pt/sch'+schC+'hodens64.dat',i,jackt=20)
	return True

	
def pixl2odens_see(res=256,t=.74,seec=100,skyc=100,redc=100,stm=10000):
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	ave = 0
	nt = 0
	ng = 0
	np = 0
	for line in f:
		ln = line.split()
		p = float(ln[3])
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			ng += see*p
			nt += 1.*p
	ave = ng/nt
	#print np,ng,ave,aver
	if redc == 100:
		ro = ''
	else:
		ro = 'red'+str(redc)
	if seec == 100:
		so = ''
	else:
		so = 'see'+str(seec)
	if skyc == 100:
		sko = ''
	else:
		sko = 'sky'+str(skyc)
	if stm == 10000:
		sto = ''
	else:
		sto = 'st'+str(stm)
		
	fo = open('pc2pt/see'+ro+so+sko+sto+'hodens64.dat','w')
	f = open('pc2pt/LRGcasBwredab0_01sghpix256.dat')
	seef = open('pc2pt/healsee256.dat')
	skyf = open('pc2pt/healsky256.dat')
	redf = open('pc2pt/healred256.dat')
	stf = open('pc2pt/allstars17.519.9n_mHealhpix256.dat')
	sump = 0
	sumn = 0
	for line in f:
		see = float(seef.readline())
		sky = float(skyf.readline())
		red = float(redf.readline())
		ln = line.split()
		p = float(ln[3])
		st = float(stf.readline().split()[2])/p
		if p > t and see < seec and sky < skyc and red < redc and st < stm:
			#od = float(ln[2])/ave - p
			od = see/ave - 1.
			fo.write(ln[0]+' '+ln[1]+' '+str(od)+' '+str(p)+'\n')
			sump += p
	fo.close()
	print sump,
 	from pcsourcezw import pcsourcezw
 	a = pcsourcezw.createSourcesLE_zw('pc2pt/seehodens64.dat')
 	for i in range(0,20):
 		a = pcsourcezw.createSourcesLE_jackzw('pc2pt/seehodens64.dat',i,jackt=20)
	return True

def filestripDR8h1(file,t=.9,md='csv'):
	ml = []
	for i in range(0,12*1024*1024):
		ml.append(0)
	f = open('/Users/ashleyr/boss/shirleyhpixNorth1024.dat')
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			if float(ln[j]) > 0:
				ml[bpix] = float(ln[j])
		i += 1
	print 'north masked'
	f = open('/Users/ashleyr/boss/shirleyhpixSouth1024.dat')
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			if float(ln[j]) > 0:
				ml[bpix] = float(ln[j])
		i += 1
	print 'south masked'
	if md == 'n':
		f = open(file+'.dat')
	if md == 'csv':
		f = open(file+'.csv')
		f.readline()
	fo = open(file+'_mHeal1.dat','w')
	h = healpix()
	for line in f:
		if line[0] != '#':
			if md == 'n':
				ln = line.split()
			if md == 'csv':
				ln = line.split(',')
			ra,dec = float(ln[0]),float(ln[1])
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(1024,th,phi))
			if float(ml[p]) > t:
				fo.write(line)
	fo.close()
	return True
	
def maskfracDR8h(file,bres=1024,ores=256,tb=.9,to=0,reg='NS'):
	mm = 0
	if ores != bres:
		mm =1
	ml = []
	for i in range(0,12*1024*1024):
		ml.append(0)
	ml2 = []
	for i in range(0,12*256*256):
		ml2.append(0)
	i = 0
	h = healpix()
	np = 0
	wsum = 0
	if reg == 'NS' or reg == 'N':
		f = open('/Users/ashleyr/boss/shirleyhpixNorth1024.dat')
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				wsum += w
				if w > tb:
					np += 1.
					if mm == 1:
						th,phi = h.pix2ang_nest(bres,bpix)
						opix = int(h.ang2pix_nest(ores,th,phi))
						ml2[opix] += w*(ores/float(bres))**2
			i += 1
		print 'north masked'
	if reg == 'NS' or reg == 'S':
		f = open('/Users/ashleyr/boss/shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				wsum += w
				if w > tb:
					np += 1.
					if mm == 1:
						th,phi = h.pix2ang_nest(bres,bpix)
						opix = int(h.ang2pix_nest(ores,th,phi))
						ml2[opix] += w*(ores/float(bres))**2
			i += 1
		print 'south masked'
	print np,wsum
	if ores == 1024:
		return np/(12.*1024*1024)
	sum = 0
	for i in range(0,len(ml2)):
		m = ml2[i]
		if m > to:
			sum += 1.
	return sum/(12.*ores*ores)
	
	

	
def filestripDR8h(file,bres=1024,ores=256,t=.9,cpc=0,mm=0,md='n',rc=0,reg='NS'):
	
	ml = []
	for i in range(0,12*1024*1024):
		ml.append(0)
	ml2 = []
	for i in range(0,12*ores*ores):
		ml2.append(0)
	i = 0
	h = healpix()
	np = 0
	if reg == 'NS' or reg == 'N':
		f = open('/Users/ashleyr/boss/shirleyhpixNorth1024.dat')
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > t:
					np += w
					ml[bpix] = 1.
					if mm == 1:
						th,phi = h.pix2ang_nest(bres,bpix)
						opix = int(h.ang2pix_nest(ores,th,phi))
						ml2[opix] += 1.*(ores/float(bres))**2
			i += 1
		print 'north masked'
	if reg == 'NS' or reg == 'S':
		f = open('/Users/ashleyr/boss/shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > t:
					np += w
					ml[bpix] = 1.
					if mm == 1:
						th,phi = h.pix2ang_nest(bres,bpix)
						opix = int(h.ang2pix_nest(ores,th,phi))
						ml2[opix] += 1.*(ores/float(bres))**2
			i += 1
		print 'south masked'
	if mm == 1:
		fm = open('hmaskDR8'+str(ores)+'.dat','w')
		for i in range(0,len(ml2)):
			fm.write(str(ml2[i])+'\n')
		fm.close()
	if md == 'n':	
		f = open(file+'.dat')
	if md == 'csv':
		f = open(file+'.csv')
	fo = open(file+'_mHeal'+reg+'.dat','w')
	for line in f:
		try:
			if md == 'n':
				ln = line.split()
			if md == 'csv':
				ln = line.split(',')
			ra,dec = float(ln[rc]),float(ln[rc+1])
			th,phi = radec2thphi(ra,dec)
			#cp = float(ln[cpcol])
			p = int(h.ang2pix_nest(bres,th,phi))
			if p > len(ml):
				print len(ml),p,th,phi,h.ang2pix_nest(bres,th,phi)
			if ml[p] > t:# and cp > cpc:
				for k in range(0,len(ln)-1):
					fo.write(ln[k]+' ')
				fo.write(ln[-1].strip('\n')+'\n')
		except:
			pass
	fo.close()
	return np

def mkrandomsDR8(N,tol):
	#this module put N objects randomly distributed in lamda/eta, without any masking
	from random import random
	ml = []
	for i in range(0,12*1024*1024):
		ml.append(0)
	f = open('shirleyhpixNorth1024.dat')
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			if float(ln[j]) > 0:
				ml[bpix] = float(ln[j])
		i += 1
	print 'north masked'
	f = open('shirleyhpixSouth1024.dat')
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			if float(ln[j]) > 0:
				ml[bpix] = float(ln[j])
		i += 1
	print 'south masked'

	fo = open('randomshDR8.dat','w')
	Nrandoms = 0
	h = healpix()
	while Nrandoms < N:
		ranx = -1.+2.*random()
		rany = -1.+2.*random()
		ranz = -1.+2.*random()
		r = sqrt(ranx**2.+rany**2.+ranz**2.)
		if r < 1:
			theta = acos(ranz/r)
			if rany > 0:
				phi = acos(ranx/(r*sin(theta)))				
			else:
				phi = acos(ranx/(r*sin(theta)))+pi
			pix = int(h.ang2pix_nest(1024,theta,phi))
			if ml[pix] > .9:
				ra,dec = thphi2radec(theta,phi)
				fo.write(str(ra)+' '+str(dec)+'\n')
				Nrandoms += 1.
	fo.close()
	return True

def mkrandomsDR83D(zmin,zmax,N,tol=.9):
	#this module put N objects randomly distributed in lamda/eta, without any masking
	from random import random
	from Cosmo import *
	d = distance(.25,.75)
	ml = []
	for i in range(0,12*1024*1024):
		ml.append(0)
	f = open('/Users/ashleyr/BOSS/shirleyhpixNorth1024.dat')
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			if float(ln[j]) > tol:
				ml[bpix] = 1
		i += 1
	print 'north masked'
	f = open('/Users/ashleyr/BOSS/shirleyhpixSouth1024.dat')
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			if float(ln[j]) > tol:
				ml[bpix] = 1
		i += 1
	print 'south masked'

	fo = open('randomshDR8z'+str(zmin)+str(zmax)+'.dat','w')
	Nrandoms = 0
	h = healpix()
	fz = open('/Users/ashleyr/BOSS/nzLRGcasBwredab1.6_'+str(zmin)+str(zmax)+'sg.dat').readlines()
	step = float(fz[1].split()[0])-float(fz[0].split()[0])
	zmax = 0
	for i in range(0,len(fz)):
		nz = float(fz[i].split()[1])
		if nz > zmax:
			zmax = nz
	print zmax
	sl = []
	for i in range(0,12*256*256):
		sl.append(0)
	fs = open('/Users/ashleyr/BOSS/pc2pt/allstars17.519.9n_mHealhpix256.dat')
	for line in fs:
		ln = line.split()
		lam,eta = float(ln[0]),float(ln[1])
		th,phi = le2thetaphi(lam,eta)
		pix = int(h.ang2pix_nest(256,th,phi))
		sl[pix] = float(ln[2])/float(ln[-1])
	#areaarc = 679827.85262104508 #for a full pixel at 256
	#r1 = 7.56
	#r2 = 10.6
	#r3 = 12.0
	#r4 = 10.9
	arc2 = 169956.96315526127
	rs = 5.25
	nb = 0
	while Nrandoms < N:
		ranx = -1.+2.*random()
		rany = -1.+2.*random()
		ranz = -1.+2.*random()
		r = sqrt(ranx**2.+rany**2.+ranz**2.)
		if r < 1:
			theta = acos(ranz/r)
			if rany > 0:
				phi = acos(ranx/(r*sin(theta)))				
			else:
				phi = acos(ranx/(r*sin(theta)))+pi
			pix = int(h.ang2pix_nest(1024,theta,phi))
			if ml[pix] == 1:
				pix2 = int(h.ang2pix_nest(256,theta,phi))
				ns = sl[pix2]
				if ns == 0:
					nb += 1

				ra,dec = thphi2radec(theta,phi)
				z = random()*.2+zmin
				zind = int(z/step)
				#rs = 0
				#if z >= .45 and z < 0.5:
				#	rs = r1
				#if z >=.5 and z < 0.55:
				#	rs = r2
				#if z >= 0.55 and z < .6:
				#	rs = r3
				#if z >= 0.6 and z < 0.65:
				#	rs = r4
				#if rs == 0:
				#	print z
				wast = 1.-pi*rs**2./arc2*ns
				if wast < 0:
					wast = 0
				zw = float(fz[zind].split()[-1])/zmax*wast
				dp = d.dc(z)
				fo.write(str(ra)+' '+str(dec)+' '+str(dp)+' '+str(zw)+'\n')
				Nrandoms += 1.
	print nb
	fo.close()
	return True


def maskup(bres=1024,ores=256):
	f = open('MegaZDR7/megaz_dr7.'+str(bres)+'.mask.dat')
	h = healpix()
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		ml.append(0)
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			pxb = i*bres+j
			th,phi = h.pix2ang_nest(bres,pxb)
			pxo = int(h.ang2pix_nest(ores,th,phi))
			ml[pxo] += float(ln[j])
		i += 1
	fo = open('MegaZDR7/megaz_dr7.'+str(ores)+'.mask.dat','w')
	for i in range(0,len(ml)):
		fo.write(str(ml[i])+'\n')
	fo.close()
	return True

def fpixl(file,fac=16.,res=256):
	f = open(file+'hpix'+str(res)+'.dat').readlines()
	fo = open(file+'hpix'+str(res)+'.dat','w')
	for i in range(0,len(f)):
		ln = f[i].split()
		w = float(ln[-1])/fac
		fo.write(ln[0]+' '+ln[1]+' '+ln[2]+' '+str(w)+'\n')
	fo.close()

def filestripM7h(file,res=1024,rc=0,dc=1,md='n'):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*res**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	mm = 0
	f = open('MegaZDR7/megaz_dr7.'+str(res)+'.mask.dat')
	i = 0
	for line in f:
		ln = line.split()		
		for j in range(0,len(ln)):
			px = i*1024+j
			ml[px] = float(ln[j])
		i += 1
	print 'masked'
	if md == 'n':
		f = open(file+'.dat')
	if md == 'csv':
		f = open(file+'.csv')
	fo = open(file+'_M7Heal.dat','w')	
	for line in f:
		try:
			if md == 'n':
				ln = line.split()
			if md == 'csv':
				ln = line.split(',')
			ra,dec = float(ln[rc]),float(ln[dc])
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))	
			if ml[p] > 0:
				for i in range(0,len(ln)-1):
					fo.write(ln[i]+' ')
				fo.write(ln[-1].strip('\n')+'\n')
		except:
			print line
	fo.close()
	return True

def mkM7healpixl(file,bres=1024,ores=1024,rc=0,dc=1,md='n',mm=0):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	#mm = 0
	f = open('MegaZDR7/megaz_dr7.'+str(ores)+'.mask.dat')
	i = 0
	for line in f:
		if ores == 1024:
			ln = line.split()		
			for j in range(0,len(ln)):
				px = i*1024+j
				ml[px] = float(ln[j])
		else:
			ml[i] = float(line)*(ores/float(bres))**2.
		i += 1
	if mm == 1:
		fmo = open('MZhmask'+str(ores)+'.dat','w')
		for i in range(0,len(ml)):
			fmo.write(str(ml[i])+'\n')
		fmo.close()
	
	print 'masked'
	if md == 'n':
		f = open(file+'.dat')
	if md == 'csv':
		f = open(file+'.csv')
	for line in f:
		if md == 'n':
			ln = line.split()
		if md == 'csv':
			ln = line.split(',')
		ra,dec = float(ln[rc]),float(ln[dc])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(ores,th,phi))
		pixl[p] += 1.
		
	fo = open(file+'hpix'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(ores,i)
			lam,eta = thphi2le(th,phi)
			fo.write(str(lam)+' '+str(eta)+' '+str(pixl[i])+' '+str(ml[i])+'\n')
	fo.close()
	return True

def mkM7healpixlsg(file,bres=1024,ores=1024,rc=0,dc=1,sgcol=-1):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	mm = 0
	f = open('MegaZDR7/megaz_dr7.'+str(ores)+'.mask.dat')
	i = 0
	for line in f:
		if ores == 1024:
			ln = line.split()		
			for j in range(0,len(ln)):
				px = i*1024+j
				ml[px] = float(ln[j])
		else:
			ml[i] = float(line)*(ores/float(bres))**2.
		i += 1
	print 'masked'
	f = open(file+'.dat')
	for line in f:
		ln = line.split()
		ra,dec = float(ln[rc]),float(ln[dc])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(ores,th,phi))
		pixl[p] += float(ln[sgcol])
		
	fo = open(file+'sghpix'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(ores,i)
			lam,eta = thphi2le(th,phi)
			fo.write(str(lam)+' '+str(eta)+' '+str(pixl[i])+' '+str(ml[i])+'\n')
	fo.close()
	return True

def pixl_up(file,reso,resn):
	h = healpix()
	f = open(file+str(reso)+'.dat')
	ol = []
	for line in f:
		ol.append(float(line))
	nl = []
	np = 12*resn*resn
	for i in range(0,np):
		nl.append(0)
	for i in range(0,len(ol)):
		th,phi = h.pix2ang_nest(reso,i)
		p = int(h.ang2pix_nest(resn,th,phi))
		nl[p] += ol[i]
	fo = open(file+str(resn)+'.dat','w')
	for i in range(0,np):
		fo.write(str(nl[i])+'\n')
	fo.close()
	return True

def ranpixl_up(file,reso,resn):
	h = healpix()
	f = open('ranHeal_pix'+str(reso)+file+'.dat')
	ol = []
	npo = 12*reso*reso
	for i in range(0,npo):
		ol.append(0)
	for line in f:
		ln = line.split()
		p = int(ln[0])
		ol[p] += float(ln[1])
	nl = []
	np = 12*resn*resn
	for i in range(0,np):
		nl.append(0)
	for i in range(0,len(ol)):
		th,phi = h.pix2ang_nest(reso,i)
		p = int(h.ang2pix_nest(resn,th,phi))
		nl[p] += ol[i]
	fo = open('ranHeal_pix'+str(resn)+file+'.dat','w')
	for i in range(0,np):
		if nl[i] != 0:
			fo.write(str(i)+' '+str(nl[i])+'\n')
	ol = []
	nl = []
	fo.close()
	return True
	
		

def mkhealpixl_simp(file,res=256,rc=0,dc=1,zc=2,md='csv'):
	h = healpix()
	pixl = []
	npo = 12*res**2
	for i in range(0,npo):
		pixl.append(0)
	
	f = open(file+'.'+md)
	f.readline()
	n = 0
	for line in f:
		if line[0] != '#':
			if md == 'csv':
				ln = line.split(',')
			else:
				ln = line.split()
			try:
				ra,dec = float(ln[rc]),float(ln[dc])
				th,phi = radec2thphi(ra,dec)
				p = int(h.ang2pix_nest(res,th,phi))
				pixl[p] += 1.
				n += 1
			except:
				pass
	print n
	fo = open(file+'hpixall'+str(res)+'.dat','w')
	for i in range(0,npo):
		th,phi = h.pix2ang_nest(res,i)
		fo.write(str(pixl[i])+'\n')
	fo.close()
	return True

def mkhealDR8odens_4alm(file,res=256,rc=0,dc=1,zc=2):
	h = healpix()
	pixl = []
	npo = 12*res**2
	for i in range(0,npo):
		pixl.append(0)

	f = open(file+'.dat')
	for line in f:
		if line[0] != '#':
			ln = line.split()
			ra,dec = float(ln[rc]),float(ln[dc])
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(res,th,phi))
			pixl[p] += 1.
	f = open('/Users/ashleyr/boss/hmaskDR8'+str(res)+'.dat')
	ml = []
	for line in f:
		ml.append(float(line))
	if len(ml) != len(pixl):
		return 'mismatched lists' + str(len(ml))+' '+str(len(pixl))
	ng = 0
	at = 0
	for i in range(0,len(ml)):
		ng += pixl[i]
		at += ml[i]
	nave = ng/at
	fo = open('/Users/ashleyr/healpix/'+file+'odens'+str(res)+'.dat','w')
	for i in range(0,len(ml)):
		od = 0
		if ml[i] > 0:
			od = (pixl[i]/ml[i]-nave)/nave
		fo.write(str(od)+' '+str(ml[i])+'\n')
	fo.close()
	return True


def mkDR8healmaskN(bres=1024,ores=256,rc=0,dc=1,zc=2):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	
	fo = open('hmaskDR8N'+str(ores)+'.dat','w')
	f = open('shirleyhpixNorth1024.dat')
	i = 0
	for line in f:
		ln = line.split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			w = float(ln[j])
			if w > 0:
				th,phi = h.pix2ang_nest(bres,bpix)
				opix = int(h.ang2pix_nest(ores,th,phi))
				ml[opix] += w*(ores/float(bres))**2
		i += 1
	print 'north masked'
	for i in range(0,npo):
		fo.write(str(ml[i])+'\n')
	fo.close()
	return True

def mkdpoffN(bres=1024,ores=256):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	
	fo = open('dpoffN'+str(ores)+'.dat','w')
	fri = open('schlafly_ri_f.dat')
	fgr = open('schlafly_gr_f.dat')
	i = 0
	for line in fri:
		ln = line.split()
		lngr = fgr.readline().split()
		for j in range(0,len(ln)):
			bpix = 1024*i+j
			w = float(ln[j])
			if w > 0:
				th,phi = h.pix2ang_nest(bres,bpix)
				opix = int(h.ang2pix_nest(ores,th,phi))
				gr = float(lngr[j])
				dp = w-gr/8.
				ml[opix] += 1.
				pixl[opix] += dp
		i += 1
	for i in range(0,len(ml)):
		if ml[i] > 0:
			pixl[i] = pixl[i]/ml[i]
	for i in range(0,npo):
		fo.write(str(pixl[i])+'\n')
	fo.close()
	return True

def densDR8_NS(file,NS,md='.dat',bres=1024,rc=0,dc=1,zc=2,tol=.9):
	h = healpix()
	pixl = []
	zpixl = []
	ml = []
	npo = 12*bres**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	mt = 0	
	if NS == 'N':
		f = open('/Users/ashleyross/BOSS/shirleyhpixNorth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > tol:
					ml[bpix] = 1.
					mt += 1.
			i += 1
		print 'north masked'
	if NS == 'S':
		f = open('/Users/ashleyross/BOSS/shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > tol:
					ml[bpix] += 1.
					mt += 1.
			i += 1
		print 'south masked'
	f = open(file+md)
	ng = 0
	for line in f:
		if md == '.dat':
			ln = line.split()
		if md == '.csv':
			ln = line.split(',')
		ra,dec = float(ln[rc]),float(ln[dc])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(bres,th,phi))
		if ml[p] > 0:
			ng += 1.
	return ng/mt

def mkDR8healpixl(file,bres=1024,ores=256,rc=0,dc=1,zc=2):
	h = healpix()
	pixl = []
	zpixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
		zpixl.append(0)
	mm = 0
	try:
		f = open('hmaskDR8'+str(ores)+'.dat')
		i = 0
		for line in f:
			ml[i] = float(line)
			i += 1
	except:
		mm = 1	
		f = open('shirleyhpixNorth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'north masked'
		f = open('shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'south masked'
	f = open(file+'.dat')
	for line in f:
		ln = line.split()
		ra,dec = float(ln[rc]),float(ln[dc])
		z = float(ln[zc])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(ores,th,phi))
		pixl[p] += 1.
		zpixl[p] += z
		
	fo = open(file+'hpix'+str(ores)+'.dat','w')
	if mm == 1:
		fm = open('hmaskDR8'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if mm == 1:
			fm.write(str(ml[i])+'\n')
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(ores,i)
			#ra,dec = thphi2radec(th,phi)
			lam,eta = thphi2le(th,phi)
			fo.write(str(lam)+' '+str(eta)+' '+str(pixl[i])+' '+str(ml[i])+'\n')
	fo.close()
	if mm == 1:
		fm.close()
	return True

def mkDR8healpixl_weight(file,bres=1024,ores=256,rc=0,dc=1,zc=2):
	h = healpix()
	pixl = []
	zpixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
		zpixl.append(0)
	mm = 0
	try:
		f = open('hmaskDR8'+str(ores)+'.dat')
		i = 0
		for line in f:
			ml[i] = float(line)
			i += 1
	except:
		mm = 1	
		f = open('shirleyhpixNorth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'north masked'
		f = open('shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'south masked'
	f = open(file+'.dat')
	for line in f:
		if line[0] != '#':
			ln = line.split()
			ra,dec = float(ln[rc]),float(ln[dc])
			z = float(ln[zc])
			th,phi = radec2thphi(ra,dec)
			wt = float(ln[3])
			p = int(h.ang2pix_nest(ores,th,phi))
			pixl[p] += wt
			zpixl[p] += z
		
	fo = open(file+'hpix'+str(ores)+'.dat','w')
	if mm == 1:
		fm = open('hmaskDR8'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if mm == 1:
			fm.write(str(ml[i])+'\n')
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(ores,i)
			#ra,dec = thphi2radec(th,phi)
			lam,eta = thphi2le(th,phi)
			fo.write(str(lam)+' '+str(eta)+' '+str(pixl[i])+' '+str(ml[i])+'\n')
	fo.close()
	if mm == 1:
		fm.close()
	return True

	
def mkDR8healpixlsg(file,bres=1024,ores=256,rc=0,dc=1,sgc=4,zc=2):
	h = healpix()
	pixl = []
	zpixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
		zpixl.append(0)
	mm = 0
	try:
		mf = 'hmaskDR8'+str(ores)+'.dat'
		print mf
		f = open(mf)
		i = 0
		for line in f:
			ml[i] = float(line)
			i += 1
	except:
		print 'exception',i,len(ml)
		mm = 1	
		f = open('shirleyhpixNorth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'north masked'
		f = open('shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'south masked'
	f = open(file+'.dat')
	for line in f:
		if line[0] == '#':
			pass
		else:
			ln = line.split()
			z = float(ln[zc])
			ra,dec = float(ln[rc]),float(ln[dc])
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(ores,th,phi))
			sg = float(ln[sgc])
			if sg < 0:
				sg = 0
			if sg > 1.:
				sg = 1.
			pixl[p] += sg
			zpixl[p] += z
		
	fo = open(file+'sghpix'+str(ores)+'.dat','w')
	if mm == 1:
		fm = open('hmaskDR8'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if mm == 1:
			fm.write(str(ml[i])+'\n')
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(ores,i)
			#ra,dec = thphi2radec(th,phi)
			lam,eta = thphi2le(th,phi)
			fo.write(str(lam)+' '+str(eta)+' '+str(pixl[i])+' '+str(ml[i])+'\n')
	fo.close()
	if mm == 1:
		fm.close()
	return True

def mkDR8healpixlsgwst(file,bres=1024,ores=256,rc=0,dc=1,sgc=4,zc=2,ast=4.74):
	h = healpix()
	pixl = []
	zpixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
		zpixl.append(0)
	mm = 0
	try:
		mf = '/Users/ashleyr/BOSS/hmaskDR8'+str(ores)+'.dat'
		print mf
		f = open(mf)
		i = 0
		for line in f:
			ml[i] = float(line)
			i += 1
	except:
		print 'exception',i,len(ml)
		mm = 1	
		f = open('/Users/ashleyr/BOSS/shirleyhpixNorth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'north masked'
		f = open('/Users/ashleyr/BOSS/shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'south masked'
	f = open(file+'.dat')
	for line in f:
		if line[0] == '#':
			pass
		else:
			ln = line.split()
			z = float(ln[zc])
			ra,dec = float(ln[rc]),float(ln[dc])
			th,phi = radec2thphi(ra,dec)
			p = int(h.ang2pix_nest(ores,th,phi))
			sg = float(ln[sgc])
			if sg < 0:
				sg = 0
			if sg > 1.:
				sg = 1.
			pixl[p] += sg
			zpixl[p] += z
		
	fo = open('LRGcmass_sg_Nside'+str(ores)+'.dat','w')
	if mm == 1:
		fm = open('hmaskDR8'+str(ores)+'.dat','w')
	stf = open('allstars17.519.9n_mHealhpix256.dat')
	arc2 = 169956.96315526127
	for i in range(0,npo):
		if mm == 1:
			fm.write(str(ml[i])+'\n')
		if ml[i] > 0:
			stln = stf.readline().split()
			th,phi = h.pix2ang_nest(ores,i)
			#ra,dec = thphi2radec(th,phi)
			lam,eta = thphi2le(th,phi)
			st = float(stln[2])
			p = ml[i]
			pm = st*pi*(ast)**2./arc2
			ps = p-pm
			if ps < 0:
				ps = 0
			fo.write(str(pixl[i])+' '+str(ml[i])+' '+str(ps)+'\n')
		else:
			fo.write('0 0 0\n')
	fo.close()
	if mm == 1:
		fm.close()
	return True


def mkDR8healpixlsgzp(file,zbinc,bres=1024,ores=256,rc=0,dc=1,sgc=4):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	mm = 0
	try:
		mf = 'hmaskDR8'+str(ores)+'.dat'
		print mf
		f = open(mf)
		i = 0
		for line in f:
			ml[i] = float(line)
			i += 1
	except:
		print 'exception',i,len(ml)
		mm = 1	
		f = open('shirleyhpixNorth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'north masked'
		f = open('shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'south masked'
	f = open(file+'.dat')
	for line in f:
		ln = line.split()
		ra,dec = float(ln[rc]),float(ln[dc])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(ores,th,phi))
		sg = float(ln[sgc])
		zb = float(ln[zbinc])
		if sg < 0:
			sg = 0
		if sg > 1.:
			sg = 1.
		pixl[p] += sg*zb
		
	fo = open(file+'sgzp'+str(zbinc)+'hpix'+str(ores)+'.dat','w')
	if mm == 1:
		fm = open('hmaskDR8'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if mm == 1:
			fm.write(str(ml[i])+'\n')
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(ores,i)
			#ra,dec = thphi2radec(th,phi)
			lam,eta = thphi2le(th,phi)
			fo.write(str(lam)+' '+str(eta)+' '+str(pixl[i])+' '+str(ml[i])+'\n')
	fo.close()
	if mm == 1:
		fm.close()
	return True

def mkDR8healpixlzp(file,zbinc,bres=1024,ores=256,rc=0,dc=1,sgc=4):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	mm = 0
	try:
		mf = 'hmaskDR8'+str(ores)+'.dat'
		print mf
		f = open(mf)
		i = 0
		for line in f:
			ml[i] = float(line)
			i += 1
	except:
		print 'exception',i,len(ml)
		mm = 1	
		f = open('shirleyhpixNorth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'north masked'
		f = open('shirleyhpixSouth1024.dat')
		i = 0
		for line in f:
			ln = line.split()
			for j in range(0,len(ln)):
				bpix = 1024*i+j
				w = float(ln[j])
				if w > 0:
					th,phi = h.pix2ang_nest(bres,bpix)
					opix = int(h.ang2pix_nest(ores,th,phi))
					ml[opix] += w*(ores/float(bres))**2
			i += 1
		print 'south masked'
	f = open(file+'.dat')
	for line in f:
		ln = line.split()
		ra,dec = float(ln[rc]),float(ln[dc])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(ores,th,phi))
		sg = float(ln[sgc])
		zb = float(ln[zbinc])
		if sg < 0:
			sg = 0
		if sg > 1.:
			sg = 1.
		pixl[p] += zb
		
	fo = open(file+'zp'+str(zbinc)+'hpix'+str(ores)+'.dat','w')
	if mm == 1:
		fm = open('hmaskDR8'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if mm == 1:
			fm.write(str(ml[i])+'\n')
		if ml[i] > 0:
			th,phi = h.pix2ang_nest(ores,i)
			#ra,dec = thphi2radec(th,phi)
			lam,eta = thphi2le(th,phi)
			fo.write(str(lam)+' '+str(eta)+' '+str(pixl[i])+' '+str(ml[i])+'\n')
	fo.close()
	if mm == 1:
		fm.close()
	return True


def healfits2radec(file,nside):
	f = open(file+'.dat').readlines()
	#h = healpix()
	fo = open(file+'rd.dat','w')
	for i in range(0,len(f)):
		row = f[i].split()
		for j in range(0,len(row)):
			#pix = int((j/nside)*(12*nside*nside**2./1024)+fmod(j,nside)+i*nside)
			pix = 1024*i+j
			th,phi = pix2ang_ring(nside, pix)
			#print pix,th,phi
			ra,dec = thphi2radec(th,phi)
			fo.write(str(ra)+' '+str(dec)+' '+row[j]+'\n')
	fo.close()

def mkschl():
	fsc = open('schlafly_masked2rd.dat').readlines()
	pl = open('LRGcasBwredsg-1_01sk100p.75lameta64.dat')
	fo = open('DR8LRGsm64p.75sch.dat','w')
	for line in pl:
		ln = line.split()
		lam,eta = float(ln[0]),float(ln[1])
		th,phi = le2thetaphi(lam,eta)
		pix = ang2pix_ring(128,th,phi)
		sch = fsc[pix].split()[2]
		fo.write(sch+'\n')
	fo.close()

def schl2hp(res=256):
	fsch = open('pc2pt/schlafly_masked2rd.dat').readlines()
	f = open('/Users/ashleyr/BOSS/pc2pt/LRGcasBwredab0_01sghpix256.dat')
	#h = healpix()
	fo = open('pc2pt/healschri256.dat','w')
	for line in f:
		ln = line.split()
		lam,eta = float(ln[0]),float(ln[1])
		th,phi = le2thetaphi(lam,eta)
		p = int(ang2pix_ring(128,th,phi))
		#print p
		sch = float(fsch[p].split()[2])
		if sch == 0:
			sch = -1
		fo.write(str(sch)+'\n')
	fo.close()

def thphi2radec(theta,phi):
	return 180./pi*phi,-(180./pi*theta-90)

def thphi2le(theta,phi):
	deg2Rad = pi/180.0
	rarad = phi
	decrad = -(theta-piover2)
	surveyCenterDEC = 32.5
	surveyCenterRA = 185.0
	etaPole = deg2Rad*surveyCenterDEC
	node = deg2Rad*(surveyCenterRA - 90.0)
	
	x = cos(rarad-node)*cos(decrad)
  	y = sin(rarad-node)*cos(decrad)
  	z = sin(decrad)

  	lam = -1.0*asin(x)/deg2Rad
  	eta = (atan2(z,y) - etaPole)/deg2Rad
  	if eta < -180.0:
  		eta += 360.0
  	if eta > 180.0:
  		eta -= 360.0
  	
  	return (lam,eta)

def le2thetaphi(lam,eta):
	deg2Rad = pi/180.0
	surveyCenterDEC = 32.5
	surveyCenterRA = 185.0
	etaPole = deg2Rad*surveyCenterDEC
	node = deg2Rad*(surveyCenterRA - 90.0)
	x = -1.0*sin(lam*deg2Rad)
	y = cos(lam*deg2Rad)*cos(eta*deg2Rad+etaPole)
	z = cos(lam*deg2Rad)*sin(eta*deg2Rad+etaPole)
	ra = (atan2(y,x) + node)
	if ra < 0.0:
		ra +=twopi
	dec = asin(z)
	return -dec+piover2,ra


def mkhealpixl_nm(file,ores=256,rc=0,dc=1):
	h = healpix()
	pixl = []
	ml = []
	npo = 12*ores**2
	for i in range(0,npo):
		pixl.append(0)
		ml.append(0)
	f = open(file+'.dat')
	for line in f:
		ln = line.split()
		ra,dec = float(ln[rc]),float(ln[dc])
		th,phi = radec2thphi(ra,dec)
		p = int(h.ang2pix_nest(ores,th,phi))
		pixl[p] += 1.
		
	fo = open(file+'hpix_nm'+str(ores)+'.dat','w')
	for i in range(0,npo):
		if pixl[i] > 0:
			fo.write(str(i)+' '+str(pixl[i])+'\n')
	fo.close()
	return True


def radec2thphi(ra,dec):
	return (-dec+90.)*pi/180.,ra*pi/180.

def mkthphi_nest(res):
	f = open('/Users/ashleyr/healpix/thphi_nest'+str(res)+'.dat','w')
	h = healpix()
	n = 12*res*res
	for i in range(0,n):
		th,phi = h.pix2ang_nest(res,i)
		f.write(str(th)+' '+str(phi)+'\n')
	f.close()
	return True

def mkYlm(res,l,m):
	f = open('/Users/ashleyr/healpix/thphi_nest'+str(res)+'.dat').readlines()
	fo = open('/Users/ashleyr/healpix/Y'+str(l)+str(m)+'_nest'+str(res)+'.dat','w')
	for i in range(0,len(f)):
		pixl = f[i].split()
		th = float(pixl[0])
		phi = float(pixl[1])
		ylm = sqrt((2.*l+1)*factorial(l-m)/(4.*pi*factorial(l+m)))*legendre(l,cos(th),m)*ee**(phi*m*1.j)
		fo.write(str(ylm.real)+' '+str(ylm.imag)+'\n')
	return True

def alm(file,l,m,res=256):
	strp = 4.*pi/(12.*res*res)
	pixl = open('/Users/ashleyr/healpix/'+file+'odens'+str(res)+'.dat').readlines()
	Ylm = open('/Users/ashleyr/healpix/Y'+str(l)+str(m)+'_nest'+str(res)+'.dat').readlines()
	if len(Ylm) != len(pixl):
		return 'mismatched lists'
	alm = 0
	tot = 0
	for i in range(0,len(pixl)):
		Ylmcon = float(Ylm[i].split()[1])#double check this is correct
		w = float(pixl[i].split()[1])
		od = float(pixl[i].split()[0])
		alm += od*w*Ylmcon
		tot += Ylmcon*w
	alm = alm*strp
	tot = tot*strp
	print tot
	return alm


class healpix:
	def __init__(self):
		self.pix2x,self.pix2y = self.mk_pix2xy()
		self.x2pix,self.y2pix = self.mk_xy2pix()

	def ang2pix_nest(self,nside,theta,phi):
		#x2pix,y2pix = mk_xy2pix()
		z  = cos(theta)
		za = fabs(z)
		z0 = 2./3.
		if phi>=twopi:
			phi = phi - twopi
		if phi<0.:
			phi = phi + twopi
		tt = phi / piover2
		if za<=z0:  #{ /* equatorial region */
		
		#/* (the index of edge lines increase when the longitude=phi goes up) */
			jp = int(floor(ns_max*(0.5 + tt - z*0.75)))# /* ascending edge line index */
			jm = int(floor(ns_max*(0.5 + tt + z*0.75)))#; /* descending edge line index */
			
			#/* finds the face */
			ifp = jp / ns_max#; /* in {0,4} */
			ifm = jm / ns_max
			
			if ifp==ifm:
				face_num = int(fmod(ifp,4)) + 4#; /* faces 4 to 7 */
			else:
				if ifp<ifm:
					face_num = int(fmod(ifp,4)) #/* (half-)faces 0 to 3 */
				else:
					face_num = int(fmod(ifm,4)) + 8#;           /* (half-)faces 8 to 11 */
			
			ix = int(fmod(jm, ns_max))
			iy = ns_max - int(fmod(jp, ns_max)) - 1
	
		else:# { /* polar region, za > 2/3 */
		
			ntt = int(floor(tt))
			if ntt>=4:
				ntt = 3
			tp = tt - ntt
			tmp = sqrt( 3.*(1. - za) )#; /* in ]0,1] */
		
		#/* (the index of edge lines increase when distance from the closest pole
		# * goes up)
		# */
		#/* line going toward the pole as phi increases */
			jp = int(floor( ns_max * tp * tmp )) 
		
		#/* that one goes away of the closest pole */
			jm = int(floor( ns_max * (1. - tp) * tmp ))
			if jp >= ns_max:
				jp = ns_max-1
			if jm >= ns_max:
				jm = ns_max-1
			#jp = int((jp < ns_max-1 ? jp : ns_max-1))
			#jm = (int)(jm < ns_max-1 ? jm : ns_max-1);
		
		#/* finds the face and pixel's (x,y) */
			if z>=0:# ) {
				face_num = ntt#; /* in {0,3} */
				ix = ns_max - jm - 1
				iy = ns_max - jp - 1
			else:
				face_num = ntt + 8#; /* in {8,11} */
				ix =  jp
				iy =  jm	
	
		ix_low = int(fmod(ix,128))
		ix_hi  = ix/128
		iy_low = int(fmod(iy,128))
		iy_hi  = iy/128
	
		ipf = (self.x2pix[ix_hi]+self.y2pix[iy_hi]) * (128 * 128)+ (self.x2pix[ix_low]+self.y2pix[iy_low]);
		ipf = (long)(ipf / pow(ns_max/nside,2))#;     /* in {0, nside**2 - 1} */
		return int( ipf + face_num*pow(nside,2))#; /* in {0, 12*nside**2 - 1} */
	
	def mk_xy2pix(self):
	#   /* =======================================================================
	#    * subroutine mk_xy2pix
	#    * =======================================================================
	#    * sets the array giving the number of the pixel lying in (x,y)
	#    * x and y are in {1,128}
	#    * the pixel number is in {0,128**2-1}
	#    *
	#    * if  i-1 = sum_p=0  b_p * 2^p
	#    * then ix = sum_p=0  b_p * 4^p
	#    * iy = 2*ix
	#    * ix + iy in {0, 128**2 -1}
	#    * =======================================================================
	#    */
	#  int i, K,IP,I,J,ID;
		x2pix = []
		y2pix = []
		for i in range(0,128):#(i = 0; i < 127; i++) x2pix[i] = 0;
			x2pix.append(0)
			y2pix.append(0)
		for I in range(1,129):#( I=1;I<=128;I++ ) {
			J  = I-1#;//            !pixel numbers
			K  = 0#;//
			IP = 1#;//
			while J!=0:
	#    truc : if( J==0 ) {
	#     x2pix[I-1] = K;
	#      y2pix[I-1] = 2*K;
	#    }
	#    else {
				ID = int(fmod(J,2))
				J  = J/2
				K  = IP*ID+K
				IP = IP*4
	#      goto truc;
			x2pix[I-1] = K
			y2pix[I-1] = 2*K
		return x2pix,y2pix
	
	def mk_pix2xy(self): 
	
	#   /* =======================================================================
	#    * subroutine mk_pix2xy
	#    * =======================================================================
	#    * constructs the array giving x and y in the face from pixel number
	#    * for the nested (quad-cube like) ordering of pixels
	#    *
	#    * the bits corresponding to x and y are interleaved in the pixel number
	#    * one breaks up the pixel number by even and odd bits
	#    * =======================================================================
	#    */
	
	#  int i, kpix, jpix, IX, IY, IP, ID;
		pix2x = []
		pix2y = []
		for i in range(0,1024):
			pix2x.append(0)
			pix2y.append(0)
	#  for (i = 0; i < 1023; i++) pix2x[i]=0;
	  
	#  for( kpix=0;kpix<1024;kpix++ ) {
		for kpix in range(0,1024):
			jpix = kpix
			IX = 0
			IY = 0
			IP = 1# ;//              ! bit position (in x and y)
			while jpix!=0:# ){// ! go through all the bits
				ID = int(fmod(jpix,2))#;//  ! bit value (in kpix), goes in ix
				jpix = jpix/2
				IX = ID*IP+IX
				
				ID = int(fmod(jpix,2))#;//  ! bit value (in kpix), goes in iy
				jpix = jpix/2
				IY = ID*IP+IY
				
				IP = 2*IP#;//         ! next bit (in x and y)
			pix2x[kpix] = IX#;//     ! in 0,31
			pix2y[kpix] = IY#;//     ! in 0,31
	  
		return pix2x,pix2y
	
	def pix2ang_nest(self,nside, ipix):
	
	#   /*
	#     c=======================================================================
	#     subroutine pix2ang_nest(nside, ipix, theta, phi)
	#     c=======================================================================
	#     c     gives theta and phi corresponding to pixel ipix (NESTED) 
	#     c     for a parameter nside
	#     c=======================================================================
	#   */
		
		#pix2x,pix2y = mk_pix2xy()
		jrll = []
		jpll = []
		for i in range(0,12):
			jrll.append(0)
			jpll.append(0)
		jrll[0]=2
		jrll[1]=2
		jrll[2]=2
		jrll[3]=2
		jrll[4]=3
		jrll[5]=3
		jrll[6]=3
		jrll[7]=3
		jrll[8]=4
		jrll[9]=4
		jrll[10]=4
		jrll[11]=4
		jpll[0]=1
		jpll[1]=3
		jpll[2]=5
		jpll[3]=7
		jpll[4]=0
		jpll[5]=2
		jpll[6]=4
		jpll[7]=6
		jpll[8]=1
		jpll[9]=3
		jpll[10]=5
		jpll[11]=7
		  
		  
		npix = 12 * nside*nside
		if ipix < 0 or ipix > npix-1:
			return 'ipix out of range'
	
	#      /* initiates the array for the pixel number -> (x,y) mapping */
	
		fn = 1.*nside
		fact1 = 1./(3.*fn*fn)
		fact2 = 2./(3.*fn)
		nl4   = 4*nside
	
	#      //c     finds the face, and the number in the face
		npface = nside*nside
		
		face_num = ipix/npface#//  ! face number in {0,11}
		ipf = int(fmod(ipix,npface))#//  ! pixel number in the face {0,npface-1}
		
	#	//c     finds the x,y on the face (starting from the lowest corner)
	#	//c     from the pixel number
		ip_low = int(fmod(ipf,1024))#;//       ! content of the last 10 bits
		ip_trunc =   ipf/1024# ;//       ! truncation of the last 10 bits
		ip_med = int(fmod(ip_trunc,1024))#;//  ! content of the next 10 bits
		ip_hi  =     ip_trunc/1024   #;//! content of the high weight 10 bits
		
		ix = 1024*self.pix2x[ip_hi] + 32*self.pix2x[ip_med] + self.pix2x[ip_low]
		iy = 1024*self.pix2y[ip_hi] + 32*self.pix2y[ip_med] + self.pix2y[ip_low]
		
	#	//c     transforms this in (horizontal, vertical) coordinates
		jrt = ix + iy#;//  ! 'vertical' in {0,2*(nside-1)}
		jpt = ix - iy#;//  ! 'horizontal' in {-nside+1,nside-1}
		
	#	//c     computes the z coordinate on the sphere
	#	//      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
		jr =  jrll[face_num]*nside - jrt - 1
	#	//      cout << "face_num=" << face_num << endl;
	#	//      cout << "jr = " << jr << endl;
	#	//      cout << "jrll(face_num)=" << jrll[face_num] << endl;
	#	//      cout << "----------------------------------------------------" << endl;
		nr = nside#;//                  ! equatorial region (the most frequent)
		z  = (2*nside-jr)*fact2
		kshift = int(fmod(jr - nside, 2))
		if jr<nside:#  { //then     ! north pole region
			nr = jr
			z = 1. - nr*nr*fact1
			kshift = 0
	
		else:# {
			if jr>3*nside:# {// then ! south pole region
				 nr = nl4 - jr
				 z = - 1. + nr*nr*fact1
				 kshift = 0
		theta = acos(z)
		
	#	//c     computes the phi coordinate on the sphere, in [0,2Pi]
	#	//      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
		jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2
		if jp>nl4:
			jp = jp - nl4
		if jp<1:
			jp = jp + nl4
		
		phi = (jp - (kshift+1)*0.5) * (piover2 / nr)
		return theta,phi

	def ring2nest(self,nside,p_ring): 
#	"""  /*
#		c=======================================================================
#		subroutine ring2nest(nside, ipring, ipnest)
#		c=======================================================================
#		c     conversion from RING to NESTED pixel number
#		c=======================================================================
#	  */
#	"""  
		ns_max=8192
	  
	#  static int x2pix[128], y2pix[128];
	#  //      common    /xy2pix/ x2pix,y2pix
	
		jrll = []
		jpll = []#;// ! coordinate of the lowest corner of each face
		for i in range(0,12):
			jrll.append(0)
			jpll.append(0)
		jrll[0]=2
		jrll[1]=2
		jrll[2]=2
		jrll[3]=2
		jrll[4]=3
		jrll[5]=3
		jrll[6]=3
		jrll[7]=3
		jrll[8]=4
		jrll[9]=4
		jrll[10]=4
		jrll[11]=4
		jpll[0]=1
		jpll[1]=3
		jpll[2]=5
		jpll[3]=7
		jpll[4]=0
		jpll[5]=2
		jpll[6]=4
		jpll[7]=6
		jpll[8]=1
		jpll[9]=3
		jpll[10]=5
		jpll[11]=7
	  
		npix = 12 * nside*nside
	#  if( ipring<0 || ipring>npix-1 ) {
	#    fprintf(stderr, "ipring out of range\n");
	#    exit(0);
	#  }
		if x2pix[127]<=0:
			self.mk_xy2pix()
	  
		nl2 = 2*nside
		nl4 = 4*nside
		npix = 12*nside*nside#;//      ! total number of points
		ncap = 2*nside*(nside-1)#;// ! points in each polar cap, =0 for nside =1
		ipring1 = p_ring + 1
	  
	#  //c     finds the ring number, the position of the ring and the face number
		if ipring1<=ncap: #//then
		
			hip   = ipring1/2.
			fihip = int(floor ( hip ))
			irn   = int(floor( sqrt( hip - sqrt(fihip) ) ) + 1)#;// ! counted from North pole
			iphi  = ipring1 - 2*irn*(irn - 1);
			
			kshift = 0
			nr = irn  # ;//               ! 1/4 of the number of points on the current ring
			face_num = (iphi-1) / irn#;// ! in {0,3}
	
		else:
			if ipring1<=nl2*(5*nside+1):# {//then
		
				ip    = ipring1 - ncap - 1
				irn   = int(floor( ip / nl4 ) + nside)#;//               ! counted from North pole
				iphi  = int(fmod(ip,nl4) + 1)
				
				kshift  = int(fmod(irn+nside,2))#;//  ! 1 if irn+nside is odd, 0 otherwise
				nr = nside
				ire =  irn - nside + 1#;// ! in {1, 2*nside +1}
				irm =  nl2 + 2 - ire
				ifm = (iphi - ire/2 + nside -1) / nside#;// ! face boundary
				ifp = (iphi - irm/2 + nside -1) / nside
				if ifp==ifm:# {//then          ! faces 4 to 7
					face_num = int(fmod(ifp,4) + 4)
				
				else:
					if ifp + 1==ifm:#  {//then ! (half-)faces 0 to 3
						face_num = ifp
				
					else:
						if ifp - 1==ifm:# {//then ! (half-)faces 8 to 11
							face_num = ifp + 7
	 
	
			else:
		
				ip    = npix - ipring1 + 1
				hip   = ip/2.
				fihip = floor ( hip )
				irs   = int(floor( sqrt( hip - sqrt(fihip) ) ) + 1)#;//  ! counted from South pole
				iphi  = 4*irs + 1 - (ip - 2*irs*(irs-1))
				
				kshift = 0
				nr = irs
				irn   = nl4 - irs
				face_num = (iphi-1) / irs + 8#;// ! in {8,11}
	  
	#  //c     finds the (x,y) on the face
	#  //  irt =   irn  - jrll[face_num+1]*nside + 1;//       ! in {-nside+1,0}
	#  //  ipt = 2*iphi - jpll[face_num+1]*nr - kshift - 1;// ! in {-nside+1,nside-1}
		irt =   irn  - jrll[face_num]*nside + 1#;//       ! in {-nside+1,0}
		ipt = 2*iphi - jpll[face_num]*nr - kshift - 1
	
	
		if ipt>=nl2:
			ipt = ipt - 8*nside#;// ! for the face #4
	  
		ix =  (ipt - irt ) / 2
		iy = -(ipt + irt ) / 2
		
		ix_low = int(fmod(ix,128))
		ix_hi  = ix/128
		iy_low = int(fmod(iy,128))
		iy_hi  = iy/128
	#  //  cout << "ix_low = " << ix_low << " ix_hi = " << ix_hi << endl;
	#  //  cout << "iy_low = " << iy_low << " iy_hi = " << iy_hi << endl;
	#  //  ipf =  (x2pix[ix_hi +1]+y2pix[iy_hi +1]) * (128 * 128)
	#  //    + (x2pix[ix_low+1]+y2pix[iy_low+1]);//        ! in {0, nside**2 - 1}
		ipf =  (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low])
	
	
	#  //  cout << "ipf = " << ipf << endl;
	#  //  for( int i(0);i<128;i++ ) cout << x2pix[i] << " || " << y2pix[i] << endl;
		return ipf + face_num* nside *nside#;//   ! in {0, 12*nside**2 - 1}
	  




