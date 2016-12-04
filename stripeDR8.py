from math import *
deg2Rad = pi/180.0
strad2Deg = 360.0*360.0/(4.0*pi*pi)


def le2radec(lam, eta):
	from math import pi,cos,sin,atan2,asin
	deg2Rad = pi/180.0
	surveyCenterDEC = 32.5
	surveyCenterRA = 185.0
	etaPole = deg2Rad*surveyCenterDEC
	node = deg2Rad*(surveyCenterRA - 90.0)
	x = -1.0*sin(lam*deg2Rad)
	y = cos(lam*deg2Rad)*cos(eta*deg2Rad+etaPole)
	z = cos(lam*deg2Rad)*sin(eta*deg2Rad+etaPole)
	ra = (atan2(y,x) + node)/deg2Rad
	if ra < 0.0:
		ra +=360
	dec = asin(z)/deg2Rad
	return (ra,dec)


def radec2le(ra,dec):
	from math import pi,cos,sin,atan2,asin
	deg2Rad = pi/180.0
	surveyCenterDEC = 32.5
	surveyCenterRA = 185.0
	etaPole = deg2Rad*surveyCenterDEC
	node = deg2Rad*(surveyCenterRA - 90.0)
	
	x = cos(deg2Rad*ra-node)*cos(deg2Rad*dec)
  	y = sin(deg2Rad*ra-node)*cos(deg2Rad*dec)
  	z = sin(deg2Rad*dec)

  	lam = -1.0*asin(x)/deg2Rad
  	eta = (atan2(z,y) - etaPole)/deg2Rad
  	if eta < -180.0:
  		eta += 360.0
  	if eta > 180.0:
  		eta -= 360.0
  	
  	return (lam,eta)

def radec2le_nyu(ra,dec):
	from math import pi,cos,sin,atan2,asin
	deg2Rad = pi/180.0
	realRa0 = 275.
	realDec0 = 0.
	surveyRa0 = 57.5
	x = sin(realDec0*deg2Rad)*sin(dec*deg2Rad)
	y = cos(realDec0*deg2Rad)*cos(dec*deg2Rad)*cos((ra-realRa0)*deg2Rad)
	z = cos(dec*deg2Rad)*sin((ra-realRa0)*deg2Rad)
	w = cos(realDec0*deg2Rad)*sin(dec*deg2Rad)-sin(realDec0*deg2Rad)*cos(dec*deg2Rad)*cos((ra-realRa0)*deg2Rad)
	lam = asin(x+y)/deg2Rad
	eta = surveyRa0+atan2(z,w)/deg2Rad
	return lam,eta


def pix2etalam(resolution,pixnum):
	from math import pi,acos
	nx0 = 36
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	strad2Deg = 360.0*360.0/(4.0*pi*pi)
	etaOffSet = 91.25
	surveyCenterRA = 185.0
	surveyCenterDEC = 32.5
	node = deg2Rad*(surveyCenterRA - 90.0)
	etaPole = deg2Rad*surveyCenterDEC
	nx = nx0*resolution
	ny = ny0*resolution
	j = int(pixnum/nx)
	i = pixnum - nx*j

	eta = rad2Deg*(2.0*pi*(i+0.5))/nx + etaOffSet
	if eta > 180.0:
		eta += -360.0
	lam = 90.0 - rad2Deg*acos(1.0-2.0*(j+0.5)/ny)
	
	return(eta,lam)	

def pix2etalama(resolution,pixnum):
	from math import pi,acos
	nx0 = 36
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	strad2Deg = 360.0*360.0/(4.0*pi*pi)
	etaOffSet = 91.25
	surveyCenterRA = 185.0
	surveyCenterDEC = 32.5
	node = deg2Rad*(surveyCenterRA - 90.0)
	etaPole = deg2Rad*surveyCenterDEC
	nx = nx0*resolution
	ny = ny0*resolution
	i = int(pixnum/ny)
	j = pixnum - ny*i

	eta = rad2Deg*(2.0*pi*(i+0.5))/nx + etaOffSet
	if eta > 180.0:
		eta += -360.0
	lam = 90.0 - rad2Deg*acos(1.0-2.0*(j+0.5)/ny)
	
	return(eta,lam)	


def pix2ang_radec(resolution,pixnum):
	lameta = pix2etalam(resolution,pixnum)
	radec = le2radec(lameta[1],lameta[0])
	return (radec)

def ang2pix( resolution, lam, eta):
	#x is eta, y is lambda
	from math import pi,acos,cos
	nx0 = 36
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	strad2Deg = 360.0*360.0/(4.0*pi*pi)
	etaOffSet = 91.25
	surveyCenterRA = 185.0
	surveyCenterDEC = 32.5
	node = deg2Rad*(surveyCenterRA - 90.0)
	etaPole = deg2Rad*surveyCenterDEC
	
	nx = nx0*resolution
  	ny = ny0*resolution
  	
  	eta -= etaOffSet
  	eta *= deg2Rad
  	
  	if eta >= 0.0:
  		eta2 = eta
  	else:
  		eta2 = eta + 2.0*pi
  

  	i = int(nx*eta2/(2.0*pi))

  	lam = (90.0 - lam)*deg2Rad

  	if lam >= pi:
 		j = ny - 1
  	else:
  		j = int((ny*((1.0 - cos(lam)))/2.0))
  

  	pixnum = int(nx*j + i*1.)
  	return (pixnum)

def ang2pixa( resolution, lam, eta):
	#x is eta, y is lambda
	from math import pi,acos,cos
	nx0 = 36
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	strad2Deg = 360.0*360.0/(4.0*pi*pi)
	etaOffSet = 91.25
	surveyCenterRA = 185.0
	surveyCenterDEC = 32.5
	node = deg2Rad*(surveyCenterRA - 90.0)
	etaPole = deg2Rad*surveyCenterDEC
	
	nx = nx0*resolution
  	ny = ny0*resolution
  	
  	eta -= etaOffSet
  	eta *= deg2Rad
  	
  	if eta >= 0.0:
  		eta2 = eta
  	else:
  		eta2 = eta + 2.0*pi
  

  	i = int(nx*eta2/(2.0*pi))

  	lam = (90.0 - lam)*deg2Rad

  	if lam >= pi:
 		j = ny - 1
  	else:
  		j = int((ny*((1.0 - cos(lam)))/2.0))
  

  	pixnum = int(1.*j + i*ny)
  	return (pixnum)


def ang2pix2( resolution, lam, eta):
	from math import pi,acos,cos
	nx0 = .25
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	strad2Deg = 360.0*360.0/(4.0*pi*pi)
	etaOffSet = 91.25
	surveyCenterRA = 185.0
	surveyCenterDEC = 32.5
	node = deg2Rad*(surveyCenterRA - 90.0)
	etaPole = deg2Rad*surveyCenterDEC
	nx = nx0*resolution
  	ny = ny0*resolution
  	stripe = ang2stripe(eta)
  	etamin = eta_bound(stripe)[0]
  	i = int(nx*(eta-etamin)/2.5)
  	lamr = lamdic(stripe)
  	lammax = lamr[1]*deg2Rad
  	lammin = lamr[0]*deg2Rad
  	if lamr[1]==1000:
  		return (-1)
  	#if stripe == 12 and lam > .2 and lam < 2.3:
  	#	return(-1)
   	#if stripe == 14 and lam > 22.25 and lam < -20.1:
  	#	return(-1) 	
  	#if stripe == 36 and lam > -10 and lam < -8.5:
  	#	return(-1)
  	#if stripe == 36 and lam > 23.5 and lam < 24.75:
  	#	return(-1)
  	#if stripe == 37 and lam > 32.75 and lam < 34.25:
  	#	return(-1)
  	lam = (90.0 - lam)*deg2Rad

  	if lam >= pi:
 		j = ny - 1
  	else:
  		j = ((ny*((cos(pi/2-lammax) - cos(lam)))/2.0))
  	
  	h = int((ny*((cos(pi/2-lammax) - cos(pi/2-lammin)))/2.0)) - 1
  	
  	if int(j) > h or j < 0.0:
  		pixnum = -1
  	else:
  		pixnum = int(nx*int(j) + i*1.)
  	
  	if pixnum < 0:
  		pixnum = -1
  	
  	return (pixnum)


def pix2pixa(res,pixnum):
	nx = 36*res
	ny = 13*res
	j = int(pixnum/nx)
	i = (pixnum - j*nx)
	pixa = int(1.*j + i*ny)
	return pixa

def pix2ang2(resolution,stripe,pixnum):
	from math import pi,acos,cos
	nx0 = .25
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	nx = nx0*resolution
  	ny = ny0*resolution
  	etamin = eta_bound(stripe)[0]
  	j = int(pixnum/nx)
	i = pixnum - nx*j
  	eta = (2.5*(i+0.5))/nx + etamin
	if eta > 180.0:
		eta += -360.0
	lamr = lamdic(stripe)
  	lammax = lamr[1]*deg2Rad
	lam = 90 - rad2Deg*acos(1.0*cos(pi/2-lammax)-2.0*(j*1.+.5)/ny)
  	return(lam,eta)
  	
def pixup(res1,res2,stripe,pixnum):
	lameta = pix2ang2(res1,stripe,pixnum)
	lam = lameta[0]
	eta = lameta[1]
	pix2 = ang2pix2(res2,lam,eta)	
	return pix2

def pixupa(res1,res2,pixnum):
	lameta = pix2etalama(res1,pixnum)
	lam = lameta[1]
	eta = lameta[0]
	pix2 = ang2pixa(res2,lam,eta)	
	return pix2



def pixupS(res1,res2,pixnum):
	etalam = pix2etalam(res1,pixnum)
	lam = etalam[1]
	eta = etalam[0]
	pix2 = ang2pix(res2,lam,eta)
	return pix2

def pixdownS(res1,res2,pixnum):
	lameta = pix2etalam(res1,pixnum)
	pixb = pix_bound(res1,pixnum)
	lamr = abs(float(pixb[1])-float(pixb[0]))
	etar = abs(float(pixb[3])-float(pixb[2]))
	lam = lameta[1]
	eta = lameta[0]
	div = res2/res1
	pixl = []
	for i in range(0,pow(div,2)):
		row = i/div
		col = (float(i)/float(div)-row)*div
		lamW = lam +((row/float(div)-(.5-1./(2.*float(div)))))*lamr
		etaW = eta +((col/float(div)-(.5-1./(2.*float(div)))))*etar
		pix = ang2pix(res2,lam,eta)
		pixl.append(pix)
	return pixl

def pixdown(res1,res2,stripe,pixnum):
	lameta = pix2ang2(res1,stripe,pixnum)
	lam = lameta[0]
	eta = lameta[1]
	pix1 = ang2pix2(res2,lam-0.0002,eta-0.0002)
	pix2 = ang2pix2(res2,lam+0.0002,eta-0.0002)
	pix3 = ang2pix2(res2,lam-0.0002,eta+0.0002)
	pix4 = ang2pix2(res2,lam+0.0002,eta+0.0002)
	return (pix1,pix2,pix3,pix4)
	
def pixmax(resolution,stripe):
	from math import pi,acos,cos
	nx0 = .25
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	nx = nx0*resolution
  	ny = ny0*resolution
  	lamr = lamdic(stripe)
  	lammax = lamr[1]*deg2Rad
  	lammin = lamr[0]*deg2Rad
  	h = int((ny*((cos(pi/2-lammax) - cos(pi/2-lammin)))/2.0)) - 1
  	pixm = int(nx*h+nx-1)
  	return (pixm)

def ang2pix_radeca(resolution,ra,dec):
	lameta = radec2le(ra,dec)
	lam = lameta[0]
	eta = lameta[1]
	
	pixnum = ang2pixa(resolution, lam, eta)
	
	return (pixnum)

def pix2gal(res,stripe,pix):
	import cocoJ2gal
	lameta = pix2ang2(res,stripe,pix)
	radec = le2radec(lameta[0],lameta[1])
	lb = cocoJ2gal.cocoJ2gal(radec[0],radec[1])
	return lb


def pix2galALL(res,pix):
	import cocoJ2gal
	radec = pix2ang_radec(res,pix)
	lb = cocoJ2gal.cocoJ2gal(radec[0],radec[1])
	return lb


def ang2stripe(eta):
	inc = eta + 32.5
	stripe = int(inc/2.5 + 10)
	a = eta_bound(stripe)
	if eta < a[0] or eta > a[1]:
		stripe = stripe + 1
	return stripe

def ang2perpstripe(lam):
	pix = ang2pixa(64,lam,0)
	st = (pix-1431095)/19
	if st < 0:
		return 0
	if st > 37:
		return 37
	return st
# 	if lam > 55:
# 		return 0
# 	if lam > 50.6:
# 		return 1
# 	if lam > 46.65
# 		return 2
# 	if lam > 42.96:
# 		return 3
# 	if lam > 39.48:
# 		return 4
# 	if lam > pix2etalama(64,1431115+19*5):
# 		return 5
# 	if lam > pix2etalama(64,1431115+19*6):
# 		return 6
# 	if lam > pix2etalama(64,1431115+19*7):
# 		return 7
# 	if lam > pix2etalama(64,1431115+19*8):
# 		return 8
# 	if lam > pix2etalama(64,1431115+19*9):
# 		return 9
# 	if lam > pix2etalama(64,1431115+19*10):
# 		return 10
	
def pix_bound(resolution, pixnum):
	from math import pi,acos,sin,atan2,asin
	nx0 = 36
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	surveyCenterDEC = 32.5
	surveyCenterRA = 185.0
	etaOffSet = 91.25
	etaPole = deg2Rad*surveyCenterDEC
	node = deg2Rad*(surveyCenterRA - 90.0)
	nx = nx0*resolution
	ny = ny0*resolution
	j = pixnum/nx
	i = pixnum - nx*j
	lammin = 90.0 - rad2Deg*acos(1.0 - 2.0*(j+1)/ny)
	lammax = 90.0 - rad2Deg*acos(1.0 - 2.0*j/ny)
	etamin = rad2Deg*2.0*pi*(i+0.0)/nx + etaOffSet
	if etamin >= 180.0:
		etamin = etamin - 360.0
	etamax = rad2Deg*2.0*pi*(i+1.0)/nx + etaOffSet
	if etamax >= 180.0:
		etamax = etamax - 360.0
	return((lammin,lammax,etamin,etamax))

def pix_bounda(resolution, pixnum):
	from math import pi,acos,sin,atan2,asin
	nx0 = 36
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	surveyCenterDEC = 32.5
	surveyCenterRA = 185.0
	etaOffSet = 91.25
	etaPole = deg2Rad*surveyCenterDEC
	node = deg2Rad*(surveyCenterRA - 90.0)
	nx = nx0*resolution
	ny = ny0*resolution
	i = int(pixnum/ny)
	j = pixnum - ny*i
	lammin = 90.0 - rad2Deg*acos(1.0 - 2.0*(j+1)/ny)
	lammax = 90.0 - rad2Deg*acos(1.0 - 2.0*j/ny)
	etamin = rad2Deg*2.0*pi*(i+0.0)/nx + etaOffSet
	if etamin >= 180.0:
		etamin = etamin - 360.0
	etamax = rad2Deg*2.0*pi*(i+1.0)/nx + etaOffSet
	if etamax >= 180.0:
		etamax = etamax - 360.0
	return((lammin,lammax,etamin,etamax))


def pix_bound2(resolution,stripe,pixnum):
	from math import pi,acos,cos
	nx0 = .25
	ny0 = 13
	deg2Rad = pi/180.0
	rad2Deg = 180.0/pi
	nx = nx0*resolution
  	ny = ny0*resolution
  	etamin = eta_bound(stripe)[0]
  	j = int(pixnum/nx)
	i = pixnum - nx*j
  	eta1 = (2.5*(i))/nx + etamin
	if eta1 > 180.0:
		eta1 += -360.0
	eta2 = (2.5*(i+1))/nx + etamin
	if eta2 > 180.0:
		eta2 += -360.0
	lamr = lamdic(stripe)
  	lammax = lamr[1]*deg2Rad
	lam1 = 90 - rad2Deg*acos(1.0*cos(pi/2-lammax)-2.0*(j*1.)/ny)	
	lam2 = 90 - rad2Deg*acos(1.0*cos(pi/2-lammax)-2.0*(j*1.+1)/ny)
	return((lam2,lam1,eta1,eta2))

def pix_area2(resolution,stripe,pixnum):	
	from math import sin, pi
	deg2Rad = pi/180.0
	strad2Deg = 360.0*360.0/(4.0*pi*pi)
	
	a = pix_bound2(resolution,stripe,pixnum)
	b = strad2Deg*(deg2Rad*(a[3]-a[2]))*(sin(deg2Rad*a[1])-sin(deg2Rad*a[0]))
	if b < 0:
		b = 0
	return b
  
def eta_bound(stripe):
	inc = stripe_inclination(stripe)
	etamin = inc - 32.5 - 1.25 + 0.0000001
	etamax = inc - 32.5 + 1.25 - 0.0000001
	if etamin > 180:
		etamin = etamin - 360.
	if etamax > 180:
		etamax = etamax - 360.
	return(etamin,etamax)


def stripe_inclination( stripe): 
	return 2.5*(stripe-10)


def pix_area( resolution,  pixnum):
	from math import sin, pi
	deg2Rad = pi/180.0
	strad2Deg = 360.0*360.0/(4.0*pi*pi)
	
	a = pix_bounda(resolution,pixnum)
	b = strad2Deg*(deg2Rad*(a[3]-a[2]))*(sin(deg2Rad*a[1])-sin(deg2Rad*a[0]))
	if b < 0:
		b = 0
	return b

def pixdis(res,pix1,pix2):
	#No guarentee this correct, check this at some point
	etalam1 = pix2etalama(res,pix1)
	etalam2 = pix2etalama(res,pix2)
	etadif = deg2Rad*(etalam1[0]-etalam2[0])
	lamdif = sin(deg2Rad*etalam1[1])-sin(deg2Rad*etalam2[1])
	#print lamdif/deg2Rad
	dif = sqrt(etadif**2.+lamdif**2)/deg2Rad
	return dif


def pixboundtest(res,pixnum):
	a = pix_bound(res,pixnum)
	b = pix2etalam(res,pixnum)
	stripe = ang2stripe(b[0])
	stripe2 = ang2stripe(a[2])
	stripe3 = ang2stripe(a[3])
	etab = eta_bound(stripe)
	lr = lrange(stripe)
	lr1 = lrange(stripe2)
	lr2 = lrange(stripe3)
	pixb = pix_bound(res,pixnum)
	if lr[0] == -1000 or lr1[0] == -1000 or lr2[0] == -1000 or a[0] < lr[0] or a[1] > lr[1]:
		return(-1)
	else:
		return(0)

def lamdic(stripe):
	lammin = -1000
	lammax = 1000
	lamr = { 10 : (-63,64.5), 9 : (-55,53.5), 11 : (-59.5,51.5), 12 : (-64.1,56), 13 : (-62.15,57.8),
	14 : (-57.143, 58.9), 15 : (-64.95, 59.8), 16 : (-65,62), 17: ( -65.5,61.5), 18: (-65.,62.8), 19: (-66,59.5), 20: (-65.,63.8), 21: (-65,61.8), 22: (-65,64), 23: (-66,60.), 24: (-65,59), 25: (-64.5,64.), 26: (-65,64.5), 27: (-64.75,63.7), 28: (-65.55,63.6), 29: (-62.8,63.4), 30: (-63.0,63.1), 31: (-60.4,62.8), 32: (-56.64,63.35), 33: (-60.0,53.94), 34: (-59.25,61.8), 35: (-53.5,60.95), 36: (-48.281,61.48), 37: (-52.3,48.312), 38: (-51,50.), 39: (-49.,57.), 68: (-25, 37.), 69: (-25, 37.5), 70: (-24.5, 36), 71: (-25.5, 36.5), 72: (-55., 55.), 73: (-26., 40.), 74: (-28., 38.5), 75: (-27.5, 37.5), 76: (-27.95, 56.), 77: (-26., 49.), 78: (-25.5, 49.), 79: (-52, 56.), 80: (-40.5, 54.), 81: (-43.5, 52.), 82: (-70.0, 60.0), 83: (-40.5, 53.0), 84: (-40.5, 49.0), 85: (-40.5, 47.0), 86: (-61.8,55.7)
    	 }
	if stripe >= 9 and stripe <= 39:
		return lamr[stripe]
	if stripe >= 68 and stripe <=86:
		return lamr[stripe]
	
	return (lammin,lammax)

def pixboundtest2(res,pixnum):
	a = 0
	while a != -1:
		b = pix2etalam(res,pixnum)
		c = pix_bound(res,pixnum)
		
		if b[1] > 65 or b[1] < -65:
			a = -1
			#print ('+')
			return (a)
			#print ('+')
			break
		stripe1 = ang2stripe(c[2])
		#print (stripe1)
		if lamdic(stripe1)[0] == -1000:
			a = -1
			#print ('++')
			return (a)
			#print ('+')
			break
		stripe2 = ang2stripe(c[3])
		#print (stripe2)
		if lamdic(stripe2)[0] == -1000:
			a = -1
			#print ('+++')
			return (a)
			#print ('+')
			break
		if c[1] > lamdic(stripe1)[1] or c[1] > lamdic(stripe2)[1] or c[0] < lamdic(stripe2)[0] or c[0] <  lamdic(stripe1)[0]:
			a = -1
			#print('++++')
			return (a)
			#print ('+')
			break
		a = 0
		return (a)
		break
		


def pixfind(angdis,res,pixnum):
#find the pixels around a pixel that are a given angular distance away.  To start, simply returns 10 pixels
	from math import pi, cos, sin
	b = 0
	try:
		import biggles
	except:
		b = 1
	
	(eta,lam) = pix2etalam(res,pixnum)
	out = []
	etal = []
	laml = []
	count = 0
	r = int(160.*angdis*res/256.)
	for i in range(0,r/2):
		ang = i*2.*pi/float(r)
		eta1 = eta +angdis*cos(ang)
		etal.append(eta1)
		lam1 = lam + angdis*sin(ang)
		laml.append(lam1)
		pix = ang2pix(res,lam1,eta1)
		if count == 0:
			out.append(pix)
			count += 1
		else:
			if pix != out[count-1]:
				out.append(pix)
				count += 1
			else:
				print('too many divisions')
				print(pix)
	peta = []
	plam = []
	for i in range(0,len(out)):
		(eta1,lam1) = pix2etalam(res,out[i])
		peta.append(eta1)
		plam.append(lam1)
	plot = biggles.FramedPlot()
	crv1 = biggles.Curve(laml, etal)
	crv2 = biggles.Curve(plam,peta)		
	plot.add(crv1)
	plot.add(crv2)
	
	plot.show()
	return out


def pixfinda(angdis,res,pixnum):
#find the pixels around a pixel that are a given angular distance away.  To start, simply returns 10 pixels
	from math import pi, cos, sin
	b = 0
	try:
		import biggles
	except:
		b = 1
	
	(eta,lam) = pix2etalama(res,pixnum)
	out = []
	etal = []
	laml = []
	count = 0
	r = int(160.*angdis*res/256.)
	for i in range(0,r/2):
		ang = i*2.*pi/float(r)
		eta1 = eta +angdis*cos(ang)
		etal.append(eta1)
		lam1 = lam + angdis*sin(ang)
		laml.append(lam1)
		pix = ang2pixa(res,lam1,eta1)
		if count == 0:
			out.append(pix)
			count += 1
		else:
			if pix != out[count-1]:
				out.append(pix)
				count += 1
			else:
				print('too many divisions')
				print(pix)
	
# 	peta = []
# 	plam = []
# 	for i in range(0,len(out)):
# 		(eta1,lam1) = pix2etalam(res,out[i])
# 		peta.append(eta1)
# 		plam.append(lam1)
# 	plot = biggles.FramedPlot()
# 	crv1 = biggles.Curve(laml, etal)
# 	crv2 = biggles.Curve(plam,peta)		
# 	plot.add(crv1)
# 	plot.add(crv2)
# 	
# 	plot.show()
	return out



def lrange(stripe):
	lammin = -1000
	lammax = 1000
	
	if stripe == 10:
		lammin = -59.6
		lammax = 53.6
	
	if stripe == 9:
		lammin = -14.2
    	lammax = 16.22
  
	if stripe == 11:
		lammin = -56.6
    	lammax = 43.5
	
	if stripe == 12:
		lammin = -61.2
    	lammax = 56.0
  
	if stripe == 13:
		lammin = -62.15
    	lammax = 57.8
  
	if stripe == 14:
		lammin = -62.4
    	lammax = 58.9
  	
  	if stripe == 15:
  		lammin = -64.95
    	lammax = 59.8
  
  	if stripe == 16:
  		lammin = -65.4
    	lammax = 62.0
  
  	if stripe == 17:
  		lammin = -63.4
    	lammax = 61.2
  
  	if stripe == 18:
  		lammin = -63.6
    	lammax = 61.8
  
  	if stripe == 19:
  		lammin = -63.7
    	lammax = 62.3
  	
  	if stripe == 20:
  		lammin = -63.8
    	lammax = 62.7
  
  	if stripe == 21:
  		lammin = -63.7
    	lammax = 63.1
  
  	if stripe == 22:
  		lammin = -63.7
    	lammax = 63.3
  
  	if stripe == 23:
  		lammin = -63.5
    	lammax = 63.5
  
  	if stripe == 24:
  		lammin = -63.3
    	lammax = 63.7
  
  	if stripe == 25:
  		lammin = -63.1
    	lammax = 63.7
  
  	if stripe == 26:
  		lammin = -62.7
    	lammax = 63.8
  
  	if stripe == 27:
  		lammin = -64.75
    	lammax = 63.7
  
  	if stripe == 28:
  		lammin = -65.55
    	lammax = 63.6
  
  	if stripe == 29:
  		lammin = -62.8
    	lammax = 63.4
    
    
    	
  	if stripe == 30:
  		lammin = -63.0
    	lammax = 63.1
  
  	if stripe == 31:
  		lammin = -60.4
    	lammax = 62.8
  
  	if stripe == 32:
  		lammin = -60.0
    	lammax = 63.35
  
  	if stripe == 33:
  		lammin = -60.0
    	lammax = 61.9
  
  	if stripe == 34:
  		lammin = -59.25
    	lammax = 61.8
  
  	if stripe == 35:
  		lammin = -53.89
    	lammax = 60.95
  
  	if stripe == 36:
  		lammin = -48.28
    	lammax = 59.6
  
  	if stripe == 37:
  		lammin = -51.7
    	lammax = 48.5
  
  	if stripe == 76:
  		lammin = -27.95
    	lammax = -1.95
  
  	if stripe == 82:
  		lammin = -54.7
    	lammax = 19.6
  
  	if stripe == 86:
  		lammin = -57.8
    	lammax = 55.7
    	
    
  
	return (lammin,lammax)
	
def pixlist3(pixnum,res1,res2):
#In new scheme designed for larger angles, find the list of base pixels from original scheme
#that go into large pixel in new scheme
	nx = 36*res1
	ny = 13*res1
	pixels = res1/res2
	pix2 = pixels*pixels
	pixl = []
	for i in range(0,pixels):
		for j in range(0,pixels):
			pixl.append(pixnum+i+nx*j)
	return pixl
	
def pixlistNEW(pixnum,res1,res2):
	nx = 36*res1
	ny = 13*res1
	pixels = (float(res1)/float(res2))/2.
	pixs = []
	pixnum = int(pixnum)
	if res2 == res1/3.:
		return(pixnum,pixnum+1,pixnum+2,pixnum+nx,pixnum+nx+1,pixnum+nx+2, pixnum+2*nx,pixnum+2*nx+1,pixnum+2*nx+2)
	if res2 == res1/5.:
		return(pixnum,pixnum+1,pixnum+2,pixnum+3,pixnum+4,pixnum+nx,pixnum+nx+1,pixnum+nx+2,pixnum+nx+3, pixnum+nx+4,pixnum+2*nx,pixnum+2*nx+1,pixnum+2*nx+2,pixnum+2*nx+3,pixnum+2*nx+4,pixnum+3*nx,pixnum+3*nx+1,pixnum+3*nx+2,pixnum+3*nx+3,pixnum+3*nx+4,pixnum+4*nx,pixnum+4*nx+1,pixnum+4*nx+2,pixnum+4*nx+3,pixnum+4*nx+4)
	if res2 == res1/7. or res2 == res1/11. or res2 == res1/13. or res2 == res1/15. or res2 == res1/17. or res2 == res1/9.:
		ind = int(res1/res2)
		for i in range(0,ind):
			for j in range(0,ind):
				pixs.append(pixnum+j+i*nx)
		return(pixs)
	#print(pixels)
	pix1 = pixnum
	pix2 = pixnum + pixels
	pix3 = pixnum +pixels*nx
	pix4 = pixnum + pixels*nx + pixels
	return (int(pix1), int(pix2), int(pix3), int(pix4))
	
def pixlistNEWa(pixnum,res1,res2):
	#nx is replaced by ny since this is designed to work with ang2pixa
	nx = 36*res1
	ny = 13*res1
	pixels = (float(res1)/float(res2))/2.
	pixs = []
	pixnum = int(pixnum)
	if res2 == res1/3.:
		return(pixnum,pixnum+1,pixnum+2,pixnum+ny,pixnum+ny+1,pixnum+ny+2, pixnum+2*ny,pixnum+2*ny+1,pixnum+2*ny+2)
	if res2 == res1/5.:
		return(pixnum,pixnum+1,pixnum+2,pixnum+3,pixnum+4,pixnum+ny,pixnum+ny+1,pixnum+ny+2,pixnum+ny+3, pixnum+ny+4,pixnum+2*ny,pixnum+2*ny+1,pixnum+2*ny+2,pixnum+2*ny+3,pixnum+2*ny+4,pixnum+3*ny,pixnum+3*ny+1,pixnum+3*ny+2,pixnum+3*ny+3,pixnum+3*ny+4,pixnum+4*ny,pixnum+4*ny+1,pixnum+4*ny+2,pixnum+4*ny+3,pixnum+4*ny+4)
	if res2 == res1/7. or res2 == res1/11. or res2 == res1/13. or res2 == res1/15. or res2 == res1/17. or res2 == res1/9.:
		ind = int(res1/res2)
		for i in range(0,ind):
			for j in range(0,ind):
				pixs.append(pixnum+j+i*ny)
		return(pixs)
	#print(pixels)
	pix1 = pixnum
	pix2 = pixnum + pixels
	pix3 = pixnum +pixels*ny
	pix4 = pixnum + pixels*ny + pixels
	return (int(pix1), int(pix2), int(pix3), int(pix4))






def pixlistNEW2(pixnum,res1,res2,stripe):
	
	nx = 0.25*res1
	ny = (pixmax(res1,stripe)+1)/nx
	#print(nx)
	pixels = (float(res1)/float(res2))/2.
	pixs = []
	pixnum = int(pixnum)
	row = pixnum/nx
	pos = (row+1.)*nx-pixnum
	print(pos)
	if pixels == 1.5:
		if pos >= 2:
			#print(pix2ang2(res1,stripe,pixnum),pix2ang2(res1,stripe,pixnum+1),pix2ang2(res1,stripe,pixnum+2),pix2ang2(res1,stripe,pixnum+nx),pix2ang2(res1,stripe,pixnum+nx+1),pix2ang2(res1,stripe,pixnum+nx+2), pix2ang2(res1,stripe,pixnum+2*nx),pix2ang2(res1,stripe,pixnum+2*nx+1),pix2ang2(res1,stripe,pixnum+2*nx+2))

			return(pixnum,pixnum+1,pixnum+2,pixnum+nx,pixnum+nx+1,pixnum+nx+2, pixnum+2*nx,pixnum+2*nx+1,pixnum+2*nx+2)
		else:
			return(-1)
	#print(pixels)
	pix1 = pixnum
	pix2 = pixnum + pixels
	pix3 = pixnum +pixels*nx
	pix4 = pixnum + pixels*nx + pixels
	if pos >= 2*pixels:
		return (int(pix1), int(pix2), int(pix3), int(pix4))
	else:
		return(-1)

def radec2xyz(ra,dec):
	from math import pi, sin, cos
	theta = (ra)*pi/180.
	phi = (90-dec)*pi/180.
	x = sin(phi)*cos(theta)
	y = sin(phi)*sin(theta)
	z = cos(phi)
	#print(x,y,z)
	return(x,y,z)


def VolPix(scale,ra,dec,Z):
	from math import pow
	from cosmo import *
	(x,y,z) = radec2xyz(ra,dec)
	dd = distance()
	sideDIM = int(2312.66/scale)
	pixm = int(pow(sideDIM,3))
	r = dd.dc(Z)
	x = r*x
	y = r*y
	z = r*z
	pix = pow(sideDIM,2.)*int(z/scale)+sideDIM*int(y/scale)+int(x/scale)
	return(pix)


if __name__=='__main__':
	#import psycotwofour
	#psycotwofour.full()
	f = open('pixdistemp.dat','w')
	start = 600000
	for i in range(0,10):
		for j in range(start,10+start):
			pix = 104*8*i + j
			dis = pixdis(64,start,pix)
			if j == 10+start-1:
				f.write(str(dis)+'\n')
			else:
				f.write(str(dis)+' ')
	f.close()
	f = open('pixdistemp.dat','r').readlines()
	#print(f)
	for i in range(0,len(f)):
		line = f[i].split()
		for j in range(0,len(line)):
			if float(line[j]) > .18*4. and float(line[j]) < .22*4.:
				print(i,j)
	

# 	from time import time
# 	t1 = time()
# 	
# 	for i in range(0,1000000):
# 		ra = 20
# 		dec = -30
# 		Z =.3
# 		pix = VolPix(1,ra,dec,Z)
# 	t2 = time()
# 	print(t2-t1)


# 	f = open('stripe10radec16.dat','w')
# 	stripe = 10
# 	res = 16
# 	pm = pixmax(res,stripe)
# 	for i in range(0,pm):
# 		le = pix_bound2(res,stripe,i)
# 		radec1 = le2radec(le[0],le[2])
# 		radec2 = le2radec(le[0],le[3])
# 		radec3 = le2radec(le[1],le[2])
# 		radec4 = le2radec(le[1],le[3])
# 		f.write(str(radec1[0])+' '+str(radec1[1])+' '+str(radec2[0])+' '+str(radec2[1])+' '+str(radec3[0])+' '+str(radec3[1])+' ' +str(radec3[0])+' '+str(radec3[1])+'\n')
# 	f = open('stripe11radec16.dat','w')
# 	stripe = 11
# 	res = 16
# 	pm = pixmax(res,stripe)
# 	for i in range(0,pm):
# 		le = pix_bound2(res,stripe,i)
# 		radec1 = le2radec(le[0],le[2])
# 		radec2 = le2radec(le[0],le[3])
# 		radec3 = le2radec(le[1],le[2])
# 		radec4 = le2radec(le[1],le[3])
# 		f.write(str(radec1[0])+' '+str(radec1[1])+' '+str(radec2[0])+' '+str(radec2[1])+' '+str(radec3[0])+' '+str(radec3[1])+' ' +str(radec3[0])+' '+str(radec3[1])+'\n')
# 	f = open('stripe11radec8.dat','w')
# 	stripe = 11
# 	res = 8
# 	pm = pixmax(res,stripe)
# 	for i in range(0,pm):
# 		le = pix_bound2(res,stripe,i)
# 		radec1 = le2radec(le[0],le[2])
# 		radec2 = le2radec(le[0],le[3])
# 		radec3 = le2radec(le[1],le[2])
# 		radec4 = le2radec(le[1],le[3])
# 		f.write(str(radec1[0])+' '+str(radec1[1])+' '+str(radec2[0])+' '+str(radec2[1])+' '+str(radec3[0])+' '+str(radec3[1])+' ' +str(radec3[0])+' '+str(radec3[1])+'\n')
# 	
# 	f = open('stripe10radec8.dat','w')
# 	
# 	stripe = 10
# 	res = 8
# 	pm = pixmax(res,stripe)
# 	for i in range(0,pm):
# 		le = pix_bound2(res,stripe,i)
# 		radec1 = le2radec(le[0],le[2])
# 		radec2 = le2radec(le[0],le[3])
# 		radec3 = le2radec(le[1],le[2])
# 		radec4 = le2radec(le[1],le[3])
# 		f.write(str(radec1[0])+' '+str(radec1[1])+' '+str(radec2[0])+' '+str(radec2[1])+' '+str(radec3[0])+' '+str(radec3[1])+' ' +str(radec3[0])+' '+str(radec3[1])+'\n')
# 	

