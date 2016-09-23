import healpy as hp
from math import *
dir = '/cosma5/data/dp002/manera/fitting/testmaps/'

def mkgalodens(samp,zbin,res=256,ores=4096,tol=.2,space='rsd'):
	np = 12*res*res
	npo = 12*4096*4096
	frac = (res/float(ores))**2.
	gal = hp.read_map(dir+'map_'+space+'_nside4096_'+samp+zbin+'.fits')
	mask = hp.read_map(dir+'mask_'+space+'_nside4096_'+samp+'.fits')
	ml = []
	gl = []
	for i in range(0,np):
		ml.append(0)
		gl.append(0)
	for i in range(0,len(mask)):
		if mask[i] > 0:
			th,phi = hp.pix2ang(ores,i)
			pix = hp.ang2pix(res,th,phi)
			ml[pix] += frac*mask[i]
			gl[pix] += gal[i]

	ng = 0
	np = 0 
	for i in range(0,len(gl)):
		if ml[i] > tol:
			ng += gl[i]
			np += ml[i]
	print ng,np
	ave = ng/np
	print ave
	cw = ''
	fo = open('mock'+space+samp+zbin+cw+str(res)+'odenspczw.dat','w')
	ft = open('mock'+space+samp+zbin+cw+str(res)+'rdodens.dat','w')
	no = 0		
	for i in range(0,len(ml)):
		if ml[i] > tol:
			th,phi = hp.pix2ang(res,i)
			sra = sin(phi)
			cra = cos(phi)
			sdec = sin(-1.*(th-pi/2.))
			cdec = cos(-1.*(th-pi/2.))
			od = gl[i]/(ave*ml[i]) -1.
			fo.write(str(sra)+' '+str(cra)+' '+str(sdec)+' '+str(cdec)+' '+str(od)+' '+str(ml[i])+'\n')
			ft.write(str(th)+' '+str(phi)+' '+str(od)+' '+str(ml[i])+'\n')
	print no
	ft.close()
	fo.close()
	return True
	
def sqdind(mockf,dataf,ind=-8,zmin=0.6,dz=.05,zmax=1.0001,res=512,pz='mean_z_bpz',w='wst'):
	sqd = 0
	nz = int((zmax-zmin)/dz)
	for i in range(0,nz):
		zm = zmin+i*dz
		zx = zmin+(i+1.)*dz
		mf = numpy.loadtxt(mockf+'0'+str(i)+str(res)+'2ptPixclb.dat').transpose()
		df = numpy.loadtxt(dataff+str(zmin)+str(zmax)+pz+w+str(res)+'2ptPixclb.dat').transpose()
		sqd += (df[1][ind]-mf[1][ind])**2.
	return sqd
