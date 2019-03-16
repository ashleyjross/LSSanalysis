#what was used to run DR12 2D fits
from baofit_pub2D import *
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
		print 'MISMATCHED data and cov matrix!'
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
	print len(dv)
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
