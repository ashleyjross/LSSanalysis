def xibaoSM(mn,c,rmin=50,rmax=200.,bs=5,md=1.,m=1.,mb='',v='n',p='',nbar=.002,fsky=.035,Bp=.4,mom='muw',nm='muwSA',covmd='s',N=63,mumin=0,mumax=1):
	from baofit_pubtest import doxi_isolike
	from Cosmo import distance
	d = distance(.31,.69)
	vol = fsky*4.*pi/3.*(d.dc(1.)**3.-d.dc(.6)**3.)
	muw = ''
	if mumin != 0:
		muw += 'mumin'+str(mumin)
	if mumax != 1.:
		muw += 'mumax'+str(mumax)	

	if nm == 'avemuwSM':
		dw = load('xiave'+mom+muw+'SM.dat').transpose()
	if nm == 'muwSA':	
		dw = load('ximuwSM/xi'+mom+muw+'SA'+str(mn)+'c'+c+'.dat').transpose()
	if nm == 'lampave':
		dw = load('xiSM'+nm+'mu'+str(mumin)+str(mumax)+str(bs)+'.dat').transpose()
	if nm == 'Y1':
		munm = ''
		if mumin != 0:
			munm += 'mum'+str(mumin)
		if mumax != 1:
			munm += 'mux'+str(mumax)
		dw = load('xi0gY1redBAO_mean_z_bpz_VFmz0.6xz1.0fcos'+munm+'5st0.dat').transpose()		
	dd = dw[1]
	rl = dw[0]
	mn = str(mn)
	if covmd == 'th':
		cov = load('covxiSM5.dat')
	if covmd == 's':
		cov = load('covxi'+mom+muw+'SM.dat')
	if covmd == 'lamp':
		cov = load('covSMlamp1.4.1'+'mu'+str(mumin)+str(mumax)+str(bs)+'.dat')	
	cov = cov*m	
	if md != 1:
		for i in range(0,len(cov)):
			cov[i][i] = cov[i][i]*md
	if p == 'dp':
		for i in range(0,len(cov)):
			np = nbar*vol*nbar*5.*4.*pi*rl[i]**2.
			cov[i][i] = cov[i][i]+2./np
				
	mdi = 1
	if mb == 'nobao':
		mdi = 2
	#mod = load('xizconv'+mom+'MICE_matterpower'+muw+'0.406.010.00.029sp1.0.dat').transpose()[mdi]
	mod = load('xizconvcMICE_matterpower'+muw+'0.4010.010.0combzsiglsp1.0.dat').transpose()[mdi]
	#mod = load('xizconvmuwPkMICE0.43.06.010.00.029sp1.0.dat').transpose()[mdi]
	Nmock = 8*N
	Nbin = (rmax-rmin)/5.
	chifac = (Nmock-Nbin-1)/(float(Nmock)-2.)
	chil = doxi_isolike(dd,cov,mod,rl,rmin=rmin,rmax=rmax,v=v,wo=mn+c+mb,Bp=Bp,chi2fac=chifac)
	fo = open('ximuwSM/BAOxichilSM'+nm+mn+c+mb+p+str(md)+str(m)+covmd+str(Bp)+muw+'.dat','w')
	for i in range(0,len(chil)):
		a = .8+.001*i+.0005
		fo.write(str(a)+' '+str(chil[i])+'\n')
	fo.close()
	from baofit import sigreg_c12
	a = sigreg_c12('ximuwSM/BAOxichilSM'+nm+mn+c+mb+p+str(md)+str(m)+covmd+str(Bp)+muw)
	#print a
	return a
