import numpy as np
from math import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from matplotlib import rcParams
from pylab import *
from numpy import loadtxt as load
import matplotlib.cm as cm

rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)



def BAOrelPlanckCurrent(wo='Current',xmax=2.5,Lya=True,BOSS=False,BOSSDR12=True,MGS=True,wz=True,sdss=False,df6=True,QSODR14=True,LRGDR14=True,des=False,desy1=True,eboss=False,desi=False):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.axes as ax
	import numpy as np
	from Cosmo import *
	from numpy import loadtxt as load

	pp = PdfPages('BAOrelPlanck'+wo+'.pdf')
	plt.clf()
	fig=plt.figure()
	ax0=fig.add_subplot(1,1,1)
	ax0.yaxis.set_ticks(np.arange(0.95,1.1,.05))
	plt.minorticks_on()
	ax.Axes.axis(ax0,'auto')
	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.ylim ( 0.911, 1.09 )
	#pe = np.loadtxt('/Users/ashleyross/DR7VAC/DVordPlanck.txt').transpose()
	pe = []
	f = open('Dvrserr.dat').readlines()
	fo = open('DVout.dat','w')
	zlf = []
	dvlf = []
	rmd = 147.78
	for i in range(0,len(f)):
		ln= f[i].split()
		pe.append(float(ln[2])/float(ln[1]))
		zlf.append(.05*(1.+i))
		dvlf.append(float(ln[1])*rmd)
	#x = pe[0]
	x=np.arange(0.05,3.,0.05)
	y=np.ones(len(x))
	yu = np.zeros(len(x))
	yd = np.zeros(len(x))
	yu2 = np.zeros(len(x))
	yd2 = np.zeros(len(x))
	omwig = .27
	hwig = .71
	obp = .02205
	obw = .0225
	omp = 0.315
	hp = 0.673
	hpu = 0.679
	delh = 0.012
	delm = .017
	hpd = 0.6665
 	df = distance(omp,1.-omp,hp,obhh=.022)
 	dh = distance(omp,1.-omp,hp+.01,obhh=.022)
 	dm = distance(omp+.01,1.-omp-.01,hp,obhh=.022)
 	wigf = (alph(.44,omp,hp,obp,omwig,hwig,obw)[0],alph(.6,omp,hp,obp,omwig,hwig,obw)[0],alph(.73,omp,hp,obp,omwig,hwig,obw)[0])
 	print wigf
 	#ru = du.rs/df.rs
 	#rd = dd.rs/df.rs
 	#re = (ru-rd)/2.
 	#fo = open('planckdvrserr_om.dat','w')
 	for i in range(0,len(x)):
 		#dvrsf = df.dV(x[i])/df.rs
 		#dd = delm*(dm.dV(x[i])/dm.rs/dvrsf-1.)/.01
		#fo.write(str(x[i])+' '+str(dd)+'\n')
		yu[i] = y[i]+pe[i]#pe[1][i]#+dd#sqrt(pe[1][i]**2.)#+re**2.)
		yd[i] = y[i]-pe[i]#pe[1][i]#dd#sqrt(pe[1][i]**2.)#+re**2.)
	plt.fill_between(x,yd,yu,color='0.75')
	#plt.fill_between(x,yd2,yu2,color='0.75')
	plt.plot(x,1.0*y,'k-')
	#plt.plot(x,1.0*yu,'k-')
	#plt.plot(x,yd,'k-')
	if df6:
		plt.text(0.06,.925,'6dFGS',fontsize=18,color='green')
		x6 = [0.1]
		y6 = [0.987]
		e6 = [0.045]
		plt.errorbar(x6, y6,e6,fmt='s',markersize=7,elinewidth=1.75,color='green',markeredgecolor='k')
		dv = np.interp(.1,zlf,dvlf)
		fo.write(str(0.1)+' '+str(dv*.987)+' '+str(dv*0.045)+'\n')
	if BOSS:
		plt.text(0.3,.945,'BOSS DR11',fontsize=18)
		xl = [0.32,0.57]
		yl = [0.974,0.983]
		el = [0.020,0.01]
		plt.errorbar(xl, yl,el,fmt='o',markersize=7,elinewidth=1.75,color='k')
	if BOSSDR12:
		plt.text(0.38,.96,'BOSS',fontsize=18,color='purple')
		xl = [0.35,0.61]
		yl = [0.999*.995,0.985*.996]
		el = [0.01,0.01]
		plt.errorbar(xl, yl,el,fmt='D',markersize=7,elinewidth=1.75,color='purple',markeredgecolor='k')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	#plt.text(0.26,.9,'LOWZ',fontsize=18)
	#plt.text(0.51,.92,'BOSS',fontsize=18)
	#plt.text(0.5,.9,'CMASS',fontsize=18)
	if MGS:
		plt.text(0.1,1.085,'SDSS',fontsize=18,color='r')
		plt.text(0.1,1.074,'MGS',fontsize=18,color='r')
		xl = [0.15]
		yl = [1.04*.9957]
		el = [0.037]
		plt.errorbar(xl,yl,el,fmt='^',markeredgecolor='k',markersize=7,elinewidth=1.75,color='r')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if wz:
		plt.text(0.45,1.072,'WiggleZ',fontsize=18,color='.5')
		ywl = [1.061*wigf[0],1.065*wigf[1],1.039*wigf[2]]
		xwl = [.44,.6,.73] 
		ewl = [0.048,0.045,0.034]
		plt.errorbar(xwl, ywl,ewl,fmt='v',markersize=6,elinewidth=1.,color='.5',markeredgecolor='k')
		for i in range(0,len(xwl)):
			dv = np.interp(xwl[i],zlf,dvlf)
			fo.write(str(xwl[i])+' '+str(dv*ywl[i])+' '+str(dv*ewl[i])+'\n')

	if sdss:		
		xsl = [0.275,0.35]
		ysl = [0.972,0.969]
		esl = [0.027,0.018]
		plt.errorbar(xsl, ysl,esl,fmt='s',markersize=6,elinewidth=1.,mfc='w',color='k')
	#plt.text(0.95,1.07,r'$\bar{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_avg)+'\pm'+'{:.3f}'.format(delta_alpha_std)+'$')
	#plt.text(0.95,1.04,r'$\tilde{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_med)+'^{'+'{:+.3f}'.format(delta_alpha_upp)+'}_{'+'{:+.3f}'.format(delta_alpha_low)+'}$')
	if eboss:
		plt.text(1.,.97,'eBOSS projected',fontsize=18,color='r')
		xl = [0.75,0.8,1.5]
		yl = [.99,.99,.99]
		el = [0.013,0.026,0.021]
		plt.errorbar(xl,yl,el,fmt='D',markeredgecolor='r',markersize=7,elinewidth=1.75,color='r')
	if QSODR14:
		plt.text(1.2,1.045,'eBOSS quasars',fontsize=18,color='orange')
		#plt.text(1.2,.94,'(predicted)',fontsize=18,color='b')
		xl = [1.5]
		yl = [0.996]
		el = [0.044]
		plt.errorbar(xl,yl,el,color='orange',markeredgecolor='k',markersize=8,elinewidth=2.25,markeredgewidth=0,linewidth=0)
		plt.plot(xl,yl,'*',color='orange',markeredgecolor='k',markersize=18,markeredgewidth=1.5)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if LRGDR14:
		plt.text(.68,.93,'eBOSS LRGs',fontsize=18,color='firebrick')
		xl = [.71]
		yl = [0.968]
		el = [0.026]
		plt.errorbar(xl,yl,el,fmt='x',markeredgecolor='k',markersize=7,elinewidth=1.75,color='firebrick')

	if des:
		plt.text(.85,.98,'DES Y3',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.8]
		yld = [1.0]
		eld = [0.025]
		#yu = [1.025,1.025]
		#yd = [0.975,0.975]
		plt.errorbar(xld,yld,eld,fmt='h',markeredgecolor='k',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
	if desy1:
		plt.text(.85,1.015,'DES',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.8]
		yld = [0.991]
		eld = [0.041]
		yu = [1.019,1.019]
		yd = [0.981,0.981]
		plt.errorbar(xld,yld,eld,fmt='-s',markeredgecolor='steelblue',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
		plt.ylim(.9,1.1)
	if desi:
		plt.text(1.4,1.03,'DESI',fontsize=18,color='purple')		
		d = load('/Users/ashleyross/BAOpredictionscode/BAO_BigBOSS.dat').transpose()
		xl = d[0]
		yl = np.ones((len(d[0])))
		el = d[-1]
		plt.errorbar(xl,yl,el,fmt='-*',markeredgecolor='purple',markersize=7,elinewidth=1.75,color='purple')
	if Lya:
 		lyf = [alph(2.33,omp,hp,obp,0.1426/.6731/.6731,.6731,0.02222)[0],alph(2.35,omp,hp,obp,.27,.7,.0227)[0]]
 		print lyf
 		lyb = [1.026/lyf[0],0.988/lyf[1]]
 		lyt = (lyb[0]/.026**2.+lyb[1]/.022**2.)/(1/.026**2.+1/.022**2.)
 		lyst = sqrt(1./(1/.026**2.+1/.022**2.))
 		print lyt,lyst
 		xl = [2.34]
 		yl = [lyt]
 		el = [lyst]
 		plt.errorbar(xl,yl,el,fmt='o',color='b',markersize=8,markeredgecolor='k')
		plt.text(1.85,.97,r'BOSS Ly$\alpha$',color='b',fontsize=18)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	plt.xlim ( 0.0, xmax )
	
	plt.xlabel ('Redshift', fontsize=18)
	#plt.ylabel (r'$(D_{\rm V}/r_{\rm d})/(D_{\rm V}/r_{\rm d})_{\rm Planck}$', fontsize=18)
	plt.ylabel (r'Distance/Distance(Planck$\Lambda$CDM)', fontsize=18)
	pp.savefig()
	pp.close()
	fo.close()
	return True

def BAOrelPlanckNear(wo='Near',xmax=2.5,Lya=True,BOSS=False,BOSSDR12e=True,MGS=True,wz=True,sdss=False,df6=True,QSO=True,QSODR14=False,ELG=True,LRG=True,LRGDR14=False,des=True,desy1=False,eboss=False,desi=False):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.axes as ax
	import numpy as np
	from Cosmo import *
	from numpy import loadtxt as load

	pp = PdfPages('BAOrelPlanck'+wo+'.pdf')
	plt.clf()
	fig=plt.figure()
	ax0=fig.add_subplot(1,1,1)
	ax0.yaxis.set_ticks(np.arange(0.95,1.1,.05))
	plt.minorticks_on()
	ax.Axes.axis(ax0,'auto')
	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.ylim ( 0.911, 1.09 )
	#pe = np.loadtxt('/Users/ashleyross/DR7VAC/DVordPlanck.txt').transpose()
	pe = []
	f = open('Dvrserr.dat').readlines()
	fo = open('DVout.dat','w')
	zlf = []
	dvlf = []
	rmd = 147.78
	for i in range(0,len(f)):
		ln= f[i].split()
		pe.append(float(ln[2])/float(ln[1]))
		zlf.append(.05*(1.+i))
		dvlf.append(float(ln[1])*rmd)
	#x = pe[0]
	x=np.arange(0.05,3.,0.05)
	y=np.ones(len(x))
	yu = np.zeros(len(x))
	yd = np.zeros(len(x))
	yu2 = np.zeros(len(x))
	yd2 = np.zeros(len(x))
	omwig = .27
	hwig = .71
	obp = .02205
	obw = .0225
	omp = 0.315
	hp = 0.673
	hpu = 0.679
	delh = 0.012
	delm = .017
	hpd = 0.6665
 	df = distance(omp,1.-omp,hp,obhh=.022)
 	dh = distance(omp,1.-omp,hp+.01,obhh=.022)
 	dm = distance(omp+.01,1.-omp-.01,hp,obhh=.022)
 	wigf = (alph(.44,omp,hp,obp,omwig,hwig,obw)[0],alph(.6,omp,hp,obp,omwig,hwig,obw)[0],alph(.73,omp,hp,obp,omwig,hwig,obw)[0])
 	print wigf
 	#ru = du.rs/df.rs
 	#rd = dd.rs/df.rs
 	#re = (ru-rd)/2.
 	#fo = open('planckdvrserr_om.dat','w')
 	for i in range(0,len(x)):
 		#dvrsf = df.dV(x[i])/df.rs
 		#dd = delm*(dm.dV(x[i])/dm.rs/dvrsf-1.)/.01
		#fo.write(str(x[i])+' '+str(dd)+'\n')
		yu[i] = y[i]+pe[i]#pe[1][i]#+dd#sqrt(pe[1][i]**2.)#+re**2.)
		yd[i] = y[i]-pe[i]#pe[1][i]#dd#sqrt(pe[1][i]**2.)#+re**2.)
	plt.fill_between(x,yd,yu,color='0.75')
	#plt.fill_between(x,yd2,yu2,color='0.75')
	plt.plot(x,1.0*y,'k-')
	#plt.plot(x,1.0*yu,'k-')
	#plt.plot(x,yd,'k-')
	if df6:
		plt.text(0.06,.925,'6dFGS',fontsize=18,color='green')
		x6 = [0.1]
		y6 = [0.987]
		e6 = [0.045]
		plt.errorbar(x6, y6,e6,fmt='s',markersize=7,elinewidth=1.75,color='green',markeredgecolor='k')
		dv = np.interp(.1,zlf,dvlf)
		fo.write(str(0.1)+' '+str(dv*.987)+' '+str(dv*0.045)+'\n')
	if BOSS:
		plt.text(0.3,.945,'BOSS DR11',fontsize=18)
		xl = [0.32,0.57]
		yl = [0.974,0.983]
		el = [0.020,0.01]
		plt.errorbar(xl, yl,el,fmt='o',markersize=7,elinewidth=1.75,color='k')
	if BOSSDR12e:
		plt.text(0.35,.967,'BOSS',fontsize=18,color='purple')
		xl = [0.35,0.52]
		yl = [0.999*.995,0.99*.996]
		el = [0.01,0.01]
		plt.errorbar(xl, yl,el,fmt='D',markersize=7,elinewidth=1.75,color='purple',markeredgecolor='k')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	#plt.text(0.26,.9,'LOWZ',fontsize=18)
	#plt.text(0.51,.92,'BOSS',fontsize=18)
	#plt.text(0.5,.9,'CMASS',fontsize=18)
	if MGS:
		plt.text(0.13,1.08,'SDSS',fontsize=18,color='r')
		plt.text(0.13,1.069,'MGS',fontsize=18,color='r')
		xl = [0.15]
		yl = [1.04*.9957]
		el = [0.037]
		plt.errorbar(xl,yl,el,fmt='^',markeredgecolor='k',markersize=7,elinewidth=1.75,color='r')
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if wz:
		plt.text(0.45,1.072,'WiggleZ',fontsize=18,color='.5')
		ywl = [1.061*wigf[0],1.065*wigf[1],1.039*wigf[2]]
		xwl = [.44,.6,.73] 
		ewl = [0.048,0.045,0.034]
		plt.errorbar(xwl, ywl,ewl,fmt='v',markersize=6,elinewidth=1.,color='.5',markeredgecolor='k')
		for i in range(0,len(xwl)):
			dv = np.interp(xwl[i],zlf,dvlf)
			fo.write(str(xwl[i])+' '+str(dv*ywl[i])+' '+str(dv*ewl[i])+'\n')

	if sdss:		
		xsl = [0.275,0.35]
		ysl = [0.972,0.969]
		esl = [0.027,0.018]
		plt.errorbar(xsl, ysl,esl,fmt='s',markersize=6,elinewidth=1.,mfc='w',color='k')
	#plt.text(0.95,1.07,r'$\bar{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_avg)+'\pm'+'{:.3f}'.format(delta_alpha_std)+'$')
	#plt.text(0.95,1.04,r'$\tilde{\Delta\alpha} ='+'{:+.3f}'.format(delta_alpha_med)+'^{'+'{:+.3f}'.format(delta_alpha_upp)+'}_{'+'{:+.3f}'.format(delta_alpha_low)+'}$')
	if eboss:
		plt.text(1.,.97,'eBOSS projected',fontsize=18,color='r')
		xl = [0.75,0.8,1.5]
		yl = [.99,.99,.99]
		el = [0.013,0.026,0.021]
		plt.errorbar(xl,yl,el,fmt='D',markeredgecolor='r',markersize=7,elinewidth=1.75,color='r')
	if QSODR14:
		plt.text(1.2,1.045,'eBOSS quasars',fontsize=18,color='orange')
		#plt.text(1.2,.94,'(predicted)',fontsize=18,color='b')
		xl = [1.5]
		yl = [0.996]
		el = [0.044]
		plt.errorbar(xl,yl,el,color='orange',markeredgecolor='k',markersize=8,elinewidth=2.25,markeredgewidth=0,linewidth=0)
		plt.plot(xl,yl,'*',color='orange',markeredgecolor='k',markersize=18,markeredgewidth=1.5)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')
	if QSO:
		plt.text(1.2,1.02,'eBOSS quasars',fontsize=18,color='orange')
		#plt.text(1.2,.94,'(predicted)',fontsize=18,color='b')
		xl = [1.5]
		yl = [0.996]
		el = [0.02]
		plt.errorbar(xl,yl,el,color='orange',markeredgecolor='k',markersize=8,elinewidth=2.25,markeredgewidth=0,linewidth=0)
		plt.plot(xl,yl,'*',color='orange',markeredgecolor='k',markersize=18,markeredgewidth=1.5)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	if LRGDR14:
		plt.text(.68,.93,'eBOSS LRGs',fontsize=18,color='firebrick')
		xl = [.71]
		yl = [0.968]
		el = [0.026]
		plt.errorbar(xl,yl,el,fmt='x',markeredgecolor='k',markersize=7,elinewidth=1.75,color='firebrick')

	if LRG:
		plt.text(.68,.95,'eBOSS+BOSS LRGs',fontsize=18,color='firebrick')
		xl = [.71]
		yl = [0.968]
		el = [0.01]
		plt.errorbar(xl,yl,el,fmt='x',markeredgecolor='k',markersize=7,elinewidth=1.75,color='firebrick')

	if ELG:
		plt.text(.68,1.03,'eBOSS ELGs',fontsize=18,color='k')
		xl = [.81]
		yl = [1.0]
		el = [0.03]
		plt.errorbar(xl,yl,el,fmt='o',markeredgecolor='k',markersize=7,elinewidth=1.75,color='k')
		plt.plot(xl,yl,'wo')

	if des:
		plt.text(.85,.98,'DES',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.85]
		yld = [0.991]
		eld = [0.02]
		#yu = [1.025,1.025]
		#yd = [0.975,0.975]
		plt.errorbar(xld,yld,eld,fmt='h',markeredgecolor='k',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
	if desy1:
		plt.text(.85,1.015,'DES',fontsize=18,color='steelblue')
		#xld = [0.7,0.9,1.1]
		#yld = [1.01,1.01,1.01]
		#eld = [0.031,0.026,0.034]
		xld = [0.8]
		yld = [0.991]
		eld = [0.041]
		yu = [1.019,1.019]
		yd = [0.981,0.981]
		plt.errorbar(xld,yld,eld,fmt='-s',markeredgecolor='steelblue',markersize=7,elinewidth=1.75,color='steelblue')
		#plt.fill_between(xld,yd,yu,color='b')
		plt.ylim(.9,1.1)
	if desi:
		plt.text(1.4,1.03,'DESI',fontsize=18,color='purple')		
		d = load('/Users/ashleyross/BAOpredictionscode/BAO_BigBOSS.dat').transpose()
		xl = d[0]
		yl = np.ones((len(d[0])))
		el = d[-1]
		plt.errorbar(xl,yl,el,fmt='-*',markeredgecolor='purple',markersize=7,elinewidth=1.75,color='purple')
	if Lya:
 		lyf = [alph(2.33,omp,hp,obp,0.1426/.6731/.6731,.6731,0.02222)[0],alph(2.35,omp,hp,obp,.27,.7,.0227)[0]]
 		print lyf
 		lyb = [1.026/lyf[0],0.988/lyf[1]]
 		lyt = (lyb[0]/.026**2.+lyb[1]/.022**2.)/(1/.026**2.+1/.022**2.)
 		lyst = sqrt(1./(1/.026**2.+1/.022**2.))
 		print lyt,lyst
 		xl = [2.34]
 		yl = [lyt]
 		el = [lyst]
 		plt.errorbar(xl,yl,el,fmt='o',color='b',markersize=8,markeredgecolor='k')
		plt.text(1.85,.97,r'BOSS Ly$\alpha$',color='b',fontsize=18)
		for i in range(0,len(xl)):
			dv = np.interp(xl[i],zlf,dvlf)
			fo.write(str(xl[i])+' '+str(dv*yl[i])+' '+str(dv*el[i])+'\n')

	plt.xlim ( 0.0, xmax )
	
	plt.xlabel ('Redshift', fontsize=18)
	#plt.ylabel (r'$(D_{\rm V}/r_{\rm d})/(D_{\rm V}/r_{\rm d})_{\rm Planck}$', fontsize=18)
	plt.ylabel (r'Distance/Distance(Planck$\Lambda$CDM)', fontsize=18)
	pp.savefig()
	pp.close()
	fo.close()
	return True

	
if __name__ == "__main__":
	 BAOrelPlanckNear()		
