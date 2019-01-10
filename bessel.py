from math import *
ACC = 100.
BIGNO = 1.e10
BIGNI = 1.e-10


def bessj0(x):
	ax = abs(x)
	if ax < 8.:
		y = x*x
		ans1 = 57568490574.+y*(-13362590354.+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))))
		ans2 = 57568490411.+y*(1029532985.+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))))
		ans = ans1/ans2
	else:
		z = 8./ax
		y = z*z
		xx = ax -.785398164
		ans1 = 1.0+y*(-.1098628627e-2+y*(.2734510407e-4+y*(-.2073370639e-5+y*0.2093887211e-6)))
		ans2 = -.1562499995e-1+y*(.1430488765e-3+y*(-.6911147651e-5+y*(.7621095161e-6-y*0.934935152e-7)))
		ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2)
	return ans
	
def bessj1(x):
	ax = abs(x)
	if ax < 8.:
		y = x*x
		ans1 = x*(72362614232.+y*(-7895059235.+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))))
		ans2 = 144725228442.+y*(2300535178.+y*(18583304.43394+y*(99447.43394+y*(376.9991397+y*1.0))))
		ans = ans1/ans2
	else:
		z = 8./ax
		y = z*z
		xx = ax -2.356194491
		ans1 = 1.0+y*(.183105e-2+y*(-.3516396496e-4+y*(.2457520174e-5+y*-0.240337019e-6)))
		ans2 = .04687499995+y*(-.2002690873e-3+y*(.8449199096e-5+y*(-.88228987e-6+y*0.105787412e-6)))
		ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2)
	return ans


def bessj(n,x):
	if n == 0:
		return bessj0(x)
	if n == 1:
		return bessj1(x)
		
	ax = fabs(x)
	if ax == 0.:
		return 0.
	else:
		if ax > float(n):
			tox = 2./ax
			bjm = bessj0(ax)
			bj = bessj1(ax)
			for j in range(1,n):
				bjp=j*tox*bj-bjm
				bjm=bj
				bj=bjp
			ans = bj
		
		else:
			tox = 2./ax
			m = 2*((n+int(sqrt(ACC*n)))/2)
			jsum = 0
			bjp=ans=sum=0.
			bj=1.
			for i in range(0,m):
				j = m-i
				bjm =j*tox*bj-bjp
				bjp = bj
				bj = bjm
				if fabs(bj) > BIGNO:
					bj *= BIGNI
					bjp *= BIGNI
					ans *= BIGNI
					sum *= BIGNI
				if jsum == 1:
					sum += bj
				if jsum == 1:
					jsum = 0
				else:
					jsum = 1
				if j == n:
					ans=bjp
					
			sum = 2.*sum-bj
			ans /= sum
	
	return ans

def gamdic(k):
	#dictionary lookup for integer/half integer gamma function arguments from -3/2 to 4
	lamr = { -1.5 : 4.*sqrt(pi)/3., -.5 : -2.*sqrt(pi), .5 : sqrt(pi), 1. : 1., 1.5 : sqrt(pi)/2.,
	2. : 1., 2.5 : 3.*sqrt(pi)/4., 3. : 2, 3.5 : 15.*sqrt(pi)/8., 4.:  6.}
	return lamr[k]
	
def chi_pdf(x,k=1):
	return (.5)**(k/2.)/gamdic(k/2.)*x**(k/2. - 1.)*exp(-x/2.)