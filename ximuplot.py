from matplotlib import pyplot as plt
import numpy as np
import sys

def P2(mu):
	return .5*(3.*mu**2.-1.)
	
def P4(mu):
	return .125*(35.*mu**4.-30.*mu**2.+3.)


print('constructing xi(mu) from xi0,xi2,xi4 for '+sys.argv[2])

mu = float(sys.argv[1])
d = np.loadtxt(sys.argv[2]).transpose()
xi0 = d[1]
xi2 =  d[2]
xi4 = d[3]
ximu = xi0 + P2(mu)*xi2 + P4(mu)*xi4


plt.plot(d[0],d[0]**2.*ximu)
plt.show()

