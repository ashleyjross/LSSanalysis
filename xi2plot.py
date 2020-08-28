from matplotlib import pyplot as plt
import numpy as np
import sys

print(len(sys.argv))

fac = 1.
mm = 0
try:
	fac = float(sys.argv[-1])
	mm = -1
	print('fac ='+str(fac)+' gets applied to first file')
except:
	print('fac = 1')


for i in range(1,len(sys.argv)+mm):
	d = np.loadtxt(sys.argv[i]).transpose()
	if i != 1:
		fac = 1.
	plt.plot(d[0],d[0]*d[1]*fac)
plt.show()