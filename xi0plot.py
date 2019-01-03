from matplotlib import pyplot as plt
import numpy as np
import sys

print(len(sys.argv))

for i in range(1,len(sys.argv)):
	d = np.loadtxt(sys.argv[i]).transpose()
	plt.plot(d[0],d[0]**2.*d[1])
plt.show()