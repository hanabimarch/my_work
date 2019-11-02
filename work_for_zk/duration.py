
import numpy as np
import matplotlib.pyplot as plt
import os
import zzh_py3_file as zhf

datadir = './data/lightcurve1.txt'

savedir = './results/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)
print('xxx')
t,v = zhf.readcol(datadir)
da = np.array([t,v]).T
#t,v = np.loadtxt(datadir,delimiter=' ')
print('xxx')



np.savetxt('./data/lightcurve1.txt',da,fmt='%0.4f')
plt.plot(t,v,color = 'k')
plt.savefig(savedir + 'lightcurve.png')
plt.close()















