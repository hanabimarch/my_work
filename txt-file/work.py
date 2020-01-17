import numpy as np
import pandas as pd
import os

datadir = './old_data/'
savedir = './new_data/'

if os.path.exists(savedir)==False:
	os.makedirs(savedir)

namelist = os.listdir(datadir)
print('name list :\n',namelist)
print('------------------------------------')
for name in namelist:
	print(name)
	x1,x2 = np.loadtxt(datadir+name,dtype=np.float32).T
	x3 = np.zeros(x1.size)
	new = np.vstack([x1,x2,x3]).T
	df = pd.DataFrame(new,columns = ['wavelength','flux','error'])
	#print(df)
	df.to_csv(savedir+'new_'+name,sep = '\t',index = False,float_format = '%.6f')

print('-----------------------------------')















