import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import zzh_py3_file as zhf
import os
from Calculate_duration import Plot
from Calculate_duration import get_txx
topdir = '/home/laojin/trigdata/2017/'

bn,ni = zhf.readcol('./samples.txt')

indexbn = np.argsort(bn)

samples = np.array([bn,ni]).T
samples = samples[indexbn]
savetop = './results/'

if os.path.exists(savetop) == False:
	os.makedirs(savetop)


time_start = -100
time_stop = 400
ch_top = 91
ch_buttum = 7
binsize = 0.064
for index,value in enumerate(samples):
	print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
	print(value)
	filedir = topdir+value[0]+'/'
	filename = zhf.findfile(filedir,'glg_tte_'+value[1]+'_'+value[0]+'*')[0]
	data = fits.open(filedir + filename)
	trigtime = data[0].header['TRIGTIME']
	time = data[2].data.field(0)
	ch = data[2].data.field(1)

	ch_n = data[1].data.field(0)
	e1 = data[1].data.field(1)
	e2 = data[1].data.field(2)

	t = time - trigtime
	t_index = np.where((t>=time_start) & (t<=time_stop))
	t = t[t_index]
	ch = ch[t_index]
	ch_index = np.where((ch>=ch_buttum) & (ch<=ch_top))

	t = t[ch_index]
	ch = ch[ch_index]

	results = get_txx(t,binsize=binsize,block_time = 40,time_unified=False,hardness = 200*1/binsize,sigma = 1)
	myplt = Plot(results)
	myplt.plot_light_curve()
	savedir = savetop+value[0]+'/'
	if(os.path.exists(savedir) == False):
		os.makedirs(savedir)
	plt.savefig(savedir + 'A_light_curve.png')
	plt.savefig(savetop+'A_'+value[0]+'_light_curve.png')
	plt.close()
	if results['good']:
		for i in range(len(results['txx'])):
			myplt.plot_distribution('90',num = i)
			plt.savefig(savedir + 'D_distribution_'+str(i)+'.png')
			plt.close()
	plt.figure(figsize = (10,10))
	plt.subplot(2,1,1)
	myplt.plot_Txx1('90')
	plt.subplot(2,1,2)
	myplt.plot_Txx2('90')
	plt.savefig(savedir + 'B_txx.png')
	plt.savefig(savetop+'A_'+value[0]+'_txx.png')
	plt.close()

	myplt.plot_normallization()
	plt.savefig(savedir+'C_normall.png')
	plt.close()
