from tools import *
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os 
from astropy.io import fits
from astropy.stats import bayesian_blocks
#数据库链接---------
data_topdir = '/media/laojin/Elements/trigdata/'
#-------------------
#输出结果地址--------
result_topdir = '/home/laojin/my_work/bxa_work/results/'
#------------------------------------------------------------

#bnname,dete,max_ = readcol('./Ep_HTS_16_samples.txt')
bnname,dete,max_ = readcol('./Ep_HTS_16_samples1.txt')

binsize = 1                            #积分光变曲线用的时间切片的大小。
#----------------------------------------------------------------
bnname = np.array(bnname)
dete = np.array(dete)
max_ = np.array(max_)
#print(bnname,dete)

name_index = np.argsort(bnname)

bnname = bnname[name_index]
dete = dete[name_index]
max_ = max_[name_index]

name0 = 0

name_list = []
sample = {}

for index,name in enumerate(bnname):
	if name != name0:
		name0 = name
		name_list.append(name)
		sample[name] = ([dete[index]],[max_[index]])
		
		
	else:
		sample[name][0].append(dete[index])
		sample[name][1].append(max_[index])

for name in name_list:
	save_dir = result_topdir + name +'/'
	if(os.path.exists(save_dir) == False):
		os.makedirs(save_dir)
	year = '20'+name[2:4]
	sampledir = data_topdir + year + '/' #指向样本数据所在文件路径
	dete,max_ = sample[name]
	dete = np.array(dete)
	max_ = np.array(max_)
	for_block = dete[np.where(max_ == 0)[0][0]]#找到信号最强的
	print(type(for_block))
	datalink = sampledir +name + '/'+ 'glg_tte_'+for_block+'_'+name+'_v*.fit'
	datalink = glob(datalink)[0]
	
	data = fits.open(datalink)
	
	trigtime = data[0].header['TRIGTIME']

	time = data[2].data.field(0)
	ch = data[2].data.field(1)
	t = time-trigtime
	ch_n = data[1].data.field(0)
	e1 = data[1].data.field(1)
	e2 = data[1].data.field(2)
	#t = t[np.where((ch>=3) & (ch <=109))]
	#ch = ch[np.where((ch>=3) & (ch <=109))]
	
	start_edges,stop_edges = get_pulse_duration(t)
	dt = 1
	edges0 = np.arange(start_edges-2,stop_edges+2+dt,dt)
	bin_n,edges0 = np.histogram(t,bins = edges0)
	bin_c = (edges0[1:]+edges0[:-1])*0.5
	edges = bayesian_blocks(bin_c,np.rint(bin_n),fitness = 'events',p0 = 0.001)[1:-1]  #从第二到倒数第二个。
	
	slic_start = edges[:-1]
	slic_top = edges[1:]
			
	block_time = []
	block_time_err = []
	Ep = []
	Ep_err1 = []
	Ep_err2 = []		
	for dete_n in dete:
		make_phaI(name,dete_n,sampledir,save_dir,slic_start,slic_top,binsize = binsize,time_start= -20)
		
	for i in range(len(slic_start)):
		filelist = []
		for dete_n in dete[np.argsort(-1*max_)]:
			filelist.append('A_'+name+'_'+dete_n+'_'+str(i)+'.pha')
	
		value,value_arr1,value_arr2 = xspec_fit_kernel(filelist,save_dir,save_dir,num = i)
		plt.savefig(save_dir+'Z_foldedspec_'+str(i)+'.png')
		plt.close()
		block_time.append((slic_top[i]+slic_start[i])*0.5)
		block_time_err.append((slic_top[i]-slic_start[i])*0.5)
		Ep.append(value)
		Ep_err1.append(value -value_arr1)
		Ep_err2.append(value_arr2 - value)
	printdatatofile(save_dir + 'Y_'+name+'_time_Ep.txt',data = [block_time,block_time_err,Ep,Ep_err1,Ep_err2])
	fig = plt.figure(figsize = (10,10))
	ax1 = fig.add_subplot(1,1,1)
	
	plt.plot(bin_c,bin_n,color = 'g')

	ax2 = ax1.twinx()
	plt.errorbar(block_time,Ep,xerr = block_time_err,yerr = [Ep_err1,Ep_err2],fmt = 'o',ecolor = 'k',color = 'K')
	plt.xlabel('time (s)')
	plt.ylabel(r'Epec ${KeV}$')
	plt.savefig(save_dir + 'Z_plot_'+name+'.png')
	plt.close()
	











