from tools import *
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os 
from astropy.io import fits
from astropy.stats import bayesian_blocks
#数据库链接---------
data_topdir = ''
#-------------------
#输出结果地址--------
result_topdir = '/home/laojin/my_work/spectrum_all/pyxpec_spectral_analysis/results/'
#---------------------

bnname,dete,max_ = readcol('./Ep_HTS_16_samples.txt')

#print(bnname,dete)

name_index = np.argsort(bnname)

bnname = bnname[name_index]
dete = dete[name_index]
max_ = max_[name_index]

name0 = 0

name_list = []
sample = {}

for index,name in enumerate(bnnname):
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
	sampledir = data_topdir + year + '/'+name + '/' #指向样本数据所在文件路径
	dete,max_ = sample[name]
	for_block = dete[np.where(max_ == 0)[0][0]]#找到信号最强的
	
	datalink = sampledir + 'glg_tte_'+for_block+'_'+name+'_v*.fit'
	datalink = glob(datalink)[0]
	
	data = fits.open(datalink)
	
	trigtime = data[0].header['TRIGTIME']

	time = data[2].data.field(0)
	ch = data[2].data.field(1)
	t = time-trigtime
	ch_n = data[1].data.field(0)
	e1 = data[1].data.field(1)
	e2 = data[1].data.field(2)
	t = t[np.where((ch>=3) & (ch <=109))]
	ch = ch[np.where((ch>=3) & (ch <=109))]
	
	start_edges,stop_edges = get_pulse_duration(t)
	dt = 0.064
	edges0 = np.arange(start_edges-2,stop_edges+2+dt,dt)
	bin_n,edges0 = np.histogram(t,bins = edges0)
	bin_c = (edges0[1:]+edges0[:-1])*0.5
	edges = bayesian_blocks(bin_c,bin_n,fitness = 'events')[1:-1]#从第一到倒数第二个。
	
	slic_start = edges[:-1]
	slic_top = edges[1:]
	
	for dete_n in dete:
		make_phaI(name,dete_n,sampledir,save_dir,slic_start,slic_top,binsize = 0.02)
	
	
	
	












