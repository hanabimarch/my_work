import os
import zzh_py3_file as zhf
import numpy as np
import matplotlib.pyplot as plt
from zjh_read import Read_tte_file
from zjh_data_analysis import *
from astropy.io import fits



localdir = '/media/laojin/Elements/trigdata/'
resultlocaldir = '/home/laojin/results/results_for_bb/'
year = '2014' 
name = 'bn140506880'
plot_light_curve = True


detector = 'b0'


time_start = -100
time_stop = 300
ch_top = 128
ch_buttum = 0


topdir = localdir + year + '/'
resultdir = resultlocaldir+year + '/'
filedir = topdir + name + '/'
if(plot_light_curve):
	filem = 'glg_tte_' + detector + '_' + name + '_'
	filename = zhf.findfile(filedir,filem + '*')[0]
	data = fits.open(filedir + filename)
	trigtime = data[0].header['TRIGTIME']
	time = data[2].data.field(0)
	ch = data[2].data.field(1)
	ch_n = data[1].data.field(0)
	e1 = data[1].data.field(1)
	e2 = data[1].data.field(2)
	t = time - trigtime
	t_index = np.where((t >= time_start) & (t <= time_stop))
	t = t[t_index]
	ch = ch[t_index]
	ch_index = np.where((ch >= ch_buttum) & (ch <= ch_top))
	tk = t[ch_index]
	t_h, rate = easy_histogram(tk, binsize=0.01)
	t_r,cs,bs = r_baseline(t_h,rate)
	new_time, new_ch = r_remv_photo(t, ch, ch_n, 1)
	t_, e_ = get_all_energy(new_time, new_ch, ch_n, e1, e2)
	found_duration_dir = resultdir + name + '/found_duration_check_2/'
	start_edges, stop_edges = get_pulse_duration(tk,savedir=found_duration_dir,gamma = np.exp(-13),r=3)


	plt.figure(figsize=(20, 20))
	plt.subplot(2, 2, 1)
	plt.tick_params(labelsize=20)
	plt.title('light curve', size=20)
	plt.plot(t_h, rate, linestyle='steps', color='k')
	plt.xlabel('time (s)', size=20)
	plt.ylabel('counts rate /s', size=20)
	plt.xlim([time_start, time_stop])

	plt.subplot(2, 2, 2)
	plt.tick_params(labelsize=20)
	plt.title('point map', size=20)
	plt.plot(t_, e_, ',', color='k')
	plt.xlabel('time (s)', size=20)
	plt.ylabel('energy (kev)', size=20)
	plt.yscale('log')
	plt.xlim([time_start, time_stop])
	#plt.ylim([5, 9.1e2])
	plt.ylim([120,3e4])

	print('estimating the edges .....')

	if(len(start_edges)>0):
		savedir = resultdir + name + '/'
		if(os.path.exists(savedir) == False):
			os.makedirs(savedir)
		edge_list = [start_edges,stop_edges]
		zhf.printdatatofile(savedir + 'pulse_edges_xxx.txt', edge_list)
		in_list.append(edge_list)
		zhf.printdatatofile(resultlocaldir+'A_add.txt',[[name],[detector]])

	plt.subplot(2,2,3)
	plt.tick_params(labelsize = 20)
	plt.plot(t_h,rate,linestyle = 'steps',color = 'k')
	if(len(start_edges)>0):
		for i in range(len(start_edges)):
			plt.axvline(x = start_edges[i],color = 'r')
			plt.axvline(x = stop_edges[i],color = 'g')
	plt.xlabel('time (s)',size = 20)
	plt.ylabel('counts rate /s',size = 20)
	plt.xlim([time_start,time_stop])
	plt.subplot(2,2,4)
	plt.tick_params(labelsize = 20)
	plt.plot(t_,e_ ,',',color = 'k')
	if(len(start_edges)>0):
		for i in range(len(start_edges)):
			plt.axvline(x = start_edges[i],color = 'r')
			plt.axvline(x = stop_edges[i],color = 'g')
	plt.xlabel('time (s)',size = 20)
	plt.ylabel('energy (kev)',size = 20)
	plt.yscale('log')
	plt.xlim([time_start,time_stop])
	#plt.ylim([5,9.1e2])
	plt.ylim([120,3e4])
	savedir = resultdir + name + '/'
	if(os.path.exists(savedir) == False):
		os.makedirs(savedir)
	plt.savefig(savedir + 'Z_plot_1_2.png')
	#savealldir =resultdir + 'A_plot_1_2_2/'
	#if(os.path.exists(savealldir) == False):
	#	os.makedirs(savealldir)
	#plt.savefig(savealldir +name +'_plot_1_2.png')
	plt.close()
	




