import os
import zzh_py3_file as zhf
import numpy as np
import matplotlib.pyplot as plt
from zjh_read import Read_tte_file
from zjh_data_analysis import *
from astropy.io import fits
from multiprocessing import Pool


#yearlist = ['2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018']
yearlist = ['2019']

def run_f(year):
	localdir = '/media/laojin/Elements/trigdata/'
	resultlocaldir = '/home/laojin/results/results_GBM/'

	NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
	BGO = ['b0','b1']

	time_start = -100
	time_stop = 300
	ch_top = 110
	ch_buttum = 8

	g_ch_top = 110
	g_ch_buttum = 7

	plot_light_curve = True

	topdir = localdir + year + '/'
	resultdir = resultlocaldir+year + '/'
	dirlist1 = os.listdir(topdir)
	print(year+' have :',len(dirlist1))
	dirlist = []
	for dirl in dirlist1:
		if(os.path.isdir(topdir+dirl)):
			dirlist.append(dirl)
	print('The number of sample is ',len(dirlist))
	dirlist = np.sort(dirlist)
	one_pulse = []
	two_pulse = []
	three_pulse = []
	complex_pulse = []
	for name in dirlist:
		filedir = topdir + name + '/'

		if(plot_light_curve):
			data_c = {}
			ni_max = []

			ni_good = []

			print('function analysis doing---')
			print(name)
			for ni in NaI:
				print(ni,end = '\r')
				filem = 'glg_tte_' + ni + '_' + name + '_'
				filename = zhf.findfile(filedir,filem + '*')
				if filename == False:
					continue
				else:
					filename = filename[0]
				data = fits.open(filedir + filename)
				data_c[ni] = data
				trigtime = data[0].header['TRIGTIME']
				time = data[2].data.field(0)
				ch = data[2].data.field(1)
				t = time - trigtime

				ch_index = np.where((ch >= ch_buttum) & (ch <= ch_top))
				t = t[ch_index]
				t_index = np.where((t >= time_start) & (t <= time_stop))
				t = t[t_index]
				tif, rif = easy_histogram(t, binsize=1)
				t_r,cs,bs = r_baseline(tif,rif)
				ev = cs+bs.sum()/bs.size
				if(rif[np.where(rif == 0)].size > 1000):
					ni_good.append(False)
				else:
					ni_good.append(True)
				ni_max.append(ev.max())

			print('next.....')
			index_ni = np.argsort(ni_max)
			nai_sort = np.array(NaI)[index_ni]
			max_ni_name = [nai_sort[-1],nai_sort[-2],nai_sort[-3],nai_sort[-4]]
			for index,ni in enumerate(max_ni_name):
				print(ni,end = '\r')
				data = data_c[ni]
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
				in_list = [t_r,rate,cs,bs]

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
				plt.ylim([5, 9.1e2])


				print('estimating the edges .....')

				if(len(start_edges)>0):
					savedir = resultdir + name + '/'
					if(os.path.exists(savedir) == False):
						os.makedirs(savedir)
					edge_list = [start_edges,stop_edges]
					zhf.printdatatofile(savedir + 'pulse_edges_1.txt', edge_list)
					in_list.append(edge_list)

					if (len(start_edges) == 1):
						one_pulse.append([name, ni])
					elif (len(start_edges) == 2):
						two_pulse.append([name, ni])
					elif (len(start_edges) == 3):
						three_pulse.append([name, ni])
					else:
						complex_pulse.append([name, ni])

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
				plt.ylim([5,9.1e2])
				savedir = resultdir + name + '/'
				if(os.path.exists(savedir) == False):
					os.makedirs(savedir)
				plt.savefig(savedir + 'Z_plot_1_2.png')
				savealldir =resultdir + 'A_plot_1_2_2/'
				if(os.path.exists(savealldir) == False):
					os.makedirs(savealldir)
				plt.savefig(savealldir +name +'_plot_1_2.png')
				plt.close()
				if(len(start_edges)>0):
					break


	if(len(dirlist)>=10):
		savedir_all = resultdir+'A_plot_1_2_2/'
		if(os.path.exists(savedir_all) == False):
			os.makedirs(savedir_all)
		zhf.printdatatofile(savedir_all+'Z_one_pulse_sample_list.txt',zhf.transpose(one_pulse))
		zhf.printdatatofile(savedir_all+'Z_two_pulse_sample_list.txt',zhf.transpose(two_pulse))
		zhf.printdatatofile(savedir_all+'Z_three_pulse_sample_list.txt',zhf.transpose(three_pulse))
		zhf.printdatatofile(savedir_all+'Z_complex_pulse_sample_list.txt',zhf.transpose(complex_pulse))

pool = Pool(4)
pool.map(run_f,yearlist)

