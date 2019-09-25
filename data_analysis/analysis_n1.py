import os
import zzh_py3_file as zhf
import numpy as np
import matplotlib.pyplot as plt
from zjh_read import Read_tte_file
from zjh_data_analysis import *
import PIL.Image as Image
from astropy.stats import bayesian_blocks
from astropy.io import fits
import astropy.units as u
from zjh_gbm import GBM
from mpl_toolkits.basemap import Basemap



yearlist = ['2008']
localdir = '/media/laojin/Elements/trigdata/'
resultlocaldir = '/home/laojin/results/results_GBM/'


time_start = -100
time_stop = 300


ch_top = 110
ch_buttum = 8

g_ch_top = 110
g_ch_buttum = 7

plot_light_curve = True

plot_single_detector = True
plot_count_map = False

plot_location = True


NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
BGO = ['b0','b1']

for year in yearlist:
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

		if(plot_light_curve or plot_count_map):
			data_c = {}
			ni_max = []
			gbo_max = []

			ni_good = []
			gbo_good = []
			print('function analysis doing---')
			print(name)
			for ni in NaI:
				print(ni)
				filem = 'glg_tte_' + ni + '_' + name + '_'
				filename = zhf.findfile(filedir,filem + '*')[0]
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
				tif, rif = easy_histogram(t, binsize=0.5)
				t_r,cs,bs = r_baseline(tif,rif)
				ev = cs+bs.sum()/bs.size
				if(rif[np.where(rif == 0)].size > 1000):
					ni_good.append(False)
				else:
					ni_good.append(True)
				ni_max.append(ev.max())

			for bgoi in BGO:
				print(bgoi)
				filem = 'glg_tte_' + bgoi + '_' + name + '_'
				filename = zhf.findfile(filedir,filem + '*')[0]
				data = fits.open(filedir + filename)
				data_c[bgoi] = data
				trigtime = data[0].header['TRIGTIME']
				time = data[2].data.field(0)
				ch = data[2].data.field(1)
				t = time - trigtime
				ch_index = np.where((ch >= g_ch_buttum) & (ch <= g_ch_top))
				t = t[ch_index]
				t_index = np.where((t >= time_start) & (t <= time_stop))
				t = t[t_index]
				tif, rif = easy_histogram(t, binsize=0.5)
				t_r,cs,bs = r_baseline(tif,rif)
				ev = cs+bs.sum()/bs.size
				if(rif[np.where(rif == 0)].size > 2000):
					gbo_good.append(False)
				else:
					gbo_good.append(True)
				gbo_max.append(ev.max())
			print('next.....')
			index_ni = np.argsort(ni_max)
			nai_sort = np.array(NaI)[index_ni]
			max_ni_name = [nai_sort[-1],nai_sort[-2],nai_sort[-3],nai_sort[-4]]
			max_ni_all = np.max(ni_max)

			index_bgo = np.argsort(gbo_max)
			bgo_sort = np.array(BGO)[index_bgo]
			max_bgo = bgo_sort[-1]
			max_bgo_all = np.max(gbo_max)

			ni_list = []
			gbo_list = []
			in_list_num = []
			in_list_name = []
			for index,ni in enumerate(NaI):
				print(ni)
				data = data_c[ni]
				trigtime = data[0].header['TRIGTIME']
				time = data[2].data.field(0)
				ch = data[2].data.field(1)
				t = time - trigtime
				ch_index = np.where((ch >= ch_buttum) & (ch <= ch_top))
				t = t[ch_index]
				t_index = np.where((t >= time_start) & (t <= time_stop))
				t = t[t_index]
				t_h, rate = easy_histogram(t, binsize=0.01)
				t_r,cs,bs = r_baseline(t_h,rate)
				in_list = [t_r,rate,cs,bs]

				for num,max_name in enumerate(max_ni_name):
					if((ni == max_name)&(ni_good[index])):
						print('estimating the edges .....')
						found_duration_dir = resultdir + name + '/NaI_found_duration_check/'
						start_edges,stop_edges = get_pulse_duration(t,savedir=found_duration_dir)
						if(len(start_edges)>0):
							edge_list = [start_edges,stop_edges]
							savedir = resultdir + name + '/'
							if(os.path.exists(savedir) == False):
								os.makedirs(savedir)
							zhf.printdatatofile(savedir+str(num)+'_pulse_edges.txt',edge_list)
							in_list.append(edge_list)
							in_list_num.append(num)
							in_list_name.append(max_name)
							if(ni == max_ni_name[0]):
								if(len(start_edges) == 1):
									one_pulse.append([name,ni])
								elif(len(start_edges) == 2):
									two_pulse.append([name,ni])
								elif(len(start_edges) == 3):
									three_pulse.append([name,ni])
								else:
									complex_pulse.append([name,ni])


						print('finish estimating .')
				ni_list.append(in_list)

			in_list_num = np.array(in_list_num)
			in_list_name = np.array(in_list_name)
			sort_num_index = np.argsort(in_list_num)
			in_list_num = in_list_num[sort_num_index]
			in_list_name = in_list_name[sort_num_index]
			savedir = resultdir + name + '/'
			if (os.path.exists(savedir) == False):
				os.makedirs(savedir)
			zhf.printdatatofile(savedir + 'good_NaI_detectors.txt', [in_list_num, in_list_name])



			for index,bgoi in enumerate(BGO):
				print(bgoi)
				data = data_c[bgoi]
				trigtime = data[0].header['TRIGTIME']
				time = data[2].data.field(0)
				ch = data[2].data.field(1)
				t = time - trigtime
				ch_index = np.where((ch >= g_ch_buttum) & (ch <= g_ch_top))
				t = t[ch_index]
				t_index = np.where((t >= time_start) & (t <= time_stop))
				t = t[t_index]
				t_h, rate = easy_histogram(t, binsize=0.01)
				t_r,cs,bs = r_baseline(t_h,rate)
				gbo_list.append([t_h,rate,cs,bs])
				if(bgoi == max_bgo and gbo_good[index] == True ):
					found_duration_dir = resultdir + name + '/BGO_found_duration_check/'
					start_edges,stop_edges = get_pulse_duration(t,savedir=found_duration_dir)
					if(len(start_edges)>0):
						savedir = resultdir + name + '/'
						if (os.path.exists(savedir) == False):
							os.makedirs(savedir)
						zhf.printdatatofile(savedir + 'good_BGO_detectors.txt', [[bgoi]])



			#到这里初步的处理就搞完了。
			if(plot_light_curve):
				print('plot light curve ....')
				plt.figure(figsize = (30,60))
				plt.subplots_adjust(left = 0.1,right = 0.9,top = 0.95,bottom = 0.05)
				for index,value in enumerate(NaI):
					print(value)
					inlist = ni_list[index]

					plt.subplot(14,1,index+1)
					plt.plot(inlist[0],inlist[1],color = 'k',linestyle = 'steps',label = value)
					plt.plot(inlist[0],inlist[3],color = 'r',label = 'back')
					if(len(inlist) == 5):
						edge_list = inlist[4]
						start_edges = edge_list[0]
						stop_edges = edge_list[1]
						for index,value1 in enumerate(start_edges):
							plt.axvline(x = value1,color = 'r')
							plt.axvline(x = stop_edges[index],color = 'g')
					plt.ylabel('the count rate (N/s)')
					plt.xlim([time_start,time_stop])
					plt.ylim([0,max_ni_all*2.])
					plt.legend(loc='upper left')
				for index,value in enumerate(BGO):
					print(value)
					inlist = gbo_list[index]
					plt.subplot(14,1,index+13)
					plt.plot(inlist[0],inlist[1],color = 'k',linestyle = 'steps',label = value)
					plt.plot(inlist[0],inlist[3],color = 'r',label = 'back')
					plt.ylabel('the count rate (N/s)')
					plt.xlim([time_start,time_stop])
					plt.ylim([0,max_bgo_all*2.5])
					plt.legend(loc='upper left')
				savedir = resultdir + name + '/'
				if(os.path.exists(savedir) == False):
					os.makedirs(savedir)
				plt.savefig(savedir+'ligth_curves_of_all_detectors.png')
				savedir_all = resultdir+'A_plot_ligth_curves/'
				if(os.path.exists(savedir_all) == False):
					os.makedirs(savedir_all)
				plt.savefig(savedir_all+name+'_ligth_curves_of_all_detectors.png')
				plt.close()
				print('finish !')

			if (plot_count_map):
				count_map = {}
				for index,ni in enumerate(NaI):
					data = data_c[ni]
					time = data[2].data.field(0)
					ch = data[2].data.field(1)
					ch_n = data[1].data.field(0)
					e1 = data[1].data.field(1)
					e2 = data[1].data.field(2)
					trigtime = data[0].header['TRIGTIME']
					t = time - trigtime
					t_index = np.where((t>=time_start)&(t<=time_stop))
					t = t[t_index]
					ch = ch[t_index]
					new_time,new_ch = r_remv_photo(t,ch,ch_n,1)
					time,energy = get_all_energy(new_time,new_ch,ch_n,e1,e2)
					count_map[ni] = [time,energy]

				for index,bgoi in enumerate(BGO):
					data = data_c[bgoi]
					time = data[2].data.field(0)
					ch = data[2].data.field(1)
					ch_n = data[1].data.field(0)
					e1 = data[1].data.field(1)
					e2 = data[1].data.field(2)
					trigtime = data[0].header['TRIGTIME']
					t = time - trigtime
					t_index = np.where((t>=time_start)&(t<=time_stop))
					t = t[t_index]
					ch = ch[t_index]
					new_time,new_ch = r_remv_photo(t,ch,ch_n,1)
					time,energy = get_all_energy(new_time,new_ch,ch_n,e1,e2)
					count_map[bgoi] = [time,energy]


				plt.figure(figsize = (40,40))
				for index,ni in enumerate(NaI):
					time,energy = count_map[ni]
					plt.subplot(4,4,index+1)
					plt.title('the count map for '+ni)
					plt.plot(time,energy,',',color = 'k')
					plt.yscale('log')
					plt.xlabel('time (s)')
					plt.ylabel('energy (KeV)')
					plt.xlim([time_start,time_stop])
					plt.ylim([5,9.1e2])
				k = len(NaI)
				for index,bgo in enumerate(BGO):
					time,energy = count_map[bgo]
					plt.subplot(4,4,index+1+k)
					plt.title('the count map for '+bgo)
					plt.plot(time,energy,',',color = 'k')
					plt.yscale('log')
					plt.xlabel('time (s)')
					plt.ylabel('energy (KeV)')
					plt.xlim([time_start,time_stop])
					plt.ylim([120,3e4])
				savedir = resultdir + name + '/'
				if(os.path.exists(savedir) == False):
					os.makedirs(savedir)
				plt.savefig(savedir + name + '_substract_baseline_count_map.eps')
				plt.close()
				png_savedir = resultdir+'A_plot_substract_baseline_count_map/'
				if(os.path.exists(png_savedir) == False):
					os.makedirs(png_savedir)
				l = Image.open(savedir + name + '_substract_baseline_count_map.eps')
				l.save(png_savedir + name + '_substract_baseline_count_map.png')
				print('I finish the plotting of '+name)
				print('------------------------------------')

	if(len(dirlist)>=10):
		savedir_all = resultdir+'A_plot_ligth_curves/'
		if(os.path.exists(savedir_all) == False):
			os.makedirs(savedir_all)
		zhf.printdatatofile(savedir_all+'Z_one_pulse_sample_list.txt',zhf.transpose(one_pulse))
		zhf.printdatatofile(savedir_all+'Z_two_pulse_sample_list.txt',zhf.transpose(two_pulse))
		zhf.printdatatofile(savedir_all+'Z_three_pulse_sample_list.txt',zhf.transpose(three_pulse))
		zhf.printdatatofile(savedir_all+'Z_complex_pulse_sample_list.txt',zhf.transpose(complex_pulse))














