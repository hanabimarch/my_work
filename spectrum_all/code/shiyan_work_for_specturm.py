from zjh_make_phaI import *
import zzh_py3_file as zhf
import os
from zjh_xspec_fit_kernel import xspec_fit_kernel
import matplotlib.pyplot as plt


data_topdir = '/home/laojin/trigdata/2017/'
result_topdir = '/home/laojin/results/results_GBM/2017/'

namelist = zhf.readcol('/home/laojin/results/results_GBM/Y_all_sample_name.txt')[0]


yearnamelist = []

for i in namelist:
	if i[:4]=='bn17' :
		yearnamelist.append(i)

print('year name list:\n',yearnamelist)


for bnname in yearnamelist:
	file_1_link = result_topdir + bnname + '/good_NaI_detectors.txt'
	file_2_link = result_topdir + bnname + '/good_BGO_detectors.txt'
	file_3_link = result_topdir + bnname + '/0_pulse_edges.txt'
	save_all = result_topdir + 'A_spectrum_1/'
	save_all2 = result_topdir + 'A_spectrum_2/'
	if(os.path.exists(save_all) == False):
		os.makedirs(save_all)
	if (os.path.exists(save_all2) == False):
		os.makedirs(save_all2)

	if(os.path.exists(file_1_link) and os.path.exists(file_2_link) and os.path.exists(file_3_link)):

		ni_s,ni_list = zhf.readcol(file_1_link)
		bgo = zhf.readcol(file_2_link)[0][0]

		if(len(ni_s)>=2):
			try:
				ni_list = ni_list[:2]

				makedirs = result_topdir + bnname + '/spectrum/'
				if(os.path.exists(makedirs) == False):
					os.makedirs(makedirs)

				edges_start,edges_stop = zhf.readcol(file_3_link)
				for i in range(len(edges_start)):
					print('--------------------------------------')
					print('time slice ',i)
					savedir = makedirs + 'time_' + str(i) + '/'
					filelist = ['A_'+bnname+'_'+bgo+'.pha']
					for ni in ni_list:
						print(ni)
						filelist.append('A_'+bnname+'_'+ni+'.pha')
						make_phaI_with_events(bnname,ni,data_topdir,savedir,edges_start[i],edges_stop[i])
						plt.savefig(save_all+'A_'+bnname+'_'+ni+'_'+str(i)+'.png')
						plt.close()
					print(bgo)
					make_phaI_with_events(bnname,bgo,data_topdir,savedir,edges_start[i],edges_stop[i])
					plt.savefig(save_all+'A_'+bnname+'_'+bgo+'_'+str(i)+'.png')
					plt.close()
					print('---------------------------------------')
					try:
						value,value_arr1,value_arr2,flux_list = xspec_fit_kernel(filelist,savedir,makedirs+'Z_'+str(i))
						#plt.savefig(save_all2 + 'A_'+bnname+'_'+str(i)+'.png')
						#plt.close()
						copy_rspI(savedir+'foldedspec.png',save_all2 + 'A_'+bnname+'_'+str(i)+'.png')
						zhf.printdatatofile(makedirs+'Z_'+str(i)+'_Flux.txt',data = flux_list)
						zhf.printdatatofile(makedirs+'Z_'+str(i)+'_Epec.txt',data = [[value],[value_arr1],[value_arr2]])
					except:
						print(bnname+' '+str(i)+'出现错误，跳过拟合！')
						continue
			except(IndexError):
				print(bnname+'文件内容错误，跳过拟合！')
				continue
		else:
			print(bnname+'探头数量不足')
	else:
		print(bnname+'信息不完整。')






























