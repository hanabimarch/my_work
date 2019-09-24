from zjh_make_phaI import make_phaI
import zzh_py3_file as zhf
import os
from zjh_xspec_fit_kernel import xspec_fit_kernel



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
	file_3_link = result_topdir + bnname + '/pulse_edges.txt'

	if(os.path.exists(file_1_link) and os.path.exists(file_2_link) and os.path.exists(file_3_link)):
		if(os.path.getsize(file_1_link) == 0 or os.path.getsize(file_2_link) == 0 or os.path.getsize(file_1_link) == 0 or os.path.getsize(file_3_link) == 0):
			print(bnname+'探头数量不足')
			continue
			
		print('---------------------------------------')
		print(bnname)
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
					filelist = ['A_slice_'+bnname+'_'+bgo+'.pha']
					for ni in ni_list:
						print(ni)
						filelist.append('A_slice_'+bnname+'_'+ni+'.pha')
						make_phaI(bnname,ni,data_topdir,savedir,edges_start[i],edges_stop[i])
					print(bgo)
					make_phaI(bnname,bgo,data_topdir,savedir,edges_start[i],edges_stop[i])
					print('---------------------------------------')
					try:
						value,value_arr1,value_arr2,flux_list = xspec_fit_kernel(filelist,savedir,makedirs+'Z_'+str(i))
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






























