import numpy as np
import matplotlib.pyplot as plt
import zzh_py3_file as zhf
import os
from glob import glob

namelist = zhf.readcol('/home/laojin/results/results_GBM/Y_all_sample_name.txt')[0]
data_topdir = '/home/laojin/results/results_GBM/2017/'

savetop = '/home/laojin/results/spectrum_all/'

if(os.path.exists(savetop)==False):
	os.makedirs(savetop)

yearnamelist = []

for i in namelist:
	if i[:4]=='bn17' :
		yearnamelist.append(i)

print('year name list:\n',yearnamelist)

Epec = []
Epec_err = []
Flux = []
Flux_err = []


for bnname in yearnamelist:

	spectrum_dir = data_topdir + bnname + '/spectrum/'
	if(os.path.exists(spectrum_dir)):
		print('----------------------------------------')
		print(bnname)
		epec_file_list = glob(spectrum_dir+'Z_*_Epec.txt')
		flux_file_list = glob(spectrum_dir+'Z_*_Flux.txt')

		if(len(epec_file_list) != len(flux_file_list) or len(epec_file_list)==0 or len(flux_file_list)==0):
			print(bnname+'文件数不相等！')
			print('----------------------------------------')
			continue
		epec_file_list = np.sort(epec_file_list)
		flux_file_list = np.sort(flux_file_list)
		for index,epec_file_name in enumerate(epec_file_list):
			print(epec_file_name)
			epec,epec_er1,epec_er2 = zhf.readcol(epec_file_name)
			print(flux_file_list[index])
			flux,flux_p,flux_a,p_n,p_n_p,p_n_a = zhf.readcol(flux_file_list[index])
			epec_value = epec[0]
			epec_err = (epec_er2[0] - epec_er1[0])*0.5
			flux_value = np.mean(flux)
			flux_err = (np.mean(flux_a)-np.mean(flux_p))*0.5
			Epec.append(epec_value)
			Epec_err.append(epec_err)
			Flux.append(flux_value)
			Flux_err.append(flux_err)
		print('----------------------------------------')

zhf.printdatatofile(savetop+'epec_flux.txt',data = [Epec,Epec_err,Flux,Flux_err])

plt.figure(figsize = (10,10))
plt.subplot(1,1,1)
plt.title('The correlation of Epec and Flux of single pulse')
plt.errorbar(Epec,Flux,xerr=Epec_err,yerr=Flux_err,zorder=1,fmt = 'o',mfc = 'r',ecolor='g')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Epec ${KeV}$')
plt.ylabel(r'Flux ${ergs/cm^{2}/s}$')
plt.ylim([1e-7,1e-3])
plt.savefig(savetop+'epec_flux.png')
plt.close()





