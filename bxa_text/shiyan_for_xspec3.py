import numpy as np
import matplotlib.pyplot as plt
from xspec import *
import os
import bxa.xspec as bxa
import scipy.stats
savedir = './'
if(os.path.exists(savedir) == False):
	os.makedirs(savedir)

#brightdet = ['A_slice_bn190114873_b1.pha','A_slice_bn190114873_n7.pha','A_slice_bn190114873_n4.pha']
#brightdet = ['A_slice_bn190114873_b1.pha','A_slice_bn190114873_n7.pha']

#rsplist = ['glg_cspec_n7_bn190114873_v02.rsp','glg_cspec_b1_bn190114873_v02.rsp','glg_cspec_n4_bn190114873_v02.rsp']
#baselist = ['A_slice_bn190114873_n7.bkg','A_slice_bn190114873_b1.bkg','A_slice_bn190114873_n4.bkg']

#alldatastr = ' '.join(brightdet)
alldatastr = 'example-file.fak'
print(alldatastr)
#print(brightdet)


print('$$$$')
'''
for index,files in enumerate(brightdet):
	obj = Spectrum(files)
	obj.background = baselist[index]
	obj.response = rsplist[index]
'''
#AllData.show()
Fit.statMethod='cstat'#设置拟合形式，这里是泊松信号加高斯背景。
Plot.xAxis = 'keV'
#s = Spectrum('example-file.fak')
#s.ignore("**"); s.notice("0.2-8.0")

AllData(alldatastr)
AllData.ignore("**")
AllData.notice("0.2-8.0")
print('AllData.notice')
print(AllData.notice)
print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

#Model('cutoffpl + bbody')#模型设置
m=Model("wabs*pow + gauss")
m.wabs.nH.values = ",,0.01,0.01,1000,1000"
m.powerlaw.norm.values = ",,1e-10,1e-10,1e1,1e1"
m.gaussian.norm.values = ",,1e-10,1e-10,1e1,1e1"
m.gaussian.LineE.values = ",,6,6,8,8"
m.gaussian.Sigma.values = ",,0.001,0.001,1,1"

def my_custom_prior(u):
	# prior distributions transform from 0:1 to the parameter range
	# here: a gaussian prior distribution, cut below 1/above 3
	x = scipy.stats.norm(1.9, 0.15).ppf(u)
	if x < 1.:
		x = 1
	if x > 3:
		x = 3
	return x

transformations = [
	# jeffreys prior for nH (but see below)
	bxa.create_jeffreys_prior_for(m, m.wabs.nH),
	# jeffreys prior for scale variable
	bxa.create_jeffreys_prior_for(m, m.powerlaw.norm),
	# custom gaussian prior function for photon index
	bxa.create_custom_prior_for(  m, m.powerlaw.PhoIndex, my_custom_prior),
	# and possibly many more
	# jeffreys prior for scale variable
	bxa.create_jeffreys_prior_for(m, m.gaussian.norm),
	# uniform prior for location variable
	bxa.create_uniform_prior_for(m, m.gaussian.LineE),
	# jeffreys prior for scale variable
	bxa.create_jeffreys_prior_for(m, m.gaussian.Sigma),
]
outputfiles_basename = 'simplest3-'
bxa.standard_analysis(transformations,outputfiles_basename = outputfiles_basename,verbose=True, resume=False)
#AllModels.calcFlux("8. 40000.0 err") #参数需要一个能量的范围，之前在拟合光谱时我们设置了一个范围，我们暂时用它。
#AllModels.calcFlux('1:**-200.,40000.0,err-** 2-3:**-8.,800.,err-**')#不能这么用
#par3=AllModels(1)(3)#第一个模型的第三个参数
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
#print(par3.values)
#print(par3.error)
#value = par3.values[0]
#value_arr1,value_arr2,ffff = par3.error
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')

#AllModels.calcFlux(".1 10.0 err")
#AllModels.calcLumin(".1 10. .05 err")

#flux = AllData(1).flux              #计算流量，‘1’代表第一个光谱的流量
#lumin = AllModels(1).lumin
#print('1 flux:',list(flux))               #流量，流量下界，流量上界，光子数，光子数下，光子数上
#print('lumin:',lumin)
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')

'''
Plot.device='/xs'
Plot.xAxis='keV'
Plot.yLog=True
'''
Plot('eeufspec')

#for i in range(len(brightdet)):
for i in [0]:
	print(i)
	energies=Plot.x(i+1)
	rates=Plot.y(i+1)
	folded=Plot.model(i+1)
	xErrs=Plot.xErr(i+1)
	yErrs=Plot.yErr(i+1)
	plt.errorbar(energies,rates,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
	#plt.plot(energies,folded,color='black',zorder=2)
	plt.plot(energies,folded,zorder=2)
#plt.axvline(x = value,color = 'r')
#plt.axvline(x = value_arr1,color = 'g')
#plt.axvline(x = value_arr2,color = 'g')
plt.xlabel('Energy KeV')
plt.ylabel(r'${KeV^{2} (Photons cm^{-2}s^{-1}keV^{-1})}$')
plt.xscale('log')
plt.yscale('log')
plt.savefig(savedir + 'foldedspec.png')
plt.close()

Plot('eeufspec')
'''
for i in range(len(brightdet)):
	energies=Plot.x(i+1)
	ufspec=Plot.y(i+1)
	folded=Plot.model(i+1)
	xErrs=Plot.xErr(i+1)
	yErrs=Plot.yErr(i+1)
	plt.errorbar(energies,ufspec,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
	plt.plot(energies,folded,color='black',zorder=2)
#plt.axvline(x = value,color = 'r')
#plt.axvline(x = value_arr1,color = 'g')
#plt.axvline(x = value_arr2,color = 'g')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy KeV')
plt.ylabel(r'${\/KeV\/^{2} (\/Photons \/cm^{-2}\/s^{-1}\/keV^{-1})}$')
plt.savefig(savedir +'eeufspec.png')
plt.close()

'''




