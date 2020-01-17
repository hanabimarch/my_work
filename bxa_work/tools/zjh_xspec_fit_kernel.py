from xspec import *
import os
import matplotlib.pyplot as plt
import numpy as np
import bxa.xspec as bxa



def xspec_fit_kernel(filelist,datadir,savedir,num = 0,flux = False):
	os.chdir(datadir)
	Plot.xAxis='keV'
	Plot.yLog=True
	alldatastr = ' '.join(filelist)
	AllData.clear()
	AllModels.clear()
	AllData(alldatastr)
	#AllData.show()
	Fit.statMethod='cstat'
	
	AllData.ignore('1:**-200.0,40000.0-** 2-4:**-8.0,800.0-**')#我觉得这里可以留一个接口。
	m = Model('grbm')
	

	#--------------------------------------------
	'''
	设置grbm模型的三个参数的范围。
	'''
	#m.grbm.alpha.values = ',,-10,-10,5,5'
	#m.grbm.beta.values = ',,-10,-10,10,10'
	m.grbm.tem.values = ',,1e-10,1e-10,1e6,1e6'
	#--------------------------------------------
	transformations = [bxa.create_uniform_prior_for(m,m.grbm.alpha),
			   bxa.create_uniform_prior_for(m,m.grbm.beta),
			   #bxa.create_uniform_prior_for(m,m.grbm.tem)]
			   bxa.create_jeffreys_prior_for(m,m.grbm.tem)]
	outputdir = savedir+'results'+str(num)+'/'
	if os.path.exists(outputdir) == False:
		os.makedirs(outputdir)
	os.chdir(outputdir)
	outputfiles_basename = 'simpest_'
	
	bxa.standard_analysis(transformations,outputfiles_basename = outputfiles_basename,verbose =True,resume=False,skipsteps = ['convolved','qq'])
	par3=AllModels(1)(3)#第一个模型的第三个参数
	value = par3.values[0]
	value_arr1,value_arr2,ffff = par3.error
	Plot('eeufspec')
	flux_list = []
	for i in range(len(filelist)):
		print(i)
		if(flux):
			flux = AllData(i+1).flux
			flux_list.append(flux)
		energies=Plot.x(i+1)
		rates=Plot.y(i+1)
		folded=Plot.model(i+1)
		xErrs=Plot.xErr(i+1)
		yErrs=Plot.yErr(i+1)
		plt.errorbar(energies,rates,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
		plt.plot(energies,folded,color='black',zorder=2)
	plt.axvline(x = value,color = 'r')
	plt.axvline(x = value_arr1,color = 'g')
	plt.axvline(x = value_arr2,color = 'g')
	plt.xlabel('Energy KeV')
	plt.ylabel(r'${KeV^{2} (Photons cm^{-2}s^{-1}keV^{-1})}$')
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig(savedir + 'foldedspec.png')
	#plt.close()
	if flux :

		return value,value_arr1,value_arr2,np.array(flux_list).T
	else:
		return value,value_arr1,value_arr2















