
from xspec import *
import bxa.xspec as bxa


Plot.xAxis='keV'
Plot.yLog=True

alldatastr='n4.pha'
AllData(alldatastr)
AllData.ignore('1:**-8.0,800.0-**  2:**-200.0,40000.0-**')
Fit.statMethod='cstat'
m = Model("grbm")
m.grbm.alpha.values = ',,-10,-10,5,10'
m.grbm.beta.values = ',,-10,-10,10,20'
m.grbm.tem.values = ',,1e-10,1e-5,1e5,1e6'


print(m.grbm.tem.values)

transformations = [bxa.create_uniform_prior_for( m, m.grbm.alpha),
                   bxa.create_uniform_prior_for( m, m.grbm.beta),
                   bxa.create_jeffreys_prior_for( m, m.grbm.tem)]

outputfiles_basename = 'simplest3-'
bxa.standard_analysis(transformations,outputfiles_basename = outputfiles_basename,verbose=True, resume=False,skipsteps = ['convolved'])

par3=AllModels(1)(3)
print(par3.values[0],par3.error[0])


