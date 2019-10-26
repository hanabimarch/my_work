#!/usr/bin/env python3
import numpy as np
from scipy import stats
from astropy.io import fits
import astropy.stats 
import matplotlib.pyplot as plt 
import time
import scipy.signal



image_file='./data/evt_T0.fits'#fits file for Bayesian block analysis
#data3 = np.loadtxt('summary_t_obs_net_back-3.txt')#analysised data from BBZhang
#fits.info(image_file)
#image_data = fits.getdata(image_file)
hdul = fits.open(image_file)
data = hdul[2].data
cols = hdul[2].columns
#print(hdul[1].data)
#print(data)
#print(cols)
t11,t22=-10,10
t_total=t22-t11

Tellus=np.loadtxt('SNR.txt')
bb=Tellus[:,0].tolist()
MAXINDEX=bb.index(max(bb))
print(bb[MAXINDEX],Tellus[MAXINDEX,0])
data1 = data['TminusT0'][:]#events recorded with time T - T0 
data11 = data['PHA'][:] #PHA channel 
Channel_min=Tellus[MAXINDEX,1]
Channel_max=Tellus[MAXINDEX,2]
#print(data11)
bhere=(np.where(np.logical_and(data1>=t11,data1<=t22)))[0]#serial number of array within time -50s<=T<=20s
#print(bhere)
#bhere0=(np.where(np.logical_and(np.logical_and(data1>=t11,data1<=t22),data11<20)))[0]
#print(bhere0)
data2 =data1[bhere[0]:bhere[int(len(bhere)-1)]]#events recorded in time -50s<=T<=20s


     
hhh=astropy.stats.bayesian_blocks(data2, fitness='events', p0=0.1)#bayesian block: return the time sequence of the new bin size
print(hhh)

Blocksequnce=open('Block.txt', 'w')
for i in range(len(hhh)):
    Blocksequnce.write(str(hhh[i]))
    Blocksequnce.write('\n')
vhh=[0 for i in range(len(hhh)-1)]
xhh=[0 for i in range(len(hhh)-1)]
whh=[0 for i in range(len(hhh)-1)]
BB=open('BB.txt', 'w')
for i in range(len(hhh)-1):
    #tem_a=(np.where(np.logical_and(a>=hhh[i],a<hhh[i+1])))[0]#times sequence serial number from BBZhang within the bayesian bin
    #print(tem_a,hhh[i],hhh[i+1],a,len(tem_a))
    xhh[i]=(hhh[i+1]+hhh[i])/2.0
    whh[i]=hhh[i+1]-hhh[i]
    countn=len(np.where(np.logical_and(data11>=Channel_min,np.logical_and(np.logical_and(data1>=hhh[i],data1<hhh[i+1]),data11<Channel_max)))[0])
    vhh[i]=(countn)/whh[i]
    BB.write(str(xhh[i]))
    BB.write(' ')
    BB.write(str(vhh[i]))
    BB.write(' ')
    BB.write(str(whh[i]))
    BB.write('\n')
BB.close()
'''

BB=np.loadtxt('BB.txt') 
 
peaks=scipy.signal.find_peaks(BB[:,1])#vhh)
UPBINSIZE=BB[peaks[0][0],2]*0.5
binwidth=Tellus[MAXINDEX,3]
AC=open('AC_all.txt', 'w')

a1=[t22+5.-(ia+0.5)*binwidth for ia in range(int((t_total+5.)/binwidth))]
c1=[0 for ic11 in range(int((t_total+5.)/binwidth))]
for ic1 in range(len(a1)):
    c1[ic1]=len(np.where(np.logical_and(data11>=Channel_min,np.logical_and(np.logical_and(data1>=a1[ic1]-0.5*binwidth,data1<a1[ic1]+0.5*binwidth),data11<Channel_max)))[0])/binwidth
    AC.write(str(a1[ic1]))
    AC.write(' ')
    AC.write(str(c1[ic1]))
    AC.write('\n')
AC.close()

'''

