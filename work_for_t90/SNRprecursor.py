#!/usr/bin/env python3
import numpy as np
from scipy import stats
from astropy.io import fits
import astropy.stats 
import matplotlib.pyplot as plt 
import time
import scipy.signal
import os

image_file='./data/evt_T0.fits'
hdul = fits.open(image_file)
data = hdul[2].data
cols = hdul[2].columns
E1_min,E1_max=10,50#8.09,
E1range=np.linspace(E1_min,E1_max,5)
E2_min,E2_max=80,120#20,200,16(22.51),18(25.4613)
E2range=np.linspace(E2_min,E2_max,5)

BB=np.loadtxt('BB.txt')
Block=np.loadtxt('Block.txt')
#print(Block[1],Block[len(Block)-2])
#print(hdul[1].data)
data1 = data['TminusT0'][:]#events recorded with time T - T0 
#print(data1,len(data1))
data11 = data['PHA'][:] #PHA channel 
#print(data11,len(data11))
peaks=scipy.signal.find_peaks(BB[:,1])

UPBINSIZE=BB[peaks[0][0],2]
t11,t22=-10.,10.
t_total=t22-t11
print(E1range,E2range,UPBINSIZE)
bsize=[UPBINSIZE/2.,UPBINSIZE/3.]#,UPBINSIZE/4.,UPBINSIZE/5.]#np.linspace(UPBINSIZE/5.,UPBINSIZE,4)[::-1]
#print(bsize)

SNR=open('SNR.txt', 'a')
N_total_step=0
for Channel_min in E1range:
    for Channel_max in E2range:
        #bhere=(np.where(np.logical_and(data11>=Channel_min,np.logical_and(np.logical_and(data1>=t11,data1<=t22),data11<=Channel_max))))[0]
        #data2 =data1[bhere[0]:bhere[int(len(bhere)-1)]]
        ibin=0
        #select=[ ]
        #iselect=0
        for ibin in range(len(bsize)):
            #print(len(bsize),ibin)
            #i=peaks[0][0]
            #whh=hhh[i+1]-hhh[i]
            #countn=len(np.where(np.logical_and(data11>=int(Channel_min),np.logical_and(np.logical_and(data1>=hhh[i],data1<hhh[i+1]),data11<int(Channel_max))))[0])
            #vhh=countn/whh
            binwidth=bsize[ibin]
            ibin=ibin+1
            #print(Channel_min,Channel_max,binwidth)                
            a1=[t22+5.-(ia+0.5)*binwidth for ia in range(int((t_total+10.)/binwidth))]
            #a10=[BB[peaks[0][0],0]+(ia+0.5)*binwidth for ia in range(int((t22-BB[peaks[0][0],0])/binwidth))]
            ahere=np.where(np.logical_and(a1>=Block[1],a1<=Block[2]))[0] #np.logical_and(a1<=(BB[peaks[0][0],0]+UPBINSIZE/2.),a1>=(BB[peaks[0][0],0]-UPBINSIZE/2.)))[0]
            backhere=np.where(np.logical_or(a1<Block[1]-5.,a1>Block[2]+5.))[0]
            #print(backhere)
            #backhere=np.where(np.logical_or(a1<Block[1],a1>Block[len(Block)-2]))
            #print(backhere)
            #print(ahere)
            if len(ahere)==0:
                continue
            #print(ahere)  
            c1=[0 for ib in range(int((t_total+10.)/binwidth))]
            for ic in range(len(a1)):
                c1[ic]=len(np.where(np.logical_and(data11>=int(Channel_min),np.logical_and(np.logical_and(data1>=a1[ic]-0.5*binwidth,data1<a1[ic]+0.5*binwidth),data11<int(Channel_max))))[0])/binwidth
            #print(c1)
            c2=[c1[iccccc] for iccccc in backhere]
            #print(c2)
            sigma=np.sqrt(np.var(c2))#(c1[:backhere[len(backhere)-1]]))            #这里就是标准差与np.std()是一样的。
            ava=np.mean(c2)#(c1[:backhere[len(backhere)-1]])
            if len(ahere)<=1:
                pe=c1[ahere[0]]
            else:
                pe=max(c1[ahere[0]:ahere[len(ahere)-1]])
            nsi=(pe-ava)/sigma#(vhh-ava)/sigma                                     #这里是信噪比，（信号最高处 - 背景平均值）/ 背景标差 
            #select.append(nsi)
            #iselect=iselect+1
            #if len(select)>=2:
            #    if select[len(select)-1]<select[len(select)-2] and nsi<3:
            #        break
            print(Channel_min,Channel_max,binwidth,nsi)
            SNR.write(str(nsi))
            SNR.write(' ')
            SNR.write(str(Channel_min))
            SNR.write(' ')
            SNR.write(str(Channel_max))
            SNR.write(' ')
            SNR.write(str(binwidth))
            SNR.write('\n')
SNR.close()



