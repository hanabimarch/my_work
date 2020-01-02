
from astropy.io import fits
from Separate_background import Separate_background
from zjh_data_analysis import *
import os
from astropy.stats import bayesian_blocks
import numpy as np
import matplotlib.pyplot as plt
from Baseline import TD_baseline

datadir = '/home/laojin/trigdata/2008/bn081224887/'
filename = 'glg_tte_n6_bn081224887_v01.fit'
savedir = './results/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)

htl = fits.open(datadir+filename)
trigtime = htl[0].header['TRIGTIME']

time = htl[2].data.field(0)
ch = htl[2].data.field(1)
t = time - trigtime

ch_n = htl[1].data.field(0)
e1 = htl[1].data.field(1)
e2 = htl[1].data.field(2)

t = t[np.where((ch>=3) & (ch <=109))]
ch = ch[np.where((ch>=3) & (ch <=109))]

dt = 0.064
edges = np.arange(-15,60+dt,dt)

bin_n,bin_edges = np.histogram(t,bins = edges)
bin_c = (bin_edges[1:] + bin_edges[:-1])*0.5
bin_rate  = bin_n/dt


t_c,cs,bs = TD_baseline(bin_c,bin_rate)
baye_edges = bayesian_blocks(t_c,np.round((cs+bs.mean())*dt),fitness='events',gamma = np.exp(-20))



#首先画出光变曲线和bayesian blocks的边界。
ch_standed = 25

t_h = t[np.where(ch>=ch_standed)[0]]
t_s_ = t[np.where(ch<ch_standed)[0]]

bin_n_h,bin_h_edges = np.histogram(t_h,bins = baye_edges)
bin_h_size = bin_h_edges[1:]-bin_h_edges[:-1]
bin_h_c = (bin_h_edges[1:]+bin_h_edges[:-1])*0.5
bin_h_rate = bin_n_h/bin_h_size

bin_n_s,bin_s_edges = np.histogram(t_s_,bins = baye_edges)
bin_s_size = bin_s_edges[1:] - bin_s_edges[:-1]
bin_s_c = (bin_s_edges[1:]+bin_s_edges[:-1])*0.5
bin_s_rate = bin_n_s/bin_h_size

HDS = bin_h_rate/bin_s_rate

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(t_c,cs,color = 'k',label = 'light curve')
for i in baye_edges:
	plt.axvline(x = i,color = 'g')
ax2 = ax1.twinx()
ax2.plot(bin_s_c,HDS,'o',color = 'r',)

plt.xlim([-15,60])
plt.legend()
plt.savefig(savedir+'A_light_curve.png')
plt.close()


#没有扣背景的情况下的硬度比变化：

Sample = Separate_background(t,ch,ch_n,time_range=[-15,60])#工具包里的内容

t_s,s_e = get_all_energy(Sample.s,Sample.s_ch,ch_n,e1,e2)#来自工具
t_b,b_e = get_all_energy(Sample.b,Sample.b_ch,ch_n,e1,e2)#来自工具
ch_standed2 = 49
e_standed2 = np.sqrt(e1[ch_standed2]*e2[ch_standed2])

plt.figure(figsize = (10,15))
plt.subplot(3,1,1)
plt.plot(t_s,s_e,',',color = 'k')
plt.plot(t_b,b_e,',',color = 'k')
plt.yscale('log')
plt.ylabel('Energy kev')
plt.xlim([-15,60])


plt.subplot(3,1,2)
plt.plot(t_s,s_e,',',color = 'k')
plt.axhline(y = e_standed2,color = 'r')
plt.yscale('log')
plt.ylabel('Energy kev')
plt.xlim([-15,60])

plt.subplot(3,1,3)
plt.plot(t_b,b_e,',',color = 'k')
plt.yscale('log')
plt.ylabel('Energy kev')
plt.xlim([-15,60])
plt.xlabel('time s')
plt.savefig(savedir + 'B_count_map.png')
plt.close()

bin_n_b,bin_b_edges = np.histogram(Sample.b,bins = edges)
bin_rate_b = bin_n_b/dt
bin_b_c = (bin_b_edges[1:]+bin_b_edges[:-1])*0.5

plt.figure()
plt.subplot(1,1,1)
plt.plot(bin_c,bin_rate)
plt.plot(bin_b_c,bin_rate_b)
plt.plot(t_c,bs)
plt.xlim([-15,60])
plt.xlabel('time s')
plt.ylabel('rate')
plt.savefig(savedir + 'C_light_curve.png')
plt.close()


t_s_h = Sample.s[np.where(Sample.s_ch>= ch_standed2)[0]]
t_s_s = Sample.s[np.where(Sample.s_ch< ch_standed2)[0]]


bin_n_h,bin_h_edges = np.histogram(t_s_h,bins = baye_edges)
bin_h_size = bin_h_edges[1:]-bin_h_edges[:-1]
bin_h_c = (bin_h_edges[1:]+bin_h_edges[:-1])*0.5
bin_h_rate = bin_n_h/bin_h_size

bin_n_s,bin_s_edges = np.histogram(t_s_s,bins = baye_edges)
bin_s_size = bin_s_edges[1:] - bin_s_edges[:-1]
bin_s_c = (bin_s_edges[1:]+bin_s_edges[:-1])*0.5
bin_s_rate = bin_n_s/bin_h_size

HDS = bin_h_rate/bin_s_rate
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(t_c,cs,color = 'k')
for i in baye_edges:
	plt.axvline(x = i,color = 'g')
ax2 = ax1.twinx()
ax2.plot(bin_s_c,HDS,'o',color = 'r',)

plt.xlim([-15,60])
plt.savefig(savedir+'D_light_curve.png')
plt.close()


ch_standed3 = 23

t_b_h = Sample.b[np.where(Sample.b_ch>= ch_standed3)[0]]
t_b_s = Sample.b[np.where(Sample.b_ch< ch_standed3)[0]]


bin_n_h,bin_h_edges = np.histogram(t_b_h,bins = baye_edges)
bin_h_size = bin_h_edges[1:]-bin_h_edges[:-1]
bin_h_c = (bin_h_edges[1:]+bin_h_edges[:-1])*0.5
bin_h_rate = bin_n_h/bin_h_size

bin_n_s,bin_s_edges = np.histogram(t_b_s,bins = baye_edges)
bin_s_size = bin_s_edges[1:] - bin_s_edges[:-1]
bin_s_c = (bin_s_edges[1:]+bin_s_edges[:-1])*0.5
bin_s_rate = bin_n_s/bin_h_size

HDS = bin_h_rate/bin_s_rate
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.plot(bin_b_c,bin_rate_b,color = 'k')
for i in baye_edges:
	plt.axvline(x = i,color = 'g')
ax2 = ax1.twinx()
ax2.plot(bin_s_c,HDS,'o',color = 'r',)

plt.xlim([-15,60])
plt.savefig(savedir+'E_light_curve.png')
plt.close()








