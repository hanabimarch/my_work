#!/user/bin/python3

import matplotlib.pyplot as plt
import numpy as np
from zjh_get_daily_data import get_daily_data
import zzh_py3_file as zhf
from zjh_gbm_time import GBMtime

def finding(t,valu):
	d_time = []
	n = len(t)
	for i in range(1,n-1):
		if((valu[i] == 0)&(valu[i-1] != 0)&(valu[i+1] == 0)):			
			d_time.append(t[i])			
		if((valu[i] == 0)&(valu[i-1] == 0)&(valu[i+1] != 0)):
			d_time.append(t[i])

	return d_time
def zhifang(data,start,stop,bin):
	met_start = start
	met_stop = stop
	num = round((met_stop - met_start)/bin+1)
	bin_array =met_start + np.arange(num)*bin
	bin_n,bin_edges = np.histogram(data,bin_array)
	bin_c = (bin_edges[1:]+bin_edges[:-1])/2
	return bin_c,bin_n
def preice(time,detector):

	timestart = GBMtime.met_to_utc(time-60).fits

	timestop = GBMtime.met_to_utc(time+60).fits

	t = [timestart,timestop]
	t0,ch,ch_n,e1,e2 = get_daily_data(t,detector)

	bin_c,bin_n = zhifang(t0,time-60,time+60,1)

	d_time = finding(bin_c,bin_n)
	return d_time
def read():
	time = ['2018-07-25T00:00:00','2018-07-25T23:59:59']
	t,ch,ch_n,e1,e2 = get_daily_data(time,'n0')
	bin = 10
	bin_c,bin_n = zhifang(t,t.min(),t.max(),bin)
	return bin_c,bin_n
detectors = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
bin_c,bin_n = read()

timelist = finding(bin_c,bin_n)
#print(timelist)
resulttime = []
for i in timelist:
	sum0 = 0
	#print(i)
	for j in detectors:
		d_time = preice(i,j)

		d_time = d_time[len(d_time)-1]
		sum0 = sum0+d_time
	resulttime.append(GBMtime.met_to_utc(sum0/len(detectors)).fits)








fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.plot(bin_c,bin_n,',')
ax.set_xlabel('Time')
ax.set_ylabel('Number of Photo')
fig.savefig('/home/laojin/shiyan/for_all_time.png')

zhf.printdatatofile('/home/laojin/shiyan/result_time.txt',[resulttime])








