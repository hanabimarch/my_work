from GRB_utils import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from zjh_data_analysis import *
from Separate_background.background_kernel import *
from scipy.interpolate import interp1d
from scipy import sparse

#datalink = '/home/laojin/trigdata/2019/bn190114873/glg_tte_n7_bn190114873_v00.fit'
datalink = '/home/laojin/trigdata/2017/bn171010792/glg_tte_n8_bn171010792_v00.fit'
#datalink = '/home/laojin/trigdata/2017/bn170127634/glg_tte_n0_bn170127634_v00.fit'
#datalink = '/home/laojin/trigdata/2017/bn170921168/glg_tte_n2_bn170921168_v00.fit'
#datalink = '/home/laojin/trigdata/2017/bn170429799/glg_tte_nb_bn170429799_v00.fit'
savedir = '/home/laojin/shiyan/LB3/'

if(os.path.exists(savedir) == False):
	os.makedirs(savedir)
data = fits.open(datalink)
trigtime = data[0].header['TRIGTIME']

time = data[2].data.field(0)
ch = data[2].data.field(1)

t = time-trigtime


ch_top = 90
ch_buttum = 8

t_index = np.where((ch>=ch_buttum) & (ch<=ch_top))

t_ch = t[t_index]
t_ch = t_ch[np.where((t_ch>=-150)&(t_ch<=400))]

binsize = 0.5

edges = np.arange(t_ch[0],t_ch[-1]+binsize,binsize)

bin_n,bin_edges = np.histogram(t_ch,bins=edges)
bin_n = bin_n/binsize

bin_c = (bin_edges[1:]+bin_edges[:-1])*0.5

bin_dn = bin_n[1:]-bin_n[:-1]
bin_dn = np.concatenate((bin_dn[:1],bin_dn))

bin_ddn = bin_dn[1:] - bin_dn[:-1]
bin_ddn = np.concatenate((bin_ddn[:1],bin_ddn))

dn_std = np.std(bin_dn)
ddn_std = np.std(bin_ddn)

plt.figure(figsize = (10,10))
plt.subplot(3,1,1)
plt.title('light curve')
plt.plot(bin_c,bin_n)


plt.subplot(3,1,2)
plt.title('dn')
plt.plot(bin_c,bin_dn)
plt.axhline(y = -dn_std,color = 'g')
plt.axhline(y = dn_std,color = 'g')

plt.subplot(3,1,3)
plt.title('ddn')
plt.plot(bin_c,bin_ddn)
plt.axhline(y = -ddn_std,color = 'g')
plt.axhline(y = ddn_std,color = 'g')


plt.savefig(savedir+ 'A_LB.png')
plt.close()


index_block = np.arange(0,50,1)

block_para = []
block_at = []
block_at1 = []
block_w = np.array([0,1.  ,1.  ,1.  ,1  ,1  ,1.  ])
block_w1 = np.array([1,1,0])

blocks_index = []

step_size =1
for i in range(int(len(bin_c)/step_size)):
	if(index_block[-1]+step_size*i >= len(bin_c)):
		break
	one_block_index = index_block+step_size*i
	blocks_index.append(one_block_index)
	para_t = np.mean(bin_c[index_block+step_size*i])
	para_n = np.mean(bin_n[index_block+step_size*i])
	para_n_std = np.var(bin_n[index_block+step_size*i])
	para_dn = np.mean(bin_dn[index_block+step_size*i])
	para_dn_std = np.std(bin_dn[index_block+step_size*i])
	para_ddn = np.mean(bin_ddn[index_block + step_size*i])
	para_ddn_std = np.std(bin_ddn[index_block + step_size*i])
	lists = [
		para_t,para_n,para_n_std,para_dn,para_dn_std,para_ddn,para_ddn_std
		]
	block_para.append(lists)
	lists = np.array(lists)

	lists1 = np.array([para_n_std*para_n  ,np.abs(para_ddn_std*para_dn_std)  ,np.abs(para_ddn*para_ddn_std)])

	ww2 = np.sum(lists1*block_w1)
	ww = np.sum(lists*block_w)
	block_at.append(ww)
	block_at1.append(ww2)
block_para = np.array(block_para).T

www2 = 1*block_para[1]*block_para[2]/np.sqrt((block_para[1]**2).sum()*(block_para[2]**2).sum()) + 0*block_para[4]*block_para[6]/np.sqrt((block_para[4]**2).sum()*(block_para[6]**2).sum())
ppppp = np.sqrt((block_para[1]**2).sum()*(block_para[2]**2).sum())

plt.figure(figsize = (10,10))
plt.subplot(3,1,1)
plt.title('light curve')
plt.plot(bin_c,bin_n,alpha = 0.5)
plt.plot(block_para[0],block_para[1])
plt.plot(block_para[0],block_para[2])

plt.subplot(3,1,2)
plt.title('dn')
plt.plot(bin_c,bin_dn,alpha = 0.5)
#plt.axhline(y = -dn_std,color = 'g')
#plt.axhline(y = dn_std,color = 'g')
plt.plot(block_para[0],block_para[3])
plt.plot(block_para[0],block_para[4])
plt.subplot(3,1,3)
plt.title('ddn')
plt.plot(bin_c,bin_ddn,alpha = 0.5)
#plt.axhline(y = -ddn_std,color = 'g')
#plt.axhline(y = ddn_std,color = 'g')
plt.plot(block_para[0],block_para[5])
plt.plot(block_para[0],block_para[6])

plt.savefig(savedir+ 'A_LB_block.png')
plt.close()

plt.figure(figsize=(10,10))
plt.subplot(3,2,1)
plt.plot(block_para[1],block_para[3],'o')
plt.xlabel('n')
plt.ylabel('dn')
plt.subplot(3,2,2)
plt.plot(block_para[3],block_para[5],'o')
plt.xlabel('dn')
plt.ylabel('ddn')
plt.subplot(3,2,3)
plt.plot(block_para[5],block_para[1],'o')
plt.xlabel('ddn')
plt.ylabel('n')

plt.subplot(3,2,4)
plt.plot(block_para[1+1],block_para[3+1],'o')
plt.xlabel('n std')
plt.ylabel('dn std')
plt.subplot(3,2,5)
plt.plot(block_para[3+1],block_para[5+1],'o')
plt.xlabel('dn std')
plt.ylabel('ddn std')
plt.subplot(3,2,6)
plt.plot(block_para[5+1],block_para[1+1],'o')
plt.xlabel('ddn std')
plt.ylabel('n std')

plt.savefig(savedir + 'A_LB_s.png')
plt.close()


block_b = AirPLS(block_at1)

block_bs = block_b.bottom_airPLS()
w = np.ones(len(block_at1))
block_at1_w = WhittakerSmooth(np.array(block_at1),w,100)
plt.figure(figsize=(10,10))
#plt.axhline(y = 220,color = 'g')
#plt.axhline(y = 300,color = 'r')
plt.plot(block_para[0],block_at)
plt.plot(block_para[0],block_at1)
plt.plot(block_para[0],block_at1_w)
plt.plot(block_para[0],block_bs)
plt.ylim([-8,80])
plt.savefig(savedir + 'A_at.png')
plt.close()
block_at1_b = np.array(block_at1)-block_bs

plt.figure(figsize = (10,10))
plt.plot(block_para[0],block_at1_b)
#plt.ylim([-8,40])
plt.savefig(savedir + 'A_block_at1_b.png')
plt.close()

max_block_at1_b = np.max(block_at1_b)
#if(max_block_at1_b > 1e5):
#	max_block_at1_b = 100000
#block_at1_b_edges = np.arange(0,max_block_at1_b,1)
block_at1_b_edges = np.linspace(0,max_block_at1_b,2000)
pn,pe = np.histogram(block_at1_b,bins = block_at1_b_edges)
pn = np.concatenate((pn[:1],pn))

zzzz = np.percentile(block_at1_b,50)
print(zzzz)
#P. Eilers和H. Boelens在2005年有一种名为“Asymmetric Least Squares Smoothing”的算法
def baseline_als(y,lam,p,niter=100,dt = 1):
	L = len(y)
	D = sparse.csc_matrix(np.diff(np.eye(L), 2))
	w = np.ones(L)
	z_list = []
	dz_max = []
	dz_it = []
	dz_max_old = 0
	for i in range(niter):
		W = sparse.spdiags(w, 0, L, L)
		Z = W + lam * D.dot(D.transpose())
		z = spsolve(Z, w * y)
		dz = z[1:] - z[:-1]
		dz_max1 = np.max(np.abs(dz))/dt
		dz_iti = dz_max1 -dz_max_old
		dz_max_old = dz_max1
		dz_it.append(dz_iti)
		print(dz_max1,dz_iti)
		z_list.append(z)
		dz_max.append(dz_max1)
		if(dz_max1<=0.03):
			print('dz good')
			break
		w = p * (y > z) + (1 - p) * (y < z)
		#w[0] = (1 - p)
		#w[-1] = (1 - p)
		#print('w',w)
		if(len(w[w>0.5]) < 3):
			print('good')
			break
		if(dz_iti == 0):
			print('dz_iti = 0')
			z_list = z_list[:-2]
			dz_max = dz_max[:-2]
			break
	z_list = np.array(z_list)
	dz_max = np.array(dz_max)
	if(len(dz_max[dz_max<=0.03]) != 0):
		z_list = z_list[dz_max<=0.03]
		dz_max = dz_max[dz_max<=0.03]
	elif(len(dz_max[dz_max<=0.08]) != 0):
		z_list = z_list[dz_max <= 0.08]
		dz_max = dz_max[dz_max <= 0.08]
	index_ = np.argmin(dz_max)
	print(dz_max[index_])
	return z_list[index_]

plt.plot(pe,pn,linestyle = 'steps')
#plt.xlim([0,100000])
plt.savefig(savedir + 'B_block_at1_b.png')
plt.close()

#block_bw = Baseline_in_time(block_para[0],www2,fitness='double')
#bs_w = block_bw.bs
#cs_w = www2 - bs_w
#bs_w = baseline_als(www2*ppppp,4,1e-99,dt=step_size)/ppppp

tt,cs_w1,bs_w1 = r_baseline(block_para[0],www2,lam=1,it = 8)
#bs_w = baseline_als(cs_w11*ppppp,4,0,dt=step_size*binsize)/ppppp
#cs_w = cs_w11-bs_w
print('$$$$')
#block_bw1 = Baseline_in_time(block_para[0],cs_w,fitness='bottom')
#bs_w1 = block_bw1.bs
print(bs_w1)
#cs_w1 = cs_w-bs_w1
#print(cs_w1)

plt.plot(block_para[0],www2,label = 'www2')
#plt.plot(block_para[0],bs_w,label = 'bs_w')
#plt.plot(block_para[0],bs_w11,label = 'bs_w11')
#plt.plot(block_para[0],cs_w,label = 'cs_w')
plt.plot(block_para[0],bs_w1,label ='bs_w1')
#plt.plot(block_para[0],cs_w11,label = 'cs_w11')
plt.plot(block_para[0],cs_w1,label = 'cs_w1')
#plt.axhline(y = ee)
plt.legend()
#plt.ylim([-0.00001,0.0005])
plt.savefig(savedir + 'B_block_www2_a.png')
plt.close()
zzzzw = np.percentile(cs_w1,30)

max_cs_w1 = np.max(cs_w1)
print('zzzzw:',zzzzw)
#if(max_block_at1_b > 1e5):
#	max_block_at1_b = 100000
block_at1_b_edges = np.arange(np.min(cs_w1),max_cs_w1+zzzzw,zzzzw)
print(block_at1_b_edges)
#block_at1_b_edges = np.linspace(0,max_block_at1_b,2000)
pn1,pe = np.histogram(cs_w1,bins = block_at1_b_edges)
pe_c = (pe[1:]+pe[:-1])*0.5
pe_c = np.concatenate((pe[:1],pe_c))

pn = np.concatenate((pn1[:1],pn1))
if(len(pe[pn==0])>0):
	es = pe[pn==0][0]
else:
	es = pe_c[-1]
new_pe = np.linspace(0,es,1000)

#new_pn = np.interp(new_pe,pe_c,pn)

new_pn = interp1d(pe_c,pn,kind='linear')(new_pe)#'cubic','linear'


plt.plot(pe,pn,linestyle = 'steps')
plt.plot(new_pe,new_pn)
#plt.xlim([0,0.0005])
plt.savefig(savedir + 'C_block_at1_b.png')
plt.close()


ee = new_pe[new_pn<np.max(new_pn)*0.5][0]
print('ee:',ee)

nornallization = np.zeros(len(bin_c))
for index,value in enumerate(cs_w1):
	if(value > ee):

		nornallization[blocks_index[index]] = nornallization[blocks_index[index]] + 1



plt.plot(block_para[0],www2,label = 'www2')
#plt.plot(block_para[0],bs_w,label = 'bs_w')
#plt.plot(block_para[0],bs_w11,label = 'bs_w11')
#plt.plot(block_para[0],cs_w,label = 'cs_w')
plt.plot(block_para[0],bs_w1,label ='bs_w1')
#plt.plot(block_para[0],cs_w11,label = 'cs_w11')
plt.plot(block_para[0],cs_w1,label = 'cs_w1')
plt.axhline(y = ee)
plt.legend()
#plt.ylim([-0.00001,0.000002])
plt.savefig(savedir + 'B_block_www2_b.png')
plt.close()
#print(pe[2])

plt.plot(bin_c,nornallization)
plt.savefig(savedir + 'D_nornallization.png')
plt.close()






