import numpy as np
from scipy.sparse import csc_matrix,eye,diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import datetime


def WhittakerSmooth(x,w,lambda_):
	'''

	:param x: array 输入数值，数组
	:param w: array 与数值对应的权重数组
	:param lambda_: 平滑参数
	:return: array 平滑结果
	'''
	X=np.matrix(x)#这里将数组转化为矩阵。矩阵之后就不可以用索引进行引用了。
	m=X.size
	#i=np.arange(0,m)
	E=eye(m,format='csc')
	D=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
	W=diags(w,0,shape=(m,m))
	A=csc_matrix(W+(lambda_*D.T*D))
	B=csc_matrix(W*X.T)
	background=spsolve(A,B)   #求解矩阵方程

	return np.array(background)


def baseline(spectra,lambda_,hwi,it,int_):
	spectra = np.array(spectra)
	spectra_index = np.arange(0,spectra.size,1)
	wl = np.ones(spectra.shape[0])
	spectra = WhittakerSmooth(spectra,wl,lambda_)

	if it != 1 :
		d1 = np.log10(hwi)
		d2 = 0
		w = np.ceil(np.concatenate((10**(d1+np.arange(0,it-1,1)*(d2-d1)/(np.floor(it)-1)),[d2])))
		w = np.array(w,dtype = int)
	else:
		w = np.array([hwi],dtype = int)
	#print(w)

	lims = np.linspace(0,spectra.size -1,int_+1)
	lefts = np.array(np.ceil(lims[:-1]),dtype = int)#这里指的是索引值
	rights = np.array(np.floor(lims[1:]),dtype = int)#同上
	minip = (lefts+rights)*0.5#索引
	xx = np.zeros(int_)
	for i in range(int_):#这里是一个rebin的过程,这里可以提速
		if rights[i]+1 < spectra.size -1:
			xx[i] = np.mean(spectra[lefts[i]:rights[i]+1])
		else:
			xx[i] = np.mean(spectra[lefts[i]:])
	'''
	mreg_n = 10 #int(spectra.size/int_)
	cutoff = int(spectra.size/mreg_n)*mreg_n
	spectra_A = spectra[:cutoff]
	spectra_index_A = spectra_index[:cutoff]
	xx_reshape = spectra_A.reshape(-1,mreg_n)
	xx = xx_reshape.mean(axis = 1)
	xx_index_reshape = spectra_index_A.reshape(-1,mreg_n)
	minip = xx_index_reshape.mean(axis = 1)
	if(spectra.size % mreg_n !=0):
		xx_b = spectra[cutoff:]
		xx_index_b = spectra_index[cutoff:]
		xx = np.concatenate((xx,[xx_b.mean()]))
		minip = np.concatenate((minip,[xx_index_b.mean()]))
	int_ = xx.size
	'''
	for i in range(it):
		# Current window width
		w0 = w[i]
		# Point-wise iteration to the right
		for j in range(1,int_-1):
			# Interval cut-off close to edges
			v = np.min([j,w0,int_-j-1])
			# Baseline suppression
			if j+v+1 < int_-1 :
				a = np.mean(xx[j-v:j+v+1])
			else:
				a = np.mean(xx[j-v:])
			xx[j] = np.min([a,xx[j]])
		for j in range(1,int_-1):
			k = int_-j-1
			v = np.min([j,w0,int_-j-1])
			if k + v+1 < int_-1:
				a = np.mean(xx[k-v:k+v+1])
			else:
				a = np.mean(xx[k-v:])
			xx[k] = np.min([a,xx[k]])
	minip = np.concatenate(([0],minip,[spectra.size-1]))
	xx = np.concatenate((xx[:1],xx,xx[-1:]))
	index = np.arange(0,spectra.size,1)
	xxx = np.interp(index,minip,xx)
	return xxx

def r_baseline(time,rate,lam = None,hwi = None,it = None,inti = None):
	dt = time[1]-time[0]

	if(lam is None):
		lam = 100/dt
	if(hwi is None):
		hwi = int(40/dt)
	if(it is None):
		it = 9
	if(inti is None):

		fillpeak_int = int(len(rate)/10)

	else:
		fillpeak_int =inti
	if(lam < 1):
		lam = 1
	bs = baseline(rate,lambda_=lam,hwi=hwi,it = it,int_ = fillpeak_int)
	return time,rate-bs,bs


topdir = '/home/laojin/trigdata/2017/'
#datalink = '/home/laojin/trigdata/2019/bn190114873/glg_tte_n7_bn190114873_v00.fit'
datalink = '/home/laojin/trigdata/2017/bn171010792/glg_tte_n8_bn171010792_v00.fit'
savedir = '/home/laojin/my_work/baseline/'


if(os.path.exists(savedir) == False):
	os.makedirs(savedir)
data = fits.open(datalink)
trigtime = data[0].header['TRIGTIME']


time = data[2].data.field(0)
ch = data[2].data.field(1)

t = time-trigtime

ch_n = 50
#t_index = np.where(ch == ch_n)
t_index = np.where((t>=-200)&(t<=500))[0]
t_ch = t[t_index]

binsize = 0.5
edges = np.arange(t_ch[0],t_ch[-1]+binsize,binsize)
bin_n,bin_edges = np.histogram(t_ch,bins=edges)

bin_c = (bin_edges[1:]+bin_edges[:-1])*0.5

rate = bin_n/(bin_edges[1:]-bin_edges[:-1])

print('---------------------------------------')
start = datetime.datetime.now()
t_c,cs,bs = r_baseline(bin_c,rate)
end = datetime.datetime.now()
print('---------------------------------------')

print('耗时：',end-start)
plt.figure()
plt.plot(bin_c,rate)
plt.plot(t_c,bs)
plt.ylim([0,2500])
plt.savefig(savedir + 'A_light_curve2.png')
plt.close()








