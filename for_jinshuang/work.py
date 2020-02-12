import numpy as np
import os 
from astropy.io import fits
from astropy.table import Table
#-------------------
#路径
topdir = './data/'
hb_link = topdir + 'hb.fits'         #数据文件名
oiii_link = topdir + 'oiii.fits'     #数据文件名
savename = 'cccc.fits'               #输出文件名（保存文件名）
#--------------------------------

hb = fits.open(hb_link)
hb_namelist = hb[1].columns.names    #获取数据列的名称
#print(hb_namelist)
oiii = fits.open(oiii_link)
oiii_namelist = oiii[1].columns.names 
#-------------------------------------------
'''
筛选1：条件筛选
hb -> LINEAREA*2.355<=1000,LINEAREA>=0
oiii -> LINEAREA>=0
'''
hb_linearea = hb[1].data['LINEAREA']
hb_index = np.where((hb_linearea*2.355 <= 1000)&(hb_linearea>=0))[0] #筛选条件
hb_data = hb[1].data[hb_index]

oiii_linearea = oiii[1].data['LINEAREA']
oiii_index = np.where(oiii_linearea>=0)[0]      #筛选条件
oiii_data = oiii[1].data[oiii_index]
print('筛选1 完毕!!')
#--------------------------------------------
'''
三个条件相等
'''
print('筛选2 ...')
n_hb = len(hb_data)
n_oiii = len(oiii_data)
print('文件1 个数：',n_hb)
print('文件2 个数：',n_oiii)
#n_hb = 10
table  = []
hb_plate = hb_data['PLATE']
hb_mjd = hb_data['MJD']
hb_fiberid = hb_data['FIBERID']

for i in range(n_hb):
	#print(i,end = '\r')
	print('完成率 %.5f %%'%((i+1)/n_hb*100),end = '\r')
	index_x = np.where((oiii_data['PLATE'] == hb_plate[i])&(oiii_data['MJD'] == hb_mjd[i])&(oiii_data['FIBERID']==hb_fiberid[i]))[0]    #筛选条件
	good_append = oiii_data[index_x]
	for j in good_append:
		table.append(np.array(j))			
	#oiii_data = np.delete(oiii_data,index_x,axis=0)      #删除已筛选项，以降低查询次数。
	#print('len:',len(oiii_data[0]))
'''
写入文件
'''
table = np.vstack(table)
t = Table(table,names = oiii_namelist)
if os.path.exists(topdir+savename):
	os.remove(topdir+savename)
t.write(topdir+savename,format = 'fits')


