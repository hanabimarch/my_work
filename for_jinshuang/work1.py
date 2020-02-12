import numpy as np
import os 
from astropy.io import fits
from astropy.table import Table
import pandas as pd

#-------------------
#路径
topdir = './data/'
hb_link = topdir + 'hb.fits'         #数据文件名
oiii_link = topdir + 'oiii.fits'     #数据文件名
savename = 'cccc.fits'               #输出文件名（保存文件名）
#--------------------------------
'''
筛选1：条件筛选
hb -> LINEAREA*2.355<=1000,LINEAREA>=0
oiii -> LINEAREA>=0
'''

hb = fits.open(hb_link)
hb_formats = hb[1].columns.formats
hb_tbl = Table.read(hb_link,hdu = 1)
hb_df = hb_tbl.to_pandas()
hb_df = hb_df[(hb_df['LINEAREA']>=0)&(hb_df['LINEAREA']*2.355<=1000)]
print(hb_df)

oiii = fits.open(oiii_link)
oiii_formats = oiii[1].columns.formats
oiii_tbl = Table.read(oiii_link,hdu = 1)
oiii_df = oiii_tbl.to_pandas()
oiii_df = oiii_df[oiii_df['LINEAREA']>=0]
print(oiii_df)
new_formats = []
for i in hb_formats:
	new_formats.append(i)
for i in oiii_formats[3:]:
	new_formats.append(i)
#-------------------------------------------
print('筛选1 完毕!!')
#--------------------------------------------
'''
三个条件相等
'''
print('筛选2 ...')

new_data = pd.merge(hb_df,oiii_df,on = ['PLATE','MJD','FIBERID'])#取交集
print(new_data)
new_data_name = new_data.columns.values

#写入文件
#data_in_array = new_data.values
#print(data_in_array)
#data_in_array = data_in_array.T

table = []
for i in range(len(new_data_name)):
	print(new_data[new_data_name[i]].values)
	table.append(fits.Column(name =new_data_name[i],array = np.array(new_data[new_data_name[i]].values),format =new_formats[i]))
hdulist = fits.BinTableHDU.from_columns(table)
if os.path.exists(topdir+savename):
	os.remove(topdir+savename)
hdulist.writeto(topdir+savename)
'''
t = Table(new_data,names = new_data_name)
print(t)

t.write(topdir+savename)
'''
print('完成！！')

