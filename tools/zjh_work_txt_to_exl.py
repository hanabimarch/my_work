import xlwt
import numpy as np
import zzh_py3_file as zhf
file = xlwt.Workbook()
table = file.add_sheet('sheet sample')
'''
topdir = '/home/laojin/results/results_GBM/'
file_all = 'Y_all_sample_name.txt'
file_1 = 'Y_sample_name_1.txt'
file_2 = 'Y_sample_name_2.txt'
file_3 = 'Y_sample_name_3.txt'
file_x = 'Y_sample_name_x.txt'

sample_all = zhf.readcol(topdir+file_all)[0]
sample_1 = zhf.readcol(topdir+file_1)[0]
sample_2 = zhf.readcol(topdir+file_2)[0]
sample_3 = zhf.readcol(topdir+file_3)[0]
sample_x = zhf.readcol(topdir+file_x)[0]


num = len(sample_all)
#输出表头
table.write(0,0,'trigger name')
table.write(0,1,'sign')
table.write(0,2,'trigger name')
table.write(0,3,'sign')
table.write(0,4,'trigger name')
table.write(0,5,'sign')
for index,name in enumerate(sample_all):
	l = int(index/3)
	h = index%3
	print(l,h)
	table.write(l+1,2*h,name)
	if(name in sample_1):
		table.write(l+1,2*h+1,'1')
	elif(name in sample_2):
		table.write(l+1,2*h+1,'2')
	elif(name in sample_3):
		table.write(l+1,2*h+1,'3')
	elif(name in sample_x):
		table.write(l+1,2*h+1,'x')
	else:
		table.write(l+1,2*h+1,'0')

file.save(topdir +'sample_all.xls')
'''
topdir = '/home/laojin/results/results_GBM/background/'
file_1 = 'b_sample.txt'
file_2 = 'b_sample2.txt'

sample_1,detector_1 = zhf.readcol(topdir+file_1)
sample_2,detector_2 = zhf.readcol(topdir+file_2)

sample = sample_1+sample_2
detector = detector_1+detector_2
table.write(0,0,'trigger name')
table.write(0,1,'detector')
table.write(0,2,'trigger name')
table.write(0,3,'detector')
table.write(0,4,'trigger name')
table.write(0,5,'detector')
for index,name in enumerate(sample):
	l = int(index/3)
	h = index%3
	table.write(l+1,2*h,name)
	table.write(l+1,2*h+1,detector[index])
file.save(topdir +'b_sample_all.xls')

