import zzh_py3_file as zhf
from zjh_fermi_daily_geometry import Fermi_daily_geometry
import numpy as np

timelist = zhf.readcol('/home/laojin/shiyan/result_time.txt')

posistion = []

for i in timelist[0]:
	print(i)
	qsj,pos,sc = Fermi_daily_geometry(i).find_right_list()
	if(sc[1]>180):
		sc[1] = sc[1]-360
	posistion.append(sc)
posistion = np.array(posistion).T
der = posistion[0]
ra = posistion[1]

#index = np.argsort(ra)
#ra = ra[index]
#der = der[index]

zhf.printdatatofile('/home/laojin/shiyan/SAA_point.txt',[ra,der])








