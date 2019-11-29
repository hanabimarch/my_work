import numpy as np
import matplotlib.pyplot as plt
import os


datalink = './data.txt'

x,y = np.loadtxt(datalink,dtype = np.float32).T

print(x,y)

savedir = './results/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)


#------------------------------------------
'''
绘制数据散点图
'''
plt.plot(x,y,'o',color = 'k')
plt.xlabel('production')
plt.ylabel('energy consumption')
plt.savefig(savedir+'A_data_scatter.png')
plt.close()
#------------------------------------------
#------------------------------------------
'''
最小二乘法拟合数据
公式推到见书 P120
拟合公式 y` = k*x + b
'''
k = (np.sum(x*y)-x.sum()*y.mean())/(np.sum(x**2) - x.size*x.mean()**2)
b = y.mean()-k*x.mean()
print('------------------------------------------------')
print('拟合公式：y` = k*x + b')
print('拟合的参数为：')
print('K = ',k)
print('b = ',b)
print('------------------------------------------------')


x_ = np.linspace(x.min(),x.max(),100)
y_ = k*x_ + b

plt.figure()
plt.subplot(1,1,1)
plt.plot(x,y,'o',color = 'k',label = 'samples')
plt.plot(x_,y_,'--',color = 'r',label = 'fitting function')
plt.xlabel('production')
plt.ylabel('energy consumption')
plt.legend()
plt.savefig(savedir + 'B_fitting_function.png')
plt.close()

plt.title('dispersion')
plt.plot(x,y-k*x-b,'o',color = 'k')
plt.axhline(y = 0,color = 'r')
plt.xlabel('production')
plt.ylabel('dispersion')
plt.ylim([-1,1])
plt.savefig(savedir + 'C_dispersion.png')
plt.close()

#-----------------------------------------------
'''
利用改进后的拟合曲线估算生成100吨甲产品需要消耗多少吨煤。在用改进前的消耗记录进行比较从而得到降低值。
零赶紧后消耗的量为p1
'''
p1 = 100*k+b
dp = 90-p1
print('-------------------------------------------------------')
if(dp >= 0):
	print('改进后消耗的煤炭量比改进前减少了'+str(dp)+'吨。')
else:
	print('改进后消耗的煤炭量比改进前并没有减少！')
print('-------------------------------------------------------')



















