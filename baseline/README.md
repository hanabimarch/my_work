# R baseline python 化

## baseline算法的基本步骤

1). 数据平滑：数据平滑的作用是减少由于统计因素造成的数据局域过低现象。平滑线反应数据的主体走势，忽略了统计过程中产生的信号波动。这样有助于使采样样本反应数据主体走势，降低了统计条件对算法结果的影响。平滑参数为：lambda

2). 数据采样：对平滑之后的数据进行均匀的二次采样，方式为采样区间中的最小值（导致结果偏低的原因）。采样数由参数 int 控制。之后的处理均采用该采样。

3). 抑制迭代：迭代次数由参数 it 设置。迭代方式为可重叠窗口的平均值与原值取最小的方法。窗口宽度随迭代次数依次以指数降低。窗口初始宽度用参数 hwi 设置。

4). 数据还原：通过插值法将采样结果还原。

## baseline应用于时域分析的适应性改进

### zjh_baseline.py ：python 化 baseline 算法-改进 baseline 算法拟合偏低

zjh_baseline.py 是参考了baseline的R代码编写的相同功能的代码，主要用于实现baseline算法的基本功能。但是不同的是我们采用了whitttaker平滑器取代了R代码中使用的平滑器。该平滑器的更改不影响最终效果。由于更改了平滑器，平滑参数的设置有所更改，平滑参数越大，数据越平滑。

为解决baseline算法估计背景基线偏低的问题，在zjh_baseline.py代码的第31行，将最小采样变为平均值采样。这样可以有效减少背景基线偏低的问题。

### zjh_baseline2.py ：python 代码优化--提升计算速度

降低第一版本算法在迭代过程中创建索引数组的次数，通过切片方式进行采样。numpy数组求平均值时利用数组自带的求平均值函数（arr.mean() 要比 numpy.mean(arr) 快），对列表求最小值时利用python自带的函数（min(list) 要比 numpy.min(list) 快）。相同数据相同迭代次数，比zjh_baseline.py提速在160%左右。比R提速70% 。

下图是python代码和R代码运行耗时的比较：

![figure 0](https://github.com/zoujinhang/my_work/blob/master/baseline/1-1.png)

### zjh_baseline3.py ：解决时域数据边缘效应对 baseline 算法的影响。

利用whitttaker平滑器的性质，通过一些统计策略忽略数据边缘下降的情况。具体操作集合在get_smooth中。集体方法通过3次迭代，依次排除 1 sigma，2 sigma，3 sigma 的数据。虽然平滑的最终结果会忽视掉一些数据比较极端的凸起或凹陷，但是将贴合背景主体，在时域背景的估算中，这是合理且必要的。

同时我们对参数进行了优化，通过多次尝试，我们发现迭代次数为5时算法的时域适应性最好，大于5时容易引起背景基线偏低（但是并不会偏得太多），小于5时容易导致背景基线不平坦。

## baseline 例子

1). zjh_baseline.py代码运行结果。蓝色线为光变曲线，黄色线为拟合的背景。

![figure 1](https://github.com/zoujinhang/my_work/blob/master/baseline/A_light_curve.png)


2). zjh_baseline2.py代码运行结果。内容同上。可以看到蓝色光变曲线在结尾处存在数据的边缘效应，边缘效应导致了光变曲线的突然下降。而这影响到了结尾处背景的拟合，使背景基线偏低。

![figure 2](https://github.com/zoujinhang/my_work/blob/master/baseline/A_light_curve2.png)


3). zjh_baseline3.py代码运行结果。内容同上。可以看到蓝色的光变曲线在结尾处由于边缘效应导致了一次急速的下降。在算法改进之前，该效应会拉低背景基线从而导致背景基线脱离背景。经过改进后算法可以有效的忽略该效应的影响。

![figure 3](https://github.com/zoujinhang/my_work/blob/master/baseline/A_light_curve3.png)

























