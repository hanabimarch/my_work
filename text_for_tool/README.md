# Calculate_duration

最新例子脚本为:

```
zjh_exmple_Calculate_duration.py 
```

##旧版弃用

之前版本的工具包更名如下：

```
Calculate_duration_old
```

由于旧版的核心并不稳定，因此暂时弃用旧版。


##以Bayesian blocks算法的新核心

这里我们采用比较简单的逻辑来进行背景区域的识别。

首先通过Bayesian blocks和baseline算法联合计算出信号的信噪比，之后通过设定阈值来判定blocks是否为背景区域。block的信噪比大于阈值则判定为爆发区间，信噪比小于阈值的区域判定为背景区域。之后采用平滑器拟合选定的背景区域。平滑器拟合背景的好处是不需要背景模型，它比多项式拟合背景的适应性更好。平滑器拟合背景会引入一个参数：平滑硬度。平滑硬度越大，拟合的背景越平滑。

Bayesian blocks和baseline算法联合计算信噪比时引入了一个描述背景浮动程度的先验参数，它描述背景blocks的变动程度。算法会根据这个参数寻找信号中符合参数描述信号区域，之后通过该区域计算信噪比。

算法的灵敏度由Bayesian blocks算法的灵敏度决定。


