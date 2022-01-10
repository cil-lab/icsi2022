# ICSI'2022 测试集

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

**中文** | [English](README.md)

这是 ICSI'2022 测试集的官方仓库。我们鼓励所有研究人员在 ICSI'2022 测试集上测试算法。该测试集包括 10 个黑盒测试函数。同时，我们以此测试集为基础举办了一个算法竞赛，欢迎所有研究群体智能算法以及进化算法的研究人员参加这次比赛。我们会为这次比赛优秀的前几名提供丰厚的奖金作为奖励。

在下面的文件中我们会告诉大家该测试集如何安装与使用。

<!-- README 文件是人们通常最先看到的第一个东西。它应该告诉人们为什么要使用、如何安装、以及如何使用你的代码。README 文件标准化能够使得创建和维护 README 文件更加简单。毕竟，要写好一个文档不是那么容易的。 -->


本仓库包含以下内容：

1. 适用于python语言的测试集函数。
2. 适用于MATLAB语言的测试集函数。
3. 此测试集内各函数的定义。[ICSI_2022_Benchmark_functions.pdf](ICSI_2022_Benchmark_functions.pdf)

## 内容列表

- [ICSI'2022 测试集](#icsi2022-测试集)
	- [内容列表](#内容列表)
	- [背景](#背景)
	- [安装](#安装)
		- [MATLAB 支持](#matlab-支持)
	- [使用说明](#使用说明)
	- [函数定义](#函数定义)
	- [相关仓库](#相关仓库)
	- [维护者](#维护者)
	- [使用许可](#使用许可)

## 背景

ICSI 是国际群体智能会议的缩写，到目前为止已经成功举办了12届会议。该会议是群体智能领域内顶尖的会议之一。2022年第十三届国际群体智能会议（ICSI’2022）将于2022年7月于中国西安举行。 [网站](http://iasei.org/icsi2022/)

为了让各位研究群体智能优化算法与进化算法的研究者更好的测试算法，我们设计了一套测试集函数，为广大研究者提供一个更好的测试平台。同时，我们以此测试集为基础举办了一个算法竞赛，欢迎所有研究群体智能算法以及进化算法的研究人员参加这次比赛。我们会为这次比赛优秀的前几名提供丰厚的奖金作为奖励。

<!-- > 请记住：是文档而非代码，定义了一个模块的功能。 -->

<!-- —— [Ken Williams, Perl Hackers](http://mathforum.org/ken/perl_modules.html#document) -->

<!-- 写 README 从某种程度上来说相当不易，一直维护下去更是难能可贵。如果可以减少这个过程，则可以让写代码与修改代码更容易，使得是否在说明中指明一处需改有无必要更加清楚，你可以花费更少的时间来考虑是否你最初的文档是否需要更新，你可以分配更多的时间来写代码而非维护文档。 -->

<!-- 同时，标准化在某些别的地方也有好处。有了标准化，用户就可以花费更少的时间来搜索他们需要的信息，他们同时可以做一个工具来从描述中搜集信息，自动跑示例代码，检查授权协议等等。 -->

<!-- 这个仓库的目标是： -->


## 安装

这个项目提供 [python](https://www.python.org/) 和 [MATLAB](https://www.mathworks.com/products/matlab.html)支持。请确保你本地安装了它们。

### MATLAB 支持

对于使用MATLAB的用户，请使用MATLAB编译文件[icsi2022.cpp](MATLAB\icsi2022.cpp)，在MATLAB中运行

```MATLAB
mex icsi2022.cpp -DWINDOWS
```

如果可以正常运行[test.m](MATLAB/test.m)说明编译成功。
<!-- 然后就可以使用下面的命令使用该数据集。

```MATLAB
f = icsi2022(x, i)
``` -->

**注意：请先确保系统中已经安装MATLAB支持的编译器。关于如何在MATLAB中安装使用编译器请参考该[网站](https://ww2.mathworks.cn/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler)。**

### python 支持

首先确保系统中有以下依赖

```txt
python>=3.6
Cython==0.29.23
matplotlib==3.3.4
numpy==1.20.1
```

对于使用Linux的用户，可以在python文件夹中运行下面的命令来完成编译

```sh
$ bash c.sh
```

对于windows用户，请确保系统中有合适的编译器，推荐使用MinGW，或者Microsoft Visual C++ 14.0或更高版本，安装完成后首先在 python/lib 文件夹下运行下面的命令编译。

```sh
gcc -c icsi2022.c icsi2022.h
ar rcs libicsi2022.a icsi2022.o
```

然后在python文件夹下运行以下的命令进行编译。

```sh
python setup.py build_ext --inplace
```

如果可以正常运行[check.py](python/check.py)文件说明编译成功。


## 使用说明

安装完成后就可以调用简单的命令来使用本测试集函数。

MATLAB用户可以使用

```MATLAB
y = icsi2022(x, i)
% x is a m*D matrix
% i is from 1 to 10 which represent the serial number of 10 functions in icsi'2022 test suite.
% 
```

python用户可以使用

```python
import numpy as np
import icsi2022
# x is a m*D matrix

y=icsi2022.eval(x,i)
# i is from 0 to 9 which represent the serial number of 10 functions in icsi'2022 test suite.

```

## 函数定义

关于测试集中的函数具体定义，请参看文档[icsi2022.pdf](ICSI_2022_Benchmark_functions.pdf)

## 相关仓库

- [Fireworks Algorithm](https://github.com/cil-lab/fwaopt/tree/master/mpopt) — 烟花算法官方仓库。


## 维护者

[@moshangsang24](https://github.com/moshangsang24)

如果使用中有任何的问题，可以发送邮件到 [icsi2022@iasei.org](icsi2022@iasei.org)或者直接在github上提交问题。


## 使用许可

[MIT](LICENSE) © CIL-Lab
