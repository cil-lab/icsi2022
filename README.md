# ICSI'2022 Benchmark Functions

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

This is the official repository of the icsi'2022 benchmark functions. We encourage all researchers to test their algorithms on the ICSI'2022 test suite.  At the same time, we are holding an algorithm competition based on this test set. All researchers who study swarm intelligence algorithms and evolutionary algorithms are welcome to participate in this competition. We will provide generous bonuses as rewards for the excellent contestants in this competition.

Next we will tell you how to install and use the test set.

This repository contains:

1. Test functions for python.
2. Test functions for MATLAB.
3. The definition of each function in this test suite. [ICSI_2022_Benchmark_functions.pdf](ICSI_2022_Benchmark_functions.pdf)


## Table of Contents

- [ICSI'2022 Benchmark Functions](#icsi2022-benchmark-functions)
	- [Table of Contents](#table-of-contents)
	- [Background](#background)
	- [Install](#install)
		- [Installation for MATLAB](#installation-for-matlab)
		- [Installation for python](#installation-for-python)
	- [Usage](#usage)
	- [Definition of test functions](#definition-of-test-functions)
	- [Related repository](#related-repository)
	- [Maintainers](#maintainers)
	- [License](#license)

## Background

ICSI is the abbreviation of International  Conference on Swarm Intelligence, which has successfully held 12 sessions so far. This conference is one of the top conferences in the field of swarm intelligence. The 13th International Swarm Intelligence Conference (ICSI’2022) in 2022 will be held in Xi’an, China in July 2022. [website of ICSI'2022](http://iasei.org/icsi2022/)

In order to make researchers who study swarm intelligence optimization algorithms and evolutionary algorithms better test their algorithms, we have designed a set of test functions to provide a better test platform. At the same time, we held an algorithm competition based on this test set. All researchers who swarm intelligence algorithms and evolutionary algorithms are welcome to participate in this competition. We will provide generous bonuses as rewards for the excellent contestants in this competition.

## Install

This test suite provides [python](https://www.python.org/) and [MATLAB](https://www.mathworks.com/products/matlab.html) support. Please make sure you install them.

### Installation for MATLAB

For MATLAB users, please use MATLAB and run following command to compile file [icsi2022.cpp](MATLAB\icsi2022.cpp).

```MATLAB
mex icsi2022.cpp -DWINDOWS
```

**Note: Please make sure that the compiler supported by MATLAB has been installed in the system.**


### Installation for python

First make sure that the following dependencies are installed in the system

```txt
Cython==0.29.23
matplotlib==3.3.4
numpy==1.20.1
```

For Linux users，please run the following command to compile and install

```sh
$ bash c.sh
```

For Windows users, please make sure you have a suitable compiler in your system.  we recommend to use Microsoft Visual C++ 14.0 or higher. You can run the following command to compile after the compiler has been installed in your system.

```sh
python setup.py build_ext --inplace
```


## Usage

After successful installation, you can run a simple command to use this test set function.

Here is a simple example for MATLAB users

```MATLAB
y = icsi2022(x, i)
% x is a m*D matrix
% i is from 1 to 10 which represent the serial number of 10 functions in icsi'2022 test suite.
% 
```

Here is a simple example for python users

```python
import numpy as np
import icsi2022
# x is a m*D matrix

y=icsi2022.eval(x,i)
# i is from 0 to 9 which represent the serial number of 10 functions in icsi'2022 test suite.

```

## Definition of test functions

Please refer to the documentation [icsi2022.pdf](ICSI_2022_Benchmark_functions.pdf)
to obtain the information of test functions.

## Related repository

- [Fireworks Algorithm](https://github.com/cil-lab/fwaopt/tree/master/mpopt) — The official repository for Fireworks Algorithm.

## Maintainers

[@moshangsang24](https://github.com/moshangsang24)

If you have any problems during use, you can contact me via email [icsi2022@iasei.org] or directly submit issue on github.

## License

[MIT](LICENSE) © CIL-Lab