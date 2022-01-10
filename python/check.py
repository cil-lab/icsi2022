import numpy as np
import icsi2022
from matplotlib import pyplot as plt





def main():
    a = 50*np.ones((1,50))
    for i in range(10):
        print(icsi2022.eval(a,i+1))




if __name__ == '__main__':
    main()
    