import numpy as np
import icsi2022
from matplotlib import pyplot as plt
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource





# sys.path.append("/home/liuyifan/code/fwaopt")
def main():
    a = 50*np.ones((1,50))
    b=np.zeros((1,10))
    # for 10 dimension
    c=np.array([[6.1770164235388009e+00,-8.2403550782428852e+00,-2.1511735494021202e+00, 1.0662928429740914e+01,-2.8853127679253522e+00, 3.1242073274795931e+01, 7.1837517429067532e+00,-6.7539761135140324e+01,-7.9194126903635507e+00, 2.7532766213422995e+01]])
    d=np.array([[ 6.1770164235388009e+00,-8.2403550782428852e+00,-2.1511735494021202e+00, 1.0662928429740914e+00,-2.8853127679253522e+00, 3.1242073274795931e+00, 7.1837517429067532e+00,-6.7539761135140324e+00,-7.9194126903635507e+00, 2.7532766213422995e+01]])
    # for 20 dimension
    # c=np.array([[6.1770164235388009e+01,-8.2403550782428852e+00,-2.1511735494021202e+01, 1.0662928429740914e+01,-2.8853127679253522e+01, 3.1242073274795931e+01, 7.1837517429067532e+01,-6.7539761135140324e+01,-7.9194126903635507e+01, 2.7532766213422995e+01,-4.6393239169564673e+01,-1.9189660498495055e+01, 4.2161826637648900e+01,-2.4089553091716745e+01,-6.4419553958074914e+01, 5.9118560349158273e+01, 5.8197810408751394e+01,-4.5515033236478935e+01,-3.6031000228613152e+01,-5.2275229828560370e+01]])

    for i in range(10):
        print(icsi2022.eval(d+1,i+1))


def main3():
    dim=10
    x=y=np.arange(-100,101,1)
    sp=x.shape[0]
    X,Y=np.meshgrid(x,y)
    # add random number for caculate 10D benchmark

    # z=np.random.uniform(-100,100,(dim-2))
    z=np.zeros((dim-2))
    x1=X.reshape(sp*sp,1)
    y1=Y.reshape(sp*sp,1)
    d=np.hstack((x1,y1,np.tile(z,(sp*sp,1))))
    for func_id in range(6,7):
        print(func_id)
        # calculate the function value
        ans=icsi2022.eval(d,func_id).reshape((sp,sp))



        fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
        ls = LightSource(270, 45)
        # To use a custom hillshading mode, override the built-in shading and pass
        # in the rgb colors of the shaded surface calculated from "shade".
        rgb = ls.shade(ans, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
        surf = ax.plot_surface(X, Y, ans, rstride=1, cstride=1, facecolors=rgb,
                            linewidth=0, antialiased=False, shade=False)



        fig.savefig("./fig/3d_{}.png".format(func_id),pad_inches=0.0)
        plt.clf()



if __name__ == '__main__':
    main3()
    