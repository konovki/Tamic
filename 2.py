import numpy as np
import matplotlib.pyplot as plt
import math as m
import pickle
def Epsilon(Coordinates, ray_freq=0):
    x, y = Coordinates[0], Coordinates[1]
    X_c = 0.5
    Y_c = 0.5
    R = 0.4
    dx = (x - X_c)
    dy = (y - Y_c)
    r = np.sqrt(dx ** 2 +dy ** 2)
    if r <= R:
        return m.sqrt(2 - (r / R) ** 2) ** 2
    else:
        return 1

    # with np.errstate(invalid='ignore'):
    #     # par = 0.3 # used for magnetic calculations
    #     par = 1
    #     return np.where(r <= R, par * np.sqrt(2 - (r / R) ** 2) ** 2, par)
def check_eps():
    x = np.linspace(0,1,100)
    y = np.linspace(0,1,100)
    X,Y = np.meshgrid(x,y)
    Eps = Epsilon([X,Y])
    plt.contourf(X,Y,Eps)
    plt.show()
def set_(name):
    f = open(f'./{name}','ab')
    f.write(bytes('\r','ascii'))
    f.close()
    f = open(f'./{name}','a')
    f.write('\n')
    f.close()
length = 10
n = 123456789

