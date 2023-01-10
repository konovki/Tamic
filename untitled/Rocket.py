
import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle

import pandas as pd


def Epsilon(Coordinates,ray_freq=0.4,Norm = True):
    # ang = m.radians(45)
    # x =
    x,z = Coordinates[0],-Coordinates[1]
    koeff = 22.1119   # last 2.21119 * (10 ** 19). included factor (10^9)^2 from freq
    r_source = -0.938854
    x2_y2 = x ** 2
    teta = np.degrees(np.where(z > 0,np.arctan2(np.sqrt(x2_y2),z),0))
    f_teta = np.where((teta >= 0) & (teta <= 85),2,
                      np.where((teta > 85) & (teta < 90),-0.1 * teta + 10.5,1))
    n_teta = np.where((teta >= 0) & (teta <= 20),-0.7 * teta + 32,
                      np.where((teta > 20) & (teta <= 40),0.2 * teta + 14,
                               np.where((teta > 40) & (teta <= 60),0.5 * teta + 2,
                                        np.where((teta > 60) & (teta <= 90),0.23 * teta + 46,1))))
    p = np.where(z <= 0,0,((np.cos(np.radians(teta / f_teta)) ** n_teta) / ((np.sqrt(x2_y2 + z ** 2) + r_source) ** 2)))
    wp = koeff * p
    wp_w =  wp / (ray_freq ** 2)
    # exp = 1 - wp_w
    exp = - wp_w # для Tamic
    return exp
    # if exp > 0:
    #     return exp
    # else:
    #     return 0  #* dielectric_constant * 10 ** (-12)
    #

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

delta = 0.02
X_boundary = 10#z
Y_boundary = 0#x
Xmin = -10#z
Ymin = -10#x
print(X_boundary/delta)
x = np.linspace(Xmin,X_boundary,int(X_boundary/delta))
y = np.linspace(Ymin,X_boundary,int(X_boundary/delta))
node = 0
nX = len(x)
nY = len(y)
nPoint = nX*nY
name = 'fakel.eps'
f = open(f'./{name}','w')
f.write('#TamicRTH_planar_DistributionDielectricPermeability_File_V2.00 2000')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write('#TopologyPrimitiv RECT_STAT ')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write(f'#dDelta {delta}')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write(f'#Xmin {Xmin}')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write(f'#Ymin {Ymin}')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write(f'#nX {nX}')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write(f'#nY {nY}')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write(f'#nPoint {nPoint}')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write(f'#nAccuracy 4')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write('#sNodeFormat NULL')
f.close()
set_(name)
f = open(f'./{name}','a')
f.write('#sValueFormat NULL')
f.close()
set_(name)
# f = open(f'./{name}','ab')
# f.write(bytes('A>','ascii'))
# f.close()
# f = open(f'./{name}','ab')
# a = int(0)
# f.write(a.to_bytes(4,'little')[2:])
# x = bytes('耾','utf-16')[::-1]
# print(len(x))
# f.write(x[:])
# f.close()
count = 0
N = len(x)*len(y)
f = open(f'./{name}','ab')
df = pd.DataFrame(columns=y,index=x)

for Yi,Y in enumerate(y):
    stroka = ''
    for Xi,X in enumerate(x):
        a = Epsilon([X,Y])
        if (Yi == 0) and (Xi == 0):
            a = a/4
        elif (Yi == len(y)) and (Xi == len(x)):
            a = a/4
        elif (Yi == 0) and (Xi == len(x)):
            a = a/4
        elif (Yi == len(y)) and (Xi == 0):
            a = a/4
        elif (Yi == 0) or (Yi == len(y)) or (Xi == 0) or (Xi == len(x)):
            a = a/2
        f.write(bytearray(struct.pack("!f", a))[::-1])
        # f.write(a.to_bytes(4,'little'))
        # b = '?'
        # f.write(bytes(b,'ascii'))
        # f.write(bytes(f'{node}','ascii'))
        f.write(node.to_bytes(4,'little'))
        node += 1
        count += 1
        if count % 1000 == 0:
            print(count / N * 100)
        df[Y][X] = a
    # f.write(b'{stroka}')
f.close()
print(df)
