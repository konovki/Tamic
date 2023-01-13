
import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle

import pandas as pd


def Epsilon_works(Coordinates,ray_freq=0.5,Norm = True):
    x,z = Coordinates[0],Coordinates[1]
    koeff = 22.1119   # last 2.21119 * (10 ** 19). included factor (10^9)^2 from freq
    r_source = -0.938854
    x2_y2 = x ** 2
    par0 = 0
    deltaEps= -1
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
def Epsilon(Coordinates,ray_freq=0.4,Norm = True):
    RocketXLoc = 0
    RocketZLow = 0
    RocketZHigh = 2.35
    RocketConeHeight = 4.08
    RocketWidth = 1.17
    zAdd = RocketZHigh+RocketZHigh
    x,z = Coordinates[0],-Coordinates[1]
    koeff = 22.1119   # last 2.21119 * (10 ** 19). included factor (10^9)^2 from freq
    r_source = -0.938854
    x2_y2 = x ** 2
    par0 = 0
    deltaEps= -1
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
    #rocket
    z = -z
    par = -1

    k = (RocketWidth) / (RocketZHigh - (RocketZHigh + RocketConeHeight))
    b = RocketXLoc - k * (RocketZHigh + RocketConeHeight)
    sqrtX = np.sqrt((x - RocketXLoc) ** 2)
    exp = np.where(((z >= RocketZLow) & (z <= RocketZHigh) & (sqrtX <= RocketWidth)), par, exp)
    exp = np.where(((z >= RocketZHigh) & (z <= (RocketZHigh + RocketConeHeight)) & (
            np.abs(x - RocketXLoc) <= k * z + b - RocketXLoc)), par, exp)
    # soplo
    def calc_k_b(z1,z2,x1,x3):
        k = (z2-z1)/(x3-x1)
        b = z1 - k*x1
        return k,b
    zSopla=0.3
    xlSopla1=0.824
    xlSopla3=0.188
    xlSopla2=0.736
    xlSopla4=0
    kl,bl = calc_k_b(RocketZLow-zSopla,RocketZLow,RocketXLoc-xlSopla1,RocketXLoc-xlSopla3)
    kr,br = calc_k_b(RocketZLow-zSopla,RocketZLow,RocketXLoc-xlSopla2,RocketXLoc-xlSopla4)
    exp = np.where(((z<=RocketZLow)&(z>=RocketZLow-zSopla)&(z>=kr*x+br)&(z<=kl*x+bl)),par,exp)
    kl,bl = calc_k_b(RocketZLow-zSopla,RocketZLow,RocketXLoc+xlSopla2,RocketXLoc+xlSopla4)
    kr,br = calc_k_b(RocketZLow-zSopla,RocketZLow,RocketXLoc+xlSopla1,RocketXLoc+xlSopla3)
    exp = np.where(((z<=RocketZLow)&(z>=RocketZLow-zSopla)&(z<=kr*x+br)&(z>=kl*x+bl)),par,exp)

    return exp
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
X_boundary = 15#z
Y_boundary = 15#x
Xmin = -15#z
Ymin = -15#x
print(X_boundary/delta)
x = np.arange(Xmin,X_boundary+delta,delta)
y = np.arange(Ymin,Y_boundary+delta,delta)
node = 0
nX = len(x)
nY = len(y)
nPoint = nX*nY
path = 'Tamic_calc/Rocket/'
name = 'rocket.eps'
file = path+name
def write_headres(file):
    f = open(f'./{file}','w')
    f.write('#TamicRTH_planar_DistributionDielectricPermeability_File_V2.00 2000')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write('#TopologyPrimitiv RECT_STAT ')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write(f'#dDelta {delta}')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write(f'#Xmin {Xmin}')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write(f'#Ymin {Ymin}')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write(f'#nX {nX}')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write(f'#nY {nY}')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write(f'#nPoint {nPoint}')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write(f'#nAccuracy 4')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write('#sNodeFormat NULL')
    f.close()
    set_(file)
    f = open(f'./{file}','a')
    f.write('#sValueFormat NULL')
    f.close()
    set_(file)
def set_a(a,Xi,Yi):
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
    return a
write_headres(file)
count = 0
N = len(x)*len(y)
f = open(f'./{file}','ab')
df = pd.DataFrame(columns=y,index=x)
angle = np.radians(270)

X1,Z1 = x * np.cos(angle)+y*np.sin(angle),-x * np.sin(angle)+y*np.cos(angle)
for Yi,Y in enumerate(y):
    stroka = ''
    for Xi,X in enumerate(x):
        a = set_a(Epsilon([X,Y]),Xi,Yi)
        f.write(bytearray(struct.pack("!f", a))[::-1])
        f.write(node.to_bytes(4,'little'))
        node += 1
        count += 1
        if count % 1000 == 0:
            print(count / N * 100)
        print(Z1[Yi],X1[Xi])
        df[Z1[Yi]][X1[Xi]] = a
f.close()
print(df)
