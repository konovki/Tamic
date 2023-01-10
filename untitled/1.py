import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle
def Epsilon(Coordinates, ray_freq=0):
    x, y = Coordinates[0], Coordinates[1]
    X_c = 0.75
    Y_c = 0.75
    R = 0.5
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
delta = 0.001
X_boundary = 1.5
Y_boundary = 1.5
print(X_boundary/delta)
x = np.linspace(0,X_boundary,int(X_boundary/delta))
y = np.linspace(0,X_boundary,int(X_boundary/delta))
node = 0
Xmin = 0
Ymin = 0
nX = len(x)
nY = len(y)
nPoint = nX*nY
name = 'Eps.eps'
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
f = open('./Eps.eps','ab')
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

    # f.write(b'{stroka}')
f.close()

