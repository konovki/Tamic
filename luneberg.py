
import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle

import pandas as pd

def Epsilon(Coordinates,  Norm=True):
    x, y= Coordinates[0], Coordinates[1]
    X_c = 0.75
    Y_c = 0.75
    R = 0.5
    dx = (x - X_c)
    dy = (y - Y_c)
    r = np.sqrt(dx ** 2 + dy ** 2)
    with np.errstate(invalid='ignore'):
        # par = 0.3 # used for magnetic calculations
        par = 0
        return np.where(r <= R, -1+ np.sqrt(2 - (r / R) ** 2) ** 2, par)
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

delta = 0.002
X_boundary = 1.5#z
Y_boundary = 1.5#x
Xmin = 0#z
Ymin = 0#x
print(X_boundary/delta)
x = np.linspace(Xmin,X_boundary,int(X_boundary/delta))
y = np.linspace(Ymin,X_boundary,int(X_boundary/delta))
node = 0
nX = len(x)
nY = len(y)
nPoint = nX*nY
path = 'Tamic_calc/'
name = 'Luneberg.eps'
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
        df[Y][X] = a
f.close()
print(df)
