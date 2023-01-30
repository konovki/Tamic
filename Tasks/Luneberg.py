import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle
from multiprocessing import Pool
import properties as prop
import Library.Library as lib
import os.path
path = ''
for item in os.getcwd().split(prop.splitter)[:-2]:
    path += item + '/'
path += 'Results/'
lib.check_path(path)
path += 'Luneberg_Lense/'
lib.check_path(path)
print('path',path)

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

delta = 0.5#0.002
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
lib.write_headres(file,delta,Xmin,Ymin,nX,nY,nPoint)
count = 0
N = len(x)*len(y)
f = open(file,'ab')
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
f.close()
def generate_TPL(freq=10):
    name = 'lense' + '.TPL'
    eps_name = 'lense'
    file = path + name
    f = open(file, 'w')
    s = f'#TMC_RT_H\n' \
        f'#define L_waveg @ ( 2000.00)\n' \
        f'#define W_file  @ {eps_name}\n' \
        f'#define W_freq  @ (  {freq})\n' \
        f'#define W_eps @ ( 4. )\n' \
        f'#define W_time @ ( 25. )\n' \
        f'#define W_type @ MAGNETIC\n' \
        f'#define W_input @ (  W_waveg-100.00)\n' \
        f'#define a1 @ ( W_input )\n' \
        f'#define a11 @ ( 20 )\n' \
        f'#define l_rupor @ ( 10 )\n' \
        f'#define h_rupor @ ( 0 )\n' \
        f'#define W_waveg @ (  1100.00)\n' \
        f'#define W_delta @ (   2.0)\n' \
        f'#define W_r @ (   40.00)\n' \
        f'\n' \
        f'#STEP\n' \
        f'\n' \
        f'#PARM\n' \
        f' ANGLE_UNIT radian;\n' \
        f' FREQ_UNIT GHz;\n' \
        f' LONG_UNIT mm;\n' \
        f' TIME_UNIT ns;\n' \
        f' DELTA W_delta;\n' \
        f' TIME 0.; W_time;\n' \
        f' X_MIN  0.0;\n' \
        f' X_MAX L_waveg;\n' \
        f' Y_MIN  0.0;\n' \
        f' Y_MAX W_waveg;\n' \
        f' FREQ  W_freq;\n' \
        f'#END_PARM\n' \
        f'\n' \
        f'#TOPOLOGY\n' \
        f'\n' \
        f' BLOCK 7;\n' \
        f'  INPUT_X  W_type; -a1/2; a1/2; 0; W_time;\n' \
        f' END_B\n' \
        f'\n' \
        f' BLOCK 5;\n' \
        f'  POLYGON_STAT W_type;\n' \
        f'   L   0;   a1/2;\n' \
        f'   L   2*a11;   a1/2;\n' \
        f'   L   2*a11+(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*((l_rupor)/4);\n' \
        f'   L   2*a11+2*(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*(2*(l_rupor)/4);\n' \
        f'   L   2*a11+3*(l_rupor)/4;   a1/2+((h_rupor)/(l_rupor))*(3*(l_rupor)/4);\n' \
        f'   L   2*a11+l_rupor;   a1/2+h_rupor;\n' \
        f'  END_B\n' \
        f' BLOCK 6;\n' \
        f'  POLYGON_STAT W_type;\n' \
        f'   L   0;   -a1/2;\n' \
        f'   L   2*a11;   -a1/2;\n' \
        f'   L   2*a11+(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*((l_rupor)/4);\n' \
        f'   L   2*a11+2*(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*(2*(l_rupor)/4);\n' \
        f'   L   2*a11+3*(l_rupor)/4;   -a1/2-((h_rupor)/(l_rupor))*(3*(l_rupor)/4);\n' \
        f'   L   2*a11+l_rupor;   -a1/2-h_rupor;\n' \
        f'  END_B\n' \
        f'\n' \
        f' BLOCK 55;\n' \
        f'  FILE W_file; 0; {10000}; 0; {10000};\n' \
        f' END_B\n' \
        f' BLOCK 8;\n' \
        f'RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0.0; W_waveg-2*W_delta;\n' \
        f' END_B\n' \
        f' BLOCK 9;\n' \
        f'  RECT_STAT ABSORBER; 1*W_delta; L_waveg-2*W_delta; W_delta; 5*W_delta;\n' \
        f' END_B\n' \
        f' BLOCK 10;\n' \
        f'  RECT_STAT ABSORBER; 1*W_delta; L_waveg; W_waveg-5*W_delta; W_waveg;\n' \
        f' END_B\n' \
        f' BLOCK 11;\n' \
        f'  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; W_waveg/2-a1/2-4*W_delta;\n' \
        f' END_B\n' \
        f' BLOCK 12;\n' \
        f'  RECT_STAT ABSORBER; 0.0; 5*W_delta; W_waveg/2+a1/2+4*W_delta; W_waveg;\n' \
        f' END_B\n' \
        f'\n' \
        f'#END_TOPOLOGY\n' \
        f'\n' \
        f'#LINK_LIST\n' \
        f' T  5; 0.0; W_waveg/2;\n' \
        f' T  6; 0.0; W_waveg/2-W_delta;\n' \
        f' T  7; 0.0; W_waveg/2;\n' \
        f' T  8; 0.0; 0.0;\n' \
        f' T  9; 0.0; 0.0;\n' \
        f' T 10; 0.0; 0.0;\n' \
        f' T 11; 0.0; 0.0;\n' \
        f' T 12; 0.0; 0.0;\n' \
        f' T 55; 50; 50;\n' \
        f'#END_LINK\n' \
        f'\n' \
        f'#OUTPUT\n' \
        f' FILE W_file;\n' \
        f' FIELDS;\n' \
        f' TOPOLOGY;\n' \
        f' FIELD_DISTRIBUTION_M  W_file; W_freq; W_time-10; W_time-5;\n' \
        f'#END_OUTPUT\n' \
        f'\n' \
        f'#END_STEP\n' \
        f'\n' \
        f'#EOF'
    f.write(s)
    f.close()
generate_TPL(freq=10)
