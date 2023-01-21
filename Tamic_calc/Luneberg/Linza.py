import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle
from multiprocessing import Pool
import pandas as pd
path = './'


delta = 0.002
def Epsilon_works(Coordinates, ray_freq=0.5, Norm=True):
    x, z = Coordinates[0], Coordinates[1]
    koeff = 22.1119  # last 2.21119 * (10 ** 19). included factor (10^9)^2 from freq
    r_source = -0.938854
    x2_y2 = x ** 2
    par0 = 0
    deltaEps = -1
    teta = np.degrees(np.where(z > 0, np.arctan2(np.sqrt(x2_y2), z), 0))
    f_teta = np.where((teta >= 0) & (teta <= 85), 2,
                      np.where((teta > 85) & (teta < 90), -0.1 * teta + 10.5, 1))
    n_teta = np.where((teta >= 0) & (teta <= 20), -0.7 * teta + 32,
                      np.where((teta > 20) & (teta <= 40), 0.2 * teta + 14,
                               np.where((teta > 40) & (teta <= 60), 0.5 * teta + 2,
                                        np.where((teta > 60) & (teta <= 90), 0.23 * teta + 46, 1))))
    p = np.where(z <= 0, 0,
                 ((np.cos(np.radians(teta / f_teta)) ** n_teta) / ((np.sqrt(x2_y2 + z ** 2) + r_source) ** 2)))
    wp = koeff * p
    wp_w = wp / (ray_freq ** 2)
    # exp = 1 - wp_w
    exp = - wp_w  # для Tamic
    return exp
    # if exp > 0:
    #     return exp
    # else:
    #     return 0  #* dielectric_constant * 10 ** (-12)
    #


def Epsilon(Coordinates, ray_freq=0.4, Norm=True):
    x,y = Coordinates[0],Coordinates[1]
    Z_c = 0
    X_c = 0.5
    Y_c = 0.5
    R = 0.5
    dx = (x - X_c)
    dy = (y - Y_c)
    r = np.sqrt(dx ** 2 + dy ** 2)
    with np.errstate(invalid='ignore'):
        # par = 0.3 # used for magnetic calculations
        par = 1
        exp = np.where(r <= R, par * np.sqrt(2 - (r / R) ** 2) ** 2, par)
        exp -=1
        return exp


def check_eps():
    x = np.linspace(0, 1, 100)
    y = np.linspace(0, 1, 100)
    X, Y = np.meshgrid(x, y)
    Eps = Epsilon([X, Y])
    plt.contourf(X, Y, Eps)
    plt.show()

def calc_rasp(angle):
    X_boundary = 1  # z
    Y_boundary = 1  # x
    Xmin = 0  # z
    Ymin = 0  # x
    x1, x2, x3 = 0, Xmin, X_boundary
    y1, y2, y3 = Y_boundary, Ymin, Ymin
    angle = - angle
    x11, y11 = x1 * np.cos(angle) + y1 * np.sin(angle), -x1 * np.sin(angle) + y1 * np.cos(angle)
    x22, y22 = x2 * np.cos(angle) + y2 * np.sin(angle), -x2 * np.sin(angle) + y2 * np.cos(angle)
    x33, y33 = x3 * np.cos(angle) + y3 * np.sin(angle), -x3 * np.sin(angle) + y3 * np.cos(angle)
    Xmin,X_boundary = min(x11,x22,x33),max(x11,x22,x33)
    Ymin,Y_boundary = min(y11,y22,y33),max(y11,y22,y33)
    angle = - angle
    x = np.arange(Xmin, X_boundary + delta, delta)
    y = np.arange(Ymin, Y_boundary + delta, delta)
    node = 0
    nX = len(x)
    nY = len(y)
    nPoint = nX * nY
    default_name = 'Linza_'
    angle_str = str(np.round(angle,0))
    def generate_TPL(angle,freq=10):
        name = default_name+angle + '.TPL'
        eps_name = default_name+angle
        file = path + name
        f = open(f'./{file}', 'w')
        s = f'#TMC_RT_H\n' \
            '#define L_waveg @ ( 2000.00)\n' \
            f'#define W_file  @ {eps_name}\n' \
            f'#define W_freq  @ (  {freq})\n' \
            f'#define W_eps @ ( 4. )\n' \
            f'#define W_time @ ( 25. )\n' \
            f'#define W_waveg @ (  1100.00)\n' \
            f'#define W_input @ (  W_waveg-100.00)\n' \
            f'#define W_delta @ (   {delta*1000})\n' \
            f'#define W_r @ (   40.00)\n' \
            f'#define W_type @ MAGNETIC\n' \
            f'#define a1 @ ( W_input )\n' \
            f'#define a11 @ ( 20 )\n' \
            f'#define l_rupor @ ( 10 )\n' \
            f'#define h_rupor @ ( 0 )\n' \
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
            f'  FILE W_file; 0; 10000; 0; 10000;\n' \
            f' END_B\n' \
            f' BLOCK 8;\n' \
            f'  RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0.0; W_waveg-2*W_delta;\n' \
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
            f' T 55; 50;  50;\n' \
            f'#END_LINK\n' \
            f'\n' \
            f'#OUTPUT\n' \
            f' FILE W_file;\n' \
            f' FIELDS;\n' \
            f' TOPOLOGY;\n' \
            f' FIELD_DISTRIBUTION_M  W_file; W_freq; 15; 20;\n' \
            f'#END_OUTPUT\n' \
            f'\n' \
            f'#END_STEP\n' \
            f'\n' \
            f'#EOF'
        f.write(s)
        f.close()
    generate_TPL(angle_str)
    name = f'{default_name}{angle_str}.eps'
    file = path + name
    def write_headres(file):
        f = open(f'./{file}', 'w')
        f.write('#TamicRTH_planar_DistributionDielectricPermeability_File_V2.00 2000')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write('#TopologyPrimitiv RECT_STAT ')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write(f'#dDelta {delta}')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write(f'#Xmin {Xmin}')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write(f'#Ymin {Ymin}')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write(f'#nX {nX}')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write(f'#nY {nY}')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write(f'#nPoint {nPoint}')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write(f'#nAccuracy 4')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write('#sNodeFormat NULL')
        f.close()
        set_(file)
        f = open(f'./{file}', 'a')
        f.write('#sValueFormat NULL')
        f.close()
        set_(file)
    def set_a(a, Xi, Yi):
        if (Yi == 0) and (Xi == 0):
            a = a / 4
        elif (Yi == len(y)) and (Xi == len(x)):
            a = a / 4
        elif (Yi == 0) and (Xi == len(x)):
            a = a / 4
        elif (Yi == len(y)) and (Xi == 0):
            a = a / 4
        elif (Yi == 0) or (Yi == len(y)) or (Xi == 0) or (Xi == len(x)):
            a = a / 2
        return a
    write_headres(file)
    count = 0
    N = len(x) * len(y)
    f = open(f'./{file}', 'ab')

    for Yi, Y in enumerate(y):
        stroka = ''
        for Xi, X in enumerate(x):
            X1, Z1 = X, Y
            a = set_a(Epsilon([X1, Z1]), Xi, Yi)
            f.write(bytearray(struct.pack("!f", a))[::-1])
            f.write(node.to_bytes(4, 'little'))
            node += 1
            count += 1
            if count % 1000 == 0:
                print(count / N * 100)
            # print(Z1[Yi],X1[Xi])
            # df[Z1[Yi]][X1[Xi]] = a
    f.close()
    # print(df)

def set_(name):
    f = open(f'./{name}', 'ab')
    f.write(bytes('\r', 'ascii'))
    f.close()
    f = open(f'./{name}', 'a')
    f.write('\n')
    f.close()
potoki = 1
if __name__ == '__main__':
    angles = np.radians(np.linspace(0,360,13))
    angles = np.array([0])
    Pool(potoki).map(calc_rasp, angles)
    Pool(potoki).close()