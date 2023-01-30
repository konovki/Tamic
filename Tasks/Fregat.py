import properties as prop
import Library.Library as lib
import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle
import os.path
import pandas as pd


def Epsilon_works(Coordinates,ray_freq=0.6,Norm = True):
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
def Epsilon(Coordinates,ray_freq=0.6,Norm = True):
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
    par = -100

    # k = (RocketWidth) / (RocketZHigh - (RocketZHigh + RocketConeHeight))
    # b = RocketXLoc - k * (RocketZHigh + RocketConeHeight)
    # sqrtX = np.sqrt((x - RocketXLoc) ** 2)
    # exp = np.where(((z >= RocketZLow) & (z <= RocketZHigh) & (sqrtX <= RocketWidth)), par, exp)
    # exp = np.where(((z >= RocketZHigh) & (z <= (RocketZHigh + RocketConeHeight)) & (
    #         np.abs(x - RocketXLoc) <= k * z + b - RocketXLoc)), par, exp)
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

path = ''
for item in os.getcwd().split(prop.splitter)[:-2]:
    path += item + '/'
path += 'Results/'
lib.check_path(path)
path += 'Fregat/'
lib.check_path(path)
delta = 0.02
X_boundary = 15#z
Y_boundary = 7#x
Xmin = -15#z
Ymin = -20#x
print(X_boundary/delta)
x = np.arange(Xmin,X_boundary+delta,delta)
y = np.arange(Ymin,Y_boundary+delta,delta)
node = 0
nX = len(x)
nY = len(y)
nPoint = nX*nY
name = 'fakel.eps'
file = path+name

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
def generate_TPL_fregat():
    name = 'fregat' + '.TPL'
    file = path + name
    f = open(file, 'w')
    s = f'#TMC_RT_H\n' \
        f'!#define L_input @ (   1.40)\n' \
        f'#define W_input @ (   1.00)\n' \
        f'#define W_delta @ (   0.02)\n' \
        f'#define W_freq @ (  0.6 )\n' \
        f'#define xPlus @ ( 2 )\n' \
        f'#define yPlus @ ( 0 )\n' \
        f'#define a1 @ ( W_input )\n' \
        '#define a2 @ ( 0.2 )\n' \
        '#define calc_time @ ( 200 )\n'\
        '#define name @ fregat\n' \
        '#define h1 @ ( 0.880 )\n' \
        '#define h2 @ ( 6.807 )\n' \
        '#define h3 @ ( 1.493 )\n' \
        '#define h4 @ ( 2.100 )\n' \
        '#define d1  @ ( 2.800 )\n' \
        '#define d2  @ ( 3.715 )\n' \
        '#define d3  @ ( 2.900 )\n' \
        '#define d_a @ ( (d1+d2)/2 )\n' \
        '#define h_a @ ( h1*(d_a-d1)/(d2-d1) )\n' \
        '#define d_a1 @ ( d_a + 6*(W_delta) )\n' \
        '#define h_a1 @ ( h1*(d_a1-d1)/(d2-d1) )\n' \
        '#define l_a @ ( 0.354 )\n' \
        '#define x_a @ ( h1*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_a @ ( ((d2-d1)/2)*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define l_r @ ( 0.470 )\n' \
        '#define x_r @ ( h1*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_r @ ( ((d2-d1)/2)*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define  d_r @ ( 0.177 )\n' \
        '#define xd_r @ ( ((d2-d1)/2)*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define yd_r @ ( h1*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define x_shift_fakel @ 1.0\n' \
        '#define y_shift_fakel @ 0\n' \
        '#define y_glob_shift @ -5\n' \
        '#define x_shift @ ( 7.000 )\n' \
        '#define y_fakel_size @ ( 25.000 )\n' \
        '#define W_type1 @ MAGNETIC\n' \
        '#define W_type  @ METAL\n' \
        '#define W_waveg @ (  h1 + h2 + h3 + h4 )\n' \
        '#define L_waveg @ (  d2 + 10 + x_shift + xPlus )\n' \
        '#STEP\n' \
        '#PARM\n' \
        ' ANGLE_UNIT radian;\n' \
        ' FREQ_UNIT GHz;\n' \
        ' LONG_UNIT m;\n' \
        ' TIME_UNIT ns;\n' \
        ' DELTA W_delta;\n' \
        ' TIME 0.;  calc_time;\n' \
        ' X_MIN  0.0;\n' \
        ' X_MAX L_waveg;\n' \
        ' Y_MIN  0.0;\n' \
        ' Y_MAX W_waveg+y_fakel_size;\n' \
        ' FREQ  W_freq;\n' \
        '#END_PARM\n' \
        '#TOPOLOGY\n' \
        ' BLOCK 1;\n' \
        '  INPUT_X  W_type1 ; h_a-2*W_delta; h_a; 0; 720; 1; 0.;\n' \
        ' END_B\n' \
        ' BLOCK 101;\n' \
        '  INPUT_X  W_type1 ; h_a1-2*W_delta; h_a1; 0; 720; 1; 3.141592653589;\n' \
        ' END_B\n' \
        ' BLOCK 4;\n' \
        '  RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 5;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; W_waveg+y_fakel_size-5*W_delta; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 6;\n' \
        '  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; W_waveg+y_fakel_size ;\n' \
        ' END_B\n' \
        ' BLOCK 7;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; 0; 5*W_delta ;\n' \
        ' END_B\n' \
        ' BLOCK 8;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        '   L   -(d2)/2;            h1;\n' \
        '   L   -(d2)/2;         h1+h2;\n' \
        '   L   -(d3)/2;      h1+h2+h3;\n' \
        '   L       0.0;   h1+h2+h3+h4;\n' \
        '   L    (d3)/2;      h1+h2+h3;\n' \
        '   L    (d2)/2;         h1+h2;\n' \
        '   L    (d2)/2;            h1;\n' \
        '   L    (d1)/2;           0.0;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        ' END_B\n' \
        ' BLOCK 9;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_a;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_a; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        ' BLOCK 10;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_r;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_r; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        '#END_TOPOLOGY\n' \
        '#LINK_LIST\n' \
        ' T  1; (L_waveg)/2+x_shift/2-d_a/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  101; (L_waveg)/2+x_shift/2-d_a1/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  4; 0.0; 0;\n' \
        ' T  5; 0.0; 0.0;\n' \
        ' T  6; 0.0; 0;\n' \
        ' T  7; 0.0; 0;\n' \
        ' T  8; (L_waveg)/2+x_shift/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  9; 0.0; 0.0 +y_fakel_size+y_glob_shift;\n' \
        ' T 10; -xd_r; yd_r +y_fakel_size+y_glob_shift;\n' \
        '#END_LINK\n' \
        '#OUTPUT\n' \
        ' FILE name;\n' \
        ' FIELDS;\n' \
        f' TOPOLOGY;\n' \
        f' FIELD_DISTRIBUTION_M name; W_freq;  calc_time -10.0;  calc_time -5.0;\n' \
        f'#END_OUTPUT\n' \
        f'#END_STEP\n' \
        f'#EOF'
    f.write(s)
    f.close()
def generate_TPL_fregat_impulse():
    name = 'fregat_impulse' + '.TPL'
    file = path + name
    f = open(file, 'w')
    s = f'#TMC_RT_H\n' \
        f'!#define L_input @ (   1.40)\n' \
        f'#define W_input @ (   1.00)\n' \
        f'#define W_delta @ (   0.02)\n' \
        f'#define W_freq @ (  0.6 )\n' \
        f'#define xPlus @ ( 2 )\n' \
        f'#define yPlus @ ( 0 )\n' \
        f'#define a1 @ ( W_input )\n' \
        '#define a2 @ ( 0.2 )\n' \
        '#define calc_time @ ( 150 )\n' \
        '#define name @ fregat_impulse\n' \
        '#define h1 @ ( 0.880 )\n' \
        '#define h2 @ ( 6.807 )\n' \
        '#define h3 @ ( 1.493 )\n' \
        '#define h4 @ ( 2.100 )\n' \
        '#define d1  @ ( 2.800 )\n' \
        '#define d2  @ ( 3.715 )\n' \
        '#define d3  @ ( 2.900 )\n' \
        '#define d_a @ ( (d1+d2)/2 )\n' \
        '#define h_a @ ( h1*(d_a-d1)/(d2-d1) )\n' \
        '#define d_a1 @ ( d_a + 6*(W_delta) )\n' \
        '#define h_a1 @ ( h1*(d_a1-d1)/(d2-d1) )\n' \
        '#define l_a @ ( 0.354 )\n' \
        '#define x_a @ ( h1*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_a @ ( ((d2-d1)/2)*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define l_r @ ( 0.470 )\n' \
        '#define x_r @ ( h1*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_r @ ( ((d2-d1)/2)*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define  d_r @ ( 0.177 )\n' \
        '#define xd_r @ ( ((d2-d1)/2)*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define yd_r @ ( h1*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define x_shift_fakel @ 1.0\n' \
        '#define y_shift_fakel @ 0\n' \
        '#define y_glob_shift @ -5\n' \
        '#define x_shift @ ( 7.000 )\n' \
        '#define y_fakel_size @ ( 25.000 )\n' \
        '#define W_type1 @ MAGNETIC\n' \
        '#define W_type  @ METAL\n' \
        '#define W_waveg @ (  h1 + h2 + h3 + h4 )\n' \
        '#define L_waveg @ (  d2 + 10 + x_shift + xPlus )\n' \
        '#STEP\n' \
        '#PARM\n' \
        ' ANGLE_UNIT radian;\n' \
        ' FREQ_UNIT GHz;\n' \
        ' LONG_UNIT m;\n' \
        ' TIME_UNIT ns;\n' \
        ' DELTA W_delta;\n' \
        ' TIME 0.;  calc_time;\n' \
        ' X_MIN  0.0;\n' \
        ' X_MAX L_waveg;\n' \
        ' Y_MIN  0.0;\n' \
        ' Y_MAX W_waveg+y_fakel_size;\n' \
        ' FREQ  W_freq;\n' \
        '#END_PARM\n' \
        '#TOPOLOGY\n' \
        ' BLOCK 1;\n' \
        '  INPUT_X  W_type1 ; h_a-2*W_delta; h_a; 0; 1.66; 1; 0.;\n' \
        ' END_B\n' \
        ' BLOCK 101;\n' \
        '  INPUT_X  W_type1 ; h_a1-2*W_delta; h_a1; 0; 1.66; 1; 3.141592653589;\n' \
        ' END_B\n' \
        ' BLOCK 4;\n' \
        '  RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 5;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; W_waveg+y_fakel_size-5*W_delta; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 6;\n' \
        '  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; W_waveg+y_fakel_size ;\n' \
        ' END_B\n' \
        ' BLOCK 7;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; 0; 5*W_delta ;\n' \
        ' END_B\n' \
        ' BLOCK 8;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        '   L   -(d2)/2;            h1;\n' \
        '   L   -(d2)/2;         h1+h2;\n' \
        '   L   -(d3)/2;      h1+h2+h3;\n' \
        '   L       0.0;   h1+h2+h3+h4;\n' \
        '   L    (d3)/2;      h1+h2+h3;\n' \
        '   L    (d2)/2;         h1+h2;\n' \
        '   L    (d2)/2;            h1;\n' \
        '   L    (d1)/2;           0.0;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        ' END_B\n' \
        ' BLOCK 9;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_a;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_a; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        ' BLOCK 10;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_r;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_r; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        '#END_TOPOLOGY\n' \
        '#LINK_LIST\n' \
        ' T  1; (L_waveg)/2+x_shift/2-d_a/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  101; (L_waveg)/2+x_shift/2-d_a1/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  4; 0.0; 0;\n' \
        ' T  5; 0.0; 0.0;\n' \
        ' T  6; 0.0; 0;\n' \
        ' T  7; 0.0; 0;\n' \
        ' T  8; (L_waveg)/2+x_shift/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  9; 0.0; 0.0 +y_fakel_size+y_glob_shift;\n' \
        ' T 10; -xd_r; yd_r +y_fakel_size+y_glob_shift;\n' \
        '#END_LINK\n' \
        '#OUTPUT\n' \
        ' FILE name;\n' \
        ' FIELDS;\n' \
        ' TOPOLOGY;\n' \
        ' FIELD_DISTRIBUTION_M name; W_freq;  calc_time -10.0;  calc_time -5.0;\n' \
        '#END_OUTPUT\n' \
        '#END_STEP\n' \
        '#EOF'
    f.write(s)
    f.close()
def generate_TPL_fregat_fakel():
    name = 'fregat_fakel' + '.TPL'
    file = path + name
    f = open(file, 'w')
    s = f'#TMC_RT_H\n' \
        f'!#define L_input @ (   1.40)\n' \
        f'#define W_input @ (   1.00)\n' \
        f'#define W_delta @ (   0.02)\n' \
        f'#define W_freq @ (  0.6 )\n' \
        f'#define xPlus @ ( 2 )\n' \
        f'#define yPlus @ ( 0 )\n' \
        f'#define a1 @ ( W_input )\n' \
        '#define a2 @ ( 0.2 )\n' \
        '#define calc_time @ ( 200 )\n' \
        '#define name @ fregat_fakel\n' \
        '#define h1 @ ( 0.880 )\n' \
        '#define h2 @ ( 6.807 )\n' \
        '#define h3 @ ( 1.493 )\n' \
        '#define h4 @ ( 2.100 )\n' \
        '#define d1  @ ( 2.800 )\n' \
        '#define d2  @ ( 3.715 )\n' \
        '#define d3  @ ( 2.900 )\n' \
        '#define d_a @ ( (d1+d2)/2 )\n' \
        '#define h_a @ ( h1*(d_a-d1)/(d2-d1) )\n' \
        '#define d_a1 @ ( d_a + 6*(W_delta) )\n' \
        '#define h_a1 @ ( h1*(d_a1-d1)/(d2-d1) )\n' \
        '#define l_a @ ( 0.354 )\n' \
        '#define x_a @ ( h1*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_a @ ( ((d2-d1)/2)*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define l_r @ ( 0.470 )\n' \
        '#define x_r @ ( h1*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_r @ ( ((d2-d1)/2)*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define  d_r @ ( 0.177 )\n' \
        '#define xd_r @ ( ((d2-d1)/2)*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define yd_r @ ( h1*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define x_shift_fakel @ 1.0\n' \
        '#define y_shift_fakel @ 0\n' \
        '#define y_glob_shift @ -5\n' \
        '#define x_shift @ ( 7.000 )\n' \
        '#define y_fakel_size @ ( 25.000 )\n' \
        '#define W_type1 @ MAGNETIC\n' \
        '#define W_type  @ METAL\n' \
        '#define W_waveg @ (  h1 + h2 + h3 + h4 )\n' \
        '#define L_waveg @ (  d2 + 10 + x_shift + xPlus )\n' \
        '#STEP\n' \
        '#PARM\n' \
        ' ANGLE_UNIT radian;\n' \
        ' FREQ_UNIT GHz;\n' \
        ' LONG_UNIT m;\n' \
        ' TIME_UNIT ns;\n' \
        ' DELTA W_delta;\n' \
        ' TIME 0.;  calc_time;\n' \
        ' X_MIN  0.0;\n' \
        ' X_MAX L_waveg;\n' \
        ' Y_MIN  0.0;\n' \
        ' Y_MAX W_waveg+y_fakel_size;\n' \
        ' FREQ  W_freq;\n' \
        '#END_PARM\n' \
        '#TOPOLOGY\n' \
        ' BLOCK 1;\n' \
        '  INPUT_X  W_type1 ; h_a-2*W_delta; h_a; 0; 720; 1; 0.;\n' \
        ' END_B\n' \
        ' BLOCK 101;\n' \
        '  INPUT_X  W_type1 ; h_a1-2*W_delta; h_a1; 0; 720; 1; 3.141592653589;\n' \
        ' END_B\n' \
        ' BLOCK 4;\n' \
        '  RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 5;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; W_waveg+y_fakel_size-5*W_delta; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 6;\n' \
        '  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; W_waveg+y_fakel_size ;\n' \
        ' END_B\n' \
        ' BLOCK 7;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; 0; 5*W_delta ;\n' \
        ' END_B\n' \
        ' BLOCK 8;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        '   L   -(d2)/2;            h1;\n' \
        '   L   -(d2)/2;         h1+h2;\n' \
        '   L   -(d3)/2;      h1+h2+h3;\n' \
        '   L       0.0;   h1+h2+h3+h4;\n' \
        '   L    (d3)/2;      h1+h2+h3;\n' \
        '   L    (d2)/2;         h1+h2;\n' \
        '   L    (d2)/2;            h1;\n' \
        '   L    (d1)/2;           0.0;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        ' END_B\n' \
        ' BLOCK 9;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_a;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_a; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        ' BLOCK 10;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_r;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_r; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        ' BLOCK 1033;\n' \
        '  FILE fakel ; -50; 50; 0.; 100.;\n' \
        ' END_B\n' \
        '#END_TOPOLOGY\n' \
        '#LINK_LIST\n' \
        ' T  1; (L_waveg)/2+x_shift/2-d_a/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  101; (L_waveg)/2+x_shift/2-d_a1/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  4; 0.0; 0;\n' \
        ' T  5; 0.0; 0.0;\n' \
        ' T  6; 0.0; 0;\n' \
        ' T  7; 0.0; 0;\n' \
        ' T  8; (L_waveg)/2+x_shift/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  9; 0.0; 0.0 +y_fakel_size+y_glob_shift;\n' \
        ' T 10; -xd_r; yd_r +y_fakel_size+y_glob_shift;\n' \
        ' T 1033; (L_waveg)/2+x_shift/2; +y_fakel_size+y_glob_shift;\n' \
        '#END_LINK\n' \
        '#OUTPUT\n' \
        ' FILE name;\n' \
        ' FIELDS;\n' \
        f' TOPOLOGY;\n' \
        f' FIELD_DISTRIBUTION_M name; W_freq;  calc_time -10.0;  calc_time -5.0;\n' \
        f'#END_OUTPUT\n' \
        f'#END_STEP\n' \
        f'#EOF'
    f.write(s)
    f.close()
def generate_TPL_fregat_fakel_impulse():
    name = 'fregat_fakel_impulse' + '.TPL'
    file = path + name
    f = open(file, 'w')
    s = f'#TMC_RT_H\n' \
        f'!#define L_input @ (   1.40)\n' \
        f'#define W_input @ (   1.00)\n' \
        f'#define W_delta @ (   0.02)\n' \
        f'#define W_freq @ (  0.6 )\n' \
        f'#define xPlus @ ( 2 )\n' \
        f'#define yPlus @ ( 0 )\n' \
        f'#define a1 @ ( W_input )\n' \
        '#define a2 @ ( 0.2 )\n' \
        '#define calc_time @ ( 150 )\n' \
        '#define name @ fregat_impulse\n' \
        '#define h1 @ ( 0.880 )\n' \
        '#define h2 @ ( 6.807 )\n' \
        '#define h3 @ ( 1.493 )\n' \
        '#define h4 @ ( 2.100 )\n' \
        '#define d1  @ ( 2.800 )\n' \
        '#define d2  @ ( 3.715 )\n' \
        '#define d3  @ ( 2.900 )\n' \
        '#define d_a @ ( (d1+d2)/2 )\n' \
        '#define h_a @ ( h1*(d_a-d1)/(d2-d1) )\n' \
        '#define d_a1 @ ( d_a + 6*(W_delta) )\n' \
        '#define h_a1 @ ( h1*(d_a1-d1)/(d2-d1) )\n' \
        '#define l_a @ ( 0.354 )\n' \
        '#define x_a @ ( h1*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_a @ ( ((d2-d1)/2)*l_a/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define l_r @ ( 0.470 )\n' \
        '#define x_r @ ( h1*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define y_r @ ( ((d2-d1)/2)*l_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define  d_r @ ( 0.177 )\n' \
        '#define xd_r @ ( ((d2-d1)/2)*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define yd_r @ ( h1*d_r/sqrt(h1*h1+(d2-d1)*(d2-d1)/4) )\n' \
        '#define x_shift_fakel @ 1.0\n' \
        '#define y_shift_fakel @ 0\n' \
        '#define y_glob_shift @ -5\n' \
        '#define x_shift @ ( 7.000 )\n' \
        '#define y_fakel_size @ ( 25.000 )\n' \
        '#define W_type1 @ MAGNETIC\n' \
        '#define W_type  @ METAL\n' \
        '#define W_waveg @ (  h1 + h2 + h3 + h4 )\n' \
        '#define L_waveg @ (  d2 + 10 + x_shift + xPlus )\n' \
        '#STEP\n' \
        '#PARM\n' \
        ' ANGLE_UNIT radian;\n' \
        ' FREQ_UNIT GHz;\n' \
        ' LONG_UNIT m;\n' \
        ' TIME_UNIT ns;\n' \
        ' DELTA W_delta;\n' \
        ' TIME 0.;  calc_time;\n' \
        ' X_MIN  0.0;\n' \
        ' X_MAX L_waveg;\n' \
        ' Y_MIN  0.0;\n' \
        ' Y_MAX W_waveg+y_fakel_size;\n' \
        ' FREQ  W_freq;\n' \
        '#END_PARM\n' \
        '#TOPOLOGY\n' \
        ' BLOCK 1;\n' \
        '  INPUT_X  W_type1 ; h_a-2*W_delta; h_a; 0; 1.66; 1; 0.;\n' \
        ' END_B\n' \
        ' BLOCK 101;\n' \
        '  INPUT_X  W_type1 ; h_a1-2*W_delta; h_a1; 0; 1.66; 1; 3.141592653589;\n' \
        ' END_B\n' \
        ' BLOCK 4;\n' \
        '  RECT_STAT ABSORBER; L_waveg-5*W_delta; L_waveg - 1*W_delta+0.00; 0; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 5;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; W_waveg+y_fakel_size-5*W_delta; W_waveg+y_fakel_size;\n' \
        ' END_B\n' \
        ' BLOCK 6;\n' \
        '  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; W_waveg+y_fakel_size ;\n' \
        ' END_B\n' \
        ' BLOCK 7;\n' \
        '  RECT_STAT ABSORBER; 5*W_delta; L_waveg-5*W_delta; 0; 5*W_delta ;\n' \
        ' END_B\n' \
        ' BLOCK 8;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        '   L   -(d2)/2;            h1;\n' \
        '   L   -(d2)/2;         h1+h2;\n' \
        '   L   -(d3)/2;      h1+h2+h3;\n' \
        '   L       0.0;   h1+h2+h3+h4;\n' \
        '   L    (d3)/2;      h1+h2+h3;\n' \
        '   L    (d2)/2;         h1+h2;\n' \
        '   L    (d2)/2;            h1;\n' \
        '   L    (d1)/2;           0.0;\n' \
        '   L   -(d1)/2;           0.0;\n' \
        ' END_B\n' \
        ' BLOCK 9;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_a;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_a; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_a-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        ' BLOCK 10;\n' \
        '  POLYGON_STAT W_type1;\n' \
        '   L      (L_waveg)/2+x_shift/2-d_a/2;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        '   L  (L_waveg)/2+x_shift/2-d_a/2-x_r;            h_a+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r;\n' \
        '   L (L_waveg)/2+x_shift/2-d_a1/2-x_r; h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2-y_r-W_delta;\n' \
        '   L     (L_waveg)/2+x_shift/2-d_a1/2;     h_a1-3*W_delta+(W_waveg)/2-(h1+h2+h3+h4)/2;\n' \
        ' END_B\n' \
        ' BLOCK 1033;\n' \
        '  FILE fakel ; -50; 50; 0.; 100.;\n' \
        ' END_B\n' \
        '#END_TOPOLOGY\n' \
        '#LINK_LIST\n' \
        ' T  1; (L_waveg)/2+x_shift/2-d_a/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  101; (L_waveg)/2+x_shift/2-d_a1/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  4; 0.0; 0;\n' \
        ' T  5; 0.0; 0.0;\n' \
        ' T  6; 0.0; 0;\n' \
        ' T  7; 0.0; 0;\n' \
        ' T  8; (L_waveg)/2+x_shift/2; (W_waveg)/2-(h1+h2+h3+h4)/2 +y_fakel_size+y_glob_shift;\n' \
        ' T  9; 0.0; 0.0 +y_fakel_size+y_glob_shift;\n' \
        ' T 10; -xd_r; yd_r +y_fakel_size+y_glob_shift;\n' \
        ' T 1033; (L_waveg)/2+x_shift/2; +y_fakel_size+y_glob_shift;\n' \
        '#END_LINK\n' \
        '#OUTPUT\n' \
        ' FILE name;\n' \
        ' FIELDS;\n' \
        ' TOPOLOGY;\n' \
        ' FIELD_DISTRIBUTION_M name; W_freq;  calc_time -10.0;  calc_time -5.0;\n' \
        '#END_OUTPUT\n' \
        '#END_STEP\n' \
        '#EOF'
    f.write(s)
    f.close()
generate_TPL_fregat()
generate_TPL_fregat_impulse()
generate_TPL_fregat_fakel()
generate_TPL_fregat_fakel_impulse()