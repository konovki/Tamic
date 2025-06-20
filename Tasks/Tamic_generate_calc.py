import time
import numpy as np
import matplotlib.pyplot as plt
import math as m
import struct
import pickle
import random
from multiprocessing import Pool
# import properties as prop
import datetime as dt

import pandas as pd

import Library as lib
import os.path
import os
import psutil
from psutil._common import bytes2human

path = ''
for item in os.getcwd().split(lib.get_split())[:-2]:
    path += item + '/'
path += 'Results/'
lib.check_path(path)
path += 'Rocket_freq_ALL/'
lib.check_path(path)
local_path = os.getcwd()
global k
potoki = 23
angles = np.radians(np.arange(0, 360, 10))
freq_l = [[0.7, 10], [0.8, 10], [0.9, 10], [1, 10], [1.5,10],[2,10],[2.5,10],[3,6],[3.5,6],[4,6],[4.5,6]]
par = []
Gen_eps = True
Empty = False

def generate_Empty_TPL(par, file_name):
    freq, kk = par[0], par[1]
    delta = 3 * (10 ** 8) / (freq * (10 ** 9)) / kk
    delta_tamic = delta * 1000
    name = file_name + '.TPL'
    file = path + name
    f = open(file, 'w')
    s = f'#TMC_RT_H\n' \
        f'#define Xsize @ ( 80000.00)\n' \
        f'#define W_file  @ {file_name}\n' \
        f'#define W_freq  @ (  {freq})\n' \
        f'#define W_eps @ ( 4. )\n' \
        f'#define W_time @ ( 350. )\n' \
        f'#define Ysize @ (  50000.00)\n' \
        f'#define W_input @ (  Ysize-1000.00)\n' \
        f'#define W_delta @ (   {delta_tamic})\n' \
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
        f' X_MAX Xsize;\n' \
        f' Y_MIN  0.0;\n' \
        f' Y_MAX Ysize;\n' \
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
        f' BLOCK 8;\n' \
        f'  RECT_STAT ABSORBER; Xsize-5*W_delta; Xsize - 1*W_delta+0.00; 0.0; Ysize-2*W_delta;\n' \
        f' END_B\n' \
        f' BLOCK 9;\n' \
        f'  RECT_STAT ABSORBER; 1*W_delta; Xsize-2*W_delta; W_delta; 5*W_delta;\n' \
        f' END_B\n' \
        f' BLOCK 10;\n' \
        f'  RECT_STAT ABSORBER; 1*W_delta; Xsize; Ysize-5*W_delta; Ysize;\n' \
        f' END_B\n' \
        f' BLOCK 11;\n' \
        f'  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; Ysize/2-a1/2-4*W_delta;\n' \
        f' END_B\n' \
        f' BLOCK 12;\n' \
        f'  RECT_STAT ABSORBER; 0.0; 5*W_delta; Ysize/2+a1/2+4*W_delta; Ysize;\n' \
        f' END_B\n' \
        f'\n' \
        f'#END_TOPOLOGY\n' \
        f'\n' \
        f'#LINK_LIST\n' \
        f' T  5; 0.0; Ysize/2;\n' \
        f' T  6; 0.0; Ysize/2-W_delta;\n' \
        f' T  7; 0.0; Ysize/2;\n' \
        f' T  8; 0.0; 0.0;\n' \
        f' T  9; 0.0; 0.0;\n' \
        f' T 10; 0.0; 0.0;\n' \
        f' T 11; 0.0; 0.0;\n' \
        f' T 12; 0.0; 0.0;\n' \
        f'#END_LINK\n' \
        f'\n' \
        f'#OUTPUT\n' \
        f' FILE W_file;\n' \
        f' !FIELDS;\n' \
        f' TOPOLOGY;\n' \
        f' FIELD_DISTRIBUTION_M  W_file; W_freq; W_time-10; W_time-5;\n' \
        f'#END_OUTPUT\n' \
        f'\n' \
        f'#END_STEP\n' \
        f'\n' \
        f'#EOF'
    f.write(s)
    f.close()


def calc_rasp(par):
    t = dt.datetime.now()
    angle, freq, kk, Only_rocket = par[0], par[1], par[2], par[3]
    delta = 3 * (10 ** 8) / (freq * (10 ** 9)) / kk
    delta_tamic = delta * 1000
    if Only_rocket == False:
        def_name = f'rocket_fakel_{freq}_k_{kk}_ang_'.replace('.', '_')
        X_boundary = 5  # z
        Y_boundary = 7  # x
        Xmin = -5  # z
        Ymin = -15  # x
    elif Only_rocket == True:
        def_name = f'rocket_{freq}_k_{kk}_ang_'.replace('.', '_')
        X_boundary = 5  # z
        Y_boundary = 7  # x
        Xmin = -5  # z
        Ymin = -3  # x
    x1, x2, x3 = 0, Xmin, X_boundary
    y1, y2, y3 = Y_boundary, Ymin, Ymin
    angle = - angle

    x11, y11 = x1 * np.cos(angle) + y1 * np.sin(angle), -x1 * np.sin(angle) + y1 * np.cos(angle)
    x22, y22 = x2 * np.cos(angle) + y2 * np.sin(angle), -x2 * np.sin(angle) + y2 * np.cos(angle)
    x33, y33 = x3 * np.cos(angle) + y3 * np.sin(angle), -x3 * np.sin(angle) + y3 * np.cos(angle)
    Xmin, X_boundary = min(x11, x22, x33), max(x11, x22, x33)
    Ymin, Y_boundary = min(y11, y22, y33), max(y11, y22, y33)
    dY, dX = 1000 * (- Ymin), 1000 * (-Xmin)
    shiftX, shiftY = -1000 * np.round((Xmin + X_boundary) / 2, 0), -1000 * np.round((Ymin + Y_boundary) / 2, 0)
    angle = - angle
    x = np.arange(Xmin, X_boundary + delta, delta)
    y = np.arange(Ymin, Y_boundary + delta, delta)
    node = 0
    nX = len(x)
    nY = len(y)
    nPoint = nX * nY
    angle_str = str(np.round(np.degrees(angle), 0))[:-2]

    def generate_TPL(angle, freq=freq):
        name = def_name + angle + '.TPL'
        eps_name = def_name + angle
        file = path + name
        f = open(file, 'w')
        s = f'#TMC_RT_H\n' \
            f'#define Xsize @ ( 80000.00)\n' \
            f'#define W_file  @ {eps_name}\n' \
            f'#define W_freq  @ (  {freq})\n' \
            f'#define W_eps @ ( 4. )\n' \
            f'#define W_time @ ( 350. )\n' \
            f'#define Ysize @ (  50000.00)\n' \
            f'#define W_input @ (  Ysize-1000.00)\n' \
            f'#define W_delta @ (   {delta_tamic})\n' \
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
            f' X_MAX Xsize;\n' \
            f' Y_MIN  0.0;\n' \
            f' Y_MAX Ysize;\n' \
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
            f'  FILE W_file; 0; {10000}; 0; {20000};\n' \
            f' END_B\n' \
            f' BLOCK 8;\n' \
            f'  RECT_STAT ABSORBER; Xsize-5*W_delta; Xsize - 1*W_delta+0.00; 0.0; Ysize-2*W_delta;\n' \
            f' END_B\n' \
            f' BLOCK 9;\n' \
            f'  RECT_STAT ABSORBER; 1*W_delta; Xsize-2*W_delta; W_delta; 5*W_delta;\n' \
            f' END_B\n' \
            f' BLOCK 10;\n' \
            f'  RECT_STAT ABSORBER; 1*W_delta; Xsize; Ysize-5*W_delta; Ysize;\n' \
            f' END_B\n' \
            f' BLOCK 11;\n' \
            f'  RECT_STAT ABSORBER; 0.0; 5*W_delta; 0.0; Ysize/2-a1/2-4*W_delta;\n' \
            f' END_B\n' \
            f' BLOCK 12;\n' \
            f'  RECT_STAT ABSORBER; 0.0; 5*W_delta; Ysize/2+a1/2+4*W_delta; Ysize;\n' \
            f' END_B\n' \
            f'\n' \
            f'#END_TOPOLOGY\n' \
            f'\n' \
            f'#LINK_LIST\n' \
            f' T  5; 0.0; Ysize/2;\n' \
            f' T  6; 0.0; Ysize/2-W_delta;\n' \
            f' T  7; 0.0; Ysize/2;\n' \
            f' T  8; 0.0; 0.0;\n' \
            f' T  9; 0.0; 0.0;\n' \
            f' T 10; 0.0; 0.0;\n' \
            f' T 11; 0.0; 0.0;\n' \
            f' T 12; 0.0; 0.0;\n' \
            f' T 55; Xsize/2 +({shiftX}); Ysize/2 +({shiftY});\n' \
            f'#END_LINK\n' \
            f'\n' \
            f'#OUTPUT\n' \
            f' FILE W_file;\n' \
            f' !FIELDS;\n' \
            f' TOPOLOGY;\n' \
            f' FIELD_DISTRIBUTION_M  W_file; W_freq; W_time-10; W_time-5;\n' \
            f'#END_OUTPUT\n' \
            f'\n' \
            f'#END_STEP\n' \
            f'\n' \
            f'#EOF'
        f.write(s)
        f.close()

    gen_TPL_fakel = True
    if gen_TPL_fakel == True:
        generate_TPL(angle_str)
    if Gen_eps == True:
        name = f'{def_name}{angle_str}.eps'
        file = path + name

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

        lib.write_headres(file, delta, Xmin, Ymin, nX, nY, nPoint)
        count = 0
        f = open(file, 'ab')
        for Yi, Y in enumerate(y):
            for Xi, X in enumerate(x):
                X1, Z1 = X * np.cos(angle) + Y * np.sin(angle), -X * np.sin(angle) + Y * np.cos(angle)
                if Only_rocket == True:
                    a = set_a(lib.Epsilon_no_fakel([X1, Z1], ray_freq=freq), Xi, Yi)
                else:
                    a = set_a(lib.Epsilon([X1, Z1], ray_freq=freq), Xi, Yi)
                f.write(bytearray(struct.pack("!f", a))[::-1])
                f.write(node.to_bytes(4, 'little'))
                node += 1
                count += 1
        f.close()
    lib.beep()


def get_names(data):
    Tpl = []
    for par in data:
        angle, freq, kk, Only_rocket = par[0], par[1], par[2], par[3]
        if Only_rocket == False:
            def_name = f'rocket_fakel_{freq}_k_{kk}_ang_'.replace('.', '_')
        elif Only_rocket == True:
            def_name = f'rocket_{freq}_k_{kk}_ang_'.replace('.', '_')
        angle_str = str(np.round(np.degrees(angle), 0))[:-2]
        Tpl.append(def_name + angle_str + '.TPL')
    return Tpl


def make_deal3(par):
    angle, f, k, Only_rocket = par[0], par[1], par[2], par[3]
    i = str(np.round(np.degrees(angle), 0))[:-2]
    defname = lib.get_defname(Only_rocket)
    name = f'{defname}{f}_k_{k}_ang_{i}'.replace('.', '_') + f'.DAT'
    f = open(path + name, 'r')
    X = []
    for tmp in f:
        tmp = tmp.split()
        try:
            X.append(int(tmp[3]))
        except:
            pass
    f.close()
    Xmax = str(np.max(np.array(X)))
    s = []
    f = open(path + name, 'r')
    for tmp in f:
        tmp = tmp.split()
        try:
            if tmp[3] == Xmax:
                tmp[11] = '0;'
                tmp[14] = '0;'
        except:
            pass
        exp = ''
        for symbol in tmp:
            exp += symbol + ' '
        exp = exp[:-1] + '\n'
        s.append(exp)
    f.close()
    lib.write_file(path + 'dif_' + name, s)


def create_name(name, defname, dop1):
    tmp = name.replace(f'{defname}{dop1[1:]}', '*')
    return f'dif_{name[:-3]}dat {tmp[5:-4]} 1 1 1 1 0\n'


def get_dop_name(name, defname):
    tmp = name.replace(defname, '*')
    dop1 = ''
    for symbol in tmp:
        if symbol != 'a':
            dop1 += symbol
        else:
            break
    return dop1


def generate_DOP_F_k_A(names, Only_rocket):
    s0 = lib.s0
    defname = lib.get_defname(Only_rocket)
    control_names = []
    dop1 = get_dop_name(names[0], defname)
    for name in names:
        control_names.append(path + create_name(name, defname, dop1).split('.')[0] + '.dat~~tempDia1')
    control_names1 = control_names
    drow_names = names
    del_files = []
    while len(drow_names) > 0:
        print(psutil.cpu_percent())
        if (float(bytes2human(psutil.virtual_memory().available)[:-1]) > 1.5) and (psutil.cpu_percent() < 95):
            current_name = drow_names[0]
            drow_names = drow_names[1:]
            print(current_name)
            tmp_file = f'{path}tmp_dif_{current_name[:-4]}.$op'
            del_files.append(tmp_file)
            f = open(tmp_file, 'w')
            s = s0
            s += create_name(current_name, defname, dop1)
            os.chdir(path)
            os.system(f'start TMC_DN.exe  {tmp_file}')
            f.write(s)
            f.close()

    while len(control_names1) > 0:
        if (os.path.exists(control_names1[0]) == False):
            print('looking for', control_names1[0])
            time.sleep(t_sleep)
        else:
            control_names1 = control_names1[1:]
    if len(names) > 20:
        file = f'{path}dif_{defname}{dop1[1:]}.$op'
        f = open(file, 'w')
        s = s0
        for name in names[:20]:
            s += create_name(name, defname, dop1)
        f.write(s)
        f.close()
        file2 = f'{path}dif_{dop1[1:]}2.$op'
        f = open(file2, 'w')
        s = s0
        for name in names[20:]:
            s += create_name(name, defname, dop1)
        f.write(s)
        f.close()
    else:
        tmp = names[0].replace(defname, '*')
        dop1 = ''
        for symbol in tmp:
            if symbol != 'a':
                dop1 += symbol
            else:
                break
        file = f'{path}dif_{dop1[1:]}.$op'
        f = open(file, 'w')
        s = s0
        for name in names:
            s += create_name(name, defname, dop1)
        f.write(s)
        f.close()
        os.chdir(path)
        os.system(f'start TMC_DN.exe  {file}')
    print('freq', current_Freq[0], 'RP created')
    time.sleep(t_wait)
    dat = lib.get_knd2(names, defname, names[-1][:-3].split('_ang')[0], path)
    os.system(f'taskkill /IM TMC_DN.exe')
    lib.del_file_list(del_files)
    time.sleep(t_wait)
    return dat


i = 0
t_sleep = 10
t_wait = 2
if __name__ == '__main__':
    LOG = f'{path}log.txt'
    Log_file = open(LOG, 'w')
    Log_file.close()
    for current_Freq in freq_l:
        t1 = dt.datetime.now()
        generate_Empty_TPL(current_Freq, f'{current_Freq[0]}_Empty')
        par = []
        for ang in angles:
            par.append([ang, current_Freq[0], current_Freq[1], True])
        for ang in angles:
            par.append([ang, current_Freq[0], current_Freq[1], False])
        names = get_names(par)
        names_MEM = names
        names_chek = names
        lib.print_log(LOG, f'generation TPL EPS freq {current_Freq} time {dt.datetime.now()}')
        Pool(potoki).map(calc_rasp, par)
        Pool(potoki).close()
        lib.print_log(LOG, f'delete files {current_Freq} time {dt.datetime.now()}')
        lib.clear_files(path, names)
        while len(names) > 0:
            print(psutil.cpu_percent())
            if (float(bytes2human(psutil.virtual_memory().available)[:-1]) > 1.5) and (psutil.cpu_percent() < 95):
                lib.print_log(LOG, f'started {names[0]}')
                os.chdir(path)
                os.system(f'start tmc_rth.exe /Ar {names[0]}')
                names = names[1:]
                time.sleep(1)
            else:
                time.sleep(t_sleep)
        while len(names_chek) > 0:
            control_name = names_chek[0]
            check_path_list = [path + control_name[:-3] + 'AMP', path + control_name[:-3] + 'FAZ',
                               path + control_name[:-3] + 'DAT']
            if (os.path.exists(check_path_list[0]) == False) or (os.path.exists(check_path_list[1]) == False) or (
                    os.path.exists(check_path_list[1]) == False):
                time.sleep(t_sleep)
            else:
                names_chek = names_chek[1:]
        lib.print_log(LOG, f'freq {current_Freq[0]} finished time {dt.datetime.now()}')
        time.sleep(t_wait)
        os.system(f'taskkill /IM tmc_rth.exe')
        time.sleep(t_wait)
        lib.print_log(LOG, f'freq {current_Freq[0]} update DAT time {dt.datetime.now()}')
        Pool(potoki).map(make_deal3, par)
        Pool(potoki).close()
        time.sleep(t_wait)
        N = int(len(names_MEM) / 2)
        namesRocket = names_MEM[0:N]
        lib.print_log(LOG, f'freq {current_Freq[0]} generate dop time {dt.datetime.now()}')
        dat1 = generate_DOP_F_k_A(namesRocket, True)
        time.sleep(t_wait)
        namesRFakel = names_MEM[N:]
        dat2 = generate_DOP_F_k_A(namesRFakel, False)
        lib.drow_pic(dat1, dat2, current_Freq, path)
        time.sleep(t_wait)
        lib.beep()
        lib.print_log(LOG, f'freq {current_Freq[0]} finished {dt.datetime.now()} time spend {dt.datetime.now() - t1}')
    # os.system(f'shutdown -s -t 60')
