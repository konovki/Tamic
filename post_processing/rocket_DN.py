import time
import numpy as np
# import pandas as pd
import pickle
import Tasks.rocket_freq as Task
import os.path
import properties as prop
from multiprocessing import Pool
import matplotlib.pyplot as plt
import Library.Library as lib
def DN_Sum_operation(A1, A2, F1, F2):
    Amp = np.sqrt((A1 * np.cos(F1) + A2 * np.cos(F2)) ** 2 + (A1 * np.sin(F1) + A2 * np.sin(F2)) ** 2)
    Fas = np.arctan2((A1 * np.sin(F1) + A2 * np.sin(F2)), (A1 * np.cos(F1) + A2 * np.cos(F2)))
    return Amp, Fas
def DN_Dif_operation(A1, A2, F1, F2):
    F1, F2 = np.radians(F1), np.radians(F2)
    Amp = np.sqrt((A1 * np.cos(F1) - A2 * np.cos(F2)) ** 2 + (A1 * np.sin(F1) - A2 * np.sin(F2)) ** 2)
    Fas = np.arctan2((A1 * np.sin(F1) - A2 * np.sin(F2)), (A1 * np.cos(F1) - A2 * np.cos(F2)))
    return Amp, np.degrees(Fas)
def take_AFR(Clear):
    A = []
    F = []
    for data in Clear[3:]:
        tmp = data.split()
        A.append(float(tmp[11][:-1]))
        F.append(float(tmp[14]))
    return A, F
def place_AFR(Clear, A, F):
    NewClear = []
    NewClear.append(Clear[0])
    NewClear.append(Clear[1])
    NewClear.append(Clear[2])
    for i, data in enumerate(Clear[3:]):
        tmp = data.split()
        tmp[11] = str(A[i]) + ';'
        tmp[14] = str(F[i])
        NewClear.append(tmp)
    return NewClear
def ConvertToStr(file):
    Str = []
    Str.append(file[0])
    Str.append(file[1])
    Str.append(file[2])
    for list in file[3:]:
        tmp = ''
        for data in list:
            tmp += data + ' '
        tmp += '\n'
        Str.append(tmp)
    return Str
def load_file(file):
    f = open(file)
    stroki = []
    for data in f:
        stroki.append(data)
    f.close()
    return stroki


def write_file(file, stroki):
    f = open(file, 'w')
    f.close()
    for data in stroki:
        f = open(file, 'a')
        f.write(data)
        f.close()
angles = Task.angles

local_path = ''
for item in os.getcwd().split(lib.get_split())[:-1]:
    local_path += item + '/'
local_path += 'Tasks/Empty_rocket_DN.DAT'
path = Task.path
defname ='rocket_fakel_'#'rocket_'#rocket_fakel_
f,k = '4',6
# for i in angles:
#     i = str(np.round(np.degrees(i),0))[:-2]
#     Rocket = load_file(path + f'rocket_{i}.DAT')
#     Empty = load_file(local_path)
#     Ar, Fr = take_AFR(Rocket)
#     Ae, Fe = take_AFR(Empty)
#     A, F = DN_Dif_operation(Ar, Ae, Fr, Fe)
#     Export = place_AFR(Rocket, A, F)
#     Export = ConvertToStr(Export)
#     write_file(path + f'dif_{i}.dat', Export)
def make_deal(angle):
    i = str(np.round(np.degrees(angle),0))[:-2]
    Rocket = load_file(path + defname+f'{i}.DAT')
    Empty = load_file(local_path)
    Ar, Fr = take_AFR(Rocket)
    Ae, Fe = take_AFR(Empty)
    A, F = DN_Dif_operation(Ar, Ae, Fr, Fe)
    Export = place_AFR(Rocket, A, F)
    Export = ConvertToStr(Export)
    write_file(path +defname+ f'dif_{i}.dat', Export)
def generate_DOP(angles):
    file = path + defname+f'dif_diag.$op'
    f = open(file, 'w')
    s = f'#TMCGROTS\n' \
        f'!Graphics parameters\n' \
        f'#PointDrawFlag 0; 1;\n' \
        f'#nXType 1\n' \
        f'#nYType 1\n' \
        f'#Xmin -180\n' \
        f'#Xmax 180\n' \
        f'#Ymin -63.2203\n' \
        f'#Ymax 12.5747\n' \
        f'#nAFlagX 1\n' \
        f'#nAFlagY 1\n' \
        f'#EditorFileName  winword.exe\n' \
        f'#sDrawTextFont Times New Roman;18;0;0;0;400;0;0;0;2;\n' \
        f'#sDrawTextColor 0;\n' \
        f'#sDrawAxiesColor 0;\n' \
        f'#sDrawGridColor 0;\n' \
        f'#sDrawPointColor 255; 16711680; 0; 45; 65280; 0; 45; 16776960; 0; 45; 16711935; 0; 45; 65535; 0; 45; 8388608; 0; 45; 32768; 0; 45; 8421376; 0; 45; 5248; 0; 45; 8388736; 0; 45; 32896; 0; 45; 8421504; 0; 45; 12582912; 0; 45; 49152; 0; 45; 12632064; 0; 45; 49344; 0; 45;\n'
    for ang in angles:
        ang = str(np.round(np.degrees(ang),0))[:-2]
        s += defname+f'dif_{ang}.dat {ang} 1 1 1 1 0\n'
    f.write(s)
    f.close()
def generate_DOP_F_k_A(F,k,A):

    if len(angles) > 20:
        file = f'{path}dif_{defname}_{F}_k_{k}_1.$op'
        f = open(file, 'w')
        s = f'#TMCGROTS\n' \
            f'!Graphics parameters\n' \
            f'#PointDrawFlag 0; 1;\n' \
            f'#nXType 1\n' \
            f'#nYType 1\n' \
            f'#Xmin -180\n' \
            f'#Xmax 180\n' \
            f'#Ymin -63.2203\n' \
            f'#Ymax 12.5747\n' \
            f'#nAFlagX 1\n' \
            f'#nAFlagY 1\n' \
            f'#EditorFileName  winword.exe\n' \
            f'#sDrawTextFont Times New Roman;18;0;0;0;400;0;0;0;2;\n' \
            f'#sDrawTextColor 0;\n' \
            f'#sDrawAxiesColor 0;\n' \
            f'#sDrawGridColor 0;\n' \
            f'#sDrawPointColor 255; 16711680; 0; 45; 65280; 0; 45; 16776960; 0; 45; 16711935; 0; 45; 65535; 0; 45; 8388608; 0; 45; 32768; 0; 45; 8421376; 0; 45; 5248; 0; 45; 8388736; 0; 45; 32896; 0; 45; 8421504; 0; 45; 12582912; 0; 45; 49152; 0; 45; 12632064; 0; 45; 49344; 0; 45;\n'
        for ang in angles[:20]:
            ang = str(np.round(np.degrees(ang),0))[:-2]
            if defname == 'rocket_':
                s += f'dif_{defname}{F}_k_{k}_ang_{ang}'.replace('.','_')+f'.dat {ang} 1 0 0 0 0\n'
            else:
                s += f'dif_{defname}freq_{F}_k_{k}_ang_{ang}'.replace('.','_')+f'.dat {ang} 1 0 0 0 0\n'

        f.write(s)
        f.close()
        file = f'{path}dif_{defname}_{F}_k_{k}_2.$op'
        print(file)
        f = open(file, 'w')
        s = f'#TMCGROTS\n' \
            f'!Graphics parameters\n' \
            f'#PointDrawFlag 0; 1;\n' \
            f'#nXType 1\n' \
            f'#nYType 1\n' \
            f'#Xmin -180\n' \
            f'#Xmax 180\n' \
            f'#Ymin -63.2203\n' \
            f'#Ymax 12.5747\n' \
            f'#nAFlagX 1\n' \
            f'#nAFlagY 1\n' \
            f'#EditorFileName  winword.exe\n' \
            f'#sDrawTextFont Times New Roman;18;0;0;0;400;0;0;0;2;\n' \
            f'#sDrawTextColor 0;\n' \
            f'#sDrawAxiesColor 0;\n' \
            f'#sDrawGridColor 0;\n' \
            f'#sDrawPointColor 255; 16711680; 0; 45; 65280; 0; 45; 16776960; 0; 45; 16711935; 0; 45; 65535; 0; 45; 8388608; 0; 45; 32768; 0; 45; 8421376; 0; 45; 5248; 0; 45; 8388736; 0; 45; 32896; 0; 45; 8421504; 0; 45; 12582912; 0; 45; 49152; 0; 45; 12632064; 0; 45; 49344; 0; 45;\n'
        for ang in angles[20:]:
            ang = str(np.round(np.degrees(ang),0))[:-2]
            if defname == 'rocket_':
                s += f'dif_{defname}{F}_k_{k}_ang_{ang}'.replace('.','_')+f'.dat {ang} 1 0 0 0 0\n'
            else:
                s += f'dif_{defname}freq_{F}_k_{k}_ang_{ang}'.replace('.','_')+f'.dat {ang} 1 0 0 0 0\n'

        f.write(s)
        f.close()
    else:
        file = f'{path}{defname}{comment}.$op'
        f = open(file, 'w')
        s = f'#TMCGROTS\n' \
            f'!Graphics parameters\n' \
            f'#PointDrawFlag 0; 1;\n' \
            f'#nXType 1\n' \
            f'#nYType 1\n' \
            f'#Xmin -180\n' \
            f'#Xmax 180\n' \
            f'#Ymin -63.2203\n' \
            f'#Ymax 12.5747\n' \
            f'#nAFlagX 1\n' \
            f'#nAFlagY 1\n' \
            f'#EditorFileName  winword.exe\n' \
            f'#sDrawTextFont Times New Roman;18;0;0;0;400;0;0;0;2;\n' \
            f'#sDrawTextColor 0;\n' \
            f'#sDrawAxiesColor 0;\n' \
            f'#sDrawGridColor 0;\n' \
            f'#sDrawPointColor 255; 16711680; 0; 45; 65280; 0; 45; 16776960; 0; 45; 16711935; 0; 45; 65535; 0; 45; 8388608; 0; 45; 32768; 0; 45; 8421376; 0; 45; 5248; 0; 45; 8388736; 0; 45; 32896; 0; 45; 8421504; 0; 45; 12582912; 0; 45; 49152; 0; 45; 12632064; 0; 45; 49344; 0; 45;\n'
        for ang in angles:
            ang = str(np.round(np.degrees(ang),0))[:-2]
            if defname == 'rocket_':
                s += f'{defname}{F}_k_{k}_ang_{ang}'.replace('.','_')+f'.dat {ang} 1 1 1 1 0\n'
            else:
                s += f'{defname}freq_{F}_k_{k}_ang_{ang}'.replace('.','_')+f'.dat {ang} 1 1 1 1 0\n'
        f.write(s)
        f.close()
def get_knd(angles):
    knd = []
    knd_dB = []
    Emitted_Power = []
    for ang in angles:
        ang = str(np.round(np.degrees(ang),0))[:-2]
        file = path + defname+f'dif_{ang}.dat~~tempDia1'
        f = open(file, 'r')
        for tmp in f:
            tmp = tmp.split()
            try:
                if (tmp[0] == 'KND'):
                    if tmp[-1] == 'dB;':
                        knd_dB.append(float(tmp[2]))
                    else:
                        knd.append(float(tmp[2][:-1]))
                if tmp[0] == 'Emitted':
                    Emitted_Power.append(float(tmp[3][:-1]))
            except: pass
    angles = np.degrees(np.array(angles))
    fig,ax = plt.subplots(3,1)
    ax[0].set_title('Ракета c двигателем')
    ax[0].plot(angles,knd,label = 'knd')
    ax[1].plot(angles,knd_dB,label = 'knd_dB')
    ax[2].plot(angles,Emitted_Power,label = 'Emitted Power')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.show()
    df = pd.DataFrame(data={'angles':angles,'knd':knd,'knd_dB':knd_dB,'EM':Emitted_Power})
    df.to_csv('./fakel.csv')
    return 0
def get_knd2(angles,fr,k,title):
    knd = []
    knd_dB = []
    Emitted_Power = []
    for ang in angles:
        ang = str(np.round(np.degrees(ang),0))[:-2]
        # name = f'dif_{defname}{f}_'+f'k_{k}_ang_{ang}'+'.DAT~~tempDia1'
        # name = f'test.DAT~~tempDia1'
        if defname == 'rocket_':
            name = f'dif_{defname}{fr}_k_{k}_'+f'ang_{ang}'+'.DAT~~tempDia0'
        else:
            name = f'dif_{defname}freq_{fr}_'+f'k_{k}_ang_{ang}'+'.DAT~~tempDia0'

        print(path + name)
        file = path + name
        # file = f'{path}{name}'
        f = open(path + name, 'r')
        for tmp in f:
            tmp = tmp.split()
            try:
                if (tmp[0] == 'KND'):
                    if tmp[-1] == 'dB;':
                        knd_dB.append(float(tmp[2]))
                    else:
                        knd.append(float(tmp[2][:-1]))
                if tmp[0] == 'Emitted':
                    Emitted_Power.append(float(tmp[3][:-1]))
            except: pass
    angles = np.degrees(np.array(angles))
    fig,ax = plt.subplots(3,1)
    ax[0].set_title(title)
    ax[0].plot(angles,knd,label = 'knd')
    ax[1].plot(angles,knd_dB,label = 'knd_dB')
    ax[2].plot(angles,Emitted_Power,label = 'Emitted Power')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.show()
    with open(path+defname+f'{fr}_k_{k}_knd.pkl','wb') as f:
        pickle.dump(knd,f)
    with open(path+defname+f'{fr}_k_{k}_knd_db.pkl','wb') as f:
        pickle.dump(knd_dB,f)
    with open(path+defname+f'{fr}_k_{k}_Emitted_Power.pkl','wb') as f:
        pickle.dump(Emitted_Power,f)
    return 0
def drow_results():
    df1 = pd.read_csv('fakel.csv')
    df2 = pd.read_csv('no_fakel.csv')
    angles = df1.angles
    fig,ax = plt.subplots(3,1)
    ax[0].set_title('knd')
    ax[0].plot(angles,df1.knd,label = 'с факелом')
    ax[0].plot(angles,df2.knd,label = 'без факела')
    ax[1].set_title('knd_dB')
    ax[1].plot(angles,df1.knd_dB,label = 'с факелом')
    ax[1].plot(angles,df2.knd_dB,label = 'без факела')
    ax[2].set_title('Emitted Power')
    ax[2].plot(angles,df1.EM,label = 'с факелом')
    ax[2].plot(angles,df2.EM,label = 'без факела')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.show()
def make_deal2(par):
    angle,k,f = par[0],par[1],par[2]
    i = str(np.round(np.degrees(angle),0))[:-2]
    if defname == 'rocket_':
        name = f'{defname}{f}_k_{k}_ang_{i}'.replace('.','_')+f'.DAT'
    else:
        name = f'{defname}freq_{f}_k_{k}_ang_{i}'.replace('.','_')+f'.DAT'
    print('open',path + name)
    Rocket = load_file(path + name)
    print(path+f'{f}_Empty.dat')
    Empty = load_file(path+f'{f}_Empty.dat')
    Ar, Fr = take_AFR(Rocket)
    Ae, Fe = take_AFR(Empty)
    # A, F = DN_Dif_operation(Ar, Ae, Fr, Fe)
    A, F = DN_Dif_operation(Ae, Ar, Fe, Fr)
    Export = place_AFR(Rocket, A, F)
    Export = ConvertToStr(Export)
    write_file(path +'dif_n_'+name, Export)
def make_deal3(par):
    angle,k,f = par[0],par[1],par[2]
    i = str(np.round(np.degrees(angle),0))[:-2]
    if defname == 'rocket_':
        name = f'{defname}{f}_k_{k}_ang_{i}'.replace('.','_')+f'.DAT'
    else:
        name = f'{defname}freq_{f}_k_{k}_ang_{i}'.replace('.','_')+f'.DAT'
    print('open',path + name)
    f = open(path + name, 'r')
    s = []
    for tmp in f:
        tmp = tmp.split()
        try:
            if tmp[3] == '79900':
                tmp[11] = '0;'
                tmp[14] = '0;'
        except: pass
        exp = ''
        for symbol in tmp:
            exp += symbol + ' '
        exp = exp[:-1]+'\n'
        s.append(exp)
    f.close()
    write_file(path +'dif_'+name, s)
def show_res(name1,name2):
    df1 = pd.read_csv(path+name1)
    print(df1)
    df2 = pd.read_csv(path+name2)
    print(df2)
    fig,ax = plt.subplots(3,1)
    ax[0].set_title('КНД')
    ax[1].set_title('КНД дБ')
    ax[2].set_title('Изл. мощ.')
    ax[0].plot(angles,df1.knd,label = '1')
    ax[0].plot(angles,df2.knd,label = '2')
    ax[1].plot(angles,df1.knd_dB,label = '1')
    ax[1].plot(angles,df2.knd_dB,label = '1')
    ax[2].plot(angles,df1.EM,label = '1')
    ax[2].plot(angles,df2.EM,label = '2')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.show()
def get_results():
    s_l = ''
    labels = [1,2,3,4]
    symbol = ['-','.','--','-.']
    k = [10,6,6,6]
    for s in path.split('/')[:-2]:
        s_l += s+'/'
    fig,ax = plt.subplots(3,1)
    ax[0].set_title('knd')
    ax[1].set_title('knd_dB')
    ax[2].set_title('Emitted Power')
    for i in range(4):
        s1 = s_l
        s1 += f'Rocket_freq_{i+1}GHz'
        with open(s1+f'/rocket_{i+1}_k_{k[i]}_knd.pkl','rb') as f:
            knd = pickle.load(f)
            ax[0].plot(angles,knd,symbol[i],label = f'без факела {i+1} GHz')
        with open(s1+f'/rocket_{i+1}_k_{k[i]}_knd_db.pkl','rb') as f:
            knd_db = pickle.load(f)
            ax[1].plot(angles,knd_db,symbol[i],label = f'без факела {i+1} GHz')
        with open(s1+f'/rocket_{i+1}_k_{k[i]}_Emitted_Power.pkl','rb') as f:
            Emitted_Power = pickle.load(f)
            ax[2].plot(angles,Emitted_Power,symbol[i],label = f'без факела {i+1} GHz')
            #с фак
        with open(s1+f'/rocket_fakel_{i+1}_k_{k[i]}_knd.pkl','rb') as f:
            knd = pickle.load(f)
            ax[0].plot(angles,knd,symbol[i],label = f'с факелом {i+1} GHz')
        with open(s1+f'/rocket_fakel_{i+1}_k_{k[i]}_knd_db.pkl','rb') as f:
            knd_db = pickle.load(f)
            ax[1].plot(angles,knd_db,symbol[i],label = f'с факелом {i+1} GHz')
        with open(s1+f'/rocket_fakel_{i+1}_k_{k[i]}_Emitted_Power.pkl','rb') as f:
            Emitted_Power = pickle.load(f)
            ax[2].plot(angles,Emitted_Power,symbol[i],label = f'с факелом {i+1} GHz')    # без фак


    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.show()
def delete_side(angles,fr,k):
    knd = []
    knd_dB = []
    Emitted_Power = []
    for ang in angles:
        ang = str(np.round(np.degrees(ang),0))[:-2]
        if defname == 'rocket_':
            name = f'{defname}{fr}_k_{k}_'+f'ang_{ang}'+'.DAT'
        else:
            name = f'{defname}freq_{fr}_'+f'k_{k}_ang_{ang}'+'.DAT'

        file = path + name
        f = open(path + name, 'r')
        s = []
        for tmp in f:
            tmp = tmp.split()
            try:
                if tmp[3] == '79900':
                    tmp[11] = '0;'
                    tmp[14] = '0;'
            except: pass
            exp = ''
            for symbol in tmp:
                exp += symbol + ' '
            exp = exp[:-1]+'\n'
            s.append(exp)
        f.close()
        write_file(path +'delete_side'+name, s)
        exit()
def generate_bat(F,k,angles):
    file = path + f'run_TMC_DN_{F}.bat'
    f = open(file, 'w')
    f.close()
    if len(angles) > 20:
        f = open(file, 'a')
        file1 = f'dif_{defname}_{F}_k_{k}_1.$op'
        s = f'start TMC_DN.exe {file1}\n'
        f.write(s)
        f.close()
        f = open(file, 'a')
        file2 = f'dif_{defname}_{F}_k_{k}_2.$op'
        s = f'start TMC_DN.exe {file2}\n'
        f.write(s)
        f.close()
        f = open(file, 'a')
        f.write(f'exit')
        f.close()
    else:
        f = open(file, 'a')
        file1 = f'dif_{defname}_{F}_k_{k}.$op'
        s = f'start TMC_DN.exe {file1}\n'
        f.write(s)
        f.close()
        f = open(file, 'a')
        f.write(f'exit')
        f.close()
potoki = 12


if __name__ == '__main__':
    #1
    # par = []
    # for ang in angles:
    #     par.append([ang,k,f])
    # Pool(potoki).map(make_deal3, par)
    # Pool(potoki).close()
    #2
    # time.sleep(5)
    # generate_DOP_F_k_A(f,k,angles)
    # time.sleep(5)
    # generate_bat(f,k,angles)
    # time.sleep(1)
    print(path+f'run_TMC_DN_{f}.bat')
    os.chdir(path)
    os.system(f'start run_TMC_DN_{f}.bat')
    os.system(f'exit')
    #3
    flag  = False
    while flag == False:
        try:
            get_knd2(angles,f,k,'')
            flag = True
        except: time.sleep(10)
    # show_res('rocket_0_5_k_10.csv','rocket_fakel_0_5_k_10.csv')
    get_results()
    lib.beep()
