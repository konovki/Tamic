import numpy as np
import matplotlib.pyplot as plt


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


operation = 1
if operation == 1:
    Rocket = load_file('./rocket_60.DAT')
    Empty = load_file('./Empty.DAT')
    Ar, Fr = take_AFR(Rocket)
    Ae, Fe = take_AFR(Empty)
    A, F = DN_Dif_operation(Ar, Ae, Fr, Fe)
    Export = place_AFR(Rocket, A, F)
    Export = ConvertToStr(Export)
    write_file('./dif60.dat', Export)
elif operation == 2:
    left = load_file('./lev.dat')
    righ = load_file('./prav.dat')
    fron = load_file('./front.dat')
    A1, F1 = take_AFR(left)
    A2, F2 = take_AFR(righ)
    A3, F3 = take_AFR(fron)
    F1, F2, F3 = np.radians(F1), np.radians(F2), np.radians(F3)
    dA, dF = DN_Sum_operation(A1, A2, F1, F2)
    A, F = DN_Sum_operation(-dA, A3, dF, F3)
    # A = np.abs(A-2)
    plt.plot(F1, label="A1")
    plt.plot(F2, label="A2")
    plt.plot(F3, label="A3")
    plt.plot(dF, label="dA")
    plt.plot(F, label="A")
    F = np.degrees(F)
    plt.legend()
    plt.show()
    Export = place_AFR(left, A, F)
    Export = ConvertToStr(Export)
    write_file('./sum.dat', Export)
