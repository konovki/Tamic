import numpy as np
import Tasks.Rocket as Task
import os.path
import properties as prop
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
for item in os.getcwd().split(prop.splitter)[:-1]:
    local_path += item + '/'
local_path += 'Tasks/Empty_rocket_DN.DAT'
path = Task.path
for i in angles:
    i = str(np.round(np.degrees(i),0))[:-2]
    Rocket = load_file(path + f'rocket_{i}.DAT')
    Empty = load_file(local_path)
    Ar, Fr = take_AFR(Rocket)
    Ae, Fe = take_AFR(Empty)
    A, F = DN_Dif_operation(Ar, Ae, Fr, Fe)
    Export = place_AFR(Rocket, A, F)
    Export = ConvertToStr(Export)
    write_file(path + f'dif_{i}.dat', Export)
