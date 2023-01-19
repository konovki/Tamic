import numpy as np
import matplotlib.pyplot as plt


def DN_operation(A1, F1, A2, F2):
    F1, F2 = np.radians(F1), np.radians(F2)
    S1, S2 = A1 * np.exp(1j * F1), A2 * np.exp(1j * F2)
    Dif = S1 - S2
    # AmpDif = np.abs(Dif)
    AmpDif = np.sqrt(np.imag(Dif)**2+np.real(Dif)**2)
    FasDif = np.degrees(np.imag(Dif))
    return AmpDif, FasDif





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


prav = load_file('./rocket_90.DAT')
Aprav, Fprav = take_AFR(prav)
# Rocket = load_file('./rocket_90.DAT')
lev = load_file('./Empty.DAT')
Alev, Flev = take_AFR(lev)
A, F = DN_operation(Aprav, Fprav, Alev, Flev)

Export = place_AFR(prav, A, F)


Export = ConvertToStr(Export)
write_file('./New1.dat', Export)
