import numpy as np
import matplotlib.pyplot as plt


def DN_operation(Aclear, Flcear, Arocket, Frocket):
    Fclear, Frocket = np.radians(Flcear), np.radians(Frocket)
    Clear, Rocket = Aclear * np.exp(1j * Fclear), Arocket * np.exp(1j * Frocket)
    Dif = Rocket + Clear
    AmpDif = np.abs(Dif)
    # AmpDif = np.sqrt(np.imag(Dif)**2+np.real(Dif)**2)
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


prav = load_file('./prav.dat')
Aprav, Fprav = take_AFR(prav)
# Rocket = load_file('./rocket_90.DAT')
lev = load_file('./lev.dat')
Alev, Flev = take_AFR(lev)
A, F = DN_operation(Aprav, Fprav, Alev, Flev)
front = load_file('./front.dat')
Afront, Ffront = take_AFR(front)
Aex, Fex = DN_operation(Afront, Ffront, A, F)
Export = place_AFR(prav, Aex, Fex)

plt.plot(Flev, label='Flev')
plt.plot(Fprav, label='Fprav')
plt.plot(Ffront, label='Ffront')
plt.plot(Fex, label='Fex')
plt.plot(F, label='F')
plt.legend()
plt.show()
Export = ConvertToStr(Export)
write_file('./New1.dat', Export)
