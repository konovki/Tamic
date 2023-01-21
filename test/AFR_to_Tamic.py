import pandas as pd
import numpy as np
df = pd.read_csv('./res.csv')
print(np.min(np.degrees(df.F)),np.max(np.degrees(df.F)))
N = len(df.F)
fr = 10
def f(x):
    s = f'1   X = {x.Rx*1000} mm; Y = {x.Ry*1000} mm; Z = {x.Rz*1000} mm; Amplitude = {x.A}; Faza = {np.degrees(x.F)} degree; FiMin = -180 degree; FiMax = 180 degree; 1.\n'
    return s
list = []
list.append(f'Frequence = {fr} Ghz;\n')
list.append(f'N = {N};\n')
# list.append(f'1+cos(f)*1\n')
list.append(f'1\n')
list.append(df.apply(lambda x:f(x),axis=1))
file = './' + 'COLSP.DAT'
f = open(f'./{file}', 'w')
f.close()
for i,data in enumerate(list):
    if i <3:
        f = open(f'./{file}', 'a')
        f.write(data)
        f.close()
    else:
        for d in data:
            f = open(f'./{file}', 'a')
            f.write(d)
            f.close()
print('done')
