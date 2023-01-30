import numpy as np
from numpy import pi
import os.path
import matplotlib.pyplot as plt
def AFR_to_RP(phi,A,F,k,d):
    N = len(A)
    n = np.arange(0,N,1)-1
    S = []
    for direction in phi:
        S.append(np.sum(A*np.exp(-1j*(F+n*k*d*np.sin(direction)))))
    S = np.array(S)
    return 20 * np.log10(np.abs(S)/np.max(np.abs(S)))
def check_path(path):
    if os.path.exists(path) == False:
        os.mkdir(path)
def write_headres(file,delta,Xmin,Ymin,nX,nY,nPoint):
    def set_(name):
        f = open(name,'ab')
        f.write(bytes('\r','ascii'))
        f.close()
        f = open(name,'a')
        f.write('\n')
        f.close()
    f = open(file,'w')
    f.write('#TamicRTH_planar_DistributionDielectricPermeability_File_V2.00 2000')
    f.close()
    s =['#TopologyPrimitiv RECT_STAT ',f'#dDelta {delta}',f'#Xmin {Xmin}',f'#Ymin {Ymin}',f'#nX {nX}',f'#nY {nY}',f'#nPoint {nPoint}',
        f'#nAccuracy 4','#sNodeFormat NULL','#sValueFormat NULL']
    for tmp in s:
        set_(file)
        f = open(file,'a')
        f.write(tmp)
        f.close()
    set_(file)