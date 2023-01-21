import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
def AFR_to_RP(phi,A,F,k,d):
    N = len(A)
    n = np.arange(0,N,1)-1
    S = []
    for direction in phi:
        S.append(np.sum(A*np.exp(-1j*(F+n*k*d*np.sin(direction)))))
    S = np.array(S)
    return 20 * np.log10(np.abs(S)/np.max(np.abs(S)))
def Sum_AFR(A1,A2,F1,F2):
    S1,S2 = A1 +1j*F1,A2 +1j*F2
    S = S1+S2
    print(S,np.imag(S),np.real(S))
    return np.sqrt(np.real(S)**2+np.imag(S)**2),np.arctan2(np.imag(S),np.real(S))
A1,A2 = 0.5,0.85
F1,F2 = pi/2,pi/4
S1,S2 = A1 +1j*F1,A2 +1j*F2
Amp = np.sqrt((A1*np.cos(F1)+A2*np.cos(F2))**2+(A1*np.sin(F1)+A2*np.sin(F2))**2)
Fas = np.arctan2((A1*np.sin(F1)+A2*np.sin(F2)),(A1*np.cos(F1)+A2*np.cos(F2)))
print('correct',Amp,Fas)
print(Sum_AFR(A1,A2,F1,F2))