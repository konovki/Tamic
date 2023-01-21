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
    Amp = np.sqrt((A1*np.cos(F1)+A2*np.cos(F2))**2+(A1*np.sin(F1)+A2*np.sin(F2))**2)
    Fas = np.arctan2((A1*np.sin(F1)+A2*np.sin(F2)),(A1*np.cos(F1)+A2*np.cos(F2)))
    return Amp,Fas

N = 40
f = 1
l = 0.3 / f
d,k = l/2,2*pi/l
A1 = np.ones(N) * 0.5
A2 = np.ones(N) * 0.5
A3 = np.ones(N)
F1 = np.array([2*pi/N*(i - 1/2)-pi for i in range(N)])
F2 = np.array([-2*pi/N*(i - 1/2)+pi for i in range(N)])
F3 = np.zeros(N)
A4,F4 = Sum_AFR(A1,A2,F1,F2)
A5,F5 = Sum_AFR(A4,A3,F4,F3)

phi = np.linspace(-pi/2,pi/2,1500)
DN1 = AFR_to_RP(phi,A1,F1,k,d)
DN2 = AFR_to_RP(phi,A2,F2,k,d)
DN3 = AFR_to_RP(phi,A3,F3,k,d)
DN5 = AFR_to_RP(phi,A5,F5,k,d)
fig,ax = plt.subplots(2,1,figsize=(10,15))
ax[0].plot(phi,DN1)
ax[0].plot(phi,DN2)
ax[0].plot(phi,DN3)
ax[0].plot(phi,DN5)
ax[0].set_ylim(-140,0)
ax[1].plot(F1,label = 'A1')
ax[1].plot(F2,label = 'A2')
ax[1].plot(F3,label = 'A3')
ax[1].plot(F4,label = 'A4')
ax[1].plot(F5,label = 'A5')
plt.legend()
plt.show()
plt.close()
S = 1 * np.exp(1j*pi/6)
print(S,np.real(S),np.imag(S))
print(np.arctan2(np.imag(S),np.real(S)),pi/6)
print(np.abs(S))


