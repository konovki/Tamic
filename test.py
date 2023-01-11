import numpy as np
import matplotlib.pyplot as plt

f = open('r7.DAT_Temp1', 'r')
angle = []
A = []
F = []
for data in f:
    tmp = data.split()
    angle.append(float(tmp[0]))
    A.append(float(tmp[1]))
    F.append(float(tmp[3]))
angle = np.array(angle)
A = np.array(A)
F = np.array(F)
fig = plt.figure(figsize=(9, 7))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2,sharex = ax1)
DN = 20*np.log10(A/np.max(A))
ax1.plot(angle,DN)
ax2.plot(angle,F)
plt.show()
plt.close()