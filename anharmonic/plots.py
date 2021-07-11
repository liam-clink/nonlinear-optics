import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import scipy
from scipy import fftpack

matplotlib.use('TkAgg')

E = np.load('E.npy')
x = np.load('x.npy')
z = np.load('z.npy')

temp = []
for i in range(len(x)):
    temp.append(x[i][0])
plt.plot(temp)
plt.ylim((-5, 5))
plt.show()

figE, axE = plt.subplots(figsize=(20, 10))
axE.set_xlabel('z')
axE.set_ylim((-5, 5))
axE.set_ylabel('Ex field')
axE.tick_params(axis='y', labelcolor='b')

axX = axE.twinx()
scale = 1. ** -1
axX.set_ylim((-1 / scale, 1 / scale))
axX.set_ylabel('Dipole Displacement')
axX.tick_params(axis='y', labelcolor='g')

plt.tight_layout()

# Number of sample points
N = len(E[0])
# sample spacing
T = z[1] - z[0]
figfft, axfft = plt.subplots(figsize=(20, 10))

plt.tight_layout()

print(len(E))
for i in range(len(E)):
    curveE, = axE.plot(z, E[i], 'b')
    curveX, = axE.plot(z, scale * x[i], 'g')
    figE.savefig('E/' + str(i).zfill(int(np.ceil(np.log10(len(E))))) + '.png')

    curveE.remove()
    del curveE
    curveX.remove()
    del curveX

    Efft = scipy.fftpack.fft(E[i])
    Xfft = scipy.fftpack.fft(x[i])
    zfft = np.linspace(0.0, 1.0 / (2.0 * T), N // 2)

    max_index = 200
    curveXfft, = axfft.plot(zfft, 2.0 / N * np.abs(Xfft[:N // 2]), 'g')
    curveEfft, = axfft.plot(zfft, 2.0 / N * np.abs(Efft[:N // 2]), 'b')
    axfft.set_yscale("log")
    figfft.savefig('E_fft/' + str(i).zfill(int(np.ceil(np.log10(len(E))))) + '.png')

    curveEfft.remove()
    del curveEfft
    curveXfft.remove()
    del curveXfft

