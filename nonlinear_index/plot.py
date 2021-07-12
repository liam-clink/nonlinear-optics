# Plot results of nonlinear_index
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('envelope_field.tsv',dtype=np.complex128)

data_norm = np.abs(data)
data_phase = np.angle(data)
data_phase = np.unwrap(data_phase)

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1,projection='3d')
ax2 = fig.add_subplot(1,2,2,projection='3d')

thickness = 10.
duration = 10.
timestep = 0.01
spacestep = 0.01
t = np.arange(0.,duration,timestep)
z = np.arange(0.,thickness,spacestep)

# Indexing='ij' makes a lot more sense
Z, T = np.meshgrid(z,t,indexing='ij')

ax1.plot_wireframe(Z, T, data_norm)
ax1.set_ylabel('t')
ax1.set_xlabel('z')

ax2.plot_wireframe(Z, T, data_phase)
ax2.set_ylabel('t')
ax2.set_xlabel('z')

plt.show()