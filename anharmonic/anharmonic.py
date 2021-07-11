from time import sleep

import matplotlib.pyplot as plt
import numpy as np

L = 5.
t = 5.
c = 1.
pulse_length = 0.1
omega = 100.
alpha = 0.  # 1. #Ne/epsilon_0 (how strongly field responds to oscillator)
beta = 300.  # e/m (how strongly oscillator responds to field)
gamma = 5.
omega_0 = 110.
a = 1.
b = 1.e4
CFL = 0.9

dz = 2 * np.pi * c / omega / 30
N = int(np.ceil(L / dz))
dt = min(CFL * dz / c, (2 * np.pi / omega_0) / 20)
Nt = int(np.ceil(t / dt))
print("CFL: ", c * dt / dz, "dt: ", dt)

z = np.linspace(0., L, N)
x = [np.zeros(N, dtype=np.float64)]
u = [np.zeros(N, dtype=np.float64)]  # starts at u[i+1/2]
E = [np.zeros(N, dtype=np.float64)]
v = [np.zeros(N, dtype=np.float64)]  # starts at v[i+1/2]

sigma = np.zeros(N)
damping_range = 0.25
damping_number = int(np.ceil(damping_range / L * N))
for i in range(damping_number):
    sigma[i - damping_number] = 1500 * (i * dz) ** 2

times = np.linspace(0., t, Nt)
E_drive = 5.e0 * np.sin(omega * times) * (np.exp(-((times - 0.2) / pulse_length) ** 2))
v[0][0] = E_drive[0].real / dt

for i in range(Nt):
    print(i)
    E.append(np.zeros(N))
    E[i + 1] = E[i] * (1 - dt * sigma) + dt * v[i]
    E[i + 1][0] = E_drive[i].real

    append = x.append(x[i] + dt * u[i])
    u_dot = -beta * E[i + 1] - gamma * u[i] - omega_0 ** 2 * x[i + 1] + a * x[i + 1] ** 2 + b * x[i + 1] ** 3
    u.append(u[i] + dt * u_dot)

    v.append(np.zeros(N))
    v[i + 1][1:N - 2] = v[i][1:N - 2] + dt * (alpha * u_dot[1:N - 2] + (c / dz) ** 2 * (
           E[i + 1][0:N - 3] - 2 * E[i + 1][1:N - 2] + E[i + 1][2:N - 1]))
    v[i + 1][-1] = E[i + 1][-1] - 2.5 * E[i + 1][-2] + 2 * E[i + 1][-3] - 0.5 * E[i + 1][-4]

np.save('E', E)
np.save('x', x)
np.save('u', u)
np.save('v', v)
np.save('z', z)
np.save('t', times)
