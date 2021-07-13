# Calculate dispersion for nonlinear index of refraction in 1D

import numpy as np
import matplotlib.pyplot as plt

### Define Physical Constants ###

# Fused Silica https://doi.org/10.1364/OE.27.011018
# Quartz and fused silica: https://doi.org/10.1103/PhysRevB.39.3337
# Quoted for 1030nm
n0 = 1.45
n2 = 2.2e-16 # cm^2/W
n2 /= 100.**2 # m^2/W
chi3 = n2*n0**2/(12*np.pi**2) # Same units as n2, Formula from Boyd Nonlinear Optics 3ed p213
c = 299792458. # m/s
mu_0 = 4e-7*np.pi # H/m
epsilon_0 = 1./(c**2*mu_0) # F/m

# Returns intensity given the amplitude of the E field wave
def I(E):
    return .5*n0*epsilon_0*c*E**2
# Returns E wave amplitude from Intensity
def E(I):
    return np.sqrt(I/(.5*n0*epsilon_0*c))

### Solving the nonlinear Schroedinger Equation ###

# Simulation parameters
thickness = 0.5e-3 # m
duration = 2.5e-12 # s
times = np.linspace(0.,duration,1000)
spaces = np.linspace(0.,thickness,1000)
timestep = times[1]-times[0]
spacestep = spaces[1]-spaces[0]
pulse_duration = 260.e-15 # s

# Self phase modulation and dispersion coefficients
wavelength = 1030.e-9 # m
angular_frequency = 2*np.pi*c/wavelength
wavenumber = 2.*np.pi/wavelength
spm_strength = 3./2.*chi3/(wavenumber*c**2)*angular_frequency**2
print(spm_strength)
print(wavenumber)
dispersion_strength = 0.

# Initialize data container
A = np.zeros((len(spaces),len(times)),dtype=np.complex128)

# Left boundary condition
# A_initial = 1.e0 * np.sin(omega * times) * (np.exp(-((times - 0.2) / pulse_duration) ** 2))
# A is only the envelope
A_initial = 1.e0 * (np.exp(-((times - duration/2.) / pulse_duration) ** 2))

A[0,:] = A_initial


### Run the simulation ###

# Right hand side of the equation
def right_hand_side(A):
    # TODO: Currently assuming no dispersion
    return -1j*spm_strength*A*A*np.conjugate(A)

def RK4(A,dz,f):
    k1 = dz*f(A)
    k2 = dz*f(A+k1/2.)
    k3 = dz*f(A+k2/2.)
    k4 = dz*f(A+k3)
    return A + k1/6. + k2/3. + k3/3. + k4/6.

for i in range(1,A.shape[0]):
    A[i,:] = RK4(A[i-1,:],spacestep,right_hand_side)

# Save results
np.savetxt("envelope_field.tsv", A, delimiter='\t')

