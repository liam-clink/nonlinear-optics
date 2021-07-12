# Calculate dispersion for nonlinear index of refraction in 1D

import numpy as np
import matplotlib.pyplot as plt

### Define Physical Constants ###

# Fused Silica https://doi.org/10.1364/OE.27.011018
# Quartz and fused silica: https://doi.org/10.1103/PhysRevB.39.3337
# Quoted for 1030nm
n0 = 1.45
n2 = 2.2e-16 # cm^2/W
chi3 = n2*n0**2/(12*np.pi**2) # Formula from Boyd Nonlinear Optics 3ed p213
c = 299792458. # m/s
mu_0 = 4e-7*np.pi # H/m
epsilon_0 = 1./(c**2*mu_0) # F/m

# Returns the polarization field in angular frequency space
def P(E):
    # The 3 is because w,w,-w can be distributed in 3 distinct ways
    # The conjugate is because E_time is real, so E_freq(-w)=E_freq(w)*
    return 3*epsilon_0*chi3*E*E*np.conjugate(E)

# Returns intensity given the amplitude of the E field wave
def I(E):
    return .5*n0*epsilon_0*c*E**2


### Solving the nonlinear Schroedinger Equation ###

# Simulation parameters
width = 10.
duration = 10.
speed = 1.
omega = 10*(2.*np.pi)
timestep = omega/20.
spacestep = 0.9*speed/timestep
times = np.arange(0.,duration,timestep)
pulse_duration = 2.

# Initialize matrix
N = int(np.ceil(width/spacestep)) # Number of interior z nodes, total z nodes add two more boundary nodes
matrix = np.zeros((N,N),dtype=np.float64)
for i in range(N-1):
    matrix[i+1,i] = -1.
    matrix[i,i+1] = 1.


# Initialize data container
A = np.zeros((len(times),N+2))

# Left boundary condition
#A_drive = 1.e0 * np.sin(omega * times) * (np.exp(-((times - 0.2) / pulse_duration) ** 2))
# A is only the envelope
A_drive = 1.e0 * (np.exp(-((times - 0.2) / pulse_duration) ** 2))

A[:,0] = A_drive

# Right boundary condition
# The right boundary is nonreflecting, which is achieved by scaling
# down the wave with increasing strength deeper into the damping zone
extinction_coefficient = .1/spacestep
damping_strength = np.zeros(N)
scaling = np.zeros(N)
damping_range = 0.25*width
damping_number = int(np.ceil(damping_range / spacestep))
for i in range(damping_number):
    damping_strength[i - damping_number] = extinction_coefficient*float(i)/float(damping_number)
for i in range(len(scaling)):
    scaling[i] = np.exp(-damping_strength[i]*timestep)


### Run the simulation ###
max_convergence_steps = 1000
max_convergence_ratio = 1.e-2

for i in range(1,A.shape[0]):
    print("shut up errors")
    # Update middle diagonal of matrix
    
    # Iterate matrix equation until satisfactory convergence

    # Damp


# Save results
np.savetxt("envelope_field.tsv", delimiter='\t')
