import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parameters
b = 7     # plant birth rate
d = 5     # plant death rate
s = 1     # plant contribution to the reservoir
r = 3     # rescue rate

# System of differential equations
def f(Z, t):
    z1, z2 = Z
    dz1_dt = -(b * z1**2 - (d + b + s) * z1 + s * z2 + d)
    dz2_dt = r * (z2 - z1)
    return [dz1_dt, dz2_dt]

# Define the grid for the quiver plot
z1 = np.linspace(0, 2.0, 30)
z2 = np.linspace(0.0, 2.0, 30)
Z1, Z2 = np.meshgrid(z1, z2)

# Time point (not used in this static phase diagram)
t = 0

# Initialize velocity field
u, v = np.zeros(Z1.shape), np.zeros(Z2.shape)

# Compute the velocity field
NI, NJ = Z1.shape
for i in range(NI):
    for j in range(NJ):
        x = Z1[i, j]
        y = Z2[i, j]
        zprime = f([x, y], t)
        u[i, j] = zprime[0]
        v[i, j] = zprime[1]

# Create the quiver plot and stream plot
plt.figure(figsize=(10, 8))

# Plot the quiver plot
Q = plt.quiver(Z1, Z2, u, v, color='r',alpha = 0.3)

# Plot the stream plot
plt.streamplot(Z1, Z2, u, v, color='g', density=1.5)
plt.scatter ( 1,1 ,color='blue')
plt.scatter ( d/b,d/b ,color='blue')
# Add labels and title
plt.xlabel('$z_1$')
plt.ylabel('$z_2$')

# # Optionally, plot trajectories from different initial conditions
# initial_conditions = [
#     [0.0, -1.0],
#     [0.0, 0.0],
#     [1.0, -1.0],
#     [-1.5, -1.5],
#     [0.2, 3],
#     [-0.5, 3.5],
#     [-4, -4]
# ]

# for y0 in initial_conditions:
#     tspan = np.linspace(0, 50, 2000)
#     ys = odeint(f, y0, tspan)
#     plt.plot(ys[:, 0], ys[:, 1], 'b-')  # path
#     plt.plot([ys[0, 0]], [ys[0, 1]], 'bo')  # start
#     plt.plot([ys[-1, 0]], [ys[-1, 1]], 'bs')  # end

# Show plot
plt.show()
