import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

k1 = 9                                              # Constants
D_base = 1e-5                                       # Diffusion coefficient for all compounds
                                                    # Grid and time step parameters
dtheta = np.pi / 180                                # Angular step size (1 degree in radians)
R = 15                                              # Radius of the semicircle (mm)
Ntheta = int(np.pi / dtheta)                        # Number of angular grid points: 180 points
dt = 1                                              # Time step size (1 second)
Nt = 50000                                          # Total number of time steps (50000 seconds)
                                                    # Initialize concentration arrays with the specified initial conditions
u = np.zeros(Ntheta + 2)                            # TDS
v = np.zeros(Ntheta + 2)                            # MCE
w = np.zeros(Ntheta + 2)                            # HEDT
x = np.zeros(Ntheta + 2)                            # DEAC
y = np.zeros(Ntheta + 2)                            # HEDS

D_theta = np.full(Ntheta + 2, D_base)               # Initialize all diffusion coefficient in to the matrix
k2 = np.zeros(Ntheta + 2)
k2[1] = 0.6
                                                    # Initial conditions
x[1] = 0.01                                         # DEAC concentration at one end
y[1] = 0.01                                         # HEDS concentration at one end
u[:] = 0.005                                        # TDS concentration in tube and the other end
v[:] = 0.01                                         # MCE concentration in tube and the other end
u[1] = 0                                            # TDS concentration at one end
v[1] = 0                                            # MCE concentration at one end
                                                    # Function to plot the DEAC concentration
def plot_deac_concentration(theta, concentration, time):
    fig, ax = plt.subplots(figsize=(10, 5))         # Output picture size
    ax.set_aspect('equal')

    for i in range(1, Ntheta + 1):
        color = 'blue' if concentration[i] > 4e-20 else 'yellow'
        arc = patches.Arc((0, 0), 2 * R, 2 * R, theta1=np.degrees(theta[i - 1]), theta2=np.degrees(theta[i]),color=color,linewidth=6)
        ax.add_patch(arc)

    ax.set_xlim(-R-.5, R+.5)                        # Decoration of picture output.
    ax.set_ylim(0, R+.5)
    plt.axis('off')
    plt.grid(False)
    plt.savefig("pic-3/temp" + str(time) + '.jpg')  # jpg files then convert to gif or video by format factory
    #plt.show()
                                                    # Time-stepping`             loop
for n in range(Nt):
    u_next = u.copy()                               # Array copied to temporarily save data
    v_next = v.copy()
    w_next = w.copy()
    x_next = x.copy()
    y_next = y.copy()
                                                    # Enforce Neumann boundary conditions(zero flux) at both ends
    u_next[0], u_next[Ntheta + 1] = u_next[1], u_next[Ntheta]
    v_next[0], v_next[Ntheta + 1] = v_next[1], v_next[Ntheta]
    w_next[0], w_next[Ntheta + 1] = w_next[1], w_next[Ntheta]
    x_next[0], x_next[Ntheta + 1] = x_next[1], x_next[Ntheta]
    y_next[0], y_next[Ntheta + 1] = y_next[1], y_next[Ntheta]

    for j in range(1, Ntheta + 1):
        if x[j] > 0.00015:
            k2[j] = 9.7

                                                    # Angular diffusion term
        u_theta_diff = D_theta[j] * (u[j + 1] - 2 * u[j] + u[j - 1]) / (R ** 2 * dtheta ** 2)
        v_theta_diff = D_theta[j] * (v[j + 1] - 2 * v[j] + v[j - 1]) / (R ** 2 * dtheta ** 2)
        w_theta_diff = D_theta[j] * (w[j + 1] - 2 * w[j] + w[j - 1]) / (R ** 2 * dtheta ** 2)
        x_theta_diff = D_theta[j] * (x[j + 1] - 2 * x[j] + x[j - 1]) / (R ** 2 * dtheta ** 2)
        y_theta_diff = D_theta[j] * (y[j + 1] - 2 * y[j] + y[j - 1]) / (R ** 2 * dtheta ** 2)
                                                    # Concentration of Last step + dt * (angular diffusion term + variation given from reaction
        u_next[j] = u[j] + dt * (u_theta_diff - k1 * u[j] * v[j])
        v_next[j] = v[j] + dt * (v_theta_diff - k1 * u[j] * v[j] - k2[j] * v[j] * w[j])
        w_next[j] = w[j] + dt * (w_theta_diff + k1 * u[j] * v[j] - k2[j] * v[j] * w[j])
        x_next[j] = x[j] + dt * (x_theta_diff+k2[j] * v[j] * w[j])
        y_next[j] = y[j] + dt * (y_theta_diff+k2[j] * v[j] * w[j])
                                                    # Update data
    u = u_next
    v = v_next
    w = w_next
    x = x_next
    y = y_next

                                                    # Plot the DEAC concentration at intervals of 5000 seconds
    if (n + 1) % 5000 == 0:
        print(x)
        plot_deac_concentration(np.linspace(0, np.pi, Ntheta + 2), x, (n + 1) * dt)
