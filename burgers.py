"""
Scrip to solve the 2D Burgers' equation using the finite difference method. Then
we plot our solutions for the time interval in a animated plot varying the value
of the viscous term changing its influence.

du/dt + u du/dt + v du/dt = nu * (d2u/dx2 + d2u/dy2)
dv/dt + u dv/dt + v dv/dt = nu * (d2v/dx2 + d2v/dy2)

where nu is the kinematic viscosity

I.C.: 2 <= x,y <= 1
      1 for everywhere else

B.C.: 1 for x = 0, 2 and y = 0, 2
"""

# Importing our favorite libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  # Colour map
import matplotlib.animation as animation


def burgers(L=2, M=2, T=0, Nx=40, Ny=40, Nt=2500, nu=0.5):
    """
    Function that solve the burgers equation and return the vectors x and y
    and the matrix of solution u.

    Parameters
    ----------
    L : float
        The maximum value of the X axis
    M : float
        The maximum value of the Y axis
    T : float
        The time interval
    Nx : int
        Number of steps for X axis
    Ny : int
        Number of steps for Y axis
    Nt : int
        Number of steps for time interval. The smaller the nu, the higher the
        number of steps for time, it has to be with the stability of the
        algorithm
    nu : float
        Value of the diffusion parameter. Kinetic viscosity if we are talking
        about momentum transport

    Returns
    -------
    list
        a list of float np.arrays (x, y, u) where x and y are the value of
        the 2D Cartesian space and u the matrix solution.
    """

    # Vectors
    x = np.linspace(0, L, Nx)
    y = np.linspace(0, M, Ny)

    # steps size
    dx = L / Nx
    dy = M / Ny
    dt = T / Nt

    # Initial conditions
    # In math notation the first index correspond to x direction and the
    # second one to the y direction. However, in arrays the first index is up
    # to down and the second one is left to right, so we have to invert them.
    u = np.ones([Ny, Nx])
    v = np.ones([Ny, Nx])
    u[int(0.5 / dx) : int(1 / dx) + 1, int(0.5 / dy) : int(1 / dy) + 1] = 2
    v[int(0.5 / dx) : int(1 / dx) + 1, int(0.5 / dy) : int(1 / dy) + 1] = 2

    for _ in range(Nt):  # The _ means we don't need the loop's index
        un = u.copy()
        vn = v.copy()

        # For u
        # Again, because in math notation the first index correspond to x
        # direction and the second one to the y direction but in arrays
        # the first index is up to down and the second one is left to right,
        # we must invert them.
        u[1:-1, 1:-1] = (
            un[1:-1, 1:-1]
            - un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, :-2])
            - vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[:-2, 1:-1])
            + nu * dt / dx ** 2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, :-2])
            + nu * dt / dy ** 2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[:-2, 1:-1])
        )

        # For v
        v[1:-1, 1:-1] = (
            vn[1:-1, 1:-1]
            - un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[:-2, 1:-1])
            - vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[1:-1, :-2])
            + nu * dt / dx ** 2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[:-2, 1:-1])
            + nu * dt / dy ** 2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, :-2])
        )

    return x, y, u


def update_data(t, nu):
    """[summary]

    Parameters
    ----------
    t : [type]
        [description]
    nu : [type]
        [description]
    """

    x, y, u = burgers(T=t, nu=nu)
    x, y = np.meshgrid(x, y)
    ax.collections = []
    ax.plot_surface(x, y, u, cmap=cm.plasma)


# Fist Frame. We need to create our fig object for the FuncAnimation()
x, y, u = burgers()
fig = plt.figure(figsize=(10, 5))
fig.suptitle("Burgers' Equation for $\\nu = 0$")
ax = fig.add_subplot(111, projection="3d")
x, y = np.meshgrid(x, y)
surf = ax.plot_surface(x, y, u, cmap=cm.plasma)
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2)
ax.view_init(60, 135)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$u$")

# Animation
t = np.linspace(0, 1, 100)
nu = 0
anim = animation.FuncAnimation(fig, update_data, t, fargs=(nu,), interval=100)
animWriter = animation.PillowWriter(fps=30) 
# Uncomment next line to save the animation
anim.save("animation.gif", animWriter)
print('Done')
plt.show()
