# Importing all the required modules
import numpy as np
import matplotlib.pyplot as plt
import argparse

###
# Defining all the functions
###
def solver(I, w, dt, T):
    """
    Solve u'' + (w**2)*u = 0 for t in (0,T], u(0) = I and u'(0) = 0
    using central finite difference method and time step dt
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    u[0] = I
    u[1] = u[0] - 0.5*(dt**2)*(w**2)*u[0]
    for n in range(1,Nt):
        u[n+1] = 2*u[n] - u[n-1] - (dt**2)*(w**2)*u[n]
    return u, t

def u_exact(t, I, w):
    """
    Returning the exact solution of u as
    u = I*cos(w*t)
    """
    return I*np.cos(w*t)

def visualize(u, t, I, w):
    """
    Visualizing the plot
    """
    plt.plot(t, u, 'r--o')
    t_fine = np.linspace(0, t[-1], 1001)
    u_e = u_exact(t_fine, I, w)

    plt.plot(t_fine, u_e, 'b-')
    plt.legend(['numerical','exact'],loc='upper left')
    plt.xlabel('t')
    plt.ylabel('u')
    dt = t[1] - t[0]
    plt.title('dt=%g' % dt)
    umin = 1.2*u.min(); umax = -umin
    plt.axis([t[0], t[-1], umin, umax])
    plt.savefig('temp1.png')
    plt.savefig('temp1.pdf')

# Allowing the parameters to be adjusted using cmdline
parser = argparse.ArgumentParser()
parser.add_argument('--I', type=float, default=1.0)
parser.add_argument('--w', type=float, default=2*np.pi)
parser.add_argument('--dt', type=float, default=0.05)
parser.add_argument('--num_periods', type=int, default=5)
a = parser.parse_args()
# Assigning all the values to the required variables
I, w, dt, num_periods = a.I, a.w, a.dt, a.num_periods
# Definition of Period
P = 2*np.pi/w

T = P*num_periods
u, t = solver(I, w, dt, T)
visualize(u, t, I, w)
