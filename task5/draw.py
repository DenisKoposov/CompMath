#!/usr/bin/python2

import numpy as np

# Drawing precise solution: u(t,x)=(x+4)**2 * exp(-2t)
# c * tau / h <= 1, max(c) = 5, tau / h <= 1/5, tau = 0.0001, h = 0.001

h = 0.01
x = np.arange(0, 1. + h/10, h)

tau = 0.01
t = np.arange(0, 1. + tau/10, tau)
u = np.ndarray([len(t), len(x)])

for i, cur_t in enumerate(t):
    u[i, ] = (x + 4) ** 2 * np.exp(-2 * cur_t)

# Drawing numerical solution
u_num = np.zeros_like(u)
# Setting initial condition
u_num[0, ] = (x + 4) ** 2
# Setting initial  condition
u_num[:, 0] = 16 * np.exp(-2 * t)
# Calculating all the rest in the cycle
for i, xc in enumerate(x[1:], 1):
    for j, tc in enumerate(t[:-1]):
        u_num[j+1, i] = ( u_num[j, i] -
                          tau / h * (xc + 4) * (u_num[j, i] - u_num[j, i-1]) )

np.savetxt('u.csv', u, delimiter=',', fmt='%.4f')
np.savetxt('u_num.csv', u_num, delimiter=',', fmt='%.4f')

print 'loss:', np.sum(u - u_num) ** 2 / (len(x) * len(t))

#animation on surface
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
def start_surface_animation():
    fig, ax = plt.subplots()
    graph, = ax.plot(x, u_num[0,:])

    def animate(i):
        graph.set_ydata(u_num[i,:])
        return graph,

    anim = animation.FuncAnimation(fig, animate, np.arange(len(t)), interval=1)
    plt.axis([-0.1, 1.1, 0, 25])
    plt.show()

#color 3d graph
def dplot(x, t, function1, function2):
    fig = plt.figure()
    X, T = np.meshgrid(x, t)
    ax = fig.gca(projection='3d')
    ax.set_title('3d graph')
    ax.set_ylabel('time')
    ax.set_xlabel('x')
    ax.set_zlabel('u(t,x)')
    surf1 = ax.plot_surface(X, T, function1, cmap='inferno')#, cmap=cm.coolwarm, linewidth=0)
    surf2 = ax.plot_surface(X, T, function2, cmap='winter')#, cmap=cm.coolwarm, linewidth=0)
    #fig.colorbar(surf1, shrink=0.5, aspect=5)
    plt.show()

dplot(x, t, u, u_num)
#start_surface_animation()
