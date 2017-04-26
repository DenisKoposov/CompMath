#!/usr/bin/python3

import matplotlib.pyplot as plt
import csv
import math

x1 = []
y1 = []
x2 = []
y2 = []
h = 0.5 / 1000
x3 = [ i*h +0.5 for i in range(0,1001) ]
y3 = [ math.log(2)*(h*i+0.5) + (h*i+0.5)*math.log(h*i+0.5) for i in range(0,1001) ]
dy_dx = []


with open("data_TDMA.csv", "r") as f1, open("data_SM.csv", "r") as f2:
    cont1 = csv.reader(f1)
    cont2 = csv.reader(f2)
    for e in cont1:
        x1.append(e[0])
        y1.append(e[1])
    for e in cont2:
        x2.append(e[0])
        y2.append(e[1])
        dy_dx.append(e[2])


plt.figure(1)
plt.subplot(121)
plt.grid(True)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Function")
plt.plot(x1, y1, 'r-')
plt.plot(x2, y2, 'b--', label="SM")
plt.plot(x3, y3, 'g--', label="AN")
plt.legend(["TDMA", "SM", "AN"])

plt.subplot(122)
plt.grid(True)
plt.xlabel("x")
plt.ylabel("y'")
plt.title("Derivative")
plt.plot(x2, dy_dx)
plt.legend(["SM"])

plt.show()
