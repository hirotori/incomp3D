import numpy as np
import matplotlib.pyplot as plt
import os

#load csv
y, u = np.loadtxt("ghia_u.csv", delimiter=",", unpack=True, skiprows=1, usecols=(0,1))
x, v = np.loadtxt("ghia_v.csv", delimiter=",", unpack=True, skiprows=1, usecols=(0,1))

#load txt
nstep = 3000
fname_u = os.path.join("cavity_result", "cavity_u{0}.txt".format(nstep))
fname_v = os.path.join("cavity_result", "cavity_v{0}.txt".format(nstep))
y2, u2 = np.loadtxt(fname_u, unpack=True, skiprows=1)
x2, v2 = np.loadtxt(fname_v, unpack=True, skiprows=1)

#plot
fig,axes = plt.subplots()
axes.set_xlabel(r"y")
axes.set_ylabel(r"u")
axes.scatter(y, u, label="Ghia (1982)")
axes.plot(y2, u2, label="present")
axes.legend()

plt.savefig("vertical_velocity")

fig,axes = plt.subplots()
axes.set_xlabel(r"x")
axes.set_ylabel(r"v")
axes.scatter(x, v, label="Ghia (1982)")
axes.plot(x2, v2, label="present")
axes.legend()

plt.savefig("horizontal_velocity")