 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, arctan2, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt

N = 10000
# generate randomly oriented rod
v = numpy.random.normal(size=(N, 3))
norm = (v**2).sum(axis=1)**0.5
v = v / norm.reshape((-1, 1))
# compute projected length
l1 = (v[:,0]**2 + v[:,1]**2)**0.5

r = 4.1
h = 1

# randomly oriented system
theta = arccos(numpy.random.uniform(size=N))
# compute length of projection
rproj = sin(theta)
# should be the same as l1
plt.hist(rproj, bins=40, histtype='step')
plt.hist(l1, bins=40, histtype='step')
plt.savefig("simcones_orientation.pdf", bbox_inches="tight")
plt.close()

Acirc = pi * r * cos(theta)
Acone = h * r * sin(theta)
Atot = Acirc + Acone

histargs = dict(bins=40, histtype='step', normed=True)
plt.hist(log10(Acirc), label="circle", **histargs)
plt.hist(log10(Acone), label="cone", **histargs)
plt.hist(log10(Atot), lw=2, label="both", **histargs)
plt.legend()
plt.savefig("simcones.pdf", bbox_inches="tight")
plt.close()


