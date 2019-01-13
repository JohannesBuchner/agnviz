 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, tan, log, log10
import astropy.constants as c
import astropy.units as u
import numpy as np

def cart2pol(x, y):
	rho = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y, x)
	return(rho, phi)

def pol2cart(rho, phi):
	x = rho * np.cos(phi)
	y = rho * np.sin(phi)
	return(x, y)

# Bolometric correction functions:
def bolcorr(Lbol, c1, k1, c2, k2):
	L10 = Lbol.to(u.L_sun) / (1e10 * u.L_sun)
	return c1 * L10**k1 + c2 * L10**k2

# Hopkins+07:
def bolcorr_hardX(Lbol):
	return Lbol / bolcorr(Lbol, 10.83, 0.28, 6.08, -0.020)

def bolcorr_softX(Lbol):
	return Lbol / bolcorr(Lbol, 17.87, 0.28, 10.03, -0.020)

def bolcorr_B_(Lbol):
	return Lbol / bolcorr(Lbol, 7.40, -0.37, 10.66, -0.014)

# Marconi+04:
def bolcorr_B(Lbol):
	L12 = log10(Lbol.to(u.L_sun) / (1e12 * u.L_sun))
	return Lbol / 10**(0.8 + 0.067*L12 + 0.017*L12**2 - 0.0023*L12**3)


# u.* are units, c.* are physical constants

L_AGN = 0.02e44 * u.erg / u.s
MBH = 0.1e8 * u.Msun

kappa = 1e3 * u.cm**2/u.g
rho_g = 300 * u.Msun / u.pc**3
rho_g0 = 1 * u.Msun / u.pc**3
dust_to_gas_ratio = 0.01

L_X = bolcorr_hardX(L_AGN)
L_UV = bolcorr_B(L_AGN)

print("MBH:", MBH)
print("Luminosities: ")
print("  L_AGN: ", L_AGN)
print("  L_X: ", L_X)
print("  L_UV: ", L_UV)

# option 1
Lambda_cool_rhog2 = 1e-22 * u.erg / u.cm**3 / u.s * (rho_g**2 / rho_g0**2)
# option 2, product is constant
Lambda_cool_rhog2 = 1e-22 * u.erg / u.cm**3 / u.s

r_0 = (3 * L_X / (4 * pi * Lambda_cool_rhog2))**(1/3.)
r_d = 1.3 * u.pc * (L_AGN / (1e46 * u.erg / u.s))**0.5
gamma_d = dust_to_gas_ratio

print("Radii:")
print("   r0:", r_0.to(u.pc))
print("   rd:", r_d.to(u.pc))

theta = numpy.linspace(0, 2*pi, 400)
L_AGN_angle = 1/2. * L_UV * numpy.abs(cos(theta))
L_AGN_dust = kappa * gamma_d * L_AGN / (8*pi*c.c)
print("  L_AGN_angle: ", L_AGN_angle)
print("  L_AGN_dust: ", L_AGN_dust)
r_max = (L_AGN_dust + c.G * MBH) / (c.G * MBH / r_0 - L_AGN_dust * (1/r_d - 2/r_0))
r_max = r_max.to(u.pc)
r_max[r_max < 0] = 0
print(r_max)

import matplotlib.pyplot as plt

plt.figure(figsize=(4,6))
#plt.plot(theta, r_max.to(u.pc))
x, y = pol2cart(r_max, theta)
plt.plot(x, y, 'o-')
plt.savefig("radfountain-angle.pdf", bbox_inches="tight")
plt.close()


