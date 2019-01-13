 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, tan, log, log10
import astropy.constants as c
import astropy.units as u


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

MBH = 1e8 * u.Msun
#L_AGN = numpy.array([1e41, 1e48]) * u.erg / u.s
L_AGN = numpy.logspace(42.5, 48, 40) * u.erg / u.s

kappa = 1e3 * u.cm**2/u.g
rho_g = 300 * u.Msun / u.pc**3
dust_to_gas_ratio = 0.01
dust_to_gas_ratio = 1./1000.
n_H = rho_g / c.m_p

L_X = bolcorr_hardX(L_AGN)
L_UV = bolcorr_B(L_AGN)

print("MBH:", MBH)
print("Luminosities: ")
print("  L_AGN: ", L_AGN)
print("  L_X: ", L_X)
print("  L_UV: ", L_UV)

Lambda_cool_rhog2 = 1e-22 * u.erg * u.cm**3 / u.s * n_H**2

r_0 = (3 * L_X / (4 * pi * Lambda_cool_rhog2))**(1/3.)
r_d = 1.3 * u.pc * (L_AGN / (1e46 * u.erg / u.s))**0.5
gamma_d = dust_to_gas_ratio

print("Radii:")
print("   r0:", r_0.to(u.pc))
print("   rd:", r_d.to(u.pc))

cos_theta_crit = (c.G * MBH / r_0 * 16 * pi * c.c / (kappa * gamma_d * L_UV) * (1 / r_d - 2 / r_0)**-1).to(1)
cos_theta_crit = numpy.where(cos_theta_crit > 1, 1, cos_theta_crit)

import matplotlib.pyplot as plt

#plt.plot(L_AGN, cos_theta_crit)
plt.figure(figsize=(4,6))
plt.subplot(3,1,1)
plt.plot(L_X, r_0.to(u.pc), label="$r_0$")
plt.plot(L_X, r_d.to(u.pc), label="$r_d$")
plt.ylabel("Radius [%s]" % u.pc)
plt.xlabel("L(AGN)")
plt.yscale("log")
plt.xscale("log")
plt.legend(loc="best")
plt.subplot(3,1,2)
plt.plot(L_X, 1/(((1 / r_d - 2 / r_0)**-1)/r_0).to(1), label="inverse radii")
plt.plot(L_X, (c.G * MBH * 16 * pi * c.c / (kappa * gamma_d * L_UV)).to(1), label="lumgrav")
plt.legend(loc="best")
plt.ylabel(r"$\cos(\theta)$ terms")
plt.xlabel("L(AGN)")
plt.yscale("log")
plt.xscale("log")

plt.subplot(3,1,3)
#plt.plot(L_AGN, cos_theta_crit)
plt.plot(L_X, numpy.arccos(cos_theta_crit)/pi*180)
plt.ylim(0, 90)
plt.ylabel(r"$\theta_{crit}$")
#plt.yscale("log")
plt.xlabel("L(AGN)")
plt.xscale("log")
plt.savefig("radfountain.pdf", bbox_inches="tight")
plt.close()


