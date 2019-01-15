 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from blrtor import pol2cart, bolcorr_hardX, bolcorr_B
from blrtor import log_lum, log_bhm, compute_dust_sublimation_radius, compute_heating_radius, compute_L_UV_at_angle, compute_L_X_at_angle, compute_rmax, compute_critical_angle, compute_eddington_luminosity

# system assumptions
dust_to_gas_ratio = 1/100.
dust_to_gas_ratio = 1/20.
kappa = 1e3 * u.cm**2/u.g
rho_g = 300 * u.Msun / u.pc**3
# ##################

plt.figure(figsize=(5,5))
theta = numpy.linspace(0, pi/2, 400)
#for L_AGN in [1e44 * u.erg/u.s, 1e45 * u.erg/u.s, 1e46 * u.erg/u.s]:
for logMBH, color in [(7.5, '#1f77b4'), (8.5, '#ff7f0e'), (9.5, '#2ca02c')]:
	MBH = 10**logMBH * u.Msun
	L_AGN_edd = compute_eddington_luminosity(MBH)
	for eddrate, ls in [(0.003, ':'), (0.03, '--'), (0.3, '-.'), (3, '-')]:
		L_AGN = eddrate * L_AGN_edd
		r_d = compute_dust_sublimation_radius(L_AGN)
		r_0 = compute_heating_radius(L_AGN, rho_g = 300 * u.Msun / u.pc**3)
		r_max = compute_rmax(MBH, L_AGN, theta, r_0, r_d, 
			dust_to_gas_ratio = dust_to_gas_ratio,
			kappa = kappa)
		y, x = pol2cart(r_max, theta)
		plt.plot(x, y, 
			ls=ls, color=color,
			label="$M_\mathrm{BH}=%d$ $\lambda=%s$" % (log_bhm(MBH), eddrate))

r_max = max(max(plt.ylim()), max(plt.xlim()))
r_max = 20
plt.ylim(0.01, r_max)
plt.xlim(0.01, r_max)
plt.xlabel("x [pc]")
plt.ylabel("y [pc]")
plt.legend(loc="best", prop=dict(size=8))
#plt.xscale("log")
#plt.yscale("log")
plt.savefig("radfountain_2d-angle.pdf", bbox_inches="tight")
plt.close()



MBH_lo, MBH_hi = 6.5, 9.5
MBH1d = numpy.logspace(MBH_lo, MBH_hi, 100) * u.Msun
L_AGN_lo, L_AGN_hi = 42, 47
L_AGN1d = numpy.logspace(L_AGN_lo, L_AGN_hi, 100) * u.erg / u.s


L_AGN, MBH = numpy.meshgrid(L_AGN1d, MBH1d)
r_0 = compute_heating_radius(L_AGN, rho_g = 300 * u.Msun / u.pc**3)
r_d = compute_dust_sublimation_radius(L_AGN)

theta_crit = compute_critical_angle(MBH, L_AGN, r_0, r_d, 
	dust_to_gas_ratio = dust_to_gas_ratio,
	kappa = kappa) / pi * 180

plt.imshow(theta_crit,
	origin="lower", aspect="auto",
	extent=[L_AGN_lo, L_AGN_hi, MBH_lo, MBH_hi], 
	cmap="Oranges", #cmap="gray_r"
)
cb = plt.colorbar()
plt.xlim(L_AGN_lo, L_AGN_hi)
plt.ylim(MBH_lo, MBH_hi)

# compute eddington rate:
L_AGN_edd = compute_eddington_luminosity(MBH1d)
plt.plot(log_lum(L_AGN_edd), log_bhm(MBH1d), "--", color="k")
plt.plot(log_lum(L_AGN_edd*0.1), log_bhm(MBH1d), "--", color="k", alpha=0.5)
plt.plot(log_lum(L_AGN_edd*0.01), log_bhm(MBH1d), "--", color="k", alpha=0.25)

plt.xlabel("Bolometric Luminosity $L_\mathrm{AGN}$  [erg/s]")
plt.ylabel("Black Hole Mass $M_\mathrm{BH}$ [$M_\odot$]")
cb.set_label(r'$\theta_\mathrm{crit}$')
plt.savefig("radfountain_2d_LAGN.pdf", bbox_inches="tight")
plt.close()






L_X_lo, L_X_hi = 41, 45
L_X1d = bolcorr_hardX(L_AGN1d)
L_AGN1d = 10**numpy.interp(numpy.linspace(L_X_lo, L_X_hi, 100), log_lum(L_X1d), log_lum(L_AGN1d)) * u.erg/u.s

L_AGN, MBH = numpy.meshgrid(L_AGN1d, MBH1d)

L_X = bolcorr_hardX(L_AGN)
L_UV = bolcorr_B(L_AGN)

r_0 = compute_heating_radius(L_AGN, rho_g = rho_g)
r_d = compute_dust_sublimation_radius(L_AGN)

theta_crit = compute_critical_angle(MBH, L_AGN, r_0, r_d, 
	dust_to_gas_ratio = dust_to_gas_ratio,
	kappa = kappa) / pi * 180

import matplotlib.pyplot as plt
plt.imshow(theta_crit,
	origin="lower", aspect="auto",
	extent=[L_X_lo, L_X_hi, MBH_lo, MBH_hi], 
	cmap="Oranges", #cmap="gray_r"
)
cb = plt.colorbar()
plt.xlim(L_X_lo, L_X_hi)
plt.ylim(MBH_lo, MBH_hi)

# compute eddington rate:
L_AGN_edd = compute_eddington_luminosity(MBH1d)
plt.plot(log_lum(L_AGN_edd), log_bhm(MBH1d), "--", color="k")
plt.plot(log_lum(L_AGN_edd*0.1), log_bhm(MBH1d), "--", color="k", alpha=0.5)
plt.plot(log_lum(L_AGN_edd*0.01), log_bhm(MBH1d), "--", color="k", alpha=0.25)

plt.xlabel("2-10keV X-ray Luminosity $L_\mathrm{X}$  [erg/s]")
plt.ylabel("Black Hole Mass $M_\mathrm{BH}$ [$M_\odot$]")
cb.set_label(r'$\theta_\mathrm{crit}$')
plt.savefig("radfountain_2d_LX.pdf", bbox_inches="tight")
plt.close()




