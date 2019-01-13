 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from blrtor import pol2cart, bolcorr_hardX, bolcorr_B
from blrtor import log_lum, log_bhm, compute_dust_sublimation_radius, compute_heating_radius, compute_L_UV_at_angle, compute_L_X_at_angle, compute_rmax, compute_critical_angle, compute_eddington_luminosity

# from Baskin & Laor (2017)
rad_efficiency = 0.1
Z = 5 # metallicity


def compute_blr_height(M_BH, M_dot, Z, R):
	M8 = M_BH / (1e8 * u.Msun)
	Mdot1 = M_dot / (u.Msun / u.yr)
	R_pc = R / u.pc
	
	# equation (47)
	H = 1.17e-4 * u.pc * M8**0.29 * Mdot1**1.29 * R_pc**-0.87 * Z
	return H

def compute_kappa(T, Z, amin=0.005 * u.um, amax=1 * u.um):
	assert 500 * u.K < T < 2500 * u.K, T
	if amin == 0.005 * u.um and amax == 1 * u.um:
		return 54 * (T / (2000*u.K))**1.16 * Z * u.cm**2 / u.g
	else:
		assert False
		return 43 * (T / (2000*u.K))**1.94 * Z * u.cm**2 / u.g


def compute_sub_height(L_AGN, a, R):
	L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
	if a == 1:
		z_sub0 = 21.7
	elif a == 0.1:
		z_sub0 = 5.7
	elif a == 0.01:
		z_sub0 = 0.6
	z_sub = z_sub0 / L_AGN_46 * (R/u.pc)**3
	return z_sub

plt.figure(figsize=(5,3))
#for L_AGN in [1e44 * u.erg/u.s, 1e45 * u.erg/u.s, 1e46 * u.erg/u.s]:
#for logMBH, color in [(7.5, '#1f77b4'), (8.5, '#ff7f0e'), (9.5, '#2ca02c')]:
for logMBH, color in [(8, '#1f77b4')]:
	MBH = 10**logMBH * u.Msun
	L_AGN_edd = compute_eddington_luminosity(MBH)
	#for eddrate, ls in [(0.003, ':'), (0.03, '--'), (0.3, '-.'), (3, '-')]:
	for eddrate, ls in [(0.023, ':'), (0.1, '--'), (0.45, '-')]:
		L_AGN = eddrate * L_AGN_edd
		L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
		L_UV = bolcorr_B(L_AGN)
		L_UV_45 = (L_UV / (1e45 * u.erg/u.s)).to(1)
		
		# just before equation (10), Richards+06, probably a bit crude
		L_opt_45 = L_AGN_46 
		M8 = MBH / (1e8 * u.Msun)
		# L = rad_efficiency * M * c**2
		M_dot = (L_AGN / c.c**2 / rad_efficiency).to(u.Msun / u.yr)
		#M_dot = 0.8 * L_AGN_46 / M8
		
		# equation (9)
		R_in  = 0.018 * u.pc * (L_opt_45)**0.5
		
		x = numpy.linspace(0, 0.05, 1000) * u.pc
		z = compute_blr_height(MBH, M_dot, Z, x)
		#HR = y / x
		#HRmax = 2**-0.5
		#y = numpy.where(HR > HRmax, 0, y)
		
		print("R_in:", R_in, M_dot)
		
		plt.plot(x, z, 
			ls=ls, color=color,
			label="$M_\mathrm{BH}=%d$ $\lambda=%.3f$ $L_{46}=%.2f$ $\dot{M}=%.1f M_\odot/yr$" % (log_bhm(MBH), eddrate, L_AGN_46, M_dot.to(u.Msun / u.yr).value))
		plt.vlines(R_in.value, 0, 0.004, linestyles=[ls], colors=['k'])
		plt.plot(x, compute_sub_height(L_AGN, 1, x), '--', color='k', lw=1)

plt.xlim(0, 0.05)
plt.ylim(0, 0.03)
plt.xlabel("R [pc]")
plt.ylabel("z [pc]")
plt.legend(loc="best", prop=dict(size=8))
plt.savefig("windyblr.pdf", bbox_inches="tight")
plt.close()


