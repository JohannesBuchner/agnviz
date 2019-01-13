 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, arctan2, tan, log, log10
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
	#assert numpy.all(500 * u.K < T), T
	#assert numpy.all(T <= 10000 * u.K), T
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
	z_sub = z_sub0 / L_AGN_46 * (R/u.pc)**3 * u.pc
	#print(z_sub0, L_AGN_46, R, z_sub)
	return z_sub

def compute_sublimation_radius(MBH, M_dot, Tsub = 2000 * u.K):
	R_in = 1. / (Tsub**4 * c.sigma_sb / (3 / (8*pi)) / (c.G * MBH * M_dot))**(1/3.)
	return R_in.to(u.pc)

def compute_disk_temperature(MBH, M_dot, R):
	# eq 21:
	F = 3 / (8*pi) * (c.G * MBH * M_dot / R**3)
	# eq 5:
	T = (F / c.sigma_sb)**(1./4)
	T = T.to(u.K)
	T[T > 10000 * u.K] = 10000 * u.K
	return T

def compute_radiation_force(MBH, M_dot, R, kappa):
	a_rad = 3 / (8*pi) * (c.G * MBH * M_dot / R**3) * kappa / c.c
	return a_rad

def compute_BHgravity(MBH, z, R):
	a_BH = c.G * MBH * (z / (R**2 + z**2)**(3./2))
	return a_BH 

def compute_grav_radius(MBH):
	return 2 * c.G * MBH/ c.c**2

def get_peak(z_sub, zmax, R):
	ipeak = numpy.where(z_sub > zmax)[0][0]
	Rpeak = R.flatten()[ipeak]
	Hpeak = zmax[ipeak]
	return Rpeak, Hpeak

def compute_covering_fraction(Rpeak, Hpeak):
	CF = cos(pi/2 - arctan2(Hpeak.value, Rpeak.value))
	return CF

def compute_minshape(z1, z2):
	return numpy.where(z1 < z2, z1, z2) * u.pc

def compute_blr_shape(MBH, L_AGN, R, Z, variant='static', dyn_factor=2):
	M_dot = (L_AGN / c.c**2 / rad_efficiency).to(u.Msun / u.yr)
	R = R.reshape((1,-1))
	z = R.reshape((-1,1)) / 2
	R_in = compute_sublimation_radius(MBH, M_dot)
	T = compute_disk_temperature(MBH, M_dot, R)
	z_sub = compute_sub_height(L_AGN, 1, R.flatten()).to(u.pc)
	
	a_rad = compute_radiation_force(MBH, M_dot, R, kappa = compute_kappa(T, Z))
	a_BH  = compute_BHgravity(MBH, z, R)
	mask = a_rad < a_BH
	Rmax = numpy.where(mask, R + 0*z, numpy.inf).min(axis=1) * u.pc
	Rmax = numpy.where(Rmax < R_in, R_in, Rmax) * u.pc

	zmax = numpy.where(mask, z + 0*R, numpy.inf).min(axis=0) * u.pc
	zmax2 = compute_minshape(zmax, z_sub)
	if variant == 'static':
		return zmax, z_sub, zmax2
	
	Rpeak, Hpeak = get_peak(z_sub, zmax2, R)
	#CF_static = compute_covering_fraction(Rpeak, Hpeak)
	# dynamic solution is a bit larger, by a factor of 2
	z1 = dyn_factor * zmax
	zmax_dyn = numpy.where(R.flatten() > Rpeak, z1, z_sub) * u.pc
	#CF_dyn = compute_covering_fraction(Rpeak, Hpeak)
	if variant == 'dynamic':
		return z1, z_sub, compute_minshape(z1, z_sub)
	Rpeak, Hpeak = get_peak(z_sub, z1, R)
	
	# linear interpolation between disk at sublimation radius and peak
	z2 = Hpeak * ((R.flatten()) - R_in) / (Rpeak - R_in)
	z2[z2 < 0] = 0
	if variant == 'dustfree':
		return z1, z2, compute_minshape(z1, z2)
	
	#zmax_dyn_dustfree = numpy.where(R.flatten() > Rpeak, z1, z2) * u.pc
	#Rpeak, Hpeak = get_peak(z2, z1, R)
	#CF_dyn_dustfree = compute_covering_fraction(Rpeak, Hpeak)
	
	assert False, variant

def get_blr_covering_factor(z1, z2, R):
	Rpeak, Hpeak = get_peak(z2, z1, R)
	CF = compute_covering_fraction(Rpeak, Hpeak)
	return CF

plot_all_lines = False

plt.figure(figsize=(15,3))
logMBH = 8
MBH = 10**logMBH * u.Msun
L_AGN_edd = compute_eddington_luminosity(MBH)
R = numpy.linspace(0, 0.1, 4000) * u.pc
#z = numpy.linspace(0, 0.05, 4001).reshape((-1,1)) * u.pc
#Rgrav = compute_grav_radius(MBH)

for eddrate, ls in [(0.023, ':'), (0.1, '--'), (0.45, '-')]:
#for eddrate, ls in [(0.45, '-')]:
	L_AGN = eddrate * L_AGN_edd
	L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
	M_dot = (L_AGN / c.c**2 / rad_efficiency).to(u.Msun / u.yr)
	R_in = compute_sublimation_radius(MBH, M_dot)
	
	"""
	R_in = compute_sublimation_radius(MBH, M_dot)
	T = compute_disk_temperature(MBH, M_dot, R)
	z_sub = compute_sub_height(L_AGN, 1, R.flatten()).to(u.pc)
	
	a_rad = compute_radiation_force(MBH, M_dot, R, kappa = compute_kappa(T, Z))
	a_BH  = compute_BHgravity(MBH, z, R)
	mask = a_rad < a_BH
	Rmax = numpy.where(mask, R + 0*z, numpy.inf).min(axis=1) * u.pc
	Rmax = numpy.where(Rmax < R_in, R_in, Rmax) * u.pc

	zmax = numpy.where(mask, z + 0*R, numpy.inf).min(axis=0) * u.pc
	zmax2 = numpy.where(zmax < z_sub, zmax, z_sub)
	
	Rpeak, Hpeak = get_peak(z_sub > zmax2)
	# dynamic solution is a bit larger, by a factor of 2
	z1 = 2 * zmax
	Rpeak, Hpeak = get_peak(z_sub > z1)
	zmax_dyn = numpy.where(R.flatten() > Rpeak, z1, z_sub) * u.pc
	z2 = Hpeak * ((R.flatten()) - R_in) / (Rpeak - R_in)
	z2[z2 < 0] = 0
	zmax_dyn_dusty = numpy.where(R.flatten() > Rpeak, z1, z2) * u.pc
	
	CF = cos(pi/2-arctan2(Hpeak.value, Rpeak.value))
	"""
	
	z1, z2, zmax            = compute_blr_shape(MBH, L_AGN, R, Z, variant='static')
	z1, z2, zmax_dyn        = compute_blr_shape(MBH, L_AGN, R, Z, variant='dynamic', dyn_factor=2)
	CF = get_blr_covering_factor(z1, z2, R)
	Rpeak, Hpeak = get_peak(z2, z1, R)
	z1, z2, zmax_dyn_nodust = compute_blr_shape(MBH, L_AGN, R, Z, variant='dustfree', dyn_factor=2)
	
	l, = plt.plot(R, zmax_dyn, 
		label="$M_\mathrm{BH}=%d$ $\lambda=%.3f$ $L_{46}=%.2f$ $\dot{M}=%.1f M_\odot/yr$ CF=%.0f%%" % (log_bhm(MBH), eddrate, L_AGN_46, M_dot.to(u.Msun / u.yr).value, CF*100)
	)
	color = l.get_color()
	plt.vlines(R_in.to(u.pc).value, 0, 0.004, linestyles=[':'], colors=[color])
	plt.plot(R.flatten(), zmax_dyn_nodust, '--', color=color)
	plt.plot(Rpeak, Hpeak, 'o', color=color)
	plt.fill_between(R.flatten(), 0, zmax_dyn, color=color, hatch='//')

plt.xlim(0, 0.1)
plt.ylim(0, 0.02)
plt.xlabel("R [pc]")
plt.ylabel("z [pc]")
plt.legend(loc="best", prop=dict(size=8))
plt.savefig("windyblr_2d.pdf", bbox_inches="tight")
plt.close()








