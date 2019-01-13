import numpy
from numpy import pi, exp, sin, cos, arccos, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
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

def log_lum(L):
	return log10(L / (u.erg / u.s))
def log_bhm(MBH):
	return log10(MBH / (u.Msun))

# u.* are units, c.* are physical constants

def compute_dust_sublimation_radius(L_AGN):
	"""
	L_AGN: bolometric luminosity
	"""
	r_d = 1.3 * u.pc * (L_AGN / (1e46 * u.erg / u.s))**0.5
	return r_d.to(u.pc)

def compute_heating_radius(L_AGN, rho_g):
	"""
	L_AGN: bolometric luminosity
	rho_g: local gas density
	"""
	L_X = bolcorr_hardX(L_AGN)
	n_H = rho_g / c.m_p
	Lambda_cool_rhog2 = 1e-22 * u.erg * u.cm**3 / u.s * n_H**2
	r_0 = (3 * L_X / (4 * pi * Lambda_cool_rhog2))**(1/3.)
	return r_0

def compute_L_UV_at_angle(L_UV, theta):
	""" L_AGN: luminosity, theta in rad """
	L_UV_angle = 1/2. * L_UV * numpy.abs(cos(theta))
	return L_UV_angle
def compute_L_X_at_angle(L_X, theta):
	""" L_AGN: luminosity, theta in rad """
	return L_X

def compute_rmax(MBH, L_AGN, theta, r_0, r_d, 
	dust_to_gas_ratio,
	kappa = 1e3 * u.cm**2/u.g):

	gamma_d = dust_to_gas_ratio
	
	#L_AGN_angle = 1/2. * L_UV * numpy.abs(cos(theta))
	L_UV = bolcorr_B(L_AGN)
	L_UV_angle = compute_L_UV_at_angle(L_UV, theta)
	L_AGN_dust = kappa * gamma_d * L_UV_angle / (8*pi*c.c)
	
	r_max = (L_AGN_dust + c.G * MBH) / (c.G * MBH / r_0 - L_AGN_dust * (1/r_d - 2/r_0))
	r_max = r_max.to(u.pc)
	r_max[r_max < 0] = numpy.nan
	return r_max

def compute_critical_angle(MBH, L_AGN, r_0, r_d, 
	dust_to_gas_ratio,
	kappa = 1e3 * u.cm**2/u.g):
	"""
	Returns degrees
	"""
	gamma_d = dust_to_gas_ratio
	L_UV = bolcorr_B(L_AGN)
	cos_theta_crit = (c.G * MBH / r_0 * 16 * pi * c.c / (kappa * gamma_d * L_UV) * (1 / r_d - 2 / r_0)**-1).to(1).value
	theta_crit = numpy.arccos(cos_theta_crit) / pi * 180
	#mask = numpy.logical_and(cos_theta_crit > 0, cos_theta_crit < 1)
	#mask = numpy.isfinite(theta_crit)
	mask = numpy.logical_and(cos_theta_crit > 0, numpy.isfinite(theta_crit))
	theta_crit = numpy.where(mask, theta_crit, 0)
	return theta_crit

def compute_eddington_luminosity(MBH):
	L_AGN_edd = 4 * pi * c.G * c.c * c.m_p * MBH / c.sigma_T
	return L_AGN_edd.to(u.erg/u.s)


