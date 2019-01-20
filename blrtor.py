import numpy
from numpy import pi, exp, sin, cos, arccos, arctan2, tan, log, log10
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
	return log10((L / (u.erg / u.s)).to(1).value)
def log_bhm(MBH):
	return log10((MBH / (u.Msun)).to(1).value)

# u.* are units, c.* are physical constants

# from Wada (2016)

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
	theta_crit = numpy.arccos(cos_theta_crit)
	#mask = numpy.logical_and(cos_theta_crit > 0, cos_theta_crit < 1)
	#mask = numpy.isfinite(theta_crit)
	mask = numpy.logical_and(cos_theta_crit > 0, numpy.isfinite(theta_crit))
	theta_crit = numpy.where(mask, theta_crit, 0)
	return theta_crit

def compute_eddington_luminosity(MBH):
	L_AGN_edd = 4 * pi * c.G * c.c * c.m_p * MBH / c.sigma_T
	return L_AGN_edd.to(u.erg/u.s)


# from Baskin & Laor (2017)


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
	if numpy.any(z_sub > zmax):
		ipeak = numpy.where(z_sub > zmax)[0][0]
	else:
		ipeak = -1
	Rpeak = R.flatten()[ipeak]
	Hpeak = zmax[ipeak]
	return Rpeak, Hpeak

def compute_covering_fraction(Rpeak, Hpeak):
	CF = cos(pi/2 - arctan2(Hpeak.value, Rpeak.value))
	return CF

def compute_blr_size(R, zmax):
	"""Compute mean radius (weighted by height)"""
	Rmean = numpy.trapz(x=R, y=zmax * R) / numpy.trapz(x=R, y=zmax_dyn)
	return Rmean


def compute_minshape(z1, z2):
	return numpy.where(z1 < z2, z1, z2) * u.pc

def compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='static', dyn_factor=2):
	R = R.reshape((1,-1))
	z = R.reshape((-1,1)) / 2
	R_in = compute_sublimation_radius(MBH, M_dot)
	T = compute_disk_temperature(MBH, M_dot, R)
	z_sub = compute_sub_height(L_AGN, 1, R.flatten()).to(u.pc)
	z_sub = numpy.where(R.flatten() < R_in, 0, z_sub) * u.pc
	
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
	
	# dynamic solution is a bit larger, by a factor of 2
	z1 = dyn_factor * zmax
	zmax_dyn = numpy.where(R.flatten() > Rpeak, z1, z_sub) * u.pc
	if variant == 'dynamic':
		return z1, z_sub, compute_minshape(z1, z_sub)
	Rpeak, Hpeak = get_peak(z_sub, z1, R)
	# linear interpolation between disk at sublimation radius and peak
	# as a coarse approximation to the numerical solution
	z2 = Hpeak * ((R.flatten()) - R_in) / (Rpeak - R_in)
	z2[z2 < 0] = 0
	if variant == 'dustfree':
		return z1, z2, compute_minshape(z1, z2)
	
	assert False, variant

def get_blr_covering_factor(z1, z2, R):
	Rpeak, Hpeak = get_peak(z2, z1, R)
	CF = compute_covering_fraction(Rpeak, Hpeak)
	return CF


def compute_alphadisk_height(MBH, M_dot, R, alpha=1):
	m = (MBH / u.Msun).to(1)
	mdot = (M_dot / (3e-8 * u.Msun / u.yr) * u.Msun / MBH).to(1)
	R_g = 2 * c.G * MBH/ c.c**2
	r = (R / (3*R_g)).to(1)
	#za = 3/(8*pi) * c.sigma_T / c.c * M_dot * (1 - r**-0.5)
	za = 3.2e6 * mdot * m * (1 - r**-0.5)
	#Ha = 1e8 * m**-0.5 * r**(-3/4.)
	maskab = r / (1 - r**-0.5)**(16/21.) <= 150 * (alpha*m)**(2/21.) * mdot**(16/21.)
	zb = 1.2e4 * alpha**(-1/10.) * mdot**(1/5.) * m**(9/10.) * r / (1 - r**-0.5)**(1/5.)
	#Hb = 1.5e9 * alpha**(1/20.) * mdot**(2/5.) * m**(-9/20.) * r**(-51/40.) * (1 - r**-0.5)**(2/5.)
	zab = (za**2 + zb**2)**0.5
	
	zc = 6.1e3 * alpha**(-1/10.) * mdot**(3/20.) * m**(9/10.) * r**(9/8.) / (1 - r**-0.5)**(3/20.)
	
	maskbc = r / (1 - r**0.5)**(2/3.) < 6.3e3 * mdot**(2/3.)
	
	# pad a bit
	idx = numpy.where(maskab)[0]
	idx = idx[numpy.logical_and(idx>0, idx<len(maskab))]
	maskab[idx-1] = True
	maskab[idx+1] = True
	idx = numpy.where(maskbc)[0]
	idx = idx[numpy.logical_and(idx>0, idx<len(maskab))]
	maskbc[idx-1] = True
	maskbc[idx+1] = True
	
	# combined in quadrature for smoothness:
	zabc = (za**2 + zb**2 + zc**2)**0.5
	# very far in (near 3Rg), the shape diverges, so we do not use it
	zabc[R < 4*R_g] = numpy.nan
	
	return (zabc * u.cm).to(R.unit)

def compute_alphadisk_temperature(MBH, M_dot, R, alpha=1):
	m = (MBH / u.Msun).to(1)
	mdot = (M_dot / (3e-8 * u.Msun / u.yr) * u.Msun / MBH).to(1)
	R_g = 2 * c.G * MBH / c.c**2
	r = (R / (3*R_g)).to(1)
	T40 = (3 * c.G * MBH * M_dot / (8 * pi * c.sigma_sb * R**3)).to(u.K**4)
	T4 = T40 * (1 - r**-0.5)
	return T4**(1/4.)

def compute_jet_edge(MBH, R):
	# from https://ui.adsabs.harvard.edu/#abs/arXiv:1611.04075
	# empirical intermediate between
	#      quasi-conical streamline (Blandford & Znajek 1977),
	#      and force-free steady jet solution (Narayan et al 2007; Tchekhovskoy et al 2008)
	
	R_g = 2 * c.G * MBH / c.c**2
	R_j = (1 + 1.3*(R/R_g)**0.7) * R_g
	return R_j.to(R.unit)

def compute_radiocore_luminosity(MBH, L_AGN):
	"""
	5GHz luminosity of radio core
	"""
	L_X = bolcorr_hardX(L_AGN)
	m = log10(MBH / u.Msun)
	# Merloni, Heinz & Di Matteo (2003)
	logLR = 0.6 * log10(L_X/(u.erg/u.s)) + 0.78 * m + 7.33
	return 10**logLR * u.erg/u.s

def get_nlr_size():
	"""
	Husemann+03 uses IFU data and finds that sizes have been underestimated previously (when using slits; Bennert+02, Greene+11)
	Size of NLR is fairly lum-independent 4.1pc, but there is 0.2dex scatter.
	The Subaru narrow filter imaging by Sun+18 confirmed this as well, finding a very weak luminosity dependence from 43.5 to 46.5, with 0.5dex scatter.
	"""
	return 4.1 * u.kpc


def compute_sphere_of_influence(MBH, sigma):
	R_infl = c.G * MBH / sigma**2
	return R_infl


def compute_outflow_rate(LAGN, Mstar, SFR):
	totLum = 1.29 * SFR / (u.Msun/u.yr) + 0.81 * LAGN / (1e43 * u.erg/u.s)
	logFlow = 1.13 * log10(totLum) - 0.37 * log10(Mstar/(1e11 * u.Msun))
	return 10**logFlow * u.Msun / u.yr



