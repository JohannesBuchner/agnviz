import numpy
from numpy import pi, exp, sin, cos, arccos, arctan2, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

def compute_height_alphadisk(MBH, M_dot, a, R):
	"""
	a: unitless spin
	"""
	Ledd = 4 * pi * c.G * c.c * c.m_p * MBH / c.sigma_T
	m = (MBH / u.Msun).to(1)
	mdot = (M_dot * c.c**2 / Ledd).to(1)
	
	r_star = R * c.c**2 / (c.G * MBH)
	y = (R / MBH)**0.5
	a_united = c.G * MBH**2 / c.c
	a_star = a_united / MBH
	print(y.si, a_star.si)
	r_G = c.G * MBH/ c.c**2
	Z1 = 1 + (1 - a_star**2)**(1/3.) * ((1 + a_star)**1/3. + (1 - a_star)**(1/3.))
	Z2 = (3 * a_star**2 + Z1**2)**0.5
	
	# ISCO
	r_ms = r_G * (3 + Z2 - ((3 - Z1) * (3 + Z1 + 2 * Z2))**0.5)
	y0 = (r_ms/MBH)**0.5
	y1 =  2 * cos((arccos(a_star) - pi)/3)
	y2 =  2 * cos((arccos(a_star) + pi)/3)
	y3 = -2 * cos((arccos(a_star))/3)
	
	A = 1 + a_star**2 * y**-4
	B = 1 + a_star * y**-3
	C = 1 - 3*y**-2 + 2*a_star * y**-3
	D = 1 - 2*y**-2 + a_star**2 * y**-4
	E = 1 + 4*a_star**2 * y**-4 - 4*a_star**2 * y**-6 + 3*a_star**4*y**-8
	Q0 = (1+a_star*y**-3) / y / (1 - 3*y**-2 + 2*a_star*y**-3)**0.5
	Q = 	Q0 * (y - y0 - 3 / 2. * a_star * log(y/y0) - 3*(y1-a_star)**2 / (y1*(y1 - y2)*(y1-y3)) * log((y - y1) / (y0 - y1))) - \
		Q0 * (3*(y2-a_star)**2 / (y2*(y2 - y1)*(y2-y3)) * log((y - y2) / (y0 - y2)) - 3*(y3-a_star)**2 / (y3*(y3 - y1)*(y3-y2)) * log((y - y3) / (y0 - y3)))
	
	# Outer region
	HO = 4e2*u.cm * alpha**(-1/10.) * m**(18/20.) * mdot**(3/20.) * r_star**(9/8.) * A**(19/20.) * B**(-11/10.) * C**(1/2.) * D**(-23/40) * E**(-19/40.) * Q**(3/20.)
	HM = 1e3*u.cm * alpha**(-1/10.) * m**(9/10.) * mdot**(1/5.) * r_star**(21/20.) * A * B**(-6/5.) * C**(1/2.) * D**(-3/5) * E**(-1/2.) * Q**(1/5.)
	HI = 1e5*u.cm * mdot * A**2 * B**-3 * C**(1/2.) * D**-1 * E**-1 * Q
	return R

def compute_height_alphadisk(MBH, M_dot, a, alpha, R):
	"""
	a: unitless spin
	"""
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

	zabc = (za**2 + zb**2 + zc**2)**0.5
	
	maskbc = r / (1 - r**0.5)**(2/3.) < 6.3e3 * mdot**(2/3.)
	
	return (zabc * u.cm).to(R.unit)
	

if __name__ == "__main__":
	MBH = 1e8 * u.Msun
	M_dot = 1 * u.Msun / u.yr
	a = 0
	Rgrav = 2 * c.G * MBH/ c.c**2
	R = numpy.logspace(0, 4.5, 4000) * Rgrav
	H = compute_height_alphadisk(MBH, M_dot, a, 1, R)
	plt.plot(R/(3*Rgrav), H/(3*Rgrav))
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig("alphadisk.pdf", bbox_inches='tight')
	plt.close()
	
