 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, arctan2, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
from blrtor import pol2cart, bolcorr_hardX, bolcorr_B
from blrtor import log_lum, log_bhm, compute_dust_sublimation_radius, compute_heating_radius, compute_L_UV_at_angle, compute_L_X_at_angle, compute_rmax, compute_critical_angle, compute_eddington_luminosity

from blrtor import get_peak, compute_blr_shape, compute_sublimation_radius, get_blr_covering_factor, compute_grav_radius, compute_jet_edge, get_nlr_size

from blrtor import compute_alphadisk_height, compute_sphere_of_influence

# from Baskin & Laor (2017)
rad_efficiency = 0.1
Z = 5 # metallicity

# Radiative Fountain system assumptions
dust_to_gas_ratio = 1/100.
dust_to_gas_ratio = 1/20.
kappa = 1e3 * u.cm**2/u.g
rho_g = 300 * u.Msun / u.pc**3
# ##################

ls = ':'
#for eddrate in [0.001, 0.03, 0.1, 2]:
for logeddrate in [-3, -1.5, -1, 0, 0.5]:
	eddrate = 10**logeddrate
	plt.figure(figsize=(10,7))
	logMBH = 8.0
	MBH = 10**logMBH * u.Msun
	L_AGN_edd = compute_eddington_luminosity(MBH)
	theta = numpy.linspace(0, pi/2, 4000)

	L_AGN = eddrate * L_AGN_edd
	L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
	L_AGN_edd = compute_eddington_luminosity(MBH)
	r_grav = compute_grav_radius(MBH).to(u.pc)

	M_dot = (L_AGN / c.c**2 / rad_efficiency).to(u.Msun / u.yr)

	R_in = compute_sublimation_radius(MBH, M_dot)
	R = numpy.logspace(-3, 3, 1000) * R_in

	plt.title("$M_\mathrm{BH}=%s$ $\lambda=%.3f$ $L_{46}=%.2f$ $\dot{M}=%.1f M_\odot/yr$" % (log_bhm(MBH), eddrate, L_AGN_46, M_dot.to(u.Msun / u.yr).value))
	
	# Accretion disk:	
	z_disk = compute_alphadisk_height(MBH, M_dot, R, alpha=1).to(u.pc)
	plt.plot(R, z_disk, '-', color='k', label="Accretion disk (SS)")
	
	idx = numpy.where(R>40*r_grav)[0][0]
	plt.text(R[idx].value, z_disk[idx].value, 'SS AD', 
		va='bottom', ha='right', color='k')
	
	# Jet:
	R_fullrange = numpy.logspace(log10(r_grav.value), log10(R.value.max()*5), 400) * u.pc
	R_fullrange = numpy.logspace(log10(r_grav.value), log10(1e8*r_grav.value), 400) * u.pc
	R_j = compute_jet_edge(MBH, R_fullrange)
	plt.plot(R_j, R_fullrange, '--', color='cyan', label='Radio jet')
	plt.text(R_j.value[-1], R_fullrange.value[-1], 'Radio jet', color='cyan', va='bottom', ha='center', bbox=dict(color='white'))
	
	# Corona:
	y, x = pol2cart(5*r_grav, theta)
	plt.plot(x, y, ls=':', color='orange', label='X-ray Corona')
	plt.text(r_grav.value, 5*r_grav.value, 'corona', va='bottom', ha='left')

	# BLR:
	z1, z2, zmax            = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='static')
	z1, z2, zmax_dyn        = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='dynamic', dyn_factor=2)
	CF = get_blr_covering_factor(z1, z2, R)
	Rpeak, Hpeak = get_peak(z2, z1, R)
	z1, z2, zmax_dyn_nodust = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='dustfree', dyn_factor=2)

	idx, = numpy.where(numpy.logical_or(zmax_dyn_nodust > z_disk, zmax_dyn > z_disk))
	lo, hi = idx.min() - 1, idx.max() + 1
	zmax_dyn_disk = numpy.where(zmax_dyn > z_disk, zmax_dyn, z_disk)
	zmax_dyn_nodust_disk = numpy.where(zmax_dyn_nodust > z_disk, zmax_dyn_nodust, z_disk)
	
	l, = plt.plot(R[lo:hi], zmax_dyn_nodust_disk[lo:hi], '--',
		label="Broad Line Region (BLCH, CF=%.0f%%)" % (CF*100)
	)
	color = l.get_color()
	plt.plot(R[lo:hi], zmax_dyn_disk[lo:hi], '-', color=color, label='Dusty BLR')
	#plt.plot(Rpeak, Hpeak, 'o', color=color)
	plt.text(Rpeak.value, Hpeak.value, 'BLCH BLR\n$\lambda=%.3f$\nCF=%.2f' % (eddrate, CF), 
		va='bottom', ha='center', color=color)
	#plt.fill_between(R.flatten(), 0, zmax_dyn, color=color, hatch='//')
	#plt.vlines(R_in.to(u.pc).value, 0, Hpeak.to(u.pc).value, linestyles=[':'], colors=[color])

	# Radiative fountain
	r_d = compute_dust_sublimation_radius(L_AGN)
	r_0 = compute_heating_radius(L_AGN, rho_g = 300 * u.Msun / u.pc**3)
	r_max = compute_rmax(MBH, L_AGN, theta, r_0, r_d, 
		dust_to_gas_ratio = dust_to_gas_ratio,
		kappa = kappa).to(u.pc)
	theta_crit = compute_critical_angle(MBH, L_AGN, r_0, r_d, 
		dust_to_gas_ratio, kappa = kappa)
	CF = cos(theta_crit)
	
	r_max[r_max > 100 * u.pc] = numpy.nan
	theta_min = arctan2(numpy.interp(r_d.value, R.value, z_disk.value), r_d.value)
	mask = theta < pi/2
	y, x = pol2cart(r_max, theta)
	color = 'darkred'
	plt.plot(x[mask], y[mask], ls='-.', color=color, label="Radiative Fountain TOR")

	mask = theta > 88./180*pi
	plt.text(x[mask][0].value, y[mask][0].value, 'RF TOR \nCF=%.2f ' % CF, 
		va='bottom', ha='right', color=color)
	
	
	# NLR
	R_NLR_hi = get_nlr_size().to(u.pc)
	# these are chosen arbitrarily: the inner edge of the NLR and the maximum angle
	# in practice they are observationally and host galaxy limited, respectively
	R_NLR_lo = 10 * u.pc
	thetai = numpy.linspace(pi/400, pi/4, 10)
	for i, (thetalo, thetahi) in enumerate(zip(thetai[:-1], thetai[1:])):
		(ylo1, ylo2), (xlo1, xlo2) = pol2cart(R_NLR_lo.value, numpy.array([thetalo, thetahi]))
		(yhi1, yhi2), (xhi1, xhi2) = pol2cart(R_NLR_hi.value, numpy.array([thetalo, thetahi]))
		label = 'NLR' if i == 0 else None
		alpha = 1. - (i*1./len(thetai))**0.8
		
		plt.fill([xlo1, xhi1, xhi2, xlo2], [ylo1, yhi1, yhi2, ylo2], 
			color='orange', 
			alpha=alpha, lw=0,
			label=label)
	
	# Sphere of influence:
	theta = numpy.linspace(0, pi/2, 400)
	# from Kormendy & Ho+13, compute sigma from BH mass
	sigma = 200 * u.km / u.s * 10**((log10((MBH / (1e9 * u.Msun)).to(1)) + 0.510) / 4.377)
	r_infl = compute_sphere_of_influence(MBH, sigma).to(u.pc)
	y, x = pol2cart(r_infl, theta)
	plt.plot(x, y, ls='--', lw=0.1, color='green', label='Sphere of influence')
	plt.text(r_infl.value, r_infl.value / 1000, 'Sphere of\ninfluence', va='bottom', ha='left')

	
	# sight-lines
	theta = numpy.array([1, 5, 15, 30, 60, 85])
	theta[theta == 0] = 5
	theta[theta == 90] = 85
	for thetai in theta:
		Rline = numpy.array([5 * r_grav.to(u.pc).value, 4])
		#Rline = numpy.linspace(5 * r_grav.to(u.pc).value, 10, 1000)
		y, x = pol2cart(Rline, thetai / 180 * pi)
		#x = R_fullrange.value
		#y = x * sin(thetai / 180 * pi)
		#mask = y > 1
		#y, x = y[mask], x[mask]
		plt.plot(x, y, ':', color='k', alpha=0.1)
		plt.text(x[-1], y[-1], '$%d^\circ$' % (thetai))

	plt.xlabel("R [pc]")
	plt.ylabel("z [pc]")
	plt.legend(loc="lower right", prop=dict(size=12))
	plt.ylim(r_grav.value / 100, 200*100)
	plt.xlim(r_grav.value, 50*100)
	plt.xscale('log')
	plt.yscale('log')
	prefix = "combinedstructure_2d_MBH%s_Mdot%s" % (logMBH, logeddrate)
	plt.savefig(prefix + '_log.pdf', bbox_inches="tight")
	plt.savefig(prefix + '_log.png', bbox_inches="tight")
	continue
	plt.ylim(0, 5)
	plt.xlim(0, 5)
	plt.xscale('linear')
	plt.yscale('linear')
	plt.savefig(prefix + ".pdf", bbox_inches="tight")
	plt.savefig(prefix + ".png", bbox_inches="tight")
	plt.close()



