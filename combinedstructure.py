 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, arctan2, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
import os
from blrtor import pol2cart, bolcorr_hardX, bolcorr_B
from blrtor import log_lum, log_bhm, compute_dust_sublimation_radius, compute_heating_radius, compute_L_UV_at_angle, compute_L_X_at_angle, compute_rmax, compute_critical_angle, compute_eddington_luminosity

from blrtor import get_peak, compute_blr_shape, compute_sublimation_radius, get_blr_covering_factor, compute_grav_radius, compute_jet_edge, compute_radiocore_luminosity, get_nlr_size

from blrtor import compute_alphadisk_height, compute_alphadisk_temperature, compute_sphere_of_influence, compute_outflow_rate

def rround(value):
	"""
	reasonable rounded value
	"""
	sigfigs = -int(numpy.floor(log10(value)))
	if sigfigs >= 0:
		fmt = "%%.%df" % (sigfigs)
	else:
		fmt = "%d"
		value = round(value, sigfigs)
	#print(value, fmt, sigfigs, '-->', fmt % value)
	return fmt % value


def get_color(T):
	# approximating functions to
	# Blackbody color datafile D58 (bbr_color_D58.html)
	# Mitchell Charity 
	# http://www.vendian.org/mncharity/dir3/blackbody/
	# Version 2016-Mar-30
	a, b, c = 1817.33, 0.922924, 483.016
	blue = numpy.where(T > 2000, c * log10(T/a)**b, 0)
	blue[blue < 0] = 0
	blue[~(blue < 255)] = 255
	a, b, c = 754.824, 0.731073, 280.19
	a2, b2, c2 = 4715.01, -0.163986, 173.712
	green = numpy.where(T < 5800, c * log10(T/a)**b, c2 * log10(T/a2)**b2)
	green[green < 0] = 0
	green[green > 255] = 255
	a, b, c = 4702.1, -0.266549, 137.742
	red = c * log10(T/a)**b
	red[red < 0] = 0
	red[~(red < 255)] = 255
	return list(zip(red/255, green/255, blue/255))

color_data_T = numpy.arange(1000,40000,100) * u.K
color_data_rgb = get_color(color_data_T.value)

def find_color_chunks(color_data_T, T_disk, color_data_rgb, R, z_disk):
	mask = T_disk < color_data_T[0]
	if mask.any():
		yield color_data_rgb[0], R[mask], z_disk[mask]
	for color, Ti in zip(color_data_rgb, color_data_T):
		mask = numpy.logical_and(T_disk > Ti - 100 * u.K, T_disk <= Ti)
		if mask.any():
			#print(Ti, R[mask], z_disk[mask])
			yield color, R[mask], z_disk[mask]
		else:
			#print(Ti, 'none', T_disk.min(), T_disk.max())
			continue
	mask = T_disk > color_data_T[-1]
	if mask.any():
		yield color_data_rgb[-1], R[mask], z_disk[mask]



def plot_log_agn_postcard(MBH, eddrate, 
	# from Baskin & Laor (2017)
	rad_efficiency = 0.1,
	Z = 5, # metallicity
	# Radiative Fountain system assumptions
	dust_to_gas_ratio = 1/20.,
	kappa = 1e3 * u.cm**2/u.g,
	rho_g = 300 * u.Msun / u.pc**3,
	show_BH = True,
	show_corona = True,
	show_disk = True,
	show_BLR = True,
	show_NLR = True,
	show_jet = True,
	show_TOR = True,
	show_SOI = True,
	show_viewing = True,
	colored_disk = True,
	show_flows = True,
):
	L_AGN_edd = compute_eddington_luminosity(MBH)

	L_AGN = eddrate * L_AGN_edd
	L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
	L_AGN_edd = compute_eddington_luminosity(MBH)
	r_grav = compute_grav_radius(MBH).to(u.pc)

	M_dot = (L_AGN / c.c**2 / rad_efficiency).to(u.Msun / u.yr)
	
	R_in = compute_sublimation_radius(MBH, M_dot)
	# from Kormendy & Ho+13, compute sigma from BH mass
	sigma = 200 * u.km / u.s * 10**((log10((MBH / (1e9 * u.Msun)).to(1)) + 0.510) / 4.377)
	r_infl = compute_sphere_of_influence(MBH, sigma).to(u.pc)

	R = numpy.logspace(log10(r_grav/R_in), 3, 1000) * R_in
	R = numpy.logspace(log10(r_grav/R_in), max(3, log10(r_infl/R_in)), 1000) * R_in

	title = "$M_\mathrm{BH}=%s$ $\lambda=%.3f$ $L_\mathrm{AGN}=%.1f$ $\dot{M}=%.1f M_\odot/yr$" % (log_bhm(MBH), eddrate, log_lum(L_AGN), M_dot.to(u.Msun / u.yr).value)
	#plt.title(title)
	
	colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
	textsize = 10
	
	xlo, xhi = 8e-8 * u.pc, 50*100 * u.pc
	ylo, yhi = xlo / 100,  200*100 * u.pc

	# Sphere of influence:
	theta = numpy.linspace(0, pi/2, 400)
	y, x = pol2cart(r_infl, theta)
	if show_SOI:
		plt.plot(x, y, ls='--', lw=0.5, color=colors[8], label='Sphere of influence')
		#plt.text(x.max().value*1.1, r_infl.value / 1000, 'Sphere of influence', 
		#	va='bottom', ha='left', rotation=90)
		plt.text(x.max().value * 1.2, ylo.value * 2, 'Sphere of influence', 
			va='bottom', ha='left', rotation=90)
		#plt.fill_between(x.value, y.value, y.value*0 + yhi.value, 
		#	color=colors[8], alpha=0.1)
		#plt.fill_between([x.max().value, xhi.value], [ylo.value]*2, [yhi.value]*2, 
		#	color=colors[8], alpha=0.1)
		#plt.fill_between(x.value, y.value, y.value*0 + ylo.value, 
		#	color=colors[8], alpha=0.06)

	# Accretion disk:
	R_disk = R
	z_disk = compute_alphadisk_height(MBH, M_dot, R_disk, alpha=1).to(u.pc)
	T_disk = compute_alphadisk_temperature(MBH, M_dot, R_disk, alpha=1).to(u.K)
	if show_disk:
		mask = R < r_infl
		if colored_disk:
			for color, R_part, z_disk_part in find_color_chunks(color_data_T, T_disk[mask], color_data_rgb, R_disk[mask], z_disk[mask]):
				plt.fill_between(R_part, z_disk_part, z_disk_part * 0 + ylo, color=color)
		else:
			plt.fill_between(R_disk[mask], z_disk[mask], z_disk[mask] * 0 + ylo, color='k', label="Accretion disk (SS)")
		
		#	'Accretion disk\n$\log \dot{M}=%.1f$' % rround(log10(M_dot.value)),
		plt.text(7 * r_grav.value, ylo.value * 2, 
			'Accretion disk\n$\log L_\mathrm{bol}=%.1f$' % log_lum(L_AGN),
			va='bottom', ha='left', color='k', size=textsize, bbox=dict(color='white'))
		#idx = numpy.where(R>20*r_grav)[0][0]
		#plt.text(R_disk[idx].value, z_disk[idx].value, 'Accretion disk', 
		#	va='top', ha='left', color='k', size=textsize)
		
		if show_flows:
			outflow = compute_outflow_rate(L_AGN, Mstar = 1e11 * u.Msun, SFR = 0*u.Msun/u.yr)
			inflow = M_dot # at accretion disk at least, or when considering steady-state
			# y = z_disk[mask][-1].value/10
			plt.text(r_infl.value*2, 1e-3, '$\\uparrow$ Outflow\n$%s M_\odot/\mathrm{yr} \cdot M_{\star,11}^{-0.4}$\n$\leftarrow$ BH Inflow\n$%s M_\odot/\mathrm{yr}$' % (rround(outflow.to(u.Msun/u.yr).value), rround(inflow.to(u.Msun/u.yr).value)), 
				va='top', ha='left', color='k', size=textsize)
		
	# Jet:
	R_fullrange = numpy.logspace(log10(r_grav.value), log10(R.value.max()*5), 400) * u.pc
	R_fullrange = numpy.logspace(log10(r_grav.value), log10(max(yhi, 1e8*r_grav).value), 400) * u.pc
	R_j = compute_jet_edge(MBH, R_fullrange)
	if show_jet:
		jetcolor = colors[4]
		plt.plot(R_j, R_fullrange, '--', color=jetcolor, label='Jet')
		plt.fill_betweenx(R_fullrange, R_j, color=jetcolor, label='Jet', alpha=0.05)
		plt.text(R_j.value[-1] / 5, R_fullrange.value[-1] / 2, 'Jet', color=jetcolor, 
			va='top', ha='right', size=textsize)

	# BH:
	theta = numpy.linspace(0, pi/2, 4000)
	y, x = pol2cart(r_grav, theta[::-1])
	if show_BH:
		plt.fill_between(x, 1e-10 + 0*y.value, y, color='k', label='Black hole')
		#plt.plot(x, y / 10, '-', color='k', label='Black hole 1')
		plt.text(r_grav.value * 0.8, r_grav.value / 3, 
			'Black Hole\n\n$\log M=%.1f$' % (log_bhm(MBH)) if MBH > 10**7.5 * u.Msun else 'BH\n\n$\log M=%.1f$' % (log_bhm(MBH)), 
			va='top', ha='right', color='white', fontweight='bold', size=textsize)

	# Corona:
	y, x = pol2cart(6*r_grav, theta)
	if show_corona:
		plt.plot(x, y, ls=':', color=colors[1], label='X-ray Corona')
		plt.text(xlo.value*1.2, 6*r_grav.value, 'Corona',
			va='bottom', ha='left', color=colors[1])

	# BLR:
	z1, z2, zmax            = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='static')
	z1, z2, zmax_dyn        = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='dynamic', dyn_factor=2)
	CF = get_blr_covering_factor(z1, z2, R)
	Rpeak, Hpeak = get_peak(z2, z1, R)
	z1, z2, zmax_dyn_nodust = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='dustfree', dyn_factor=2)

	idx, = numpy.where(numpy.logical_or(zmax_dyn_nodust > z_disk, zmax_dyn > z_disk))
	if idx.any():
		lo, hi = idx.min() - 1, idx.max() + 1
		zmax_dyn_disk = numpy.where(zmax_dyn > z_disk, zmax_dyn, z_disk)
		zmax_dyn_nodust_disk = numpy.where(zmax_dyn_nodust > z_disk, zmax_dyn_nodust, z_disk)
	else:
		lo, hi = 0, -1
		zmax_dyn_disk = zmax_dyn
		zmax_dyn_nodust_disk = zmax_dyn_nodust
	
	if show_BLR:
		color = colors[0]
		l, = plt.plot(R[lo:hi], zmax_dyn_nodust_disk[lo:hi], '--',
			label="Broad Line Region (BLCH, CF=%.0f%%)" % (CF*100),
			color=color,
		)
		
		plt.plot(R[lo:hi], zmax_dyn_disk[lo:hi], '-', color=color, label='Dusty BLR')
		plt.fill_between(R[lo:hi], z_disk[lo:hi], zmax_dyn_disk[lo:hi], color=color)
		#plt.plot(Rpeak, Hpeak, 'o', color=color)
		plt.text(Rpeak.value, Hpeak.value, 'BLCH BLR\nCF=%.2f' % (CF),
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
	
	r_max[r_max > 1000 * u.pc] = numpy.nan
	theta_min = arctan2(numpy.interp(r_d.value, R.value, z_disk.value), r_d.value)
	y, x = pol2cart(r_max, theta)
	color = colors[3]
	
	if show_TOR:
		plt.plot(x, y, ls='-.', color=color, label="Radiative Fountain TOR")

	mask = theta > 89.7/180*pi
	if show_TOR and mask.any() and numpy.isfinite(x[mask][0].value):
		plt.text(x[mask][0].value, y[mask][0].value, 'RF TOR \nCF=%.2f ' % CF, 
			va='bottom', ha='right', color=color)
	
	# NLR
	R_NLR_hi = get_nlr_size().to(u.pc)
	# these are chosen arbitrarily: the inner edge of the NLR and the maximum angle
	# in practice they are observationally and host galaxy limited, respectively
	R_NLR_lo = r_infl
	R_NLR = numpy.array([R_NLR_lo.value, R_NLR_hi.value]) * u.pc
	thetai = numpy.linspace(pi/400, pi/4, 10)
	for i, (thetalo, thetahi) in enumerate(zip(thetai[:-1], thetai[1:])):
		if not show_NLR: break
		ylo1, ylo2 = R_NLR.value
		xlo1, xlo2 = compute_jet_edge(MBH, R_NLR).value
		(yhi1, yhi2), (xhi1, xhi2) = pol2cart(R_NLR.value, thetahi)
		
		#(ylo1, ylo2), (xlo1, xlo2) = pol2cart(R_NLR_lo.value, numpy.array([thetalo, thetahi]))
		#(yhi1, yhi2), (xhi1, xhi2) = pol2cart(R_NLR_hi.value, numpy.array([thetalo, thetahi]))
		label = 'NLR' if i == 0 else None
		alpha = (1. - (i*1./len(thetai))**0.8)*0.5
		
		plt.fill([xlo1, xhi1, xhi2, xlo2], [ylo1, yhi1, yhi2, ylo2], 
			color=colors[2], 
			alpha=alpha, lw=0,
			label=label)
		plt.text(xlo2, ylo2*0.9, 'NLR',
			color='white', va='top', ha='left', size=textsize)
	
	
	# sight-lines
	theta = numpy.array([1, 5, 15, 30, 60, 85])
	theta[theta == 0] = 5
	theta[theta == 90] = 85
	for thetai in theta:
		if not show_viewing: break
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
	plt.xticks()
	#plt.legend(loc="lower right", prop=dict(size=12))
	#plt.ylim(r_grav.value / 100, 200*100)
	#plt.xlim(r_grav.value / 2, 50*100)
	plt.ylim(8e-8 / 100, 200*100)
	plt.xlim(8e-8, 50*100)
	plt.xscale('log')
	plt.yscale('log')

if __name__ == "__main__":
	fL = open("luminosities.js", 'w')
	lumlist = ['L_X', 'L_R', 'L_B', 'lambda_edd', 'M_in', 'M_out']
	fL.write("var luminosities_keys = [%s];\n" % ', '.join(['"%s"' % k for k in lumlist]));
	fL.write("var luminosities = [\n");
	
	for logMBH in [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]:
		#if logMBH != 8.0: continue
		for logeddrate in [-3, -1.5, -1, 0, 0.5]:
			eddrate = 10**logeddrate
			plt.figure(figsize=(10,7))
			MBH = 10**logMBH * u.Msun
			#print(logMBH, eddrate)
			prefix = "combinedstructure_2d_MBH%s_Mdot%s" % (logMBH, logeddrate)
			#print(prefix)
			plot_log_agn_postcard(MBH, eddrate)
			plt.savefig(prefix + '_log.pdf', bbox_inches="tight")
			plt.savefig(prefix + '_log.png', bbox_inches="tight")
			plt.close()
		fL.write("\t[\n");

		for logLAGN in numpy.arange(42, 47.4, 0.2):
			L_AGN = 10**logLAGN * (u.erg/u.s)
			L_AGN_edd = compute_eddington_luminosity(MBH)
			eddrate = L_AGN / L_AGN_edd
			plt.figure(figsize=(10,7))
			MBH = 10**logMBH * u.Msun
			#print(logMBH, eddrate)
			prefix = "combinedstructure_2d_MBH%s_LAGN%.1f" % (logMBH, logLAGN)
			plot_log_agn_postcard(MBH, eddrate)
			plt.savefig(prefix + '_log.pdf', bbox_inches="tight")
			plt.savefig(prefix + '_log.png', bbox_inches="tight")
			plt.close()
			
			Ls = dict(
				L_X = log_lum(bolcorr_hardX(L_AGN)),
				L_B = log_lum(bolcorr_B(L_AGN)),
				L_R = log_lum(compute_radiocore_luminosity(MBH, L_AGN)),
				lambda_edd = rround(eddrate.to(1).value),
				M_in = log10((L_AGN / c.c**2 / 0.1).to(u.Msun / u.yr).value),
				M_out = log10(compute_outflow_rate(L_AGN, Mstar = 1e11 * u.Msun, SFR = 0*u.Msun/u.yr).value),
			)
			fL.write("\t\t[" + ",".join(["%s" % Ls[k] for k in lumlist]) + "],\n")
			fL.flush();

			plt.figure(figsize=(10,7))
			plot_log_agn_postcard(MBH, eddrate,
				show_BH = True,
				show_corona = True,
				show_disk = True,
				show_BLR = True,
				show_NLR = False,
				show_jet = False,
				show_TOR = True,
				show_SOI = False,
				show_viewing = True,
				show_flows = False,
			)
			prefix = "combinedstructure_2d_MBH%s_LAGN%.1f" % (logMBH, logLAGN)
			plt.savefig(prefix + '_log_simple.pdf', bbox_inches="tight")
			plt.savefig(prefix + '_log_simple.png', bbox_inches="tight")
			plt.close()

		fL.write("\t],\n");
		fL.flush();

	fL.write("];");

