 # calculation code I used
import numpy
from numpy import pi, exp, sin, cos, arccos, arctan2, tan, log, log10
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
from blrtor import pol2cart, bolcorr_hardX, bolcorr_B
from blrtor import log_lum, log_bhm, compute_dust_sublimation_radius, compute_heating_radius, compute_L_UV_at_angle, compute_L_X_at_angle, compute_rmax, compute_critical_angle, compute_eddington_luminosity

from blrtor import get_peak, compute_blr_shape, compute_sublimation_radius, get_blr_covering_factor

# from Baskin & Laor (2017)
rad_efficiency = 0.1
Z = 5 # metallicity

plt.figure(figsize=(15,3))
#logMBH = 8.0
#MBH = 10**logMBH * u.Msun
#L_AGN_edd = compute_eddington_luminosity(MBH)

#for eddrate, ls in [(0.00001, ':'), (0.023, ':'), (0.1, '--'), (0.45, '-')]:
#	L_AGN = eddrate * L_AGN_edd
#	L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
L_AGN = 1e45 * u.erg/u.s

for logMBH, ls in [(7, ':'), (8, '--'), (9, '-')]:
	MBH = 10**logMBH * u.Msun
	M_dot = (L_AGN / c.c**2 / rad_efficiency).to(u.Msun / u.yr)

	L_AGN_edd = compute_eddington_luminosity(MBH)
	eddrate = L_AGN / L_AGN_edd

	L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
	R_in = compute_sublimation_radius(MBH, M_dot)
	R = numpy.logspace(-3, 1, 400) * R_in
	
	z1, z2, zmax            = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='static')
	z1, z2, zmax_dyn        = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='dynamic', dyn_factor=2)
	CF = get_blr_covering_factor(z1, z2, R)
	Rpeak, Hpeak = get_peak(z2, z1, R)
	z1, z2, zmax_dyn_nodust = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='dustfree', dyn_factor=2)
	
	l, = plt.plot(R, zmax_dyn, 
		label="$M_\mathrm{BH}=%d$ $\lambda=%.3f$ $L_{46}=%.2f$ $\dot{M}=%.1f M_\odot/yr$ CF=%.0f%%" % (log_bhm(MBH), eddrate, L_AGN_46, M_dot.to(u.Msun / u.yr).value, CF*100)
	)
	color = l.get_color()
	#plt.plot(R.flatten(), zmax, ':', color=color)
	plt.plot(R.flatten(), zmax_dyn_nodust, '--', color=color)
	plt.plot(Rpeak, Hpeak, 'o', color=color)
	plt.fill_between(R.flatten(), 0, zmax_dyn, color=color, hatch='//')
	plt.vlines(R_in.to(u.pc).value, 0, Hpeak.to(u.pc).value, linestyles=[':'], colors=[color])
	#break

#plt.xlim(0, 0.1)
#plt.xlim(0.001, R_in.value * 2)
#plt.ylim(0, 0.02)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel("R [pc]")
plt.ylabel("z [pc]")
plt.legend(loc="best", prop=dict(size=8))
plt.savefig("windyblr_2d.pdf", bbox_inches="tight")
plt.close()

cmap = plt.cm.RdBu

MBH_lo, MBH_hi = 6.5, 9.5
MBH1d = numpy.logspace(MBH_lo, MBH_hi, 41) * u.Msun
L_AGN_lo, L_AGN_hi = 40, 46
L_AGN1d = numpy.logspace(L_AGN_lo, L_AGN_hi, 41) * u.erg / u.s

f = plt.figure('BLRsize')
plt.figure()
lightyear = c.c * u.day

CFgrid = []
sizegrid = []
for L_AGN in L_AGN1d:
	CFrow = []
	sizerow = []
	print(L_AGN)
	for MBH in MBH1d:
		color = cmap((log_bhm(MBH) - MBH_lo) / (MBH_hi - MBH_lo))
		L_AGN_46 = (L_AGN / (1e46 * u.erg/u.s)).to(1)
		M_dot = (L_AGN / c.c**2 / rad_efficiency).to(u.Msun / u.yr)
		R_in = compute_sublimation_radius(MBH, M_dot)
		R = numpy.logspace(-3, 1, 400) * R_in
		
		z1, z2, zmax_dyn = compute_blr_shape(MBH, M_dot, L_AGN, R, Z, variant='dynamic', dyn_factor=2)
		CF = get_blr_covering_factor(z1, z2, R)
		#plt.plot(L_AGN_46 / (MBH/(1e8 * u.Msun))**(-0.8), CF, 'o', color=color)
		LL = L_AGN_46**0.18 * (MBH/(1e8 * u.Msun))**0.15
		plt.plot(LL, CF, 'o', color=color)
		
		CFrow.append(CF)
		Rpeak, Hpeak = get_peak(z2, z1, R)
		Rmean = numpy.trapz(x=R, y=zmax_dyn * R) / numpy.trapz(x=R, y=zmax_dyn)
		#print(Rmean, Rpeak)
		sizerow.append(Rmean)
		f.gca().plot(L_AGN_46 * 1e46, (Rmean/lightyear).to(1), 'o', color=color)
	CFgrid.append(CFrow)
	sizegrid.append(sizerow)

x = numpy.logspace(-1.5, 0.2, 40)
#plt.plot(x, 0.2 * numpy.exp(log10(x)/0.45))
#plt.xscale('log')
plt.xlabel("Lift Luminosity $M_8^{0.15}L_{46}^{0.18}$")
plt.ylabel("Covering factor")
plt.savefig("windyblr_2d_LL.pdf", bbox_inches="tight")
plt.close()

plt.figure('BLRsize')
plt.yscale("log")
plt.xscale("log")
# from "Updating quasar bolometric luminosity corrections", Runnoe+12, lambda5100~LAGN
# BLR size data from THE RELATIONSHIP BETWEEN LUMINOSITY AND BROAD-LINE REGION SIZE IN ACTIVE GALACTIC NUCLEI, Kaspi+2005
plt.plot([1e42, 1e43, 1e44, 1e45, 1e46], [4.5, 7, 30, 100, 300], 'x', color='k', ms=8, mew=3)
plt.xlabel("Luminosity $L_{46}$ [erg/s]")
plt.ylabel("BLR size [lightdays]")
plt.savefig("windyblr_2d_BLRsize.pdf", bbox_inches="tight")
plt.close()



plt.imshow(CFgrid,
	origin="lower", aspect="auto",
	extent=[L_AGN_lo, L_AGN_hi, MBH_lo, MBH_hi], 
	cmap="Oranges", #cmap="gray_r"
)
cb = plt.colorbar()
CS = plt.contour(log_lum(L_AGN1d), log_bhm(MBH1d), numpy.transpose(CFgrid), [0.02, 0.05, 0.1, 0.2, 0.3])
plt.clabel(CS, inline=True)
plt.xlim(L_AGN_lo, L_AGN_hi)
plt.ylim(MBH_lo, MBH_hi)

#from labellines import labelLine, labelLines

# compute eddington rate:
def plot_edd_rate_line(MBH1d, eddrate=1, **kwargs):
	L_AGN_edd = compute_eddington_luminosity(MBH1d)
	x, y = log_lum(L_AGN_edd*eddrate), log_bhm(MBH1d)
	plt.plot(x, y, "--", label='$\lambda=%s$' % eddrate, **kwargs)
	plt.text(x[4], y[4], '$\lambda=%s$' % eddrate, 
		va='center', ha='center', rotation=60,
		bbox=dict(color='white', pad=0), size=8)
	#x = numpy.linspace(40, 46, 40)
	#plt.plot(x, (9 - (x-42)), ":", **kwargs)
plot_edd_rate_line(MBH1d, eddrate=1, color="k")
plot_edd_rate_line(MBH1d, eddrate=0.1, color="k", alpha=0.5)
plot_edd_rate_line(MBH1d, eddrate=0.01, color="k", alpha=0.25)
#lines = plt.gca().get_lines()
#labelLines(lines[-3:],zorder=2.5)
#L_AGN_edd = compute_eddington_luminosity(MBH1d)
#plt.plot(log_lum(L_AGN_edd), log_bhm(MBH1d), "--", color="k")
#plt.plot(log_lum(L_AGN_edd*0.1), log_bhm(MBH1d), "--", color="k", alpha=0.5)
#plt.plot(log_lum(L_AGN_edd*0.01), log_bhm(MBH1d), "--", color="k", alpha=0.25)

plt.xlabel("Bolometric Luminosity $L_\mathrm{AGN}$  [erg/s]")
plt.ylabel("Black Hole Mass $M_\mathrm{BH}$ [$M_\odot$]")
cb.set_label(r'Covering factor')
plt.savefig("windyblr_2d_LAGN.pdf", bbox_inches="tight")
plt.close()






