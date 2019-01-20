============================================================
Visualisation of the Structure of Active Galactic Nuclei
============================================================

Model-driven illustrations of centres of galaxies.

============================
Ingredients and Assumptions
============================

This visualisation allows extrapolation outside the range of validity of the used models. This includes going to very high and low Eddington rates. Most reliable it should be at Eddington rates of 0.01-1.

This visualisation is an idealisation of a steady-state situation. In reality, the host galaxy morphology, chaotic accretion events leading to re- and mis-alignments and variability make the picture more complex.

The following components are implemented:

* BH: Black hole

  * Schwarzschild radius (r_g)

* AD: Accretion disk

  * Shakura & Sunyaev (1972) alpha-disk (see caveats below)
  * color according to local black body temperature, http://www.vendian.org/mncharity/dir3/blackbody/

* Corona: X-ray emitting corona

  * Geometry unknown, but compact (< 6 r_g), `Chartas et al. (2016) <https://ui.adsabs.harvard.edu/#abs/2016AN....337..356C/abstract>`_

* BLR: Broad line emission region

  * BLCH model (Baskin & Laor (2017), Czerny & Hryniewicz (2011))
  * describes where accretion disk radiation pressure lifts dusty material up (can be clumpy)
  * This model is known to produce a little too low covering factors
  * Note that at low L/M^(2/3), the BLR goes away (see Elitzur & Netzer (2016))

* RF TOR: Toroidal Obscurer

  * analytic form of the radiative fountain model from Wada (2016)
  * describes where the BH gravity and radiation pressure cancel each other, and the maximum velocity that can be reached
  
* NLR: Narrow line emission region

  * Appears to be shaped substantially by the host galaxy (see e.g. http://adsabs.harvard.edu/abs/2014MNRAS.442.2145P) as well as the nuclear obscurer

* Jet

  * from `Algaba et al. (2017) <https://ui.adsabs.harvard.edu/#abs/2017ApJ...834...65A/abstract>`_
  * Size is an empirical intermediate between

     * quasi-conical streamline (Blandford & Znajek 1977)
     * force-free steady jet solution (Narayan et al 2007; Tchekhovskoy et al 2008)
  * Base is assumed at r_g, maximum 10^8 pc.

* Sphere of influence

  * Where gravity from the BH dominates over the host galaxy
  * Computed the host gravity with M_BH-sigma from Kormendy & Ho (2013)

* Outflows & Inflows

  * Inflow is given by bolometric luminosity and accretion efficiency (assumed 10%)
  * Outflow relation from `Fluetsch et al. (2018) <https://ui.adsabs.harvard.edu/#abs/arXiv:1805.05352>`_, when setting SFR=0.

* Viewing angles

  * Are a bit difficult to interpret in log-log plots. For example, 
  * Note that views to the corona, to the BLR and to the TOR are not necessarily identical.

==========
Usage
==========

Clone this repository, and look at the end of combinedstructure.py how to use
the plot_log_agn_postcard function::

	from combinedstructure import plot_log_agn_postcard
	
	import astropy.units as u
	
	MBH = 1e8 * u.Msun
	eddrate = 0.1
	
	import matplotlib.pyplot as plt
	plt.figure(figsize=(10,7))
	
	plot_log_agn_postcard(MBH, eddrate, 
		# for Baskin & Laor (2017)
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
	)
	
	plt.show()


If useful for your work and paper, please cite this repository URL.
Code is MIT licensed.




