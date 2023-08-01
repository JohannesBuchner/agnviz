============================================================
Visualisation of the Structure of Active Galactic Nuclei
============================================================

Model-driven illustrations of centres of galaxies.

See it in action at https://johannesbuchner.github.io/agnviz/

============================
Ingredients and Assumptions
============================

This visualisation allows extrapolation outside the range of validity of the used models. This includes going to very high and low Eddington rates. Most reliable it should be at Eddington rates of 0.01-1.

This visualisation is an idealisation of a steady-state situation. In reality, the host galaxy morphology, chaotic accretion events leading to re- and mis-alignments and variability make the picture more complex.

.. image:: https://raw.githubusercontent.com/JohannesBuchner/agnviz/master/img/combinedstructure_2d_MBH8.0_LAGN44.0_log.png
	:width: 300
	:target: https://johannesbuchner.github.io/agnviz/
	:alt: Click to play with visualisation

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
  * The left edge is the dust sublimation radius of the smallest graphites (2000K).
  * This model is known to produce a little too low covering factors
  * Note that at low L/M^(2/3), the BLR goes away (see Elitzur & Netzer (2016))

* RF TOR: Toroidal Obscuring Region

  * Material is likely clumpy
  * analytic form of the radiative fountain model from `Wada (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJ...812...82W/abstract>`_
  * describes where dust is no longer accelerated by radiation pressure. Material can pile up near and to the right of the line up to >>10 pc
  * Note that at low L/M^(2/3), the region is closed (high covering fraction).
  * This region starts near the sublimation radius of most dust grains (~1.3pc sqrt(L46)), but will be stratified by dust grain size (possibly down to the BLR, see above).

* NLR: Narrow line emission region

  * Appears to be shaped substantially by the host galaxy (see e.g. `Prieto et al. 2014 <http://adsabs.harvard.edu/abs/2014MNRAS.442.2145P>`_) as well as the nuclear obscurer.

* Jet

  * Width is an empirical intermediate (from `Algaba et al. (2017) <https://ui.adsabs.harvard.edu/#abs/2017ApJ...834...65A/abstract>`_) between

     * quasi-conical streamline (Blandford & Znajek 1977)
     * force-free steady jet solution (Narayan et al 2007; Tchekhovskoy et al 2008)
  * Base is assumed at r_g, maximum 10^8 pc.

* Sphere of influence

  * Where gravity from the BH dominates over the host galaxy
  * Computed the host gravity with M_BH-sigma from Kormendy & Ho (2013)

* Outflows & Inflows

  * Inflow onto BH is given by bolometric luminosity and accretion efficiency (assumed 10%)
  * Outflow relation from `Fluetsch et al. (2018) <https://ui.adsabs.harvard.edu/#abs/arXiv:1805.05352>`_, when setting SFR=0.

* Viewing angles

  * Are a bit difficult to interpret in log-log plots.
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
Code is BSD licensed, see LICENSE file.


Conferences and talks about AGN
=================================

Conferences:
--------------

* "Massive black holes in evolving galaxies: from quasars to quiescence" Institut d'Astrophysique de Paris, June 25-29, 2018 http://www.iap.fr/vie_scientifique/colloques/Colloque_IAP/2018/i-program.php
* "BHI Conference 2018" May 9-11, 2018 https://www.youtube.com/playlist?list=PL8_xPU5epJddhcKBgbwzKYPOWIMKSV29J
* "Hidden Monsters"  8-12 August 2016 Dartmouth https://www.dartmouth.edu/~hiddenmonsters/presentations_tab.php
* "Confronting MHD Theories of Accretion Disks with Observations" Kavli Institute for Theoretical Physics (Jan 9 - Mar 30, 2017) http://online.kitp.ucsb.edu/online/disks17/
* "Massive Black Holes: Birth, Growth and Impact" Kavli Institute for Theoretical Physics (Aug 5-9, 2013)" http://online.itp.ucsb.edu/online/bholes-c13/


Talks:
-------

* Black Hole Demographics in Active and Inactive Galaxies - L. Ho https://youtu.be/rzfIZ1xGH6s?t=133
* Christine Done : https://youtu.be/2LH9xFv8uCc?t=235 
* The Monster Roars: Feedback and the Co-Evolution of Galaxies and Black Holes - P. Hopkins https://www.youtube.com/watch?v=_xk2pDxQmls

More general:
-------------

* Star Formation, Magnetic Fields, and Diffuse Matter in the Galaxy

  * https://www.youtube.com/playlist?list=PLR0icYcBhA7Rv-QymeskXF15Kd-Mnh_T2

* CfA High-energy Astrophysics devision seminars & lectures

  * https://www.youtube.com/channel/UCMFEeX24_lviXNhek5-FFLA

* Rashid Sunyaev lecture https://www.youtube.com/watch?v=8V-brhbokTQ https://www.youtube.com/watch?v=1SbwwWlDkeQ

* https://www.youtube.com/user/ICTStalks/search?query=agn
* https://www.youtube.com/user/IACvideos/videos
* https://www.youtube.com/channel/UCTuACIrLKPTlp6XMZbeipig/videos
