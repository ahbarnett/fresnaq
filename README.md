# fresnap

**FRESN**el **AP**ertures: 
Fast non-uniform Fourier areal quadrature
method for computation of scalar Fresnel diffraction from
0-1 (hard-edged) apertures and occulters. The simulation of the optics
of [starshades](http://sister.caltech.edu/) for exoplanet imaging
is one application.

*Here is an example. The intensity due to diffraction from a smooth kite-shaped occulter at Fresnel number around 20 is evaluated at one million points, to 9-digit accuracy, in 0.05 sec on a laptop. All fringes are correct and not numerical or sampling artifacts. The occulter boundary is shown in white:*

![fresnap demo image](pics/kite_grid.png "kite occulter example")

# requirements

* MATLAB
* [FINUFFT](https://github.com/flatironinstitute/finufft) preferably v2.0

# installation and testing

Install FINUFFT and compile its MATLAB interface.
Add ``finufft/matlab`` to your MATLAB path.
In MATLAB run ``startup`` then ``fresnap_grid`` to run tests which should
produce very small error outputs and take <0.1 sec to run.

The main library functions we provide are:

* ``fresnap_grid.m`` : compute diffracted amplitude on regular square centered grid
* ``fresnap_pts.m`` : compute diffracted amplitude at arbitrary target points

Both routines are documented and can be tested by calling them with no arguments.

# demos

The user should run ``demo_*.m`` to produce pretty pictures and see
how to call the library.
The above example is produced by setting ``verb=2`` in ``demo_paramcurve.m``.
Since areal quadratures are key, there are various
quadrature helper routines in ``util`` used in the tests and demos.

# comparisons with boundary line integral methods

See ``bdrymeths`` directory.

