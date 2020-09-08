# fresnap

**FRESN**el **AP**ertures: 
fast Fourier method for computation of scalar Fresnel diffraction from
0-1 (hard-edged) apertures and occulters.

<img align="left" src="pics/kite_grid.png" width="500">

Example computation of diffracted amplitude on a million-point grid
to 9-digit accuracy, in 0.05 sec on a laptop.

---

# requirements

* MATLAB
* [FINUFFT](https://github.com/flatironinstitute/finufft) preferably v2.0

# installation and testing

Install FINUFFT and compile its MATLAB interface.
Add ``finufft/matlab`` to your MATLAB path.
In MATLAB run ``startup`` then ``fresnap_grid`` to run tests which should
produce very small error outputs and take <0.1 sec to run.

The main user functions are:

* ``fresnap_grid.m`` : compute diffracted amplitude on regular square centered grid
* ``fresnap_pts.m`` : compute diffracted amplitude at arbitrary target points

Both routines can be tested by calling them with no arguments.
Also see ``demo_*.m``.
