# fresnap

**FRESN**el **AP**ertures: 
fast Fourier method for computation of scalar Fresnel diffraction from
0-1 (hard-edged) apertures and occulters.

# requirements

* MATLAB
* [FINUFFT](https://github.com/flatironinstitute/finufft) preferably v2.0

# installation

Install FINUFFT and compile its MATLAB interface.
Add ``finufft/matlab`` to your MATLAB path.
In MATLAB run ``startup`` then ``fresnap_grid`` to run tests which should
produce very small error outputs and take <0.1 sec to run.
