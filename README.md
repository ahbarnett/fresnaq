# fresnap

FRESNel APertures: 
fast Fourier method for computation of scalar diffraction from
0-1 apertures and occulters.

# requirements

* MATLAB
* [FINUFFT](https://github.com/flatironinstitute/finufft) preferably v2.0

# installation

Install FINUFFT and compile its MATLAB interface.
Add ``finufft/matlab`` to your MATLAB path.
In MATLAB run ``startup`` then ``test_fresnap`` to run tests which should
produce very small error outputs.
