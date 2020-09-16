% minimal running and math test of making a PSF basis, hacked to use FRESNAQ.
% Run time is a few seconds. See README.
% Barnett 9/15/20

clear opt;
opt=WFIRST_starshade_basis_alex;
sister_basis(opt);    % writes to files, apparently

% now test it against one of the files...
u = load('/home/alex/physics/starshade/SISTER/sister_basis/non-spinning/NI2_24_16/425_552_nm//starshade_UtotL_Nx_16_pix_0425_0552_dl_10nm_dr_0.0_mas_psi_0.0_deg_ideal.mat');
v = load('starshade_UtotL_Nx_16_pix_0425_0552_dl_10nm_dr_0.0_mas_psi_0.0_deg_ideal.mat');

uerr = norm(u.UTotL(:)-v.UTotL(:),inf);
magerr = norm(abs(u.UTotL(:))-abs(v.UTotL(:)),inf);
fprintf('testing first NI2 PSF basis file against BDWF reference version:\n\tmax E ampl error %.3g, max E complex number error %.3g\n',magerr,uerr)
