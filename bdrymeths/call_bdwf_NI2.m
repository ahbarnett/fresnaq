% document/test/demo Cady's bdwf code for Fresnel diffraction from aperture
% Barnett 9/2/20

clear

% load in a boundary curve for fun:
load /home/alex/physics/starshade/SISTER/input_scenes/locus/in/NI2_test_case_1em10.mat
figure; plot(xVals,yVals,'.'); axis equal tight;




% some params pulling from whatever calls bdwf when do:
% addpath config
% addpath software
% clear opt;opt=WFIRST_starshade_basis;sister_basis(opt);
% as per handbook.
Z = 37242256.6835035;
deltaX = 0;
deltaY = 0;
dxO = 0.15;     % ?
lambda = (425:10:550) * 1e-9;    % wavelenghts in m
nO = 16;
psi1 = 0;   % then 4.8e-9 and its multiples...
psi2 = 0;   % then pi/2   - why?


t0=tic;
E = bdwf(xVals, yVals, [], Z, lambda, dxO, nO, psi1, psi2, deltaX, deltaY);
fprintf('bdwf done in %.3g s\n',toc(t0))


