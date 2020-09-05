% driver for fresnap code for Fresnel aperture diffraction, various starshades.
% Barnett 9/5/20
clear
verb = 2;  % verbosity

%Np=24; eval_apod = @eval_apod_NI2;  % actual NI2 starshad
Np=16; eval_apod = @eval_apod_erf;  % my analytic toy model, similar r0,r1

[~,r0,r1] = eval_apod(0);   % get apodization radius domain [r0,r1]
Afunc = @(r) eval_apod(r);  % func handle

n = 40; m = 120;
[xj yj wj bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb);

if verb>1, figure(1); clf; scatter(xj,yj,10,wj); axis equal tight; colorbar;
  hold on; plot([bx;bx(1)], [by;by(1)], 'k-'); drawnow
end
