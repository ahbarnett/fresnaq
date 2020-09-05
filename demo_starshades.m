% driver for fresnap code for Fresnel aperture diffraction, various starshades.
% Barnett 9/5/20
clear
verb = 1;  % verbosity

% NI2 starshade...
[~ r0 r1] = eval_apod_NI2(0);   % get apodization radius domain [r0,r1]
Afunc = @(r) eval_apod_NI2(r);  % func handle

n = 40; m = 120;

[xj yj wj bx by] = starshadequad(Np,Afunc,r1,r0,n,m,verb);

if verb>1, figure(2); clf; scatter(xj,yj,10,wj); axis equal tight; colorbar;
  hold on; plot([bx;bx(1)], [by;by(1)], 'k-'); drawnow
end
