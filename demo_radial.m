% Demo of FRESNAQ Frensel diffraction from occulter defined as radial function,
% evaluating on grid and at arbitrary target points.
% Barnett 9/8/20.

g = @(t) 1 + 0.3*cos(3*t);           % smooth radial func on [0,2pi)
fresnum = 10.0;        % Fresnel number (if char radius were R=1).
lambdaz=1/fresnum;     % wavelength times dist, recall Fres # = R^2/(lambda.z)
n=350; m=120;          % depend on Fresnel number; convergence must be tested
[xq yq wq bx by] = polarareaquad(g,n,m);   % areal quadrature
tol = 1e-9;            % desired accuracy
verb = 1;              % text output

ximax = 1.5; ngrid = 1e3;   % million-pt grid
[u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb);
ut = u(1,1); xi = xigrid(1); eta = xi;     % math check u in SW corner of grid
fprintf('u(%.3g,%.3g) = %.12g + %.12gi\n',xi,eta,real(ut),imag(ut))
% return                     % useful for convergence testing
kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2)) .* wq);
fprintf('abs error vs direct Fresnel quadr at (%.3g,%.3g) = %.3g\n',xi,eta,abs(ut-ud))

figure; subplot(1,2,1);
u = 1-u;               % convert aperture to occulter
imagesc(xigrid,xigrid,log10(abs(u)'.^2));  % note transpose: x is fast, y slow.
colorbar; hold on; plot(bx,by,'k-'); axis xy equal tight;
xlabel('\xi'); ylabel('\eta'); title('log_{10} |u|^2, occulter, grid');

M=1e6; xi = ximax*(2*rand(M,1)-1); eta = ximax*(2*rand(M,1)-1);  % million pts
u = fresnaq_pts(xq, yq, wq, lambdaz, xi, eta, tol, verb);
subplot(1,2,2);
u = 1-u;               % again convert aperture to occulter
scatter(xi,eta,10,log10(abs(u).^2));       % plotting a bit slow
colorbar; hold on; plot(bx,by,'k-'); axis xy equal tight;
xlabel('\xi'); ylabel('\eta'); title('log_{10} |u|^2, occulter, arb pts');
