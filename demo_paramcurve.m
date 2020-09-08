% Demo of FRESNAP Frensel diffraction from occulter defined by parameterized
% curve, evaluating on grid and at arbitrary target points.
% Barnett 9/8/20.

% parameterization of smooth closed curve, for t in [0,2pi)...
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);
fresnum = 20.0;        % Fresnel number (if char radius were R=1).
lambdaz=1/fresnum;     % wavelength times dist, recall Fres # = R^2/(lambda.z)
n = 600; m = 100;      % depend on Fresnel number; convergence must be tested
t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);      % bdry points
bxp = perispecdiff(bx); byp = perispecdiff(by);  % derivatives wrt t (if smooth)
[xq yq wq] = curveareaquad(bx,by,bxp,byp,m);   % areal quadrature
% The rest of this code is same as demo_radial.m ...
tol = 1e-9;            % desired accuracy
verb = 1;              % verbosity

ximax = 1.5; ngrid = 1e3;   % million-pt grid
[u xigrid] = fresnap_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb);
ut = u(1,1); xi = xigrid(1); eta = xi;     % math check u in SW corner of grid
fprintf('u(%.3g,%.3g) = %.12g + %.12gi\n',xi,eta,real(ut),imag(ut))
% return                     % useful for convergence testing
kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2)) .* wq);
fprintf('abs error vs direct Fresnel quadr at (%.3g,%.3g) = %.3g\n',xi,eta,abs(ut-ud))

figure(1); clf; subplot(1,2,1);
u = 1-u;               % convert aperture to occulter
imagesc(xigrid,xigrid,log10(abs(u)'.^2));  % note transpose: x is fast, y slow.
colorbar; hold on; plot(bx,by,'k-'); axis xy equal tight;
xlabel('\xi'); ylabel('\eta'); title('log_{10} |u|^2, occulter, grid');

if verb>1, figure(2); clf; imagesc(xigrid,xigrid,abs(u)'.^2); % fig for repo
  colormap(hot(256)); colorbar; hold on; plot(bx,by,'w-'); axis xy equal tight;
  xlabel('\xi'); ylabel('\eta'); title('|u|^2, occulter, grid');
  v=caxis; caxis([0 v(2)]); drawnow; print -dpng -r200 kite_grid.png
end

M=1e6; xi = ximax*(2*rand(M,1)-1); eta = ximax*(2*rand(M,1)-1);  % million pts
u = fresnap_pts(xq, yq, wq, lambdaz, xi, eta, tol, verb);
figure(1); subplot(1,2,2);
u = 1-u;               % again convert aperture to occulter
scatter(xi,eta,10,log10(abs(u).^2));       % plotting a bit slow
colorbar; hold on; plot(bx,by,'k-'); axis xy equal tight;
xlabel('\xi'); ylabel('\eta'); title('log_{10} |u|^2, occulter, arb pts');


