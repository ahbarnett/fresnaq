function [u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb)
% FRESNAQ_GRID  fast scalar Fresnel aperture diffraction, square target grid
%
% [u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol)
%  returns list of complex amplitude values u on ngrid-by-ngrid grid of targets
%  a distance z downstream, for an aperture for which nodes (xq,yq) and weights
%  wq are a good areal quadrature scheme in the (x,y) aperture plane, to
%  accuracy tol, for a single-wavelength unit amplitude plane incident wave,
%  in the scalar Fresnel approximation, not including propagation prefactor.
%  The grid is set by ximax and ngrid, and is the product of 1D grids xigrid,
%  in the xi,eta plane.
%
% ... = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb)
%  also gives diagnostics.
%
% Without any arguments, does self-test.
%
% Inputs:
%  xq, yq   - quadrature nodes over aperture (real vectors length N), in meters
%  wq       - corresponding quadrature weights (real vector length N)
%  lambdaz  - product of wavelength and downstream distance (in meters^2).
%             Fresnel number is R^2/lambdaz where R is a characteristic
%             aperture radius in meters.
%  ximax    - half-width of target grid region, in meters.
%  ngrid    - number of grid points in each dimension.
%  tol      - desired approximate error in u
%  verb     - verbosity (0,1,..)
%
% Outputs:
%  u        - complex scalar diffracted amplitude at each target point
%  xigrid   - list of target grid points used, along each axis (in meters)
%
% Notes:
%  The input quadrature scheme should integrate smooth 2D functions f, with
%  oscillations up to the needed Fresnel zone spatial frequency, to accuracy
%  tol, in the sense that
%
%    sum_{j=1}^N f(xq(j),yq(j)) wq(j) - int_Omega f(x,y) dx dy
%
%  is relatively no larger than of order tol, where Omega is the aperture
%  domain in R^2.
%
%  The formula for u approximated is aperture Fresnel diffraction in the
%  scalar Kirchhoff approximation:
%
%  u(xi,eta) = (i.lambdaz)^{-1} int_Omega
%                exp { i.pi.lambdaz [(x-xi)^2+(y-eta)^2] }  dxdy
%
%  This is evaluated for each point (xi(i),eta(i)) in target list i=1,...ngrid^2
%  This is as in (4) Cady, 2012 Opt. Expr. but without
%  the A or plane-z-propagation prefactors, and without the "1-" from Babinet.
%  Simply subtract from 1 to turn an aperture into an occulter.
%
%  The grid values in each dim are:
%      ximax*(2*(0:ngrid-1)/ngrid - 1)         if ngrid is even,
%      ximax*(2*(0.5:ngrid-0.5)/ngrid - 1)     if ngrid is odd.
%  The grid is thus centered on the origin, or half a grid-point off in the
%  even case.  (This arises from the "modeord" of FINUFFT, and could easily
%  be changed by phasing.)
%  Output values go in increasing xi fast then increasing eta slow (major axis),
%  ie [xi, eta] = ndgrid(g,g) where g is the above list of 1D grid values.
%
%  The algorithm uses a 2D type 1 NUFFT; depends on FINUFFT library. It is
%  faster than fresnaq_pts for large grids.
%
%  Example: see test_fresnaq_grid
%
%  Also see: fresnaq_pts


% Barnett 9/7/20
if nargin==0, test_fresnaq_grid; return; end
if nargin<8, verb = 0; end

t0=tic;
xigrid = ximax*(2*(0:ngrid-1)/ngrid - 1);            % for user edification
dxi = 2*ximax/ngrid;                                 % target grid spacing
if mod(ngrid,2), xigrid = xigrid + dxi/2; end        % hangle the odd case
sc = 2*pi/lambdaz;                                   % scale factor to become FT
dk = sc*dxi;                                         % scaled grid spacing
maxNU = max(abs(dk*[xq(:);yq(:)]));                  % max NU coord
cq = exp((1i*pi/lambdaz)*(xq.^2+yq.^2)) .* wq;       % premult by quadratic bit
if maxNU>3*pi                            % only do if needed, since a bit slow
  xq = mod(xq,2*pi/dk); yq = mod(yq,2*pi/dk);     % wrap in case grid too coarse
end
o = []; if verb>1, o.debug = 2; end                  % FINUFFT reporting
u = finufft2d1(dk*xq, dk*yq, cq, -1, tol, ngrid, ngrid, o);    % do it
% postmult by quadr bit (outer prod since separable grid)...
u = u .* (exp((1i*pi/lambdaz)*xigrid(:).^2) * exp((1i*pi/lambdaz)*xigrid.^2));
kirchfac = 1/(1i*lambdaz);                           % Kirchhoff prefactor
u = kirchfac * u;
if verb, fprintf('fresnaq_grid: N=%d quadr, M=%d targs, %.3g s\n',numel(xq),numel(u),toc(t0)), end


%%%%
function test_fresnaq_grid
fresnum = 10.0;        % Fresnel number
lambdaz=1/fresnum;   % since we test with O(1) radius aperture
g = @(t) 1 + 0.3*cos(3*t);   % smooth radial func on [0,2pi)
n=350; m=120; [xq yq wq] = polarareaquad(g,n,m);   % areal quadrature
tol = 1e-9;

ximax = 1.5; ngrid = 100;    % small grid for math test
[u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol);
[xi,eta] = ndgrid(xigrid,xigrid);    % recreate grid
i = 1;                       % check one target (grid SW corner)
kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi(i)).^2+(yq-eta(i)).^2)) .* wq);
fprintf('abs error vs direct Fresnel quadr at (%.3g,%.3g) = %.3g\n',xi(i),eta(i),abs(u(i)-ud))

% validate against slower fresnaq_pts for all targets...
upts = fresnaq_pts(xq, yq, wq, lambdaz, xi(:), eta(:), tol);
fprintf('max abs err fresnaq_pts vs grid: %.3g\n',norm(upts-u(:),inf))

ximax = 1.5; ngrid = 1e3; verb = 1;             % big grid for speed test
[u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb);
