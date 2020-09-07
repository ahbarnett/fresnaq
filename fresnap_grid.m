function [u x1] = fresnap_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb)
% FRESNAP_GRID  fast scalar Fresnel aperture diffraction, square target grid
%
% u = fresnap_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol)
%  returns list of complex amplitude values u on ngrid-by-ngrid grid of targets
%  a distance z downstream, for an aperture for which nodes (xq,yq) and weights
%  wq are a good areal quadrature scheme in the (x,y) aperture plane, to
%  accuracy tol, for a single-wavelength unit amplitude plane incident wave,
%  in the scalar Fresnel approximation, not including propagation prefactor.
%
% u = fresnap_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb)
%  also reports text
%
% Without any arguments, does self-test.
%
% Inputs:
%  xq, yq   - quadrature nodes over aperture (real vectors length N), in meters
%  wq       - corresponding quadrature weights (real vector length N)
%  lambdaz  - Frensel parameter (in meters^2) = product of wavelength and
%             downstream distance. Fresnel number is R^2/lambdaz where R is a
%             characteristic aperture radius.
%  ximax    - width of target grid in meters.

***
(real vectors, length M) x,y coords of target points in detector

%             plane.
%  tol      - desired approximate error in u
%  verb     - verbosity (0,1,..)
%
% Outputs:
%  u        - complex scalar diffracted amplitude at each target point
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
%  This is evaluated for each point (xi(i),eta(i)) in target list i=1,...M.
%  This is as in (4) Cady, 2012 Opt. Expr. but without
%  the A or plane-z-propagation prefactors, and without the "1-" from Babinet.
%  Simply subtract from 1 to turn an aperture into an occulter.
%
%  The algorithm uses a 2D type 3 NUFFT; depends on FINUFFT library.
%
%  Example: see test_fresnap_pts

% Barnett 9/5/20
if nargin==0, test_fresnap_pts; return; end
if nargin<8, verb = 0; end

t0=tic;
% cq will be input strengths to NUFFT...
cq = exp((1i*pi/lambdaz)*(xq.^2+yq.^2)) .* wq;       % premult by quadratic bit
sc = 2*pi/lambdaz;                                   % scale factor to become FT
o = []; if verb>1, o.debug = 2; end                  % FINUFFT reporting
u = finufft2d3(xq,yq, cq, -1, tol, sc*xi,sc*eta, o); % do the work
kirchfac = 1/(1i*lambdaz);                           % Kirchhoff prefactor
u = kirchfac * (u .* exp((1i*pi/lambdaz)*(xi.^2+eta.^2)));  % postmult quadratic
if verb, fprintf('fresnap_pts: N=%d quadr, M=%d targs, %.3g s\n',numel(xq),numel(xi),toc(t0)), end

%%%%
function test_fresnap_pts
fresnum = 10.0;        % Fresnel number
lambdaz=1/fresnum;   % since we test with O(1) radius aperture
g = @(t) 1 + 0.3*cos(3*t);   % smooth radial func on [0,2pi)
n=350; m=120; [xq yq wq] = polarareaquad(g,n,m);   % areal quadrature
tol = 1e-6;

xi = 1.5; eta = -0.5;   % math test: target to test at
u = fresnap_pts(xq, yq, wq, lambdaz, xi, eta, tol)
kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2)) .* wq)
fprintf('abs error vs direct Fresnel quadr at (%.3g,%.3g) = %.3g\n\n',xi,eta,abs(u-ud))

M=1e6; xi = rand(M,1); eta = rand(M,1);   % speed test: bunch of targets
verb = 1;
u = fresnap_pts(xq, yq, wq, lambdaz, xi, eta, tol, verb);

%
% [u x1] = fresnelgrid(fbar, xj, yj, wj, xmax, M1, tol) evaluates scalar Fresnel
%  diffraction integral for aperture defined by quadrature nodes (xj,yj) and
%  weights wj.
%
%  Centered uniform target grid of same size and spacing in x and in y, for now.
%  (so that type 1 can be simply used, not type 3)
%  *** formula.
%
%  Includes 1/(i.lambda.d) prefactor but not exp(i.k.d) downstream phase factor.
%  This makes u -> 1+0i for a big aperture.
%
% Inputs:
%  fbar   - "increased" Fresnel #: fbar = 2pi.f = k_free/d = 2pi/(lambda.d)
%  xj, yj, wj - quadrature rule for the aperture (assumed converged for desired
%           fbar and target grid. Nodes (xj,yj), weights wj.
%  xmax   - half-size of target box
%  N1     - number of points
%  tol    - desired precision in u, eg 1e-6.
% 
% Outputs:
%  u      - 2D complex array (x down, y across) of amplitude answers
%  x1     - corresponding 1D grid of size N1 (holds for both dims):
%           {-xmax, -xmax+dx, ... xmax-dx}  where dx = 2*xmax/M1
%
% To do:
% * need to check dk*(xj,yj) doesn't fall outside +-3pi, and 2pi-wrap it.
%   Could fail for coarse (undersampled) output grids, but user could want such
%   a grid.

% Barnett 8/24/20
if nargin==0, test_fresnelgrid; return; end

dx = 2*xmax/M1;
x1 = dx * (-M1/2:(M1/2-1));    % target grid each dim
t0 = tic;

dk = fbar*dx;
k1 = fbar*x1;                  % freq grid is simply scaled location grid
maxNU = max(abs(dk*[xj(:);yj(:)]));
fprintf('fresnelgrid: fbar=%.3g, xmax=%.3g: kmax=%.3g, max NU input=%.3g\n',fbar,xmax,max(k1),maxNU)
cj = exp(0.5i*fbar*(xj.^2+yj.^2)) .* wj;        % premult by a quadratic bit
if maxNU>3*pi                            % only do if needed, since a bit slow
  xj = mod(xj,2*pi/dk); yj = mod(yj,2*pi/dk);     % wrap in case grid too coarse
end
u = finufft2d1(dk*xj, dk*yj, cj, -1, tol, M1, M1);   % do it: M1^2 output nodes
u = u .* (exp(0.5i*fbar*x1(:).^2) * exp(0.5i*fbar*x1.^2));  % postmult by quadr bit
kirchfac = fbar/(2i*pi);   % Kirchhoff approx prefactor = 1/i.lambda.d
u = kirchfac * u;

fprintf('fresnelgrid done in %.3g s\n',toc(t0))

%%%%%%%
function test_fresnelgrid
verb = 1;
fbar = 100;   % "increased" Fresnel #, fbar = 2pi.f = k_free/d = 2pi/(lambda.d)
g = @(t) 1 + 0.3*cos(3*t);   % radial func on [0,2pi)
n=350; m=120; [xj yj wj] = polarareaquad(g,n,m);
xmax = 2.0;       % targ params: box half-size
M1 = 1e3;         % target location grid size per dim (eg, <=50 to trigger mod)
tol = 1e-9;
xt = -xmax; yt=-xmax;   % test pt: make sure it's on the targ grid

if verb
  t=2*pi*(1:n)/n; bx=g(t).*cos(t); by=g(t).*sin(t);      % get bdry pts
  figure(1); clf; subplot(1,2,1);
  plot([bx bx(1)],[by by(1)],'-'); hold on; plot(xj,yj,'.'); axis equal tight;
  plot(xt,yt,'*'); xlabel('x_1'); ylabel('x_2'); title('Aperture + nodes');
  subplot(1,2,2); scatter(xj,yj,10,cos(0.5*fbar*((xj-xt).^2+(yj-yt).^2)));
  axis equal tight; xlabel('x_1'); ylabel('x_2'); title('Re Fresnel integrand');
end

[u x1] = fresnelgrid(fbar, xj, yj, wj, xmax, M1, tol);         % do it

if verb, figure(2); clf;
  subplot(1,2,1); imagesc(x1,x1,abs(u)'); colorbar; hold on;
  title('|u| for aperture'); plot([bx bx(1)],[by by(1)],'-'); plot(xt,yt,'*');
  axis xy equal tight; colormap(hot(256));
  subplot(1,2,2); imagesc(x1,x1,abs(1.0-u)'); colorbar; hold on;
  title('|u| for shade'); plot([bx bx(1)],[by by(1)],'-'); plot(xt,yt,'*');
  axis xy equal tight; colormap(hot(256));
end

% check it...
kirchfac = fbar/(2i*pi);   % Kirchhoff approx prefactor = 1/i.lambda.d
ut = kirchfac * sum(exp(0.5i*fbar*((xj-xt).^2+(yj-yt).^2)) .* wj)
j1 = find(x1==xt); j2 = find(x1==yt); errt = u(j1,j2) - ut;
fprintf('abs error vs slow Fresnel quad at (%.3g,%.3g) = %.3g\n\n',xt,yt,abs(errt))
