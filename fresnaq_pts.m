function u = fresnaq_pts(xq, yq, wq, lambdaz, xi, eta, tol, verb)
% FRESNAQ_PTS  fast scalar Fresnel aperture diffraction, arbitrary target points
%
% u = fresnaq_pts(xq, yq, wq, lambdaz, xi, eta, tol)
%  returns list of complex amplitude values u at (xi,eta) target points a const
%  distance z downstream, for an aperture for which nodes (xq,yq) and weights
%  wq are a good areal quadrature scheme in the (x,y) aperture plane, to
%  accuracy tol, for a single-wavelength unit amplitude plane incident wave,
%  in the scalar Fresnel approximation, not including propagation prefactor.
%
% u = fresnaq_pts(xq, yq, wq, lambdaz, xi, eta, tol, verb) also reports text.
%
% Without any arguments, does self-test.
%
% Inputs:
%  xq, yq   - quadrature nodes over aperture (real vectors length N), in meters
%  wq       - corresponding quadrature weights (real vector length N)
%  lambdaz  - product of wavelength and downstream distance (in meters^2).
%             Fresnel number is R^2/lambdaz where R is a characteristic
%             aperture radius in meters.
%  xi, eta  - (real vectors or arrays, M elements) x,y coords of target points
%             in detector plane.
%  tol      - desired approximate error in u
%  verb     - verbosity (0,1,..)
%
% Outputs:
%  u        - complex scalar diffracted amplitude at each target point
%             (same shape as xi)
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
%  This makes u -> 1 in the limit of large Fresnel number.
%  Simply subtract from 1 to turn an aperture into an occulter.
%
%  The algorithm uses a 2D type 3 NUFFT; depends on FINUFFT library.
%
%  Example: see test_fresnaq_pts

% Barnett 9/5/20
if nargin==0, test_fresnaq_pts; return; end
if nargin<8, verb = 0; end

t0=tic;
% cq will be input strengths to NUFFT...
cq = exp((1i*pi/lambdaz)*(xq.^2+yq.^2)) .* wq;       % premult by quadratic bit
sc = 2*pi/lambdaz;                                   % scale factor to become FT
o = []; if verb>1, o.debug = 2; end                  % FINUFFT reporting
u = finufft2d3(xq,yq, cq, -1, tol, sc*xi(:),sc*eta(:), o);  % do the work
kirchfac = 1/(1i*lambdaz);                           % Kirchhoff prefactor
u = kirchfac * (u .* exp((1i*pi/lambdaz)*(xi(:).^2+eta(:).^2)));  % postmult bit
u = reshape(u,size(xi));
if verb, fprintf('fresnaq_pts: N=%d quadr, M=%d targs, %.3g s\n',numel(xq),numel(xi),toc(t0)), end

%%%%
function test_fresnaq_pts
fresnum = 10.0;        % Fresnel number
lambdaz=1/fresnum;   % since we test with O(1) radius aperture
g = @(t) 1 + 0.3*cos(3*t);   % smooth radial func on [0,2pi)
n=350; m=120; [xq yq wq] = polarareaquad(g,n,m);   % areal quadrature
tol = 1e-6;

xi = 1.5; eta = -0.5;   % math test: target to test at
u = fresnaq_pts(xq, yq, wq, lambdaz, xi, eta, tol)
kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2)) .* wq)
fprintf('abs error vs direct Fresnel quadr at (%.3g,%.3g) = %.3g\n\n',xi,eta,abs(u-ud))

M=1e6; xi = rand(M,1); eta = rand(M,1);   % speed test: bunch of targets
verb = 1;
u = fresnaq_pts(xq, yq, wq, lambdaz, xi, eta, tol, verb);
