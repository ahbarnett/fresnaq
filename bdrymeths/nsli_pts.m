function u = nsli_pts(bx,by,wx,wy, lambdaz, xi, eta)
% NSLI_PTS  Fresnel aperture diffr via non-singular line-integral, arb targets
%
% u = nsli_pts(bx, by, wx, wy, lambdaz, xi, eta)
%  returns list of complex amplitude values u at (xi,eta) target points a const
%  distance z downstream, for an aperture for which nodes (bx,by) and vector
%  weights (wx,wy) are a good quadrature scheme for vector line integrals over
%  the aperture boundary. The incident wave is single-wavelength unit amplitude
%  plane wave. The scalar Fresnel approximation is used, and u does not include
%  the overall plane z propagation prefactor.
%
% Inputs:
%   bx, by - quadrature nodes (x,y) coords for aperture boundary line integral
%   wx, wy - corresponding vector quadrature weights for boundary line integral
%  lambdaz  - product of wavelength and downstream distance (in meters^2).
%             Fresnel number is R^2/lambdaz where R is a characteristic
%             aperture radius in meters.
%  xi, eta  - (real vectors, length M) x,y coords of target points in detector
%             plane.
%
% Outputs:
%  u        - complex scalar diffracted amplitude at each target point
%
% Notes:
% 1) We denote the aperture by Omega in R2. The input must be any good
%  quadrature for the 1D vector line integral on its boundary, dOmega, ie,
%  lists of nodes bx,by and vector weights bx,by, of same length n, such that
%
%    sum_{i=1}^n F(bx(i),by(i)) dot (wx(i),wy(i))  \approx
%        int_dOmega F(x,y) dot d(x,y)
%
%  for smooth (in R2 -> R2) functions F. In particular this should hold to the
%  desired accuracy for functions with oscillations down to the Fresnel zone
%  spacing for the largest source-target distances.
%  This allows various aperture boundary
%  parameterizations and unions of smooth or nonsmooth boundary components.
%
% 2) The formula for u approximated is aperture Fresnel diffraction in the
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
% 3) This is nearly a drop-in replacement for bdwf, in the flat aperture,
%  psi1=0 (on-axis), small-angle case. The difference is that a line integral
%  quadrature is demanded, allowing high-order accuracy, and that no in/out
%  target point testing is needed.

% Barnett 9/13/20
if nargin==0, test_nsli_pts; return; end

u = complex(nan*xi);                    % alloc complex outputs
for i=1:numel(xi)
  dx = bx - xi(i); dy = by - eta(i);    % displacements of nodes from ith targ
  r2 = dx.*dx + dy.*dy;
  f = (1-exp((1i*pi/lambdaz)*r2))./r2;  % unstable r2->0, but u st (barycentric)
  f(r2==0.0) = 0.0;                     % kill NaNs (targ = some node)
  u(i) = (1/2/pi) * sum((dx.*wy - dy.*wx) .* f);
end

%%%%%%%%%%%%
function test_nsli_pts
verb = 0;
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);        % smooth kite shape
fresnum = 20.0;        % Fresnel number (if char radius were R=1)
lambdaz=1/fresnum;     % wavelength times dist, recall Fres # = R^2/(lambda.z)
n = 500; m = 100;      % number of boundary, and radial, quadrature points;
                       % depend on Fresnel number & target pt; tested to 1e-12
t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
bxp = perispecdiff(bx); byp = perispecdiff(by);  % param-derivatives (if smooth)
wx = (2*pi/n)*bxp; wy = (2*pi/n)*byp;            % LI vector quadr wei

xi = 0.3; eta = -0.5;   % math test, one target, inside...
kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor for direct Fresnel integral
[xq yq wq] = curveareaquad(bx,by,wx,wy,m);       % areal quadrature
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2)) .* wq);
u = nsli_pts(bx,by,wx,wy, lambdaz, xi, eta);     % do it, w/ same bdry LI quadr
fprintf('NSLI (n=%d) 1-targ err, vs direct Fresnel quadr: %.3g\n',n,abs(u-ud));

j=400; d = [0 logspace(-17,-1,17)];  % check close targ dists from a node...
xi = bx(j) + 0*d; eta = by(j)-d;  ud = complex(nan*xi);
for i=1:numel(xi)
  ud(i) = kirchfac*sum(exp((1i*pi/lambdaz)*((xq-xi(i)).^2+(yq-eta(i)).^2)).*wq);
end
u = nsli_pts(bx,by,wx,wy, lambdaz, xi, eta);
disp('target distance from a node, resulting errors (vs areal method):')
disp([d(:) abs(u(:)-ud(:))])
if verb, figure; plot(bx,by,'.'); hold on; plot(xi,eta,'r*'); axis equal tight;
overlay_zones(lambdaz,xi(1),eta(1),'g-'); drawnow; end

lambdaz = 1/100.0; n = 2e3;   % speed test, with more bdry pts (higher F #) ...
t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
bxp = perispecdiff(bx); byp = perispecdiff(by);  % param-derivatives (if smooth)
wx = (2*pi/n)*bxp; wy = (2*pi/n)*byp;            % LI vector quadr wei
ximax = 1.0; ngrid = 1e2; g = ximax*(2*(0:ngrid-1)/ngrid - 1);   % targets
[xi eta] = ndgrid(g,g); M = numel(xi);
t0=tic;
u = nsli_pts(bx,by,wx,wy, lambdaz, xi, eta);     % do it
%u(ngrid/2+1,ngrid/2+1)       % use to check n-conv at the higher F #
t0=toc(t0);
fprintf('NSLI one lambda (n=%d, nTargs=%d) in %.3g s = %.3g targ-bdry pairs/s\n',n,M,t0,M*n/t0)
% 4e7 pairs/sec, v similar to BDWF, but is 6 lines of code!

if verb, figure; imagesc(g,g,abs(u.^2)'); colormap(hot(256)); colorbar; hold on;
plot(bx,by,'.'); axis xy equal tight; title('NSLI method: aperture'); end
