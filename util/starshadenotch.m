function [xj yj wj bx by] = starshadenotch(Np,Afunc,notch,n,m)
% STARSHADENOTCH  areal quadr for angular notch on one side of one petal.
%
% [xj yj wj bx by] = starshadenotch(Np,Afunc,notch,n,m)
%
%  returns areal nodes (xj,yj) and weights (wj) given notch struct describing
%  a notch in the upper edge of the 0th petal of a given ideal starshade with
%  Np petals, apodization Afunc, petal radius range [r0,r1].
%  For now the notch is cut in the theta direction, *not* the normal direction,
%  since the latter would require a little more calculus. (The difference
%  would be O(dT^2).)
%  The designed usage is to negate these weights so to subtract notch from
%  ideal SS.
%
% Inputs:
%  Np = # petals
%  Afunc = apodization from roughly 1 down to roughly 0, over domain [r0,r1]
%  notch = struct with fields:
%          notch.R0 = start radius
%          notch.R1 = end radius
%          notch.dT = function handle of radius, giving absolute change in
%                     apodization A(r) for the one petal, due to notch.
%                     Note that the notch depth thus scales with r, being
%                     2pi/Np * r * dT(r)
%  n = nodes covering theta of notch (convergence param)
%  m = nodes in radial range of notch (")
%
% Outputs:
% xj,yj - (areal) nodes
% wj    - corresp (areal) weights, positive
% bx,by - boundary nodes of notch area for plotting purposes only

% Barnett 12/15/20 to help Dumont/Shaklan.

if nargin==0, test_starshadenotch; return; end
R0=notch.R0; R1=notch.R1;   % notch radius range
if R0>R1, error('R0 must be <= R1'); end

[z w] = lgwt(m,R0,R1);              % note z descending order
[x v] = lgwt(n,0,1);
rr = ones(n,1)*z';    % n*m
thetaedge = (pi/Np)*Afunc(z);       % ideal top edge of 0th petal, at r nodes
tt = nan*rr; wj=tt;
for j=1:m                           % radii
  theta = (2*pi/Np)*notch.dT(z(j));  % this Theta of notch
  tt(:,j) = thetaedge(j) - theta*x;
  wj(:,j) = w(j) * z(j) * theta*v;        % weights from r dr dtheta
end
xj = rr.*cos(tt); xj = xj(:);
yj = rr.*sin(tt); yj = yj(:);
wj = wj(:);

% for plotting notch bdry only...
rr = [z; R0*ones(n,1); z(end:-1:1); R1*ones(n,1)];
tt = [thetaedge; (pi/Np)*(Afunc(R0)-2*notch.dT(R0)*x(end:-1:1)); thetaedge(end:-1:1)-(2*pi/Np)*notch.dT(z(end:-1:1)); (pi/Np)*(Afunc(R1)-2*notch.dT(R1)*x)];
bx = rr.*cos(tt); bx=bx(:);
by = rr.*sin(tt); by=by(:);


%%%%%%%%
function test_starshadenotch        % demo it w/ NI2
addpath ~/numerics/finufft/matlab/
verb = 1;
Np = 16;
cwd = fileparts(mfilename('fullpath')); file = [cwd '/../occulter/NI2'];
[~,r0,r1] = eval_sister_apod(file,0);  % get r0,r1 for full SS
A = @(r) eval_sister_apod(file,r);     % func

m = 400; n=40;                     % for ideal SS, known converged to 1e-6
[xj0 yj0 wj0 bx0 by0] = starshadequad(Np,A,r0,r1,n,m);  % ideal

notch.R0 = 7; notch.R1 = 8.3;
T0 = 0.1;   % notch depth in terms of apodization of the one petal (0: no notch)
notch.dT = @(r) T0+0*r;  % const (in theta width) notch.
%notch.dT = @(r) T0*sin(4*r).^2;  % variable-width notch
m=30; n=8; % need to check convergence wrt these two
[xj yj wj bx by] = starshadenotch(Np,A,notch,n,m);

area = sum(wj);   % test area various ways
parea = sum((circshift(bx,1)-bx).*by - (circshift(by,1)-by).*bx) / 2;  % polygon
%earea = T0*(2*pi/Np)*(notch.R1+notch.R0)/2*(notch.R1-notch.R0) % exact dT=const
fprintf('notch (m=%d,n=%d): area=%.6g (abs err vs poly: %.3g)\n',m,n,area,area-parea);

if verb, figure(1); clf; plot(bx,by,'+-'); hold on;
  plot(bx0,by0,'.-k'); scatter(xj,yj,10,wj); colorbar; axis equal tight; end

xj = [xj0; xj]; yj = [yj0; yj]; wj = [wj0; -wj];   % combined AQ: note -ve notch

tol=1e-9;                      % do diffraction...
ximax = 14.0; ngrid = 1e3;     % ngrid^2 centered target grid out to +-ximax
lambdaz = 21;       % lam.z in [16,21] m^2, low Fres#~5.
[u xigrid] = fresnaq_grid(xj, yj, wj, lambdaz, ximax, ngrid, tol, verb);
it = ngrid/2; jt = ngrid/2;    % grid pt to test, eg (1,1) is SW corner of grid
ut = u(it,jt); xi = xigrid(it); eta = xigrid(jt);     % math check u
fprintf('u(%.3g,%.3g) = %.12g + %.12gi\n',xi,eta,real(ut),imag(ut))

if verb, figure(2); clf;
  u = 1-u;               % convert aperture to occulter
  imagesc(xigrid,xigrid,log10(abs(u)'.^2));  % note transpose: x fast, y slow
  colorbar; caxis([-11 0.2]); hold on; axis xy equal tight;
  plot([bx0;bx0(1)], [by0;by0(1)], 'k+-','markersize',1);  % SS
  plot([bx;bx(1)], [by;by(1)], 'b+-','markersize',1);      % notch
  xlabel('\xi (m)'); ylabel('\eta (m)'); if verb>1, overlay_zones(lambdaz); end
  title(sprintf('NI2: log_{10} |u|^2, lambda.z=%.3g m^2',lambdaz))
end
% looks like this notch ruins the deep shadow, now only 1e-4 intensity
