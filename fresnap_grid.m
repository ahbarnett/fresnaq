function [u x1] = fresnap_grid(fbar, xj, yj, wj, xmax, M1, tol)
% FRESNELGRID   use quadr for Fresnel scalar diffraction on centered target grid
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
