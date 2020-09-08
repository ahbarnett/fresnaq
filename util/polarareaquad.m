function [xq yq wq bx by] = polarareaquad(g,n,m)
% POLARAREAQUAD   quadrature for integrals over area inside polar function
%
% [xq yq wq] = polarareaquad(g,n,m) returns lists xq,yq of coordinates of nodes,
%  and corresponding weights wq, for 2D quadrature over a polar domain, ie,
%
%    sum_{j=1}^N f(xq(j),yq(j)) wq(j)   \approx   int_Omega f(x,y) dx dy
%
%  for all smooth functions f, where Omega is the polar domain defined by
%  radial function r=g(theta).
%
% [xq yq wq bx by] = polarareaquad(g,n,m) also returns boundary points, useful
%  for plotting.
%
% Inputs:
%   g = function handle for g(theta), theta in [0,2pi)
%   n = # theta nodes
%   m = # radial nodes
%
% Outputs:
%  xq, yq = column vectors (real, length N) of x,y coordinates of nodes
%  wq = column vector (real, length N) of weights
%
%  Note: no deriv of g func handle is needed.

% Barnett 8/24/20
if nargin==0, test_polarareaquad; return; end

t = 2*pi*(1:n)/n; wt = (2*pi/n);                       % theta nodes, const weights
bx = cos(t).*g(t); by = sin(t).*g(t);                  % boundary points
[xr,wr] = lgwt(m,0,1);                                 % rule for (0,1)
xq = nan(n*m,1); yq = xq; wq = xq;
for i=1:n                                              % loop over angles
  r = g(t(i)); jj = (1:m) + (i-1)*m;                   % this radius; index list
  xq(jj) = cos(t(i))*r*xr; yq(jj) = sin(t(i))*r*xr;    % line of nodes
  wq(jj) = wt*r^2*xr.*wr;            % theta weight times rule for r.dr on (0,r)
end

%%%%%%%%%%%
function test_polarareaquad
a = 0.3;
g = @(t) 1 + a*cos(3*t);   % radial func on [0,2pi)
n=50; m=20; [xq yq wq] = polarareaquad(g,n,m);
sum(wq) - pi*(1+a^2/2)        % analytic area formula
figure; scatter(xq,yq,10,wq); axis equal tight; colormap;
