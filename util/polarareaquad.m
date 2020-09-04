function [xj yj wj] = polarareaquad(g,n,m)
% POLARAREAQUAD   quadrature for integrals over area inside polar function
%
% [x y w] = polarareaquad(g,n,m)
%  n = # theta nodes, m = # radial nodes.
%  note no deriv of g func handle needed.

% Barnett 8/24/20
if nargin==0, test_polarareaquad; return; end

t = 2*pi*(1:n)/n; wt = (2*pi/n);                       % theta nodes, const weights
bx = cos(t).*g(t); by = sin(t).*g(t);                  % boundary points
[xr,wr] = lgwt(m,0,1);                                 % rule for (0,1)
xj = nan(n*m,1); yj = xj; wj = xj;
for i=1:n                                              % loop over angles
  r = g(t(i)); jj = (1:m) + (i-1)*m;                   % this radius; index list
  xj(jj) = cos(t(i))*r*xr; yj(jj) = sin(t(i))*r*xr;    % line of nodes
  wj(jj) = wt*r^2*xr.*wr;            % theta weight times rule for r.dr on (0,r)
end

%%%%%%%%%%%
function test_polarareaquad
a = 0.3;
g = @(t) 1 + a*cos(3*t);   % radial func on [0,2pi)
n=50; m=20; [xj yj wj] = polarareaquad(g,n,m);
sum(wj) - pi*(1+a^2/2)        % analytic area formula
figure; scatter(xj,yj,10,wj); axis equal tight; colormap;
