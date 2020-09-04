function [xj yj wj bx by] = starshadequad(Np,Afunc,r1,r0,n,m,verb)
% STARSHADEQUAD  quadrature for area integral over apodized petal 0-1 occulter.
%
% [xj yj wj bx by] = starshadequad(Np,Afunc,r1,r0,n,m)
%
% Inputs:
%  Np = # petals
%  Afunc is apodization from roughly 1 down to roughly 0, over the domain [0,1].
%  [r1,r0] is apodization range of radii (note differs from Cash'11 defn a,b)
%   radial apodization is 1 for r<r1, A((r-r1)/(r0-r1)) for r1<r<r0, 0 for r>r0.
%  n = nodes per petal width
%  m = nodes over disc and over radial apodization [r1,r0]
%
% Outputs:
% xj,yj - nodes.
% wj - weights
% bx,by - boundary nodes just for plotting purposes
%
% Symmetry is not exploited, in case perturbations are desired.

% Barnett 8/24/20
if nargin==0, test_starshadequad; return; end
if nargin<7, verb = 0; end

if r1>=r0, error('r1 must be < r0'); end
[z w] = lgwt(m,0,1);

% central disc...
N = ceil(0.3*n*Np);   % PTR in theta, rough guess so dx about same
t = 2*pi*(1:N)/N; wt = (2*pi/N);          % theta nodes, const weights
xj = nan(N*m,1); yj = xj; wj = xj;
for i=1:N                                              % loop over angles
  jj = (1:m) + (i-1)*m;                   % index list
  xj(jj) = cos(t(i))*r1*z; yj(jj) = sin(t(i))*r1*z;    % line of nodes
  wj(jj) = wt*r1^2*z.*w;            % theta weight times rule for r.dr on (0,r)
end
if verb
  fprintf('disc area err: %.3g\n',sum(wj) - pi*r1^2)
  fprintf('disc: %d nodes (%d angular, %d radial)\n',numel(wj),N,m)
end

% petals... radial outer loop (using z,w scheme from above)
[zn wn] = lgwt(n,-1,1);   % [-1,1] becomes petal width in theta. col vecs
for i=1:m
  ri = r1 + (r0-r1)*z(i);                       % this radius
  wi = (r0-r1)*w(i)*ri;                         % this radial weight for r.dr
  A = Afunc(z(i));                              % this apodization
  if A>1 || A<0, warning('apodization A should be in [0,1]'); end
  t = kron(ones(Np,1), (pi/Np)*A*zn) + kron(2*pi*(1:Np)'/Np, ones(n,1)); %thetas
  ww = kron(ones(Np,1), wi*(pi/Np)*A*wn);
  xj = [xj; ri*cos(t)]; yj = [yj; ri*sin(t)];    % append new nodes
  wj = [wj; ww];                                 % their weights
end
if verb
  fprintf('petals: %d nodes (%d angular, %d radial)\n',Np*n*m,Np*n,m)
  fprintf('tot nodes: %d\n',numel(wj))
end

bx = nan(m,2*Np); by=bx;   % make bdry pt list, in order
z = z(end:-1:1);
r = r1 + (r0-r1)*z; A = Afunc(z);  % incr radial list and their apodizations
for p=1:Np
  t = (2*pi/Np)*p - pi/Np*A; bx(:,2*p-1) = r.*cos(t); by(:,2*p-1) = r.*sin(t);
  t = (2*pi/Np)*p + pi/Np*A(end:-1:1); bx(:,2*p) = r(end:-1:1).*cos(t); by(:,2*p) = r(end:-1:1).*sin(t);
end
bx = bx(:); by = by(:);


%%%%%%%%
function test_starshadequad
verb = 1;
Np = 16;

A = @(t) 1+0*t;    % first test disc radius r0
r1 = 0.7; r0 = 1.3; n = 20; m = 30;
[xj yj wj] = starshadequad(Np,A,r1,r0,n,m,1);
if verb, figure(1); clf; scatter(xj,yj,10,wj); axis equal tight; colorbar; end
sum(wj) - pi*r0^2

beta = 3.5; A = @(t) erfc(2*beta*(t-0.5))/2;        % my own apodization for now
fprintf('1-A(0)=%.3g, A(1)=%.3g\n',1-A(0),A(1))     % check its decay
if verb, figure(2); clf; r = linspace(0,1,1e3); plot(r,A(r),'-'); end
r1 = 0.6; r0 = 1.5; n = 20; m = 30;
[xj yj wj] = starshadequad(Np,A,r1,r0,n,m,1);
if verb, figure(1); clf; scatter(xj,yj,10,wj); axis equal tight; colorbar; end

% check area converged
[xj yj wj] = starshadequad(Np,A,r1,r0,20,30); disp(sum(wj))
[xj yj wj] = starshadequad(Np,A,r1,r0,25,30); disp(sum(wj))
[xj yj wj] = starshadequad(Np,A,r1,r0,25,40); disp(sum(wj))
