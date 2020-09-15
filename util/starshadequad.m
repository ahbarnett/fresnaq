function [xj yj wj bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb,Aquad)
% STARSHADEQUAD  quadrature for area integral over apodized petal 0-1 occulter.
%
% [xj yj wj bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb)
%
%  Uses Theta(r) formula of (1)-(2) in Cady '12 to build area quadrature scheme
%  over starshade, given apodization function and other geometric parameters.
%
% Inputs:
%  Np = # petals
%  Afunc = apodization from roughly 1 down to roughly 0, over domain [r0,r1].
%  r0,r1 = apodization domain of radii (note differs from Cash'11 defn a,b)
%   radial apodization is 1 for r<r0, A((r-r0)/(r1-r0)) for r0<r<r1, 0 for r>r1.
%  n = nodes per petal width
%  m = nodes over disc and over radial apodization [r0,r1]. (exact may vary)
%  verb = (optional) 0,1, etc, verbosity.
%  Aquad = (optional) apodization quadrature type: 'a' Alpert (default),
%          'g' Gauss, 'u' Uniform grid (including end points).
%
% Outputs:
% xj,yj - (areal) nodes
% wj    - corresp (areal) weights
% bx,by - boundary nodes (those lying along profile, not petal tips or gaps)
%
% Design is perfectly symmetric; no perturbations for now.

% Barnett 8/24/20; changed to Afunc(r), swapped r0,r1, 9/5/20.
if nargin==0, test_starshadequad; return; end
if nargin<7, verb = 0; end
if nargin<8, Aquad = 'a'; end

if r0>r1, error('r0 must be <= r1'); end
switch Aquad          % set up a [0,1] quadr scheme
 case 'g'
  [z w] = lgwt(m,0,1);                  % radial petal quadrature, descending z
 case 'a' % or use Alpert, which may make couple more nodes than requested...
  [z w] = QuadNodesInterval(0,1,m,[],1,1,32); m=numel(z); z=z(end:-1:1); w=w(end:-1:1);  % flip to descending order
 case 'u'   % uniform grid, including endpoints, composite trap rule...
  z = linspace(0,1,m); h=z(2)-z(1); z=z(end:-1:1)'; w=h*[.5;ones(m-2,1);.5];
 otherwise
  error('unknown Aquad!')
end

% central disc...
N = ceil(0.3*n*Np);   % PTR in theta, rough guess so dx about same
t = 2*pi*(1:N)/N; wt = (2*pi/N);          % theta nodes, const weights
xj = nan(N*m,1); yj = xj; wj = xj;
for i=1:N                                              % loop over angles
  jj = (1:m) + (i-1)*m;                   % index list
  xj(jj) = cos(t(i))*r0*z; yj(jj) = sin(t(i))*r0*z;    % line of nodes
  wj(jj) = wt*r0^2*z.*w;            % theta weight times rule for r.dr on (0,r)
end
if verb
  fprintf('disc area err: %.3g\n',sum(wj) - pi*r0^2)
  fprintf('disc: %d nodes (%d angular, %d radial)\n',numel(wj),N,m)
end

% petals... radial outer loop (using z,w scheme from above)
[zn wn] = lgwt(n,-1,1);   % [-1,1] becomes petal width in theta. col vecs
r = r0 + (r1-r0)*z;       % radius quadrature nodes
A = Afunc(r);                                    % get all apodizations at once
for i=1:m
  wi = (r1-r0)*w(i)*r(i);                        % this radial weight for r.dr
  if A(i)>1 || A(i)<0, warning('apodization should obey 0<=A<=1'); end
  t = kron(ones(Np,1), (pi/Np)*A(i)*zn) + kron(2*pi*(1:Np)'/Np, ones(n,1)); %thetas
  ww = kron(ones(Np,1), wi*(pi/Np)*A(i)*wn);
  xj = [xj; r(i)*cos(t)]; yj = [yj; r(i)*sin(t)];    % append new nodes
  wj = [wj; ww];                                     % their weights
end
if verb
  fprintf('petals (Aquad=%s): %d nodes (%d angular, %d radial)\n',Aquad,Np*n*m,Np*n,m)
  fprintf('tot nodes: %d\n',numel(wj))
end

bx = nan(m,2*Np); by=bx;   % make bdry pt list, in order
r = r(end:-1:1); A = A(end:-1:1);
for p=1:Np
  t = (2*pi/Np)*p - pi/Np*A; bx(:,2*p-1) = r.*cos(t); by(:,2*p-1) = r.*sin(t);
  t = (2*pi/Np)*p + pi/Np*A(end:-1:1); bx(:,2*p) = r(end:-1:1).*cos(t); by(:,2*p) = r(end:-1:1).*sin(t);
end
bx = bx(:); by = by(:);


%%%%%%%%
function test_starshadequad
verb = 1;
Np = 16;

A = @(t) 1+0*t;    % first test disc radius r0 (filled petals)
r0 = 0.7; r1 = 1.3;
n = 20;    % theta across each petal
m = 30;    % radial
[xj yj wj] = starshadequad(Np,A,r0,r1,n,m,1);
if verb, figure(1); clf; scatter(xj,yj,10,wj); axis equal tight; colorbar; end
fprintf('disc err: %.3g\n',sum(wj) - pi*r1^2)

% choose either...
%eval_apod = @eval_apod_erf; Np = 16; ms=40:10:60; % analytic, w/ m-conv vals
eval_apod = @eval_apod_NI2; Np = 24; ms=60:20:180;  % actual

[~,r0,r1] = eval_apod(0);  % get r0,r1
A = @(r) eval_apod(r);    % func
%profile clear; profile on;
[xj yj wj] = starshadequad(Np,A,r0,r1,20,ms(1)); %disp(sum(wj))
%profile off; profile viewer;   % mostly loading file :)
if verb, figure(2); clf; scatter(xj,yj,10,wj); axis equal tight; colorbar; end

for m=ms        % check area converged ...
  [xj yj wj] = starshadequad(Np,A,r0,r1,20,m); disp(sum(wj))
end
% NI2 case jumps around at 1e-6 rel level; sucky function, due to only C^1.
