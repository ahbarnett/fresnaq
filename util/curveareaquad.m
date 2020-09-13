function [xq yq wq] = curveareaquad(bx,by,wx,wy,m)
% CURVEAREAQUAD  areal quadrature rule given nodes & vec weights for bdry quadr
%
% function [xq yq wq] = curveareaquad(bx,by,wx,wy,m) returns lists xq,yq of
%  coordinates of nodes, and corresponding weights wq, for areal quadrature
%  over a bounded 2D domain Omega, ie such that
%
%    sum_{j=1}^N f(xq(j),yq(j)) wq(j)   \approx   int_Omega f(x,y) dx dy
%
%  for all smooth functions f.
%  It uses m-pt Gauss-Legendre in radial direction. The input must be any good
%  quadrature for the 1D vector line integral on the boundary, dOmega, ie,
%  lists of nodes bx,by and vector weights bx,by, of same length n, such that
%
%    sum_{i=1}^n F(bx(i),by(i)) dot (wx(i),wy(i))  \approx
%        int_dOmega F(x,y) dot d(x,y)
%
%  for smooth (in R2 -> R2) functions F.  This allows various
%  parameterizations and unions of smooth or nonsmooth boundary components.
%
% Inputs:
%   bx, by - quadrature nodes (x,y) coords for boundary line integral
%   wx, wy - corresponding vector quadrature weights for line integral
%   m      - desired # radial nodes
%
% Outputs:
%  xq, yq  - column vectors (real, length N) of x,y coordinates of nodes
%  wq      - column vector (real, length N) of weights
%
% Notes:
%  * curve need not enclose origin, but should be vaguely centered on it
%    since then m need not be too large, and weights do not have too many
%    cancellations. If Omega is not star-shaped about 0, the areal nodes lie
%    outside Omega (strangely), so the integrand they apply to must be smooth
%    in this larger domain formed from convex combinations of Omega and 0.

% To do:
% * make m adaptive to length and a worst-case dx requirement?
%
% Barnett 8/24/20, changed input to generic dOmega quadr
if nargin==0, test_curveareaquad; return; end

[xr,wr] = lgwt(m,0,1);                                 % rule for (0,1)
wr = wr.*xr;                                           % rule a da on (0,1)
n = numel(bx);
xq = nan(n*m,1); yq = xq; wq = xq;
for i=1:n                                              % loop over angles
  jj = (1:m) + (i-1)*m;                                % index list
  xq(jj) = bx(i)*xr; yq(jj) = by(i)*xr;                % line of nodes
  c = bx(i)*wy(i) - by(i)*wx(i);                       % area factor
  wq(jj) = c*wr;        % theta weight times rule for a da on (0,1)
end

%%%%%
function test_curveareaquad

a = 1; b = 0.5;   % ellipse
e2 = 1-(b/a)^2;   % square of eccentricity
[~,per] = ellipke(e2); per=4*a*per;  % exact perimeter
n = 100;
t = 2*pi*(1:n)/n;                    % bdry param
bx = a*cos(t); by = b*sin(t); bxp = -a*sin(t); byp = b*cos(t);
sp = sqrt(bxp.^2+byp.^2);
bw = (2*pi/n)*sp;           % scalar LI weights ("speed weights")
sum(bw)-per        % warm up by checking the perim, so bw are good bdry weights
wx = (2*pi/n)*bxp; wy = (2*pi/n)*byp;  % vector LI weights

m = 50;
[x y w] = curveareaquad(bx,by,wx,wy,m);
sum(w) - pi*a*b        % check ellipse area
figure; scatter(x,y,10,w); axis equal tight; colorbar;
title('ellipse nodes & weights');
