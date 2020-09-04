function [xj yj wj] = curveareaquad(bx,by,bxp,byp,m)
% CURVEAREAQUAD  quadrature rule for area inside curve given nodes & 1st derivs
%
% function [x y w] = curveareaquad(bx,by,bxp,byp,m)
%   uses m-pt Gauss-Legendre in radial direction, given list (bx,by) of n
%   nodes and (bxp,byp) corresponding derivative vectors (with respect to a
%   parameterization on [0,2pi)). No 2pi/n factor should be given.
%
% curve need not enclose origin, but should be somewhat centered on it.
%
% To do:
% * make m adaptive to length and a worst-case dx requirement

% Barnett 8/24/20
if nargin==0, test_curveareaquad; return; end

[xr,wr] = lgwt(m,0,1);                                 % rule for (0,1)
wr = wr.*xr;                                           % rule a da on (0,1)
n = numel(bx);
xj = nan(n*m,1); yj = xj; wj = xj;
for i=1:n                                              % loop over angles
  jj = (1:m) + (i-1)*m;                                % index list
  xj(jj) = bx(i)*xr; yj(jj) = by(i)*xr;                % line of nodes
  c = bx(i)*byp(i) - by(i)*bxp(i);                     % area factor
  wj(jj) = (2*pi/n)*c*wr;        % theta weight times rule for a da on (0,1)
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
bw = (2*pi/n)*sp;           % speed weights
sum(bw)-per        % warm up by checking the perim, so bw are good bdry weights

m = 50;
[x y w] = curveareaquad(bx,by,bxp,byp,m);
sum(w) - pi*a*b        % check ellipse area
figure; scatter(x,y,10,w); axis equal tight; colorbar;
title('ellipse nodes & weights');
