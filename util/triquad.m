function [xq yq wq] = triquad(vx,vy,p)
% TRIQUAD  areal quadrature for a triangle.
%
% [xq yq wq] = triquad(vx,vy,p) returns areal quadrature with nodes
%  xq,yq, and weights wq, for the triangle T with counter-clockwise ordered
%  vertices (vx,vy).
%
% Inputs:
%  vx, vy - 3-element lists of x and y coords of CCW ordered vertices.
%  p      - order of scheme. There will be p^2 nodes.
%
% Outputs:
%  xq, yq - node x and y coordinates of nodes in T
%  wq     - corresponding positive weights
%
% With no arguments, does self-test

% Barnett 10/1/20
if nargin==0, test_triquad; return; end

dx = vx(2:3)-vx(1);   % 1st pt is the base pt. vectors to other two verts
dy = vy(2:3)-vy(1);
[t w] = lgwt(p,0,1);   % col vecs
[aa tt] = ndgrid(t,t);
xq = vx(1) + aa.*(dx(1)*(1-tt) + dx(2)*tt);
yq = vy(1) + aa.*(dy(1)*(1-tt) + dy(2)*tt);
wq = (dx(1)*(vy(3)-vy(2))-dy(1)*(vx(3)-vx(2))) * (t.*w) * w';   % outer
xq = xq(:);
yq = yq(:);
wq = wq(:);


%%%%%%%%%%%
function test_triquad
v = exp(2i*pi/3*(1:3)); vx = real(v); vy= imag(v);
area = 0.5 * ((vx(2)-vx(1))*(vy(3)-vy(1)) - (vy(2)-vy(1))*(vx(3)-vx(1)));
p=20;
[xq yq wq] = triquad(vx,vy,p);
sum(wq)-area
sum(wq.*yq.^3)   % y-odd should be zero

figure; plot([vx vx(1)], [vy vy(1)], '+-'); axis equal
hold on; scatter(xq,yq,20,wq); colormap(jet(256)); colorbar;
