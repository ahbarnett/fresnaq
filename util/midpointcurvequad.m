function [bx by wx wy] = midpointcurvequad(bx, by)
% MIDPOINTCURVEQUAD  midpt weights % new nodes for vector closed line integral
%
% [bx by wx wy] = midpointcurvequad(bx0, by0) given coordinates of nodes
%  on a curve in the plane, which need not correspond to any smooth
%  parameterization, return simple midpoint rule (new) nodes and weights.
%  Ie, such that,
%             sum_{i=1}^n F(bx(i),by(i)) dot (wx(i),wy(i))  \approx
%                   int_dOmega F(x,y) dot d(x,y)
%
%  for all smooth (in R2 -> R2) functions F. The midpoint rule is used by BDWF.
%  However it appears that the difference between this and the trapezoid rule
%  of crudecurvequad.m is minimal, and doesn't explain BDWF's occulter
%  shadow advantage over NSLI.
%
%  With no arguments, does self test, showing expected 2nd order (1/n^2)
%  accuracy.

% Barnett 9/28/20
if nargin==0, test_midpointcurvequad; return; end

wx = (circshift(bx,-1)-bx);            % vector weights simply displacements
wy = (circshift(by,-1)-by);
bx = (circshift(bx,-1)+bx)/2;          % replace by midpoint nodes
by = (circshift(by,-1)+by)/2;

%%%%%%%%
function test_midpointcurvequad        % test weights work for known LI
a = 1; b = 0.5;                        % ellipse
ns = 100*2.^(0:5);                     % convergence study
es = nan*ns;                           % errors
for i=1:numel(ns), n = ns(i);
  t = 2*pi*(1:n)/n;                    % bdry param
  bx = a*cos(t); by = b*sin(t);
  [bx by wx wy] = midpointcurvequad(bx, by);  % changes bx, by
  es(i) = sum(bx.*wy - by.*wx)/2 - pi*a*b;  % area = (int_dOmega x cross dl)/2
end
disp([ns(end), es(end)])
figure; loglog(ns,abs(es),'+-'); xlabel('n'); hold on; plot(ns,10*ns.^-2, 'r-');
legend('area error','O(1/n^2)'); axis tight; title('midpointcurvequad test');
