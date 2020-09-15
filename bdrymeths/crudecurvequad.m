function [wx wy] = crudecurvequad(bx, by)
% CRUDECURVEQUAD  trapezoid weights for vector closed line integral, given nodes
%
% [wx wy] = crudecurvequad(bx, by) given coordinates of nodes defining a
%  curve in the plane, which need not correspond to any smooth
%  parameterization, return simple composite trapezoid rule weights. Ie, such
%  that,
%             sum_{i=1}^n F(bx(i),by(i)) dot (wx(i),wy(i))  \approx
%                   int_dOmega F(x,y) dot d(x,y)
%
%  for all smooth (in R2 -> R2) functions F. The trap. rule corresponds to
%  linear interpolation of the integrand between nodes, of course.
%
%  With no arguments, does self test, showing 2nd order (1/n^2) accuracy.

% Barnett 9/13/20
if nargin==0, test_crudecurvequad; return; end

wx = (circshift(bx,-1)-circshift(bx,1))/2;
wy = (circshift(by,-1)-circshift(by,1))/2;

%%%%%%%%
function test_crudecurvequad           % test weights work for known LI
a = 1; b = 0.5;                        % ellipse
ns = 100*2.^(0:5);                     % convergence study
es = nan*ns;                           % errors
for i=1:numel(ns), n = ns(i);
  t = 2*pi*(1:n)/n;                    % bdry param
  bx = a*cos(t); by = b*sin(t);
  [wx wy] = crudecurvequad(bx, by);
  es(i) = sum(bx.*wy - by.*wx)/2 - pi*a*b;  % area = (int_dOmega x cross dl)/2
end
disp([ns(end), es(end)])
figure; loglog(ns,abs(es),'+-'); xlabel('n'); hold on; plot(ns,10*ns.^-2, 'r-');
legend('area error','O(1/n^2)'); axis tight; title('crudecurvequad test');


