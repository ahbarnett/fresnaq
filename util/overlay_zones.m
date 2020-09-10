function h = overlay_zones(lambdaz,xi,eta,lt,nmax)
% OVERLAY_ZONES  plots Fresnel zones on current plot centered at (xi,eta)
%
% h = overlay_zones(lambdaz,xi,eta,lt) plots on current axes a set of
%  Fresnel zones (as circles), without changing the current axes, and
%  returning a graphics handle to these circles. If (xi,eta), the center
%  of the zone, is not given, it is taken as (0,0).
%
%  Recall zones are defined as integer multiples of lambda/2 in extra path
%  length, in the Fresnel approximation. Thus for a source at infinity, the
%  nth zone has radius r_n := sqrt(n.lambda.z).
%
% Inputs:
%  lambdaz   - product of wavelength and downstream distance (in meters^2)
%  xi, eta   - (optional) center in the (xi,eta) plane of the Fresnel zones
%  lt        - (optional) MATLAB plot linetype for the circles
%  nmax      - (optional) maximum n to add zones up to
%
% Outputs:
%  h         - vector of MATLAB graphics handles to the set of circles

% Barnett 9/9/20
if nargin<2, xi=0; eta=0; end
if nargin<4 || isempty(lt), lt = 'w-'; end
if nargin<5, nmax = 100; end

hold on; v = axis;
n=1e3; t=2*pi*(0:n)/n; x=cos(t); y=sin(t);    % unit circle
h = [];
for zone=1:nmax
  r = sqrt(zone*lambdaz);
  h = [h plot(xi+r*x,eta+r*y,lt)];            % append plot handles
end
axis(v);
