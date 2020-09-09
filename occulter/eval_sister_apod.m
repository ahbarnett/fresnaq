function [A r0 r1] = eval_sister_apod(filehead,r)
% EVAL_SISTER_APOD  SISTER starshade apodization profile A(r) at arbitrary radii
%
% [A r0 r1] = eval_sister_apod(filehead,r) reads the file in SISTER format
%  giving the sampled profile, then uses interpolation to approximate A(r) at
%  the given list of r values.
%  The apodization radius range [r0,r1] is also extracted from other fields.
%
% With no arguments, does self-test and plot.
%
%  Notes:
%  1) apodization function A(r) is in notation of E. Cady, Opt. Expr. 2012,
%     starting at 1 for r<r0 and becoming 0 for r>r1.
%  2) reads a multi-MB file so don't call too often!

% Barnett 9/9/20
if nargin==0, test_eval_sister_apod; return; end

o = load(filehead);
r1 = o.occulterDiameter/2;      % upper apodization radius in meters
r0 = r1 - o.petalLength;        % lower apodization radius in meters

inds = (o.r>=r0);               % kill off useless inner radii samples
x = o.r(inds); x = x(:); y = o.Profile(inds);   % x=ordinates, y=data, col vecs

meth = 'spline'; % 'pchip'      % no different to 1e-8 with their sample dr=2mm
A = 0*r;
A(r<r0) = 1;                    % pad with 1 on left and 0 on right of domain
ii = find(r>=r0 & r<=r1);
if sum(ii)>0
  A(ii) = interp1(x,y,r(ii),meth);        % interpolate off given samples to r
end
  
%%%%%%
function test_eval_sister_apod
cwd = fileparts(mfilename('fullpath')); filehead = [cwd '/NI2'];  % in repo
eps = 1e-10;    % negligible but needed to avoid hitting 0 or 1 regions
[~,r0,r1] = eval_sister_apod(filehead,0);  % get r0,r1
A = @(r) eval_sister_apod(filehead,r);     % func
fprintf('1-A(%.5g)=%.3g, A(%.5g)=%.3g\n',r0,1-A(r0+eps),r1,A(r1-eps))  % check decay
figure; r = linspace(0,r1,1e4); plot(r,A(r),'-');
% note every time A(r) used, the file is loaded :(
title('eval\_sister\_apod NI2 test'); xlabel('r (meters)'); ylabel('A(r)');
vline([r0 r1]);