function [A r0 r1] = eval_apod_NI2(r)
% EVAL_APOD_NI2  NI2 starshade apodization profile A(r) at arbitrary radii
%
% [A r0 r1] = eval_apod_NI2(r) read the file giving the sampled NI2 profile and
%  uses interpolation to evaluate A(r) at the given list of r values.
%  The apodization radius domain [r0,r1] is also returned.
%
% With no arguments, does self-test and plot.
%
%  Notes:
%  1) apodization function A(r) is in notation of E. Cady Opt. Expr. 2012,
%     starting at 1 for r<r0 and becoming 0 for r>r1.
%  2) reads a 3MB file so don't call too often!

% Barnett 9/5/20
if nargin==0, test_eval_apod_NI2; return; end

cwd = fileparts(mfilename('fullpath'));
o = load([cwd '/NI2.mat']);
r1 = o.occulterDiameter/2;      % upper apodization radius in meters
r0 = r1 - o.petalLength;        % lower apodization radius in meters

inds = (o.r>=r0);               % kill off useless inner radii samples
x = o.r(inds); x = x(:); y = o.Profile(inds);   % x=ordinates, y=data, col vecs

meth = 'spline'; % 'pchip'      % no different to 1e-8
% pad with 1 on left and 0 on right of domain...
A = 0*r;
ii = find(r<r0);
A(ii) = 1;
ii = find(r>=r0 & r<=r1);
if sum(ii)>0
  A(ii) = interp1(x,y,r(ii),meth);        % interpolate off given samples to r
end
  
%%%%%%
function test_eval_apod_NI2
eps = 1e-10;    % negligible but needed to avoid hitting 0 or 1 regions
[~,r0,r1] = eval_apod_NI2(0);  % get r0,r1
A = @(r) eval_apod_NI2(r);    % func
fprintf('1-A(%.5g)=%.3g, A(%.5g)=%.3g\n',r0,1-A(r0+eps),r1,A(r1-eps))  % check decay
figure; r = linspace(0,r1,1e4); plot(r,A(r),'-');
% note every time A(r) used, the file is loaded :(
title('eval\_apod\_NI2 test'); xlabel('r'); ylabel('A(r)');
