function [Ar r0 r1] = eval_apod_erf(r)
% EVAL_APOD_ERF  simple erf-based analytic starshade apodization profile A(r)
%
% [Ar r0 r1] = eval_apod_erf(r)
%  evaluates A(r) at the given list of r values, for simple starshade design
%  of analytic form. The apodization radius domain [r0,r1] is also returned.
%
% Without arguments, does self-test and plot.
%
% Simple design based on blending function; Barnett 8/24/20

r0 = 7; r1 = 14;
beta = 3.0;       % width param, eg 3.0 for 1e-5 decay.
a = (r0+r1)/2;    % 1/2-way radius (Cash '11 notation)
A = @(r) (r<r0) + (r<=r1 & r>=r0).*erfc(2*beta*(r-a)/(r1-r0))/2;  % pad by 1,0
%A = @(t) exp(-(t/0.6).^6);        % Cash'11 hyper-Gaussian in [0,1], to rescale

if nargin>0
  Ar = A(r);      % do it  
else              % test
  eps = 1e-10;    % negligible but needed to avoid hitting 0 or 1 regions
  fprintf('1-A(%.5g)=%.3g, A(%.5g)=%.3g\n',r0,1-A(r0+eps),r1,A(r1-eps))  % check decay
  figure; r = linspace(0,r1,1e4); plot(r,A(r),'-');
  title('eval\_apod\_erf test'); xlabel('r'); ylabel('A(r)');
end
