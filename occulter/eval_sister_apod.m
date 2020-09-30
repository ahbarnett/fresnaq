function [A r0 r1 Ap] = eval_sister_apod(filehead,r)
% EVAL_SISTER_APOD  Interpolate SISTER starshade apodization profile A(r), A'(r)
%
% [A r0 r1] = eval_sister_apod(filehead,r) reads the file in SISTER format
%  giving the sampled profile, then uses interpolation to approximate A(r) at
%  the given list of r values.
%  The apodization radius range [r0,r1] is also extracted from other fields.
%
% [A r0 r1 Ap] = eval_sister_apod(filehead,r) also returns A'(r) values.
%
% If r is empty or absent, function handles are returned for A and Ap instead.
%
% With no arguments, does self-test and plot.
%
%  Notes:
%  1) apodization function A(r) is in notation of E. Cady, Opt. Expr. 2012,
%     starting at 1 for r<r0 and becoming 0 for r>r1.
%  2) reads a multi-MB file so don't call too often!
%  3) cubic spline interp is used.

% Barnett 9/9/20, deriv & func handle out 9/28/20
if nargin==0, test_eval_sister_apod; return; end
if nargin<2, r = []; end
wantAp = nargout>3;

o = load(filehead);
r1 = o.occulterDiameter/2;      % upper apodization radius in meters
r0 = r1 - o.petalLength;        % lower apodization radius in meters

inds = (o.r>=r0);               % kill off useless inner radii samples
x = o.r(inds); x = x(:); y = o.Profile(inds);   % x=ordinates, y=data, col vecs

meth = 'spline'; % 'pchip'    % no different, to 1e-8 with their sample dr=2mm
if isempty(r)                   % if func handles ok
  A = @(r) Aintpad(r0,r1,x,y,r,meth);      % turn eval meth into handle
  if wantAp
    Ap = @(r) Apintpad(r0,r1,x,y,r,meth);  % "
  end
else
  A = Aintpad(r0,r1,x,y,r,meth);       % just evaluate on given r
  if wantAp
    Ap = Apintpad(r0,r1,x,y,r,meth);   % "
  end
end

%%%%%%%% helpers which do the eval, so can easily be wrapped as funcs...
function A = Aintpad(r0,r1,x,y,r,meth)  % A interp vals, padded outside [r0,r1]
A = 0*r;
A(r<r0) = 1;                    % pad with 1 on left and 0 on right of domain
ii = find(r>=r0 & r<=r1);
if sum(ii)>0
  A(ii) = interp1(x,y,r(ii),meth);        % interpolate off given samples to r
end
function Ap = Apintpad(r0,r1,x,y,r,meth)  % A' interp vals, pad outside [r0,r1]
Ap = 0*r;
ii = find(r>=r0 & r<=r1);
if sum(ii)>0
  h = 1e-5;                   % centered-diff estimates (expect O(h^2) err)
  Ap(ii) = (interp1(x,y,r(ii)+h,meth) - interp1(x,y,r(ii)-h,meth)) / (2*h);
end  



%%%%%%
function test_eval_sister_apod
cwd = fileparts(mfilename('fullpath')); filehead = [cwd '/NI2'];  % in repo

% test handles out...
[A,r0,r1,Ap] = eval_sister_apod(filehead);
eps = 1e-10;    % negligible but needed to avoid hitting 0 or 1 regions
a = r0+eps; b = r1-eps;
fprintf('1-A(%.5g)=%.3g, A(%.5g)=%.3g\n',r0,1-A(a),r1,A(b))  % check decay
[z w] = lgwt(500,a,b);
%[z, A(z), Ap(z)], figure; plot(z,[A(z),Ap(z)],'+-'); return
fprintf('G-L A'' FTC test: %.3g\n',sum(Ap(z).*w) + A(a)-A(b))

% test values out...
r = linspace(0,r1,1e4);
[Ar,~,~,Apr] = eval_sister_apod(filehead,r);
fprintf('crude A'' FTC test: %.3g\n',sum(Apr*(r(2)-r(1)))+A(a)-A(b))
% (since A'~0 at ends, just expect O(h^2) err due to A'' in BV class.)

figure; plot(r,[A(r);Ar;Ap(r);Apr],'-'); axis tight;       % compare
hold on; o = load(filehead); plot(o.r,o.Profile,'r.');     % original data
legend('A(r)','A vals','A''(r)','A'' vals','data');
title('eval\_sister\_apod NI2 test'); xlabel('r (meters)'); vline([r0 r1]);


