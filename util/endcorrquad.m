function [x w] = endcorrquad(N,a,b,k)
% ENDCORRQUAD  Alpert low-to-medium order end-corrected rule on regular nodes
%
% [x w] = endcorrquad(N,a,b) returns arithmetic sequence of nodes x, and weights
%  w to match, for N-node quadrature of a smooth function on the 1D interval
%  [a,b].  Both x=linspace(a,b,N) and w are returned as column vectors.
%
% [x w] = endcorrquad(N,a,b,k) uses kth-order scheme: k=2,4,6,8,10 avail.
%
% without arguments, does convergence test for various orders

% Barnett 9/7/20; coeffs from Table 3.5 of YALEU/DCS/RR-814 by B. Alpert, 1990.
if nargin==0, test_endcorrquad; return; end
if nargin<4, k=4; end

switch k
 case 2
  d = 1; c = [0];
 case 4
  d=24; c = [-3,4,-1];    % h^4 corrs as in Table 3.5 of Alpert Tech Rep, p.43.
 case 6
  d=1440; c = [-245 462 -336 146 -27];   % etc
 case 8
  d=120960; c = [-23681, 55688, -66109, 57024, -31523, 9976, -1375];
 case 10
  d=7257600; c = [-1546047, 4274870, -6996434, 9005886, -8277760, 5232322, -2161710, 526154, -57281];
 otherwise
  error('unknown order k');
end
nc = k-1;          % nodes to correct per end

x = linspace(a,b,N); x = x(:);
h = (b-a)/(N-1);
w = [0.5;ones(N-2,1);0.5];   % h^2 trap rule, w/o h prefactor, col vec

if N<2*nc, warning('N must be at least 2(k-1); not correcting TR');
  w = w*h; return; end
w(1:nc) = w(1:nc) + c'/d;
w(N-nc+1:N) = w(N-nc+1:N) + c(end:-1:1)'/d;
w = w*h;

%%%%
function test_endcorrquad
W = 30;    % freq param: roughly # half-wavelengths
f = @(x) sin(W*x); F = @(x) -cos(W*x)/W;    % f and analytic antideriv
a = 4; b = 7;
Iex = F(b)-F(a);
Ns = ceil(logspace(1.5,3,16));   % note exponent is 10
figure;
for k= [2 4 6 8 10 inf]         % last case compares to G-L, sadly much better
  errs = nan*Ns;
  for i=1:numel(Ns), N=Ns(i);
    if isfinite(k), [x w] = endcorrquad(N,a,b,k); else, [x w] = lgwt(N,a,b); end
    I=sum(f(x).*w);
    errs(i) = I-Iex;
  end
  loglog(Ns,abs(errs),'+-'); hold on; plot(Ns,(3/W)*(Ns/W).^-k,'--');  % model
end
xlabel('N'); ylabel('quadr error'); axis tight;
title('endcorrquad: Alpert unif grid quadr err test vs Gauss');
