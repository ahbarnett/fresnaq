function [fun p res] = globalinterp1d(x,y,a,b,tol,fixedp)
% GLOBALINTERP1D  high-order global 1D interpolant func from point samples
%
%  fun = globalinterp1d(x,y,a,b,tol) returns a function handle to an interpolant
%   f(x) to the samples (x_i,y_i) for i=1..numel(x) on the interval [a,b].
%   The estimated accuracy is tol. Crude Global Chebyshev expansion is used,
%   not very fast, but fine for small problems.
%
%  [fun p res] = ... also returns interpolant order used and max residual
%
%  .. = globalinterp1d(x,y,a,b,tol,fixedp) fixes a single p to use.
%
%  Without arguments does a self test.

% Barnett 9/4/20

if nargin==0, test_globalinterp1d; return; end
ysize = max(abs(y));
x = (2*x(:)-(a+b))/(b-a);     % rescale to [-1,1], col vec
y = y(:);                     % data col vec
if nargin<6                   % adaptive p
  p = 4;                      % half the starting # Cheby poly's
  res = inf;
  pmax = 1e4;                 % quit
  while res>tol*ysize & p<pmax
    p = 2*p;                  % refine
    V = cos(acos(x)*(0:p-1)); % outer prod for Vandermonde matrix for Cheby's
    c = V\y;                  % col vec of coeffs
    res = norm(V*c-y,inf);
  end
else
  p = fixedp;
  V = cos(acos(x)*(0:p-1));   % outer prod for Vandermonde matrix for Cheby's
  c = V\y;                    % col vec of coeffs
  res = norm(V*c-y,inf);
end
% set up rescale then eval Cheby...
fun = @(t) reshape(cos(acos((2*t(:)-(a+b))/(b-a))*(0:p-1))*c, size(t));

%%%%%%%
function test_globalinterp1d
f = @(x) sin(10*x);
a = 4; b = 12;
N = 1e3;                         % since ordinates random, need enough of them
x = sort(rand(1,N))*(b-a) + a;   % samples
y = f(x);                        % data
tol = 1e-10;
[fun p res] = globalinterp1d(x,y,a,b,tol);
t = sort(rand(1,N))*(b-a) + a;
%[f(t(:)), fun(t(:))]
fprintf('test_globalinterp1d: residual=%.3g, degree p=%d, test sup err = %.3g\n',res,p,norm(f(t)-fun(t),inf))
figure; plot(x,y,'+'); hold on; plot(t,fun(t),'.-');
