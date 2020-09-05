% Test convergence of quadrature on interpolated apodization profile from NI2.
% This makes sure we can interpolate from the NI2 A(r) samples on regular
% grid and use our high-order quadrature on it, at least in 1D.
% Barnett 9/4/20

cwd = fileparts(mfilename('fullpath'));
o = load([cwd '/NI2.mat']);
r1 = o.occulterDiameter/2;      % upper apodization radius in meters
r0 = r1 - o.petalLength;        % lower apodization radius in meters

inds = (r>=r0);                 % kill off useless inner radii
x = o.r(inds); x = x(:); y = o.Profile(inds);   % x=ordinates, y=data, col vecs

Ns = ceil(2.^(4:0.25:10));                 % numbers of quadr pts to test
Is = 0*Ns;                       % values of integrals
meth = 'spline'; % 'pchip'       % no different to 1e-8
for i=1:numel(Ns)
  [xj wj] = lgwt(Ns(i),r0,r1);
  yj = interp1(x,y,xj,meth);     % interpolate off given samples to xj
  yj = yj .* sin(10*xj);        % throw in oscillation to check
  Is(i) = sum(yj.*wj);
end
disp([Ns(:), Is(:)])
ii = 1:numel(Ns)-1;
figure; subplot(2,1,1); plot(xj,yj,'.-');
subplot(2,1,2);
semilogy(Ns(ii),abs(Is(ii)-Is(end)),'+-'); title('self-convergence of I');
axis tight;
% N=300 enough for 1e-6 at freq 0.
% N=400 enough for 1e-6 at freq <=100.



