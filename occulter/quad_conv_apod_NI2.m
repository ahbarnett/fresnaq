% Test convergence of quadrature on interpolated apodization profile from NI2.
% This makes sure we can interpolate from the NI2 A(r) samples on regular
% grid and use our high-order quadrature on it, at least in 1D.
% We also warm up with exploring the underlying smoothness of A: is only C^1.
% Barnett 9/4/20

clear
cwd = fileparts(mfilename('fullpath'));
o = load([cwd '/NI2.mat']);
r1 = o.occulterDiameter/2;      % upper apodization radius in meters
r0 = r1 - o.petalLength;        % lower apodization radius in meters

inds = (o.r>=r0);               % kill off useless inner radii
x = o.r(inds); x = x(:); y = o.Profile(inds);   % x=ordinates, y=data, col vecs

figure; plot(x, [y, 1e5*[0;diff(y,2);0]],'.','markersize',10); xlabel('r');
legend('samples','1e5 * D^2(samples)');  % A'' is piecewise const!
title('strangenss of sampled profile A(r) and its 2nd discrete deriv');

% first check Fourier decay for this sampled func...
% (uncomment next 2 lines to replace by samples of generic model C^1 func)
%y = kron((-1).^(1:20)',ones(100,1))/1e5; y = [0; cumsum(cumsum(y(1:end-1)))];
%x=(1:numel(y))'; %y = y + 0.1*cos(17*pi*x/numel(y));   % non-poly term
yp = [y; y(end:-1:1)]; Ahat = ifft(yp);    % flip back & DFT to check smoothness
%diff(yp,2), diff([yp(end-10:end); yp(1:10)],2)  % check cleanly C^1
i=1:numel(y)+1; Ahata = abs(Ahat(i));       % Fourier mag data to plot
figure;
subplot(1,2,1); semilogy(i, Ahata, '+-'); set(gca,'ylim',[1e-10 1]); title('A samples: DFT decay');
subplot(1,2,2); loglog(i, Ahata, '.'); axis tight; set(gca,'ylim',[1e-10 1]); title('decay and n^{-3}');
hold on; plot(i,i.^-3,'r-');  % seems to be 1/n^3, ie A is only C^1 (surprise!)

% We also notice 250 Fourier modes (125 PTR nodes over transition region) would
% be enough to capture A to 1e-6 absolute error. This would be a better
% quadr scheme than G-L by factor pi/2, asymptotically. Alpert, below, seems
% to do as predicted.

Ns = ceil(2.^(4:0.2:10));                 % numbers of quadr pts to test
Is = 0*Ns;                       % values of integrals
meth = 'spline'; % 'pchip'       % no different to 1e-8
for i=1:numel(Ns)
  %[xj wj] = lgwt(Ns(i),r0,r1);   % pick one of the this or next two quadr's...
  %ord = 10; [xj wj] = endcorrquad(Ns(i),r0,r1,ord);  % reg grid no help for osc
  ord = 32; [xj wj] = QuadNodesInterval(r0,r1,Ns(i),[],1,1,ord);  % Alpert grid
  yj = interp1(x,y,xj,meth);     % interpolate off given samples to xj
  yj = yj .* sin(10*xj);         % throw in oscillation to check
  Is(i) = sum(yj.*wj);
end
disp([Ns(:), Is(:)])
ii = 1:numel(Ns)-1;
figure; subplot(2,1,1); plot(xj,yj,'.-');
subplot(2,1,2);
loglog(Ns(ii),abs(Is(ii)-Is(end)),'+-'); title('self-convergence of I');
axis tight; hold on; plot(Ns, 20*Ns.^-3,'r-'); legend('err','N^{-3}');

% for either Gauss-L or Alpert (latter slightly better):
% N=300 enough for 1e-6 at freq 0, or 10, say.
% N=400 enough for 1e-6 at freq <=100.



