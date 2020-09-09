% tester for smooth end corrections of Alpert, 1999, SISC, in Quad*.m
% Note that the number of nodes produced may be up to 4 larger than requested.
% Also see endcorrquad.m
% Barnett 9/8/20.

W = 30;    % freq param: roughly # half-wavelengths
f = @(x) sin(W*x); F = @(x) -cos(W*x)/W;    % f and analytic antideriv
a = 4; b = 7;
Iex = F(b)-F(a);
Ns = ceil(logspace(1.5,3,16));   % note exponent is 10
figure;
ks = [3:8, 12:4:32, inf];    % last case compares to G-L
for k=ks
  errs = nan*Ns;
  for i=1:numel(Ns), N=Ns(i);
    if isfinite(k), [x w] = QuadNodesInterval(a,b,N,[],1,1,k); else, [x w] = lgwt(N,a,b); end
    I=sum(f(x).*w);
    errs(i) = I-Iex;
  end
  coi = get(gca,'ColorOrderIndex'); loglog(Ns,abs(errs),'+-'); hold on;
  set(gca,'ColorOrderIndex',coi);  % reuse color
  plot(Ns,(3/W)*(Ns/W).^-k,'--');  % model
end
xlabel('N'); ylabel('quadr error'); axis tight; set(gca,'ylim',[1e-15 1]);
title('Alpert graded aux end-corr quadr err test vs Gauss');
orders = kron(ks,[1 1]); orders = orders(1:end-1); legnum(orders);


