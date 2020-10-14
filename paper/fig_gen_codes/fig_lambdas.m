% figure of lambda sweep, max over theta as function of target radius.
% Barnett 9/30/20
clear
addpath ~/numerics/finufft/matlab
fdir = '/home/alex/physics/starshade/fresnaq';
addpath(fdir)
startup
verb = 1;

designs = {'NI2','HG'};  % choose design from list...(taken from fig_valid)
for d=1:numel(designs), design = designs{d}   % -------------------------
switch design
 case 'HG'
  Np = 16;         % appears from Fig. 4, but mentions Np=12 too.
  a = 12.5; b = a; pow=6;             % a,b in meters.  pow is his "n"
  Z=8e7;                              % 80000km (big and far, to get IWA)
  Afunc = @(r) exp(-((r-a)/b).^pow);  % "offset hyper-Gaussian" in r
  Apfunc = @(r) (-pow/b^pow)*(r-a).^(pow-1).*Afunc(r);    % A'
  % r=10; h=1e-5; Apfunc(r) - (Afunc(r+h)-Afunc(r-h))/(2*h) % test A'
  r0 = a; r1 = 31; Reff=a+b;          % r1 = outer max radius, meters.
  lambdaint = [3e-7 1e-6];            % listed p.2.
  lambdasho = [2e-7 1.1e-6];
  n = 30; m = 60;                     % keep m low so LI meths have some err
  Rmax = 6;                           % biggest rho to explore
 case 'NI2'                         % actual NI2 starshade, cubic interpolated
  file = [fdir '/occulter/NI2'];
  n = 40; m = 400;    % use quad_conv_apod_NI2 converged m (to 1e-6), lam.z>=5
  Rmax = 3;
  lambdasho = [3.7e-7 6.3e-7];
end
if ~exist('lambdaint')              % then get usage params from SISTER file
  [Afunc,r0,r1,Apfunc] = eval_sister_apod(file);  % func handle
  o = load(file);
  Z = o.Z; Np = o.numPetals;
  lambdaint = [o.lowerScienceWavelength o.upperScienceWavelength];
  i=find(o.Profile<0.5); Reff = o.r(i(1));      % Reff is where A=0.5, for now
end
FN = Reff^2./(lambdaint*Z);
mas = pi/180/60^2/1e3;            % one milliarcsecond
fprintf('design %s: Reff=%.3g m, Z=%.0f km, Fres#=[%.2g,%.2g], geo(eff)IWA=%.3g(%.3g) mas\n',design,Reff,Z/1e3,min(FN),max(FN),(r1/Z)/mas,(Reff/Z)/mas)

% targets: polar grid (rho,phi) in Vanderbei notation; radii (m) and angles...
Nr=200; rho = linspace(0,Rmax,Nr); Na = 3e2; phi = 2*pi*(1:Na)/Na;
[rr pp] = ndgrid(rho,phi); xi = rr.*cos(pp); eta = rr.*sin(pp);  % (xi,eta)'s

[xq yq wq bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb,'g'); % fill areal quadr

% wavelength range to explore...
Nl = 50; lambda = logspace(log10(lambdasho(1)),log10(lambdasho(2)),Nl);
verb = 0; tol = 1e-8;
tic; maxu2 = nan(Nr,Nl);   % alloc output, and do the calc...
for l=1:Nl
  u = fresnaq_pts(xq,yq,wq, lambda(l)*Z, xi,eta, tol, verb);
  maxu2(:,l) = max(abs(1-u).^2,[],2);   % Babinet u-1 for occulter, max over phi
end
toc

figure(d); clf;
[C,h] = contourf(lambda*1e6,rho,log10(maxu2),[-20, -11:-3]);
axis xy tight; ylabel('radius $\rho$ (m)','interpreter','latex');
xlabel('$\lambda$ ($\mu$m)','interpreter','latex');
colormap(jet);
if d==2, colorbar; end
caxis([-13 -4]); vline(lambdaint*1e6);
clabel(C,h);
title(sprintf('(%c)  %s ($N=%d$): max $\\log_{10} |u^{\\rm oc}|^2$ on target radius $\\rho$','a'+d-1,design,numel(xq)),'interpreter','latex');
%title(sprintf('%s: log_{10} (max |u|^2 over %s) vs radius and lambda',design,sets{p}))
set(gcf,'paperposition',[0 0 4+(d/2) 2.5]);
print('-depsc2',['../lambdas_' design '.eps'])

end  % -------------------------------------------------------------------