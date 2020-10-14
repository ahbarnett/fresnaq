% figure validating FRESNAQ vs BDWF etc, in starshades, incl shadow region.
% Barnett 9/28/20
clear
addpath ~/numerics/finufft/matlab
fdir = '/home/alex/physics/starshade/fresnaq';
addpath(fdir)
startup

design = 'HG';   % choose design from below list...
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
  n = 30; m = 60;                     % keep m low so LI meths have some err
  ximax = 32;
  xilab = 12;        % hline label location
  lambda = 5e-7;
  anal = true;
 case 'NI2'                         % actual NI2 starshade, cubic interpolated
  file = [fdir '/occulter/NI2'];
  n = 40; m = 400;    % use quad_conv_apod_NI2 converged m (to 1e-6), lam.z>=5
  ximax = 14;          % what to grid
  xilab = 5;        % hline label location
  lambda = 5e-7; % near longest wavelength
  anal = false;
end
if ~exist('lambdaint')              % then get usage params from SISTER file
  [Afunc,r0,r1,Apfunc] = eval_sister_apod(file);  % func handle
  o = load(file);
  %r1 = max(o.r);     % override? no, off by 2mm, messes up shadow.
  Z = o.Z; Np = o.numPetals;
  lambdaint = [o.lowerScienceWavelength o.upperScienceWavelength];
  i=find(o.Profile<0.5); Reff = o.r(i(1));      % Reff is where A=0.5, for now
end
FN = Reff^2./(lambdaint*Z);
mas = pi/180/60^2/1e3;            % one milliarcsecond
fprintf('design %s: Reff=%.3g m, Z=%.0f km, Fres#=[%.2g,%.2g], geo(eff)IWA=%.3g(%.3g) mas\n',design,Reff,Z/1e3,min(FN),max(FN),(r1/Z)/mas,(Reff/Z)/mas)
fprintf('Fres# using R and lambda = %.3g\n',r1^2/lambda/Z)
fprintf('\tapod: rel gap wid 1-A(%.5g)=%.3g, rel tip wid A(%.5g)=%.3g\n',r0,1-Afunc(r0),r1,Afunc(r1))

ngrid = 1e3;         % grid for fresnaq t1 & line will be extracted
verb = 1;
% fresnaq method
[xq yq wq bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb,'g');   % areal quadr
tol=1e-8;
tic
[u xigrid] = fresnaq_grid(xq,yq,wq, lambda*Z, ximax, ngrid, tol);
fprintf('fresnaq_grid %.3g s\n',toc)
u = 1-u;                       % convert aperture to occulter
xi = xigrid(ngrid/2+1:end); eta=0*xi;  % targs to non-neg line along (xi,0)
inds = ngrid^2/2 + ngrid/2 + (1:ngrid/2);  % corresp u grid vals
u = u(inds); u = u(:);

% LI meths
if ~anal   % NI2
  Bx=o.xVals(1:end-1)'; By=o.yVals(1:end-1)';  % raw locus
  %[wx wy] = crudecurvequad(bx, by);  % is much worse than midpt for raw NI2
  [bx by wx wy] = midpointcurvequad(Bx,By);   % low-ord midpt rule as BDWF has
else             % build our own LI quad, eg high-order
  [Bx By wx wy] = starshadeliquad(Np,Afunc,Apfunc,r0,r1,m,4,verb,'g');
  bx=Bx; by=By;   % use high-ord
  % still fails to match BDWF convergence in shadow, due to ugeom! (see shadfix)
end
tic;
ub = bdwf_pts([Bx; Bx(1)],[By; By(1)], [], Z, lambda, xi, eta, 0,0);
t = toc; fprintf('bdwf_pts %.3g s (%.3g pt-targ pairs/s)\n',t,numel(inds)*numel(Bx)/t)
tic;
shadfix = 0;    % skip, too complicated for now
un = nsli_pts(bx,by,wx,wy, lambda*Z, xi, eta, shadfix);  % uses high-ord quadr
t = toc; fprintf('nsli_pts %.3g s (%.3g pt-targ pairs/s)\n',t,numel(inds)*numel(bx)/t)
un = 1-un(:);                     % convert aperture to occulter
ifit = ngrid/2;                  % where to match phase (intermediate r better?)
ub = ub * (u(ifit)/ub(ifit)) * abs(ub(ifit)/u(ifit));  % fit overall bdwf phase
% (since overall phase shift meaningless: Z/lambda ~ 1e14)
nbdry = numel(bx);

% timing table tests (HG only; NI2 is estim)
%[xi eta] = ndgrid(xigrid,xigrid); xi=xi(:); eta=eta(:);
%tic;
%ub = bdwf_pts([Bx; Bx(1)],[By; By(1)], [], Z, lambda, xi, eta, 0,0);
%t = toc; fprintf('bdwf_pts %.3g s (%.3g pt-targ pairs/s)\n',t,numel(xi)*numel(Bx)/t)


figure(1); clf;   % ---------------------------------------------------------
semilogy(xi,abs(ub), 'r-','linewidth',1.2); hold on;
plot(xi,abs(un),'g-','linewidth',0.6);
plot(xi,abs(u),'k-', 'linewidth',0.1);    % try to let others show behind!
plot(xi,abs(u-ub),':','color',[.5 0 0],'markersize',1);
plot(xi,abs(u-un),'--','color',[0 .5 0],'markersize',1);
if strcmp(design,'NI2')
  h=legend('BDWF','NSLI',sprintf('t1 ($m{=}%d$)',m),'t1-BDWF','t1-NSLI');
else
  m = 600;      % oversamp to demo poor shadow perf cf bdwf
  [Bx By wx wy] = starshadeliquad(Np,Afunc,Apfunc,r0,r1,m,4,verb,'g');
  [bx by wx wy] = midpointcurvequad(Bx, By);   % low-ord midpt rule as BDWF
  unlo = nsli_pts(bx,by,wx,wy, lambda*Z, xi, eta, shadfix);
  unlo = 1-unlo(:);
  plot(xi,abs(u-unlo),'-.','color',[0 1 1],'markersize',1);
  h=legend('BDWF','NSLI','t1','t1-BDWF','t1-NSLI','t1-NSLIlo');
end
set(h,'location','northwest','interpreter','latex');
h.Position = h.Position + [0.04 -0.05 0 0];    % doesn't match screen!
if 0, u =abs(u); un=abs(un); ub = abs(ub);  % now discard phases, incl errors...
  plot(xi,abs(u-ub),'*:','color',[.5 0 0]);
  plot(xi,abs(u-un),'*:','color',[0 .5 0]); end
hline(1e-5,'k:'); text(xilab,1.7e-5,'intensity suppression $10^{-10}$','interpreter','latex');
hline(1,'k:'); text(xilab,0.6,'incident intensity $1$','interpreter','latex');
xlabel('$\xi$ (meters)','interpreter','latex'); ylabel('$|u^{\rm oc}(\xi,0)|$, or abs diff in $u^{\rm oc}(\xi,0)$','interpreter','latex');
%axis([0 3 3e-8 1.5e-4])
axis([0 ximax 3e-8 1.5])
title(sprintf('%s ($n=%d$ bdry nodes): $\\lambda=$ %.3g m, $Z=$ %.3g m',design,nbdry,lambda,Z),'interpreter','latex')
set(gcf,'paperposition',[0 0 5 4]);
print('-depsc2',['../valid_' design '.eps'])

% now add arrows w/ xfig



