% FRESNAQ for Fresnel diffraction from various starshades: shadow vs lambda.
% Ie, max_phi |u(rho,phi)|^2 vs rho and lambda, over design range & beyond.
% Barnett 9/14/20
clear; verb = 1;

% targets: polar grid (roh,phi) in Vanderbei notation; radii (m) and angles...
Nr=80; rho = linspace(0,4,Nr); Np = 1e2; phi = 2*pi*(1:Np)/Np;
[rr pp] = ndgrid(rho,phi); xi = rr.*cos(pp); eta = rr.*sin(pp);  % (xi,eta)'s

cwd = fileparts(mfilename('fullpath'));

design = 'erf';   % choose design from below list...
switch design
 case 'HG'
  A = @(t) exp(-(t/0.6).^6);        % a Cash'11 "hyper-Gaussian" on [0,1]
  Np = 16; r0 = 7; r1 = 14;         % # petals, inner, outer radii in meters
  Reff = (r0+r1)/2;
  Afunc = @(r) A((r-r0)/(r1-r0));   % apodization vs radius in meters
  Z = 3e7; lambdaint = [4e-7 1.1e-6];
  n = 30; m = 80;
 case 'erf'                         % my analytic toy model, similar r0,r1
  Np=16; r0 = 7; r1 = 14;           % # petals, inner, outer radii in meters
  beta = 3.0;                       % good for A or 1-A to decay to 1e-5
  Afunc = @(r) erfc(beta*(2*r-(r0+r1))/(r1-r0))/2;
  Reff = (r0+r1)/2;
  Z = 3e7; lambdaint = [4e-7 1.1e-6];
  n = 30; m = 100;       % good to lambdaz>=9; convergence must be tested
 case 'NI2'                         % actual NI2 starshade, cubic interpolated
  file = [cwd '/occulter/NI2'];
  n = 40; m = 300;    % use quad_conv_apod_NI2 converged m (to 1e-6), lam.z>=5
 case 'TV3'
  file = [cwd '/../SISTER/input_scenes/locus/in/TV3'];  % or your SISTER dir
  n = 40; m = 800;                  % converged m (to 1e-6), lam.z>23
 case 'NW2'
  file = [cwd '/../SISTER/input_scenes/locus/in/NW2'];
  n = 40; m = 800;                  % 
end
if ~exist('lambdaint')              % get usage params from SISTER file
  o = load(file);
  [~,r0,r1] = eval_sister_apod(file,0);   % get apodization range [r0,r1]
  Afunc = @(r) eval_sister_apod(file,r);  % func handle (reads file when called)
  Z = o.Z; Np = o.numPetals;
  lambdaint = [o.lowerScienceWavelength o.upperScienceWavelength];
  i=find(o.Profile<0.5); Reff = o.r(i(1));      % Reff is where A=0.5, for now
end
FN = Reff^2./(lambdaint*Z);
mas = pi/180/60^2/1e3;            % one milliarcsecond
fprintf('design %s: Reff=%.3g m, Z=%.0f km, Fres#=[%.2g,%.2g], geoIWA=%.3g mas\n',design,Reff,Z/1e3,min(FN),max(FN),(r1/Z)/mas)
fprintf('\tapod: rel gap wid 1-A(%.5g)=%.3g, rel tip wid A(%.5g)=%.3g\n',r0,1-Afunc(r0),r1,Afunc(r1))

[xq yq wq bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb);   % fill areal quadr

% wavelength range to explore...
Nl = 30; lambda = logspace(log10(lambdaint(1)*0.7),log10(lambdaint(2)*1.4),Nl);
tol = 1e-6;
tic; maxu2 = nan(Nr,Nl);   % alloc output, and do the calc...
for l=1:Nl
  u = fresnaq_pts(xq,yq,wq, lambda(l)*Z, xi,eta, tol);
  maxu2(:,l) = max(abs(1-u).^2,[],2);   % Babinet u-1 for occulter, max over phi
end
toc
figure(1); clf; Iplt=maxu2; sets = {'ring','disc'};
for p=1:2, if p==2, Iplt = cummax(maxu2,1); end   % 2nd time max over disc
  subplot(2,1,p); [C,h] = contourf(lambda*1e6,rho,log10(Iplt),[-20, -11:-3]);
  axis xy tight; ylabel('radius \rho (m)'); xlabel('\lambda (\mu m)');
  colormap(jet); colorbar; caxis([-11 -3]); vline(lambdaint*1e6); clabel(C,h);
  title(sprintf('%s: log_{10} (max |u|^2 over %s) vs radius and lambda',design,sets{p}))
  print('-dpng',['occulter/shadow_rholambda_' design '.png']);
end
