% test BDWF with generic targets on starshades as in ../demo_starshades.m
% Compares results to fresnaq, and NSLI.
% Barnett 9/10/20
clear all; verb = 1;

lambda = 4.27e-7;      % generic wavelength, blue (meters)
Z = 2.5e7;             % downstream distance (meters) ... good for Reff~10m

% we copy the starshade profile shapes (for now) from ../demo_starshades.m:
design = 'NI2';   % choose design from below list...
switch design 
 case 'disc'
  Np = 1; r0=5; r1=10; Afunc = @(r) 1+0*r; % 10m radius disc: Poisson spot u=1!
  n = 600; m = 80;                  % Afunc=1, r0, dummy hacks. Good to lz>=5 
 case 'lin'
  Np = 24; r0=5; r1=13; Afunc = @(r) (r1-r)/(r1-r0);   % simple linear apod.
  n = 600; m = 80;        % good to lambdaz>=9; convergence must be tested
 case 'HG'
  A = @(t) exp(-(t/0.6).^6);        % Cash'11 hyper-Gaussian on [0,1]
  Np = 16; r0 = 7; r1 = 14;         % # petals, inner, outer radii in meters
  Afunc = @(r) A((r-r0)/(r1-r0));   % apodization vs radius in meters
  n = 30; m = 80;
 case 'erf'                         % my analytic toy model, similar r0,r1
  Np=16; r0 = 7; r1 = 14;           % # petals, inner, outer radii in meters
  beta = 3.0;                       % good for A or 1-A to decay to 1e-5
  Afunc = @(r) erfc(beta*(2*r-(r0+r1))/(r1-r0))/2;
  %e=1e-3; Afunc = @(r) e + (1-e)*Afunc(r); % fatten the tips: get u(0,0)~e
  %e=1e-2; Afunc = @(r) e*(r1-r)/(r1-r0) + (1-e)*Afunc(r);  % sharp linear tips
  %e=1e-3; Afunc = @(r) (1-e)*Afunc(r); % fatten the gaps: get u(0,0)~e
  n = 30; m = 100;       % good to lambdaz>=9; convergence must be tested
 case 'NI2'                         % actual NI2 starshade, cubic interpolated
  cwd = fileparts(mfilename('fullpath')); file = [cwd '/../occulter/NI2'];
  Np = 24;                          % petals told cut off harshly at r1 ...bad
  [~,r0,r1] = eval_sister_apod(file,0);   % get apodization range [r0,r1]
  Afunc = @(r) eval_sister_apod(file,r);  % func handle (reads file when called)
  % try destroying the 1e-2 Poisson spot due to both petal flat ends...
  %Afunc = @(r) Afunc(r).*((r<=12) + (r>12).*cos(pi/2*(r-12)).^2); % C^1 blend
  %Afunc = @(r) 1 - (1-Afunc(r)).*((r>=6) + (r<6).*cos(pi/2*(6-r)).^2);  % "
  n = 40; m = 400;    % use quad_conv_apod_NI2 converged m (to 1e-6), lz>=5
 case 'NW2'                         % actual NI2 starshade, cubic interpolated
  file = '/home/alex/physics/starshade/SISTER/input_scenes/locus/in/TV3';
  Np = 24;                          % petals told cut off harshly at r1 ...bad
  [~,r0,r1] = eval_sister_apod(file,0);   % get apodization range [r0,r1]
  Afunc = @(r) eval_sister_apod(file,r);  % func handle (reads file when called)
  % try destroying the 1e-2 Poisson spot due to both petal flat ends...
  %Afunc = @(r) Afunc(r).*((r<=12) + (r>12).*cos(pi/2*(r-12)).^2); % C^1 blend
  %Afunc = @(r) 1 - (1-Afunc(r)).*((r>=6) + (r<6).*cos(pi/2*(6-r)).^2);  % "
  n = 50; m = 700;    % use quad_conv_apod_NI2 converged m (to 1e-6), lz>=5
end
fprintf('apod ends: 1-A(%.5g)=%.3g, A(%.5g)=%.3g\n',r0,1-Afunc(r0),r1,Afunc(r1))
Reff = (r0+r1)/2;
fprintf('lambdaZ = %.3g m^2.  approx Fresnel # = %.3g\n', lambda*Z, Reff^2/(lambda*Z))

[xq yq wq bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb);   % fill areal quadr

ximax = 1.0; ngrid = 20;       % ngrid^2 centered target grid out to +-ximax
tol = 1e-9;                    % fresnaq accuracy
[u xigrid] = fresnaq_grid(xq, yq, wq, lambda*Z, ximax, ngrid, tol, verb);
u = 1-u;                       % convert aperture to occulter
u = u * exp(2i*pi*Z/lambda);   % plane z-propagation (since fresnaq doesn't),
                               % since Z/lam~1e14 this ruins the overall phase!

it = ceil(ngrid/2+1); jt = it;      % indices of center xi=eta=0
ut = u(it,jt); xi = xigrid(it); eta = xigrid(jt);     % math check u
fprintf('|u|^2 intensity at (%.3g,%.3g) is %.6g\n',xi,eta,abs(u(it,jt)).^2)

nO = ngrid; dxO = 2*ximax/nO;  % convert our grid to bdwf grid params...
deltaX = 0; if mod(ngrid,2)==0, deltaX = +dxO/2; end   % deltaX is backwards!
deltaY = deltaX;
[xx yy] = make_grid_bdwf(dxO, nO, deltaX, deltaY);     % ugh
fprintf('grid match: %.3g\n',norm(xx(:,1)-xigrid(:)))

% compare to bdwf's and nsli's answer (w/ psi1=0) using raw & resampled locus:
if strcmp(design,'NI2'), o=load(file);
  tic
  ub = bdwf(o.xVals,o.yVals,[], Z, lambda, dxO, nO, 0,0, deltaX, deltaY);
  toc
  fprintf('|ub|^2 intensity at (0,0) is %.6g\n',abs(ub(it,jt)).^2)
  fprintf('max |ub| error on grid : %.3g\n',norm(abs(ub(:)) - abs(u(:)),inf))
  % note because of overall phase garbage, have to compare magnitudes.
  figure; subplot(2,3,1); imagesc(xigrid,xigrid,abs(u)'); colorbar
  title('u fresnaq (areal from resampled locus)'); axis xy equal tight
  subplot(2,3,2); imagesc(xigrid,xigrid,abs(abs(u)-abs(ub))'); colorbar
  title('abs diff |u|: fresnaq vs bdwf (raw locus)'); axis xy equal tight
  % even though m=400 enough for areal quadr to 1e-6, get 1e-4 BDWF, 3e4 NSLI
  % even m=1000 gets 7e-5 BDWF, 5e-5 NSLI. m=3000 needed to get sim to raw locus
  m = 2000;    % recall original was sampled at only 2462 pts in r=[5,13]
  dummyn=2; [xq yq wq bx by] = starshadequad(Np,Afunc,r0,r1,dummyn,m,verb,'u');
  tic
  ubr = bdwf(bx,by,[], Z, lambda, dxO, nO, 0,0, deltaX, deltaY);  % 1e-4, bad
  toc
  subplot(2,3,3); imagesc(xigrid,xigrid,abs(abs(u)-abs(ubr))'); colorbar
  title('abs diff |u|: fresnaq vs bdwf (resamp locus)'); axis xy equal tight

  % locus resampling seems nicely resolved rel to fresnel zones...
  subplot(2,3,4); plot(bx,by,'.'); axis equal; overlay_zones(lambda*Z,0,0,'r-');
  title(sprintf('bdry nodes (locus) from resamp A(r): m=%d',m));
  
  % compare nsli answer (w/ psi1=0) using raw or resamp locus samples...
  xi = xx(:); eta = yy(:);
  [wx wy] = crudecurvequad(o.xVals,o.yVals);          % weights for raw locus
  tic
  un = nsli_pts(o.xVals,o.yVals,wx,wy, lambda*Z, xi, eta);
  toc
  un = 1-un;     % turn aperture into occulter
  un = reshape(un,[ngrid ngrid]);
  fprintf('|un|^2 intensity at (0,0) is %.6g\n',abs(un(it,jt)).^2)
  subplot(2,3,5); imagesc(xigrid,xigrid,abs(abs(u)-abs(un))'); colorbar
  title('abs diff |u|: fresnaq vs nsli (raw locus)'); axis xy equal tight
  [wx wy] = crudecurvequad(bx,by);
  tic;
  unr = nsli_pts(bx,by,wx,wy, lambda*Z, xi, eta);      % resampled bdry
  toc
  unr = 1-unr;     % turn aperture into occulter
  unr = reshape(unr,[ngrid ngrid]);
  fprintf('|unr|^2 intensity at (0,0) is %.6g\n',abs(unr(it,jt)).^2)
  subplot(2,3,6); imagesc(xigrid,xigrid,abs(abs(u)-abs(unr))'); colorbar
  title('abs diff |u|: fresnaq vs nsli (resamp locus)'); axis xy equal tight
  % looks like NSLI is faster than BDWF. resampling affects (low-order) acc of
  % both, as expect, but BDWF suffers more with resampling.
end
