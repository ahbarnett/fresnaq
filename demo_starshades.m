% Demo of FRESNAQ computing Fresnel diffraction from various starshades.
% Barnett 9/8/20
clear; verb = 1;  % verbosity

lambdaz = 9.0;    % wavelength*dist, in m^2, recall Fres # = Reff^2/(lambda.z)
ximax = 14.0; ngrid = 1e3;     % ngrid^2 centered target grid out to +-ximax
tol = 1e-9;                    % desired accuracy in u (1e-6 sufficient)

design = 'NI2';   % choose design from below list... (may override above params)
switch design 
 case 'disc'
  Np = 1; r0=5; r1=10; Afunc = @(r) 1+0*r; % 10m radius disc: Poisson spot u=1!
  n = 600; m = 80;                  % Afunc=1, r0, dummy hacks. Good to lz>=5 
 case 'lin'
  Np = 24; r0=5; r1=13; Afunc = @(r) (r1-r)/(r1-r0);   % simple linear apod.
  n = 600; m = 80;        % good to lambdaz>=9; convergence must be tested
 case 'HG'
  A = @(t) exp(-(t/0.6).^6);        % a Cash'11 "hyper-Gaussian" on [0,1]
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
  cwd = fileparts(mfilename('fullpath')); file = [cwd '/occulter/NI2'];
  Np = 24;                          % petals told cut off harshly at r1
  [~,r0,r1] = eval_sister_apod(file,0);   % get apodization range [r0,r1]
  Afunc = @(r) eval_sister_apod(file,r);  % func handle (reads file when called)
  % try destroying the 1e-2 Poisson spot due to both petal flat ends...
  %Afunc = @(r) Afunc(r).*((r<=12) + (r>12).*cos(pi/2*(r-12)).^2); % C^1 blend
  %Afunc = @(r) 1 - (1-Afunc(r)).*((r>=6) + (r<6).*cos(pi/2*(6-r)).^2);  % "
  lambdaz = 16;       % lam.z in [16,21] m^2, low Fres#~5. Fails outside this!
  n = 40; m = 300;    % use quad_conv_apod_NI2 converged m (to 1e-6), lam.z>=5
 case 'TV3'
  cwd = fileparts(mfilename('fullpath'));    % insert your SISTER location
  file = [cwd '/../SISTER/input_scenes/locus/in/TV3'];
  Np = 24;                          % petals told cut off harshly at r1
  [~,r0,r1] = eval_sister_apod(file,0);   % get apodization range [r0,r1]
  Afunc = @(r) eval_sister_apod(file,r);  % func handle (reads file when called)
  n = 40; m = 800;               % converged m (to 1e-6), lam.z>23
  lambdaz = 30; ximax = 26;      % lam.z in [23,76] m^2. shadow only 4e-9.
end
fprintf('apod ends: 1-A(%.5g)=%.3g, A(%.5g)=%.3g\n',r0,1-Afunc(r0),r1,Afunc(r1))

tic
[xq yq wq bx by] = starshadequad(Np,Afunc,r0,r1,n,m,verb);   % fill areal quadr
fprintf('starshade quad gen time: %.3g s\n',toc)

if verb>1, figure(1); clf; scatter(xq,yq,10,wq); axis equal tight; colorbar;
  hold on; plot([bx;bx(1)], [by;by(1)], 'k-'); title('starshade quadr');
  xi = -5; eta = -10;     % see Fresnel integrand resolved for a target...
  int = exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2));  % integrand w/o prefac
  figure(2); clf; scatter(xq,yq,10,real(int)); axis equal tight; colorbar;
  hold on; plot([bx;bx(1)], [by;by(1)], 'k-'); title('Fresnel integrand');
  drawnow
end

[u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb);
it = 1; jt = 1;  % grid indices to test, eg (1,1) is SW corner of grid
ut = u(it,jt); xi = xigrid(it); eta = xigrid(jt);     % math check u
fprintf('u(%.3g,%.3g) = %.12g + %.12gi\n',xi,eta,real(ut),imag(ut))
%return                        % useful for convergence testing the above u

figure(3); clf;
u = 1-u;               % convert aperture to occulter
imagesc(xigrid,xigrid,log10(abs(u)'.^2));  % note transpose: x is fast, y slow.
colorbar; caxis([-11 0.2]); hold on; axis xy equal tight;
plot([bx;bx(1)], [by;by(1)], 'k+-','markersize',1);
xlabel('\xi (m)'); ylabel('\eta (m)'); if verb>1, overlay_zones(lambdaz); end
title(sprintf('%s: log_{10} |u|^2, lambda.z=%.3g m^2',design,lambdaz))
if strcmp(design,'NI2') % || strcmp(design,'NW2')
  show_occulter_locus(file); axis(ximax*[-1 1 -1 1]); end    % check raw locus

% check on-axis case and compare to tip/gap fractional Poisson spot prediction:
it = ceil(ngrid/2+1); jt = it;      % indices of center xi=eta=0
fprintf('|u|^2 intensity at (0,0) is %.3g\n',abs(u(it,jt)).^2)
ugap = (1-Afunc(r0))*exp((1i*pi/lambdaz)*r0^2);   % simple fraction-of-Pois-spot
utip = Afunc(r1)*exp((1i*pi/lambdaz)*r1^2);
fprintf('tip+gap only scatt pred is  %.3g\n',abs(ugap+utip).^2)   % phased add?
