% timings and acc table for param smooth curves (kite). Barnett 9/21/20
addpath ~/numerics/finufft/matlab
addpath ../../
startup

if 1 %%%%%%%%%%
clear
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);   % param, C^infty smooth
tols = [1e-6 1e-12];                % fresnaq to test
xi0 = -1.5; eta0 = -1.5;            % test pt 
ximax = abs(xi0); ngrid = 1e3;      % test grid (1,1 must be test pt, extreme)
M=ngrid^2; xi = ximax*(2*rand(M,1)-1); eta = ximax*(2*rand(M,1)-1);   % arb pts

lambdazs = [0.1 0.01]; ms = [80 560]; ns = [320 2400];   % matching pairs

for i=1:numel(lambdazs)  %......................................... main loop
  lambdaz=lambdazs(i); n=ns(i); m=ms(i);
  fprintf('\nlambda.Z=%.3g: n=%d, m=%d...\n',lambdaz, n,m)
  t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
  wx = (2*pi/n)*perispecdiff(bx); wy = (2*pi/n)*perispecdiff(by);  % quadr wei
  t = tic;
  u0 = nsli_pts(bx,by,wx,wy, lambdaz, xi,eta);      % use as ref
  fprintf('nsli:\t\t\t\t\t\t\ttime %.3g s\n',toc(t))
  Z = 1e4; lambda = lambdaz/Z;                      % BDWF needs lambda & Z
  t = tic;
  ub = bdwf_pts([bx bx(1)],[by by(1)],[],Z, lambda, xi,eta,0,0);
  u0ocprop = (1-u0)*exp(2i*pi*Z/lambda);           % what bdwf computes
  err = abs(u0ocprop(:)-ub(:));
  fprintf('bdwf:\t\t\tmederr %.2g, maxerr %.2g\ttime %.3g s\n',median(err),max(err),toc(t))
  [xq yq wq] = curveareaquad(bx,by,wx,wy,m);       % areal quadr for fresnaq
  for j=1:2
    t = tic;
    uf = fresnaq_pts(xq, yq, wq, lambdaz, xi,eta, tols(j));
    err = abs(u0(:)-uf(:));
    fprintf('faq_pts (%g):\tmederr %.2g, maxerr %.2g\ttime %.3g s\n',tols(j),median(err),max(err),toc(t))
  end
  fprintf('recomputing ref now on grid...\n')
  g = 2*ximax*(-ngrid/2:ngrid/2-1)/ngrid; [xi eta] = ndgrid(g,g);  % the grid
  u0 = nsli_pts(bx,by,wx,wy, lambdaz, xi,eta);      % use as ref
  for j=1:2
    t = tic;
    ug = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tols(j));
    err = abs(u0(:)-ug(:));
    fprintf('faq_grid (%g):\tmederr %.2g, maxerr %.2g\ttime %.3g s\n',tols(j),median(err),max(err),toc(t))
  end
end

% diagnosing bdwf max errors associated with near geom shadow...
%figure; scatter(xi,eta,10,log10(abs(u0ocprop(:)-ub(:)))); colormap

else %%%%%%%%%%


% ----------------------- final high-fresnel # example
clear
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);   % param, C^infty smooth
tol = 1e-6;                % fresnaq to test
xi0 = -1.5; eta0 = -1.5;            % test pt 
ximax = abs(xi0); ngrid = 3162;     % test grid (1,1 must be test pt, extreme)
M=ngrid^2; xi = ximax*(2*rand(M,1)-1); eta = ximax*(2*rand(M,1)-1);   % arb pts
itest = 1:1e4;
verb = 2;

lambdaz=1e-3; m = 5600; n = 24000;  % guessed, 10x the lamz=1e-2 cases.

t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
wx = (2*pi/n)*perispecdiff(bx); wy = (2*pi/n)*perispecdiff(by);  % quadr wei
u0 = nsli_pts(bx,by,wx,wy, lambdaz, xi(itest),eta(itest));      % use as ref
[xq yq wq] = curveareaquad(bx,by,wx,wy,m);       % areal quadr for fresnaq
t = tic;
uf = fresnaq_pts(xq, yq, wq, lambdaz, xi,eta, tol, verb);
err = abs(u0-uf(itest));
fprintf('faq_pts (%g):\tmederr %.2g, maxerr %.2g\ttime %.3g s\n',tol,median(err),max(err),toc(t))
fprintf('recomputing ref now on grid...\n')
g = 2*ximax*(-ngrid/2:ngrid/2-1)/ngrid; [xi eta] = ndgrid(g,g);  % the grid
itest = randi(ngrid^2,[1e4 1]);
u0 = nsli_pts(bx,by,wx,wy, lambdaz, xi(itest),eta(itest));      % use as ref
t = tic;
ug = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tol, verb);
err = abs(u0-ug(itest));
fprintf('faq_grid (%g):\tmederr %.2g, maxerr %.2g\ttime %.3g s\n',tol,median(err),max(err),toc(t))

end %%%%%%%%%%%
