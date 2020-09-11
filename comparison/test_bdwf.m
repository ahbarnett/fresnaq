% test BDWF with generic targets on simple occulter shapes that we understand.
% Barnett 9/10/20

% parameterization of smooth closed curve, for t in [0,2pi)...
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);          % kite shape
fresnum = 20.0;        % Fresnel number (if char radius were R=1)
lambdaz=1/fresnum;     % wavelength times dist, recall Fres # = R^2/(lambda.z)
n = 500; m = 100;      % number of boundary, and radial, quadrature points;
                       % depend on Fresnel number & target pt; tested to 1e-12
t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);      % bdry points

xi = 1.5; eta = -0.5;   % math test: generic target to test at

kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor, direct Fresnel integral...
bxp = perispecdiff(bx); byp = perispecdiff(by);  % derivatives wrt t (if smooth)
[xq yq wq] = curveareaquad(bx,by,bxp,byp,m);   % areal quadrature
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2)) .* wq);
ud = 1-ud              % Babinet the aperture into an occulter

% now we try BDWF with the same n bdry pts...
psi1=0; psi2=0;     % on-axis incident wave, and no Z variation
lambda = 1e-5;      % typ wavelength rel to size, makes Z huge, non-convergent!
Z = lambdaz/lambda;
ub = bdwf_pts(bx,by, [], Z, lambda, xi, eta, psi1, psi2)
fprintf('bdwf (n=%d) err vs direct Fresnel quadr: %.3g\n',n,abs(ub-ud))

% convergence study of bdwf wrt n...
ns = 100*2.^(1:9);
errs = 0*ns;
disp('                      n         err')
for i=1:numel(ns), n=ns(i);
  t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);      % bdry points
  ub = bdwf_pts(bx,by, [], Z, lambda, xi, eta, psi1, psi2);
  errs(i) = abs(ub-ud);
  disp([n errs(i)])
end
figure; loglog(ns, errs, '+-'); axis tight; title('bdwf convergence to true');
xlabel('n (bdry nodes)'); ylabel('bdwf(n) error');
hold on; plot(ns,1.0./ns,'r-');   % appears to be only 1st-order (!) in h~1/n
% also: if lambda<=1e6 then Z>=5e4, and does not converge! (round-off error?)

% check translational invariance...
ub2 = bdwf_pts(bx-xi, by-eta, [], Z, lambda, 0,0, psi1, psi2);
fprintf('bdwf (n=%d) change due to translation: %.3g\n',n,abs(ub-ub2))

% check incident angle corresp to translation of targs by (-xi,-eta)...
psi1 = sqrt(xi^2+eta^2)/Z;
psi2 = atan2(-eta,-xi);
ub3 = bdwf_pts(bx, by, [], Z, lambda, 0,0, psi1, psi2);
fprintf('bdwf (n=%d) change of incidence: %.3g, should be O(psi1^2)=%.3g\n',n,abs(ub-ub3),psi1^2)

% now test bdwf: check a (xi,eta) target grid matches fresnap_pts result...
dxO = 0.1; nO = 10; deltaX = 0.23; deltaY = -0.16;   % generic
[xi,eta] = make_grid_bdwf(dxO, nO, deltaX, deltaY);
t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);      % bdry points
t0 = tic;
ub = bdwf(bx,by,[], Z, lambda, dxO, nO, 0,0, deltaX, deltaY);  % what we test
t0 = toc(t0);
tol = 1e-9;
u = fresnap_pts(xq, yq, wq, lambdaz, xi(:), eta(:), tol);    % col vec of targs
u = 1-u;
fprintf('bdwf (n=%d) grid max err vs fresnap_pts: %.3g\n',n,norm(ub(:)-u,inf))
% is consistent w/ bdwf_pts error, fine.
figure; plot(xi,eta,'k.'); hold on; plot([bx bx(1)], [by by(1)], '-');
axis equal tight; title('test\_bdwf grid targets');

fprintf('bdwf one lambda (n=%d, nTargs=%d) in %.3g s = %.3g targ-bdry pairs/s\n',n,nO^2,t0,nO^2*n/t0)
% bdwf speed = 3e7 pairs/sec on laptop.

