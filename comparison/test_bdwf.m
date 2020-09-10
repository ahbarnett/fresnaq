% test BDWF on simple occulter shapes that we understand.
% Barnett 9/10/20

% parameterization of smooth closed curve, for t in [0,2pi)...
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);          % kite shape
fresnum = 20.0;        % Fresnel number (if char radius were R=1)
lambdaz=1/fresnum;     % wavelength times dist, recall Fres # = R^2/(lambda.z)
n = 500; m = 100;      % number of boundary, and radial, quadrature points;
                       % depend on Fresnel number & target pt; tested to 1e-12
t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);      % bdry points

xi = 1.5; eta = -0.5;   % math test: target to test at

kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor, direct Fresnel integral...
bxp = perispecdiff(bx); byp = perispecdiff(by);  % derivatives wrt t (if smooth)
[xq yq wq] = curveareaquad(bx,by,bxp,byp,m);   % areal quadrature
ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi).^2+(yq-eta).^2)) .* wq);
ud = 1-ud              % Babinet the aperture into an occulter

% now we try BDWF with the same n bdry pts...
psi1=0; psi2=0;     % on-axis incident wave, and no Z variation
ub = bdwf_pts(bx,by, [], 1.0, lambdaz, xi, eta, psi1, psi2)
fprintf('bdwf (n=%d) err vs direct Fresnel quadr: %.3g\n',n,abs(ub-ud))

% convergence study of bdwf wrt n...
ns = 100*2.^(1:7);
errs = 0*ns;
for i=1:numel(ns), n=ns(i);
  t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);      % bdry points
  ub = bdwf_pts(bx,by, [], 1.0, lambdaz, xi, eta, psi1, psi2);
  errs(i) = abs(ub-ud);
  disp([n errs(i)])
end
figure; loglog(ns, errs, '+-'); axis tight; title('bdwf convergence to true');
xlabel('n (bdry nodes)'); ylabel('bdwf(n) error');
hold on; plot(ns,1.0./ns,'r-');   % appears to be only 1st-order (!) in h~1/n

