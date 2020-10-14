% forward rounding error analysis of NSLI for target near node.
% taken from nsli_pts.m. Needs couple of fresnaq utils.
% Barnett 9/24/20
addpath ../../fresnaq/utils

x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);        % smooth kite shape
fresnum = 20.0;        % Fresnel number (if char radius were R=1)
lambdaz=1/fresnum;     % wavelength times dist, recall Fres # = R^2/(lambda.z)
n = 500; m = 100;      % number of boundary, and radial, quadrature points;
                       % depend on Fresnel number & target pt; tested to 1e-12
t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
bxp = perispecdiff(bx); byp = perispecdiff(by);  % param-derivatives (if smooth)
wx = (2*pi/n)*bxp; wy = (2*pi/n)*byp;            % LI vector quadr wei

kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor for direct Fresnel integral
[xq yq wq] = curveareaquad(bx,by,wx,wy,m);       % areal quadrature

% here 1.5e-5 causes real(exp(...))=1 exactly, but 2e-5 not so.
% thus worst-case falls between these two distances! around 1.7e-5.
j=400; d = [1.7e-5]; % 1e-8]; %[0 logspace(-17,-1,17)];  % check close targ dists from a node...
xi = bx(j) + 0*d; eta = by(j)-d;  ud = complex(nan*xi);
for i=1:numel(xi)
  ud(i) = kirchfac*sum(exp((1i*pi/lambdaz)*((xq-xi(i)).^2+(yq-eta(i)).^2)).*wq);
end

% break out the code for:   u = nsli_pts(bx,by,wx,wy, lambdaz, xi, eta);
jj = 399:403;
u = complex(nan*xi);                    % alloc complex outputs, shape of xi
for i=1:numel(xi)
  dx = bx - xi(i); dy = by - eta(i);    % displacements of nodes from ith targ
  r2 = dx.*dx + dy.*dy;
  t = 1-exp((1i*pi/lambdaz)*r2);
  %t = t + 2e-16;   % real pert
  t(jj)'
  dlp = (dx.*wy - dy.*wx) ./r2;
  dlp(jj)'
  f = t ./ r2;                % my way, or ...
  f(jj)'
  %al = 1i*pi/lambdaz; y = exp(al*r2); f = al*(1-y)./log(y); % Higham 1.14 trick
  % Higham fails due to complex log destroying the imag winding number part!!!
  % ... need to apply only for |al|*r2 < 1.   Not interested in trying this.
  %f(jj)'
  f(r2==0.0) = 0.0;                     % kill NaNs (targ = some node)
  u(i) = (1/2/pi) * sum((dx.*wy - dy.*wx) .* f);
end



disp('target distance from a node, resulting errors (vs areal method):')
disp([d(:) abs(u(:)-ud(:))])
