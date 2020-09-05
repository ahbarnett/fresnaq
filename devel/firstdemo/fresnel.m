% Demo of scalar Fresnel diffraction from polar characteristic function,
% using the NUFFT and a planar quadrature rule, in MATLAB. Barnett 8/24/20

% Key point is aperture func is not bandlimited, unlike those in Mas'99, etc.

clear; close all; verb = 1;

g = @(t) 1 + 0.3*cos(3*t);                             % boundary shape

fbar = 100;   % "increased" Fresnel #, fbar = 2pi.f = k_free/d = 2pi/(lambda.d)

kirchfac = fbar/(2i*pi);   % Kirchhoff approx prefactor = 1/i.lambda.d

% build the quadrature...
n = 350;                                               % # theta nodes
t = 2*pi*(1:n)/n; wt = (2*pi/n);                       % theta nodes, const weights
bx = cos(t).*g(t); by = sin(t).*g(t);                  % boundary points
m = 120;                                               % # r nodes
[xr,wr] = lgwt(m,0,1);                                 % rule for (0,1)
xj = nan(n*m,1); yj = xj; wj = xj;
for i=1:n                                              % loop over angles
  r = g(t(i)); jj = (1:m) + (i-1)*m;                   % this radius; index list
  xj(jj) = cos(t(i))*r*xr; yj(jj) = sin(t(i))*r*xr;    % line of nodes
  wj(jj) = wt*r^2*xr.*wr;            % theta weight times rule for r.dr on (0,r)
end

xmax = 2.0;       % targ params: box half-size
M1 = 1e3;         % target location grid size per dim
tol = 1e-9;
dx = 2*xmax/M1;
x1 = dx * (-M1/2:(M1/2-1));    % target grid each dim

%xt = 1.0; yt = 0.5;   % slow eval Fresnel integral at single target
%xt=0; yt=0;
xt = -xmax; yt=-xmax;   % make sure it's on the targ grid
ut = kirchfac * sum(exp(0.5i*fbar*((xj-xt).^2+(yj-yt).^2)) .* wj)

if verb, figure(1); clf; subplot(1,2,1);
  plot([bx bx(1)],[by by(1)],'-'); hold on; plot(xj,yj,'.'); axis equal tight;
  plot(xt,yt,'*'); xlabel('x_1'); ylabel('x_2'); title('Aperture + nodes');
  subplot(1,2,2); scatter(xj,yj,10,cos(0.5*fbar*((xj-xt).^2+(yj-yt).^2)));
  axis equal tight; xlabel('x_1'); ylabel('x_2'); title('Re Fresnel integrand');
end



tic % do it...
dk = fbar*dx;
k1 = fbar*x1;                  % freq grid is simply scaled location grid
fprintf('fbar=%.3g, xmax=%.3g: kmax = %.3g\n',fbar,xmax,max(k1))
cj = exp(0.5i*fbar*(xj.^2+yj.^2)) .* wj;        % premult by a quadratic bit
u = finufft2d1(dk*xj, dk*yj, cj, -1, tol, M1, M1);   % M1^2 output nodes
u = u .* (exp(0.5i*fbar*x1(:).^2) * exp(0.5i*fbar*x1.^2));  % postmult by quadr bit
u = kirchfac * u;
toc

% check
j1 = find(x1==xt); j2 = find(x1==yt); errt = u(j1,j2) - ut;
fprintf('abs error at (%.3g,%.3g) = %.3g\n\n',xt,yt,abs(errt))
figure(2); clf;
subplot(1,2,1); imagesc(x1,x1,abs(u)'); colorbar; hold on; title('aperture');
plot([bx bx(1)],[by by(1)],'-'); plot(xt,yt,'*'); axis xy equal tight;
colormap(hot(256));
subplot(1,2,2); imagesc(x1,x1,abs(1.0-u)'); colorbar; hold on; title('blocker');
plot([bx bx(1)],[by by(1)],'-'); plot(xt,yt,'*'); axis xy equal tight;
colormap(hot(256));
