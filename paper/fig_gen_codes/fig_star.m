% paper figs for starshade. Just the AQ geom plot for now.
clear

% ------- Cash'11 paper:
Np = 16;         % appears from Fig. 4, but mentions Np=12 too.
a = 12.5; b = a; pow=6;   % a,b in meters.  pow is his "n".
R = 31;   % outer max radius, meters.
z=8e7;   % 80000km
lambda = 5e-7;  % longest wavelength (smallest Fr), at which Cash claims good?
A = @(r) exp(-((r-a)/b).^pow);  % "offset hyper-Gaussian" in r
% -------

n=30;
m=80;
rquad = 'g';
[xj yj wj] = starshadequad(Np,A,a,R,n,m,1,rquad);
figure(1); clf; ms=0.1; scatter(xj,yj,ms,wj,'.'); axis equal;
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis([-2 31.5 -4 7]);
colormap(jet(256)); colorbar
lc = [.6 .8 1];  % line color
hold on; t=2*pi*(0:1e3)/1e3; plot(R*cos(t),R*sin(t),'-','color',lc);
plot([0 R*cos(pi/Np)],[0 R*sin(pi/Np)],'-','color',lc);
plot([0 R*cos(pi/Np)],[0 -R*sin(pi/Np)],'-','color',lc);
set(gcf,'paperposition',[0 0 8 2]); print -depsc2 aqstar.eps
