% methods overview figure
% needs: mArrow3, fresnaq
clear
addpath ~/numerics/finufft/matlab
addpath ../../fresnaq
startup

x = @(t) -0.3 + 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);  % kite
N = 100; t = 2*pi*(1:N)/N;
bx = x(t); by = y(t);      % bdry points
wx = (2*pi/N)*perispecdiff(bx); wy = (2*pi/N)*perispecdiff(by);  % LI weights
m = 20;
% shift the areal origin so star-shaped, less confusing...
xsh = 0.3;
[xq yq wq] = curveareaquad(bx-xsh,by,wx,wy,m);   % areal quadr
xq = xq+xsh;

figure(1); clf
n = 32;
s = 1.0;   % half-size
g = linspace(-s,s,n); [gx gy] = ndgrid(g,g);  % source plane grid
inside = 0*gx;  % line integral (=2D Laplace DLP) to compute grid pts in/out...
for j=1:N
  dx = bx(j)-gx; dy = by(j)-gy;
  inside = inside + (dx*wy(j)-dy*wx(j))./(dx.*dx + dy.*dy);
end
inside = inside * (1/(2*pi));   % ~0 outside, ~1 inside
siz = 1.12;

subplot(1,3,1);
%image(g,g,reshape(kron([1 1 1],inside'),[n n 3]));
surf(gx,gy,0*gx,double(inside<0.5));  % surf doesn't like boolean
c = jet(256); c(1,:) = [1 1 1]; c(end,:) = [0 0 0];   % jet w/ b&w at ends
colormap(c);
view(2); axis equal tight xy off;
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
text(-0.1,-0.3,'$\Omega$','interpreter','latex','color',[0 0 0],'fontsize',16);
title('(a) uniform 2D grid sampling');
csrc = [0,.5,0];
arrow([0 0],[0.5 0],'color',csrc);
arrow([0 0],[0 0.5],'color',csrc);
text(.5,.1,'$x$','interpreter','latex','color',csrc);
text(0,.6,'$y$','interpreter','latex','color',csrc);
axis(siz*s*[-1 1 -1 1]);   % tight fails to match other subplots

subplot(1,3,2);
plot([bx bx(1)],[by by(1)],'.-','linewidth',1.0,'markersize',12);
axis equal tight off;
axis(siz*s*[-1 1 -1 1]);   % tight fails to match other subplots
title('(b) line integral quadrature');
text(-1.1,0.7,'$\partial\Omega$','interpreter','latex','color',[0 0 0],'fontsize',14);
arrow([0 0],[0.5 0],'color',csrc);
arrow([0 0],[0 0.5],'color',csrc);
text(.5,.1,'$x$','interpreter','latex','color',csrc);
text(0,.6,'$y$','interpreter','latex','color',csrc);


subplot(1,3,3);
scatter(xq,yq,2.0,wq,'.'); caxis([0.01 0.99]*max(caxis)); % avoid the b&w ends
hold on; plot([bx bx(1)],[by by(1)],'k-','linewidth',0.1);
axis equal tight off;
title('(c) high-order areal quadrature');
colorbar
arrow([xsh 0],[0.5+xsh 0],'color',csrc);    % fake-shift origin to (xsh,0)
arrow([xsh 0],[xsh 0.5],'color',csrc);
text(xsh+.5,.1,'$x$','interpreter','latex','color',csrc);
text(xsh,.6,'$y$','interpreter','latex','color',csrc);

axis(siz*s*0.9*[-1 1 -1 1]);   % tight fails to match other subplots

set(gcf,'paperposition',[0 0 9 2.2]);
print -depsc2 -r300 ../meth.eps
%hgexport(gcf,'../meth.eps');  % ignores paperposition

