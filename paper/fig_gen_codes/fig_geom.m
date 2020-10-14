% geom overview figure
% needs: mArrow3, fresnaq
clear
addpath ~/numerics/finufft/matlab
addpath ../../fresnaq
startup

x = @(t) -0.3 + 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);  % kite
N = 200; t = 2*pi*(1:N)/N;
bx = x(t); by = y(t);      % bdry points
wx = (2*pi/N)*perispecdiff(bx); wy = (2*pi/N)*perispecdiff(by);  % LI weights

figure(1); clf

%patch(bx,by,0*t,'-');
plot3([bx bx(1)],[by by(1)],0*[t 0],'-','linewidth',1.5,'color',[0 0 0]);
hold on;
% params for Frensel integrand...
s = 1.2;  % half-size of each "plane"
z = 2.5;  % downstream
xi = 0.2; eta = 0.3; lambda = 0.07;
sx = -0.1; sy = 0.2;        % src pt
g = linspace(-s,s,300); [gx gy] = ndgrid(g,g);  % source plane grid
int = exp((1i*pi/(lambda*z))*((gx-xi).^2+(gy-eta).^2));  % Fresnel integrand
inside = 0*gx;  % line integral (=2D Laplace DLP) to compute grid pts in/out...
for j=1:N
  dx = bx(j)-gx; dy = by(j)-gy;
  inside = inside + (dx*wy(j)-dy*wx(j))./(dx.*dx + dy.*dy);
end
inside = inside * (1/(2*pi));   % ~0 outside, ~1 inside
inside = min(1,max(0,inside));  % clip
h=surf(gx,gy,0*gx,real(int),'linestyle','none','facecolor','texturemap','facealpha','texturemap','alphadata',inside);
grey=0.4; colormap(grey*[1 1 1]+(1-grey)*jet(256)); shading interp
patch([bx bx(1) s*[-1 1 1 -1 -1]], [by by(1) s*[-1 -1 1 1 -1]], zeros(1,N+6),'facecolor',0.5*[1 1 1],'linestyle','none');  % nonconvex: square minus Omega

csrc = [0 .5 0];          % source plane label colors
plot3(s*[-1 1 1 -1 -1],s*[-1 -1 1 1 -1],zeros(1,5),'-');   % src plane
plot3(sx,sy,0,'.','markersize',10,'color',csrc);
text(-s,-s-.2,0,'source plane','color',csrc,'rotation',+15,'fontsize',7);
ctrg = [.5 0 0];          % targ plane label colors
h = patch('xdata',s*[-1 1 1 -1],'ydata',s*[-1 -1 1 1],'zdata',z*ones(1,4),'edgecolor',[.5 0 0],'facecolor',[1 1 1],'facealpha',0.5);  % targ plane
plot3(xi,eta,z,'.','markersize',10,'color',ctrg);
text(-s+.2,-s+.2,z,'target plane','color',ctrg,'rotation',+15,'fontsize',8);

plot3([sx xi],[sy eta],[0 z],'k--','linewidth',0.7);  % src-targ dist rho
text(0,0.2,0.75*z,'$\rho$','interpreter','latex','color',[0 0 0]);

axis equal vis3d off;
camproj('perspective')
set(gca, 'CameraPosition', 5*[-2,1,3])
set(gca, 'CameraTarget', [0,0.2,0.6*z])
set(gca, 'CameraUpVector', [0,1,0])
set(gca,'cameraviewangle',10)

% arrows...
p1 = zeros(5,3); p1(4:5,3) = z;    % start pts as rows
p2 = [0,0,z+0.8; s+0.4,0,0; 0,s+0.4,0; s+0.4,0,z; 0,s+0.4,z]; % end pts
cc = ones(5,1)*csrc; cc(4:5,:) = [1;1]*ctrg;    % colors as rows
for a=1:size(p1,1)
  h = mArrow3(p1(a,:),p2(a,:),'color',cc(a,:),'stemWidth',.01,'tipWidth',.05);
end

csrc = csrc*0.6;          % source colors make darker
text(-0.15,-0.05,0,'$0$','interpreter','latex','color',csrc);
text(0,0,z+1,'$z$','interpreter','latex','color',csrc)
text(s+.5,0,0,'$x$','interpreter','latex','color',csrc);
text(0,s+.6,0,'$y$','interpreter','latex','color',csrc);
text(-0.2,-0.5,0,'$\Omega$','interpreter','latex','color',[0 0 0],'rotation',+10,'fontsize',12);
text(-1.1,0.7,0,'$\partial\Omega$','interpreter','latex','color',[0 0 0],'rotation',+10);
text(s+.5,0,z,'$\xi$','interpreter','latex','color',ctrg);
text(0,s+.6,z,'$\eta$','interpreter','latex','color',ctrg);
text(xi,eta,z,'$\;(\xi,\eta)$','interpreter','latex','color',ctrg,'rotation',+12)


set(gcf,'paperposition',[0 0 3 2.5]);
print -depsc2 -r300 ../geom.eps
%hgexport(gcf,'../meth.eps');  % ignores paperposition

