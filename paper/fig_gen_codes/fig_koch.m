% Fresnel from Koch snowflake. gluing general triangle to each edge.
% Barnett 10/1/20
clear
addpath ~/numerics/finufft/matlab
fdir = '/home/alex/physics/starshade/fresnaq';
addpath(fdir)
startup
verb = 0;  % 1 to make paper figs

% warm up w/ usual equilateral Koch ..............................
hmax = 4e-3;    % worst-case node sep (big tri's)
%pmin = 2;       % what p to bottom out at (small tri's)
Tb = exp(2i*pi/3*(1:3)');   % base tri, CCW
ntri = 1;
basearea = 0.5 * imag(conj(Tb(2)-Tb(1))*(Tb(3)-Tb(1)));
len = abs((Tb(2)-Tb(1)));   % its side len
p = ceil(1.5 * len/hmax);   % rule for p at large scales (resolve osc)
[xq yq wq] = triquad(real(Tb),imag(Tb),p);   % base quadr
% the generator for [0,1], and its area...
T = [1/3; 1/2-1i*sqrt(3)/6; 2/3];     % note order & flipped below x axis
Tarea = 0.5 * imag(conj(T(2)-T(1))*(T(3)-T(1)));
f = abs(T(3)-T(2));         % zoom scale of fractal per level
vb = Tb;           % base vert list
arealim = basearea + numel(vb)*abs(vb(2)-vb(1))^2*Tarea/(1 - 4*f^2);  % true

ftol = 1e-6; ximax = 1; ngrid = 1e3; lambdaz = 0.01;
atol = 1e-5;   % areal target tol
p0 = ceil( 0.8 * log(1/atol) );   % cross-over p, at scales~hmax; 0.8=converged
L=13;   % full is L=13, but takes 4 mins; for plotting can get away with less
v = vb;
tquad = 0;
for l=1:L         % can be for convergence too
  newlen = abs((v(2)-v(1))*(T(2)-T(1)));  % new tri scale, all abs(diff(v)) same
  n = numel(v);
  ntri = ntri + n;
  newarea = n*abs(v(2)-v(1))^2*Tarea;
  if newlen>=hmax
    p = max(p0, ceil(1.5 * newlen/hmax));   % rule for p at large scales
  else
    p = min(p0, ceil( 0.5 * (log(1/atol)/log(hmax/newlen) - 1) ) );  % small
  end
  fprintf('lev %d: %d tris (p=%d) of len %.3g, new area %.3g: %d nodes...\n',l,n,p,newlen,newarea,p^2*n)
  tic
  [v x y w] = addkochlayer(v,T,p);
  xq = [xq;x]; yq = [yq;y]; wq = [wq;w]; clear x y w
  tquad = tquad + toc;
  if verb & l==L-2, figure(1); clf;   % plot not at finest (can't tell): use L-2
    scatter(xq,yq,1,wq,'.'); hold on; plot([v;v(1)], '-'); axis equal tight;
    colormap(jet(256)); colorbar; xlabel('\xi'); ylabel('\eta');
    % note specific for displaying at l==L-2...
    title(sprintf('(a) level-%d Koch snowflake ($N=%d$ nodes)',L,numel(xq)+4*n+16*n),'interpreter','latex')
    cx = 0.58; cy = -0.51; cw = 0.04;    % where to zoom (cw = half-width)
    plot(cx + cw*[-1 -1 1 1 -1], cy + cw*[-1 1 1 -1 -1],'k-');  % zoom box
    plot([cx 0.8], [cy+cw 0.6],'k-');          % pointing line
    axes('position',[0.57 0.62 0.25 0.25]);
    scatter(xq,yq,10,wq,'.'); hold on; plot([v;v(1)], '-'); axis equal
    axis([cx-cw cx+cw cy-cw cy+cw]); box on; set(gca,'xtick',[],'ytick',[]);
    drawnow
    set(gcf,'paperposition',[0 0 5 4]); print -depsc2 -r200 ../kochquad.eps
    %stop
  end  
  tic; [u xigrid] = fresnaq_grid(xq,yq,wq, lambdaz, ximax, ngrid, ftol);
  fprintf('tot nodes so far N=%d (edges %d), quadr so far %.3g s, this t1 %.3g s\n',numel(xq),numel(v),tquad,toc)
  if l>=10, us{l}=u; end % save it
  if l==L, clear v; end   % save RAM
  fprintf('\trel area err = %.3g\t\t\t|u(test)|=%.8g\n',1-sum(wq)/arealim,abs(u(300,400)))
end
%fprintf('tot nodes N=%d, in %.3g s\n',numel(xq),toc)

max(abs(us{L-1}(:)-us{L}(:)))
max(abs(us{L-2}(:)-us{L}(:)))
% their ratio is close to 4/13 as expect from a sum_k (4/9)^k areal convergence

if verb>1
  figure(1); clf; plot([v;v(1)], '-'); axis equal; hold on;
  plot([vb;vb(1)],'k-');
  scatter(xq,yq,20,wq,'.'); colormap(jet(256)); colorbar;
end

if verb        % *** to do: add inset?
  figure(2); clf;
  imagesc(xigrid,xigrid,abs(u)'.^2);
  colormap(hot(256)); colorbar; hold on;
  %plot([bx,bx(1)],[by,by(1)],'w-');
  axis xy equal tight;
  xlabel('\xi'); ylabel('\eta');
  caxis([0 1.6]);   % gets most of the high end, not all
  title(sprintf('(b) $$\\lambda z=%.3g$$: t1 aperture intensity $$|u^{ap}|^2$$',lambdaz),'interpreter','latex')
  set(gcf,'paperposition',[0 0 5 4])
  print -depsc2 ../kochu.eps
end



%%%%%%%%%%%
function [v xq yq wq] = addkochlayer(v,T,p)
  % AQ for 1 Koch layer, spitting out new vertex list. T is tri to add to [0,1]
  % v is list of C-#s in CCW, and T is list of 3 C-#s in CCW. p is order
  n=numel(v);
  N = p^2*n;   % AQ nodes
  xq = nan(N,1); yq=xq; wq=xq;   % alloc
  vold = v;
  v = zeros(4*n,1);       % alloc new vert list
  [x y w] = triquad(real(T),imag(T),p);    % call only once (much faster)
  w = w*abs(vold(2)-vold(1))^2;            % size-scale w's only once
  z = x+1i*y;
  for i=1:n
    inxt = mod(i,n)+1;
    vi = vold(i);           % base pt for this edge
    d = vold(inxt)-vi;     % current edge vector
    v((1:4)+(i-1)*4) = [vi; vi + d*T(:)];    % new vert in order
    ii = (1:p^2)+(i-1)*p^2;
    wq(ii) = w;                % rot or trans no effect
    xq(ii) = real(vi) + real(d*z);
    yq(ii) = imag(vi) + imag(d*z);
  end
end
