% convergence for simple kite domain. Incl self-conv on 1e6 targs, hence 30 s.
% Barnett 9/20/20

clear
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);   % param, C^infty smooth
tols = [1e-6 1e-12];
xi0 = -1.5; eta0 = -1.5;            % test pt 
ximax = abs(xi0); ngrid = 1e3;      % test grid (1,1 must be test pt, extreme)
M=ngrid^2; xi = ximax*(2*rand(M,1)-1); eta = ximax*(2*rand(M,1)-1);   % arb pts
xi(1) = xi0; eta(1) = eta0;         % (make 1 the test pt)

figure(1); clf; let = 'a';          % subfig letter labels
lambdazs = [0.1 0.01];
for i=1:numel(lambdazs)  %......................................... main loop
  lambdaz=lambdazs(i); fprintf('lambda.Z = %.3g...\n',lambdaz)
  % establish the reference value via NSLI
  if i==1; ns = 200:40:400; else, ns = 2000:200:3000; end
  for q=1:numel(ns), n=ns(q);
    t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
    bxp = perispecdiff(bx); byp = perispecdiff(by);  % param-derivs (if smooth)
    wx = (2*pi/n)*bxp; wy = (2*pi/n)*byp;            % high-ord LI vec quadr wei
    u = nsli_pts(bx,by,wx,wy, lambdaz, xi0,eta0);
    fprintf('n=%d:\tu=%.15g\t+ %.15gi\n',n,real(u),imag(u))
  end
  uref0 = u;  % ref
  
  % fix n, do conv m...
  if i==1; n=360; ms = 25:5:85; else, n=2400; ms = 260:20:560; end
  eds = nan*ms; efs=eds; egs=eds; efs2=eds; egs2=eds; % keep various error conv
  clear uf2 ug2
  for q=1:numel(ms), m=ms(q); fprintf('m=%d...\n',m)
    [xq yq wq] = curveareaquad(bx,by,wx,wy,m);       % areal quadrature  
    kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor for direct Fresnel
    ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi0).^2+(yq-eta0).^2)) .* wq);
    eds(q) = abs(ud-uref0);                          % direct err
    u = fresnaq_pts(xq, yq, wq, lambdaz, xi,eta, tols(1), q==numel(ms)-1);
    efs(q) = abs(u(1)-uref0);                        % fresnaq t3 err (tol 1)
    u = fresnaq_pts(xq, yq, wq, lambdaz, xi,eta, tols(2), q==numel(ms)-1);
    efs2(q) = abs(u(1)-uref0);                       % fresnaq t3 err (tol 2)
    uf2{q} = u;
    [u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tols(1), q==numel(ms)-1);
    egs(q) = abs(u(1,1)-uref0);             % fresnaq grid t1 err (tol 1)
    [u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tols(2), q==numel(ms)-1);
    egs2(q) = abs(u(1,1)-uref0);            % fresnaq grid t1 err (tol 2)
    ug2{q} = u;
  end
  efss=nan*ms; egss=efss;                   % self-conv sup over targs
  for q=1:numel(ms);                        % gather self errs
    egss(q)=norm(ug2{q}(:)-ug2{end}(:),inf);
    efss(q)=norm(uf2{q}(:)-uf2{end}(:),inf);
  end
  figoff = (i-1)*3;                                  % subplot # offset
  subplot(2,3,1+figoff);
  semilogy(ms,eds,'b.-', ms,efs,'k+:', ms,efs2,'k+-', ms,efss,'r+--', ms,egs,'go:',ms,egs2,'go-', ms,egss,'ro--');
  xlabel('$m$','interpreter','latex');
  ylabel('abs. error in $u$','interpreter','latex');
  axis([ms(1) ms(end-1) 1e-15 1]);
  h=legend('direct','t3 $\varepsilon{=}10^{-6}$','t3 $\varepsilon{=}10^{-12}$','t3 self $\varepsilon{=}10^{-12}$','t1 $\varepsilon{=}10^{-6}$','t1 $\varepsilon{=}10^{-12}$','t1 self $\varepsilon{=}10^{-12}$');
  set(h,'interpreter','latex','location','southwest')
  h.Position = h.Position + [0.007 0.009 0 0];    % doesn't match screen
  title(sprintf('(%c) $$\\lambda z=%.3g$$: AQ error vs $$m$$ (fix $$n=%d$$)',let,lambdaz,n),'interpreter','latex')
  drawnow
  let = let+1;
  
  % fix m, do conv n...
  m = ms(end-1);  % need to have checked converged
  if i==1; ns = 180:20:380; else, ns = 2000:100:2600; end
  eds = nan*ns; efs=eds; egs=eds; efs2=eds; egs2=eds; % keep various error conv
  clear uf2 ug2
  for q=1:numel(ns), n=ns(q); fprintf('n=%d...\n',n)
    t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
    wx = (2*pi/n)*perispecdiff(bx); wy = (2*pi/n)*perispecdiff(by);  % quadr wei
    [xq yq wq] = curveareaquad(bx,by,wx,wy,m);       % areal quadrature 
    kirchfac = 1/(1i*lambdaz);   % Kirchhoff prefactor for direct Fresnel
    ud = kirchfac * sum(exp((1i*pi/lambdaz)*((xq-xi0).^2+(yq-eta0).^2)) .* wq);
    eds(q) = abs(ud-uref0);                          % direct err
    u = fresnaq_pts(xq, yq, wq, lambdaz, xi,eta, tols(1));
    efs(q) = abs(u(1)-uref0);                        % fresnaq t3 err (tol 1)
    u = fresnaq_pts(xq, yq, wq, lambdaz, xi,eta, tols(2));
    efs2(q) = abs(u(1)-uref0);                       % fresnaq t3 err (tol 2)
    uf2{q} = u;
    [u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tols(1));
    egs(q) = abs(u(1,1)-uref0);             % fresnaq grid t1 err (tol 1)
    [u xigrid] = fresnaq_grid(xq, yq, wq, lambdaz, ximax, ngrid, tols(2));
    egs2(q) = abs(u(1,1)-uref0);            % fresnaq grid t1 err (tol 2)
    ug2{q} = u;
  end
  efss=nan*ns; egss=efss;                   % self-conv sup over targs
  for q=1:numel(ns);                        % gather self errs
    egss(q)=norm(ug2{q}(:)-ug2{end}(:),inf);
    efss(q)=norm(uf2{q}(:)-uf2{end}(:),inf);
  end
  subplot(2,3,2+figoff);
  semilogy(ns,eds,'b.-', ns,efs,'k+:', ns,efs2,'k+-', ns,efss,'r+--', ns,egs,'go:',ns,egs2,'go-', ns,egss,'ro--');
  xlabel('$n$','interpreter','latex');
  ylabel('abs. error in $u$','interpreter','latex');
  axis([ns(1) ns(end-1) 1e-15 1]);
  h=legend('direct','t3 $\varepsilon{=}10^{-6}$','t3 $\varepsilon{=}10^{-12}$','t3 self $\varepsilon{=}10^{-12}$','t1 $\varepsilon{=}10^{-6}$','t1 $\varepsilon{=}10^{-12}$','t1 self $\varepsilon{=}10^{-12}$');
  set(h,'interpreter','latex','location','southwest')
  h.Position = h.Position + [0.007 0.009 0 0];
  title(sprintf('(%c) $$\\lambda z=%.3g$$: AQ error vs $$n$$ (fix $$m=%d$$)',let,lambdaz,m),'interpreter','latex')
  drawnow
  let = let+1;

  subplot(2,3,3+figoff);    % plot last grid output field...
  u = 1-u;
  imagesc(xigrid,xigrid,abs(u)'.^2);
  colormap(hot(256)); colorbar; hold on; plot([bx,bx(1)],[by,by(1)],'w-');
  axis xy equal tight;
  xlabel('\xi'); ylabel('\eta');
  title(sprintf('(%c) $$\\lambda z=%.3g$$: t1 occulter intensity $$|u^{oc}|^2$$',let,lambdaz),'interpreter','latex')
  hold on; plot(xi0,eta0,'g.','markersize',20);
  caxis([0 1.6]);   % gets most of the high end, not all
  drawnow
  let = let+1;

end   % ............................................... end lambdaz loop

set(gcf,'paperposition',[0 0 14 9])
%print -depsc2 ../kite.eps


% biggest source-targ dist r:
%max(sqrt((xi0-bx).^2 + (eta0-by).^2))      ans =       2.98958740566711

% R = max(sqrt((-bx).^2 + (-by).^2))   1.13292275505978


% report some NSLI timings:
ngrid = 3e2;       % less pts!
M=ngrid^2; xi = ximax*(2*rand(M,1)-1); eta = ximax*(2*rand(M,1)-1);   % arb pts
fresnums = [10 100]; ns = [360 2400];
for i=1:numel(fresnums), fr=fresnums(i); fprintf('Fresnel # = %.3g...\n',fr)
  n = ns(i);
  lambdaz = 1/fr;
  t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);        % LI bdry nodes
  wx = (2*pi/n)*perispecdiff(bx); wy = (2*pi/n)*perispecdiff(by);  % quadr wei
  tic
  u = nsli_pts(bx,by,wx,wy, lambdaz, xi,eta);
  t = toc;
  fprintf('NSLI (n=%d): %d targs in %.3g s = %.3g src-targ pts/sec \n',n,M,t,n*M/t)
end

%Fresnel # = 10...
%NSLI (n=360): 90000 targs in 2.44 s = 1.33e+07 src-targ pts/sec 
%Fresnel # = 100...
%NSLI (n=2400): 90000 targs in 6.74 s = 3.2e+07 src-targ pts/sec 
