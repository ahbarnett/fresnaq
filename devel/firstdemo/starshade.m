% test starshade occulter with 0-1 scalar diffraction Fresnel scheme.
% Barnett 8/24/20
clear; close all;

verb = 1;
fbar = 80;   % "increased" Fresnel #, fbar = 2pi.f = k_free/d = 2pi/(lambda.d)
Np = 16;        % # petals
r1 = 0.6; r0 = 1.5;
beta = 3.0; A = @(t) erfc(2*beta*(t-0.5))/2;     % my own apodization for now
%A = @(t) exp(-(t/0.6).^6);        % Cash'11 hyper-Gaussian in [0,1], also ok

fprintf('1-A(0)=%.3g, A(1)=%.3g\n',1-A(0),A(1))     % check its decay
if verb>1, figure(1); clf; r = linspace(0,1,1e3); plot(r,A(r),'-'); end

%n = 30; m = 100;   % converged for fbar = 60
n = 40; m = 120;   % converged for fbar = 80 (f=12.7)
tic; [xj yj wj bx by] = starshadequad(Np,A,r1,r0,n,m,verb); toc
if verb>1, figure(2); clf; scatter(xj,yj,10,wj); axis equal tight; colorbar;
  hold on; plot([bx;bx(1)], [by;by(1)], 'k-'); drawnow
end

xmax = 1.5;       % targ params: box half-size
M1 = 1e3;         % target location grid size per dim (eg, <=50 to trigger mod)
tol = 1e-9;
xt = -xmax; yt=-xmax;   % test pt (make sure it's on the targ grid)
% good to be far out since that's worst-case.

if verb
  figure(3); clf; subplot(1,2,1); plot(xj,yj,'.'); hold on; axis equal tight;
  plot(xt,yt,'*'); xlabel('x_1'); ylabel('x_2'); title('Aperture + nodes');
  subplot(1,2,2); scatter(xj,yj,10,cos(0.5*fbar*((xj-xt).^2+(yj-yt).^2)));
  axis equal tight; xlabel('x_1'); ylabel('x_2'); title('Re Fresnel integrand');
end

[u x1] = fresnelgrid(fbar, xj, yj, wj, xmax, M1, tol);         % diffract!

if verb, figure(4); clf;
  imagesc(x1,x1,log10(abs(1.0-u)'.^2)); caxis([-11 0.2]); colorbar; hold on;
  plot([bx;bx(1)], [by;by(1)], 'k-'); axis xy equal tight;
  title(sprintf('log_{10} |u|^2 for occulter: Fresnel # %6.3f',fbar/2/pi));
  %plot(xt,yt,'*');
end

% check it...
kirchfac = fbar/(2i*pi);   % Kirchhoff approx prefactor = 1/i.lambda.d
ut = kirchfac * sum(exp(0.5i*fbar*((xj-xt).^2+(yj-yt).^2)) .* wj)
j1 = find(x1==xt); j2 = find(x1==yt); errt = u(j1,j2) - ut;
fprintf('abs error vs slow Fresnel quad at (%.3g,%.3g) = %.3g\n\n',xt,yt,abs(errt))

if 0 % animation over fbar...
M1 = 900; tic;
fresnels = linspace(8,13,100);   % Fresnel numbers
fbars = 2*pi*fresnels;
nam = 'fres_sweep';
writerObj = VideoWriter(sprintf('%s.avi',nam),'Uncompressed AVI');
writerObj.FrameRate = 15;     % movie, FPS
open(writerObj);  % start writing to a movie
fig = figure; set(fig,'position',[700 0 1000 1000]);
set(gca,'position',[.05 0.05 0.9 0.9]);
for l=1:numel(fbars), fbar = fbars(l); fprintf('fbar=%.3g...\n',fbar)
  [u x1] = fresnelgrid(fbar, xj, yj, wj, xmax, M1, tol);         % diffract!
  imagesc(x1,x1,log10(abs(1.0-u)'.^2)); caxis([-11 0.2]); colorbar; hold on;
  plot([bx;bx(1)], [by;by(1)], 'k-');
  title(sprintf('log_{10} intensity:   Fresnel # = %6.3f',fbar/(2*pi)));
  axis xy equal tight;
  drawnow; writeVideo(writerObj,getframe(fig));
end
toc % too 26 sec
close(writerObj);          % finish movie file
system(sprintf('ffmpeg -i %s.avi -an -y -c:v libx264 -crf 20 -r 30 %s.mp4'))
% fails due to GLIBC in the env
% Have to call from another shell:
% ffmpeg -i fres_sweep.avi -an -y -c:v libx264 -crf 20 -r 30 fres_sweep.mp4
end
