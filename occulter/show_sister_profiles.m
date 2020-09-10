% plot all the SISTER profiles I can find in locus files, plus 2nd derivs,
% in their [r0,r1] radius ranges. Barnett 9/9/20

clear
S = '/home/alex/physics/starshade/SISTER/input_scenes/locus/in/';
fs = {'NI2', 'NW2', 'TV3', 'UH17'};

nfs = numel(fs);
figure; p=get(gcf,'position'); set(gcf,'position',[p(1) p(2)-600 500 1000]);
for i=1:nfs
  file = [S fs{i}];
  o = load(file);
  r1 = o.occulterDiameter/2;      % upper apodization radius in meters
  r0 = r1 - o.petalLength;        % lower apodization radius in meters
  inds = (o.r>=r0);               % kill off useless inner radii
  x = o.r(inds); x = x(:); y = o.Profile(inds); % x=ordinates, y=data, col vecs
  subplot(nfs,1,i);
  ypp = [0;diff(y,2);0];          % discrete 2nd deriv assuming unif samples
  plot(x, [y, 0.5*ypp/max(ypp)],'.','markersize',10); xlabel('r (m)');
  vline([r0 r1]);
  title(sprintf('%s: Np=%d, 1-A(%.3g)=%.3g, A(%.3g)=%.3g',fs{i},o.numPetals,r0,1-y(1),r1,y(end)))
  axis tight; legend('raw Profile samples', 'scaled discr 2nd deriv');
  % optionally overlay the interpolated func...
  %A = @(r) eval_sister_apod(file,r);
  %r = linspace(r0,r1,1e4); hold on; plot(r,A(r),'r-');
end

%print -dpng all_sister_profiles.png
