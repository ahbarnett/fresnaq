function [Afunc r0 r1] = get_NI2_Afunc
% GET_NI2_FUNC   apodization profile func for NI2 starshade from NI2.mat object
%
%  [Afunc r0 r1] = get_NI2_Afunc returns a function handle to A(r) valid in
%   the apodization domain [r0,r1].
%
%  Note: apodization function A(r) is in notation of E. Cady Opt. Expr. 2012,
%  starting at 1 for r<r0 and becoming 0 for r>r1.

% Barnett 9/4/20

cwd = fileparts(mfilename('fullpath'));
load([cwd '/../occulter/NI2.mat']);
r1 = occulterDiameter/2;      % upper apodization radius in meters
r0 = r1 - petalLength;        % lower apodization radius in meters

%tol = 3e-3;   % v sens; really shouldn't use global approx on equispaced grid
inds = (r>=r0);
x = r(inds); x = x(:); y = Profile(inds);   % x=ordinates, y=data, as col vecs
% fit function on [r0,r1]...  this is strange data since not that smooth
p = 180;   % any higher -> wildy oscillates - serves me right for Cheby on grid
[Afunc p res] = globalinterp1d(r(inds),Profile(inds),r0,r1,[],p);
fprintf('%d data pts, res=%.3g @ degree p=%d\n',sum(inds),res,p)
t = linspace(r0,r1,5e3);
figure; subplot(2,1,1);
plot(x,y,'+'); hold on; plot(t,Afunc(t),'-'); set(gca,'ylim',[-1 2])
subplot(2,1,2); plot(x,y-Afunc(x),'.-'); title('residual on samples');
% observe around 4-5 digits before global approx gets too osc and blows up.
% -> switch to local.
