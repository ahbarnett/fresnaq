function E = nsli_emulates_bdwf(xVals, yVals, zVals, Z, lambda, dxO, nO, psi1, psi2, deltaX, deltaY)
% NSLI_EMULATES_BDWF  drop-in replacement for BDWF using NSLI & inc dir transl
%
% E = nsli_emulates_bdwf(xVals, yVals, zVals, Z, lambda, dxO, nO, psi1, psi2,
%  deltaX, deltaY) sets up target grid, loops over lambda, and calls NSLI for
%  each Fresnel scalar diffraction. It uses 2nd-order accurate locus quadrature,
%  thus does not achieve the potential (high-order) accuracy of NSLI.
%  Incident wave direction (psi1,psi2) is emulated to accuracy O(psi1^2) via
%  target translations.
%  No nonzero zVals are allowed. Arguments are in SISTER format, matching BDWF.
%
% See bdwf.m for full documentation of what is computed and arguments.
%
% Without arguments, tests this routine against BDWF on smooth occulter
%
% Also see: NSLI 

% Barnett 9/14/20
if nargin==0, test_nsli_emulates_bdwf; return; end

if ~isempty(zVals) || sum(zVals~=0.0), warning('zVals~=0 not implemented!'), end
[xi,eta] = make_grid_bdwf(dxO, nO, deltaX, deltaY);    % make target grid
xi = xi - Z*psi1*cos(psi2);                % emulate off-axis wave by targ shift
eta = eta - Z*psi1*sin(psi2);              % "
[wx wy] = crudecurvequad(xVals,yVals);     % make 2nd-ord accurate locus weights
Nl = numel(lambda);
E = nan(nO,nO,Nl);
for l=1:Nl
  u = nsli_pts(xVals,yVals, wx,wy, lambda(l)*Z, xi, eta);
  E(:,:,l) = exp(2i*pi*Z/lambda(l)) * (1-u);  % plane z phase, Babinet->occulter
end

%%%%%%%%%%
function test_nsli_emulates_bdwf          % some code from test_bdwf
x = @(t) 0.5*cos(t)+0.5*cos(2*t); y = @(t) sin(t);   % smooth kite shape, Reff~1
lambda = linspace(4e-7,5e-7,10);          % some wavelengths (meters)
Z = 2.5e5;                  % downstream distance (meters) ... good for Reff~1m
psi1 = 4e-7; psi2 = pi/5;   % incident wave dir (psi1.Z must be < O(1))
dxO = 0.1; nO = 10; deltaX = 0.23; deltaY = -0.16;   % generic target grid
Reff=1; FN = Reff^2./(lambda*Z);
fprintf('Fresnel # range = [%.3g,%.3g] (Reff=1) \t psi1.Z=%.3g\n', min(FN),max(FN), psi1*Z)

ns = [1e3 3e3 1e4];        % convergence study: #s of boundary quadrature points
                       % depends on max FN & target pt; expect 1/n^2 for NSLI
for n=ns
  t = 2*pi*(0:n-1)/n; bx = x(t); by = y(t);      % make bdry points
  % compare the two methods...
  tic; Eb = bdwf(bx,by,[], Z, lambda, dxO, nO, psi1,psi2, deltaX,deltaY); %toc
  tic; En = nsli_emulates_bdwf(bx,by,[], Z, lambda, dxO, nO, psi1,psi2, deltaX,deltaY); %toc
  fprintf('n=%d\t max E diff:  %.3g (ampl only) \t (%.3g including phase)\n', n, max(abs(abs(Eb(:))-abs(En(:)))), max(abs(Eb(:)-En(:))))
end

% Notes:

% 1) We know NSLI is O(1/n^2) accurate using simple trap rule weights;
% test_bdwf shows BDWF is only O(1/n) accurate. Explains difference only O(1/n)

% 2) since kz ~ 1e13, overall phase poor, kills all but 3 digits in phase!

% 3) strangely BDWF for "smaller" n (1e3) gets v bad (err>0.1) for psi1>=1e-6
% why?
