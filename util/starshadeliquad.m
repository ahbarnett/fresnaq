function [bx by wx wy] = starshadeliquad(Np,Afunc,Apfunc,r0,r1,m,mtip,verb,Aquad)
% STARSHADELIQUAD  vector line integral quadr, apodized petal 0-1 occulter
%
% [bx by wx wy] = starshadeliquad(Np,Afunc,Apfunc,r0,r1,m,mtip,verb,Aquad)
%
%  Uses Theta(r) formula of (1)-(2) in Cady '12 to build ideal starshade then
%  vector line integral quadrature scheme, given apodization function and
%  other geometric parameters.
%  Has option to override the Afunc by given (r_j,A(r_j)) samples.
%
% Inputs:
%  Np = # petals
%  Afunc = if a function handle: apodization function from roughly 1 down to
%          roughly 0, over domain [r0,r1].
%          If a struct, overrides function with given fields r and A containing
%          point samples of A at each radius r: Apfunc, m, Aquad, ignored.
%          Just the values in [r0,r1] are used.
%  Apfunc = function handle for derivative of apodization A'(r) (optional;
%          if empty, weights not computed).
%  r0,r1 = apodization domain of radii (note differs from Cash'11 defn a,b)
%   radial apodization is 1 for r<r0, A((r-r0)/(r1-r0)) for r0<r<r1, 0 for r>r1.
%  m = requested nodes over radial apodization [r0,r1]. (exact may vary)
%  mtip = nodes per arc over each tip and gap (eg 0 or 16).
%  verb = (optional) 0,1, etc, verbosity.
%  Aquad = (optional) apodization quadrature type: 'a' Alpert (default),
%          'g' Gauss, 'u' Uniform grid (including end points).
%
% Outputs:
%  bx,by - boundary quadrature nodes for vector line integrals
%  wx,wy - (optional) boundary quadrature weights for vector line integrals,
%          if Apfunc empty or Afunc is a struct of samples.
%
% Design is perfectly symmetric; no perturbations for now.

% Barnett 9/28/20
if nargin==0, test_starshadeliquad; return; end
if nargin<8, verb = 0; end
if nargin<9, Aquad = 'a'; end
if r0>r1, error('r0 must be <= r1'); end
calcw = ~isempty(Apfunc) && ~isstruct(Afunc);      % weights?

if ~isstruct(Afunc)     % func handle case...
  switch Aquad          % set up a [0,1] quadr scheme
   case 'g'
    [z w] = lgwt(m,0,1);               % radial petal quadrature, descending z
   case 'a' % or use Alpert, which may make couple more nodes than requested...
    [z w] = QuadNodesInterval(0,1,m,[],1,1,32); m=numel(z); z=z(end:-1:1); w=w(end:-1:1);  % flip to descending order
   case 'u'   % uniform grid, including endpoints, composite trap rule...
    z = linspace(0,1,m); h=z(2)-z(1); z=z(end:-1:1)'; w=h*[.5;ones(m-2,1);.5];
   otherwise
    error('unknown Aquad!')
  end
  r = r0 + (r1-r0)*z(:);                            % r quadr nodes
  wr = (r1-r0)*w(:);
  A = Afunc(r);                                  % get all apodizations at once
  if calcw, Ap = Apfunc(r); else Ap=nan*A; end   % get all A'(r) values
  eps = r1*2e-16;
  hgap = (pi/Np)*(1-Afunc(r0+eps)); htip = (pi/Np)*Afunc(r1-eps);  % half-angs
else
  ii = Afunc.r>=r0 & Afunc.r<=r1;
  r = Afunc.r(ii); A = Afunc.A(ii); m = numel(r);     % override w/ samples
  r = r(:); A = A(:);
  hgap = (pi/Np)*(1-max(A)); htip = (pi/Np)*min(A);   % half-angs (assumes r well-sampled)
end
[r i] = sort(r,'ascend'); A = A(i); if calcw, Ap = Ap(i); wr = wr(i); end  % col
if verb>1, fprintf('starshadeliquad hgap=%.3g htip=%.3g\n',hgap,htip); end
% for scaled tip or gap theta quadr, ascending...
[t wt] = lgwt(mtip,-1,1); t = t(end:-1:1); wt = wt(end:-1:1);
% lower edge of petal, r param ascend...
th = -pi/Np*A;         % theta(r)
c = cos(th); s = sin(th); bx = r.*c; by = r.*s;
if calcw, thp = -pi/Np*Ap;                    % theta'(r)
  wx = wr.*(c - r.*thp.*s); wy = wr.*(s + r.*thp.*c); end
% tip, t param ascend...
th = htip*t;          % theta(t)
c = cos(th); s = sin(th);  bx = [bx; r1*c]; by = [by; r1*s];
if calcw, wx = [wx; -htip*r1*s.*wt]; wy = [wy; htip*r1*c.*wt]; end
% upper edge of petal, r param descend...
r = r(end:-1:1); th = pi/Np*A(end:-1:1);         % theta(r)
c = cos(th); s = sin(th); bx = [bx; r.*c]; by = [by; r.*s];
if calcw, thp = +pi/Np*Ap(end:-1:1);                    % theta'(r)
  wx = [wx; -wr.*(c - r.*thp.*s)]; wy = [wy; -wr.*(s + r.*thp.*c)]; end   % sign
% gap, t param ascend...
th = pi/Np + hgap*t;          % theta(t)
c = cos(th); s = sin(th);  bx = [bx; r0*c]; by = [by; r0*s];
if calcw, wx = [wx; -hgap*r0*s.*wt]; wy = [wy; hgap*r0*c.*wt]; end
n = numel(bx);
for p=1:Np-1                          % loop copying to remaining petals...
  t = 2*pi/Np*p; c = cos(t); s = sin(t); R = [c -s;s c];
  bb = R*[bx(1:n),by(1:n)]';         % rotate pts
  bx((1:n)+p*n) = bb(1,:)'; by((1:n)+p*n) = bb(2,:)';   % extract nodes
  if calcw
    ww = R*[wx(1:n),wy(1:n)]';       % rotate wei
    wx((1:n)+p*n) = ww(1,:)'; wy((1:n)+p*n) = ww(2,:)';  % extract wei
  end
end


%%%%%%%%
function test_starshadeliquad
verb = 1;
designs = {'erf','NI2'};
for d=1:numel(designs)
  switch designs{d}
   case 'erf'                         % my analytic toy model
    Np=16; r0 = 7; r1 = 14;           % # petals, inner, outer radii in meters
    beta = 3.0;                       % 3.0 good for A gap and 1-A tip to 1e-5
    A = @(r) erfc(beta*(2*r-(r0+r1))/(r1-r0))/2;
    erfcp = @(x) (-2/sqrt(pi))*exp(-x.^2);    % deriv of erfc
    %x=1; erfcp(x), h=1e-5;(erfc(x+h)-erfc(x-h))/(2*h)  % test erfc'
    Ap = @(r) (beta/(r1-r0)) * erfcp(beta*(2*r-(r0+r1))/(r1-r0));
    %r=10; Ap(r), h=1e-5; (A(r+h)-A(r-h))/(2*h)   % test A'
    ms = 10:10:50;                      % m-conv vals
    n = 20;                           % (conv for area check)
   case 'NI2'                         % actual NI2 starshade, cubic interpolated
    cwd = fileparts(mfilename('fullpath')); file = [cwd '/../occulter/NI2'];
    [A,r0,r1,Ap] = eval_sister_apod(file);   % r1=13, and raw xVals goes to it.
    o = load(file); Np = o.numPetals;
    ms = 200:50:500;                      % m-conv vals
    n = 40;                           % (conv for area check)
  end
  m =ms(end); [xj yj wj] = starshadequad(Np,A,r0,r1,n,m);
  area = sum(wj);         % AQ converged area
  Aquad = 'g'; mtip = 16;
  fprintf('LI sum & area err (true:%.8g) for design %s (Aquad %s)...\n',area,designs{d},Aquad)
  if exist('file','var'); o = load(file);
    bxr = o.xVals(1:end-1); byr = o.yVals(1:end-1);   % raw bdry pts
    [wxr wyr] = crudecurvequad(bxr, byr);
    liarea0 = 0.5*sum(bxr.*wyr - byr.*wxr);      % area via vec LI
    fprintf('crude wei from raw: area error vs AQ (m=%d) = %.3g\n',m,liarea0-area)
    % following not so useful since raw r,Profile only goes to 12.998 ...
    %s = []; s.r = o.r; s.A = o.Profile;         % grab raw A(r) samples.
    %[bx0 by0] = starshadeliquad(Np,s,[],r0,r1,0,0);     % no weights poss
    %[wx0 wy0] = crudecurvequad(bx0, by0);
    %liarea0 = 0.5*sum(bx0.*wy0 - by0.*wx0);      % area via vec LI
    %fprintf('crude wei from sslq: area err vs AQ (m=%d) = %.3g\n',m,liarea0-area)
  end
  
  for m=ms
    [bx by wx wy] = starshadeliquad(Np,A,Ap,r0,r1,m,mtip,verb,Aquad);
    sumw = abs(sum(wx))+abs(sum(wy));       % basic test for vec LI wei
    liarea = 0.5*sum(bx.*wy - by.*wx);      % area via vec LI
    disp([m sumw liarea-area])
  end
    figure; plot(wx,wy,'.-'); drawnow
  if verb, figure(d); clf; scatter(xj,yj,10,wj,'.'); hold on;
    plot([bx;bx(1)],[by;by(1)],'k+-');    % our output
    if exist('bxr','var'), plot([bxr bxr(1)],[byr byr(1)],'r.'); end  % raw x,y
    axis equal tight; colorbar; title(sprintf('%s',designs{d})); end
end
