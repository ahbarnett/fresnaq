function Afunc = points2apodization(x,y,r0,r1,Npetal)
% POINTS2APODIZATION  Fit apodization profile A(r) from xy points on half-petal
%
% Afunc = points2apodization(x,y,r0,r1,Npetal)
%  reverse-engineers the apodization profile A(r) over the transition radius
%  domain [r0,r1], from the points with x,y coordinates on the upper half of the
%  zeroth petal, of a symmetric starshade with Npetals.
%  Returns a function handle.
%  Note: apodization function A(r) is in notation of E. Cady Opt. Expr. 2012,
%  starting at 1 for r<r0 and becoming 0 for r>r1.

% Barnett 9/4/20

