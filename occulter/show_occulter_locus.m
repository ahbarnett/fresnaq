function [h xVals yVals] = show_occulter_locus(filehead)
% SHOW_OCCULTER_LOCUS  plot (xVals,yVals) from a occulter MAT file
%
%  [h xVals yVals] = show_occulter_locus(filehead)
%   plots on the current figure the boundary points ("locus") in occulter file
%   of the format used in SISTER, and returns the line object handle, and
%   locus coordinates from that file.

% Barnett 9/8/20
o = load(filehead);
xVals = o.xVals; yVals = o.yVals;
hold on;
h = plot([xVals xVals(1)], [yVals yVals(1)], 'b.','markersize',10); axis equal
