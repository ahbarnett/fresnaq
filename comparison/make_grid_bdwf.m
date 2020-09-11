function [xx yy] = make_grid_bdwf(dxO, nO, deltaX, deltaY)
% MAKE_GRID_BDWF  output (x,y) coords of grid that bdwf internally assumes
%
% [xx yy] = make_grid_bdwf(dxO, nO, deltaX, deltaY)
%  returns the product of the
%          1D grids:  xi = -deltaX + ((-nO+1)/2 : (nO-1)/2) * dxO and
%                    eta = -deltaY + ((-nO+1)/2 : (nO-1)/2) * dxO
%  to match the grid internal to bdwf. It's like MATLAB ndgrid in 2D, in that
%  x is fast, y is slow.
%
% Inputs:
%   dxO - target grid spacing in both eta and xi directions, in meters
%   nO - target grid size in both directions
%   deltaX, deltaY - target grid negative translation in eta, xi, in meters.
% Outputs:
%   xx, yy - column vectors (length nO^2) of x and y coords of grid points.

% Barnett 9/10/20

% we pull lines from Cady's bdwf.m:
width = nO*dxO;
xO = -width/2+dxO/2:dxO:width/2-dxO/2;     % PS this is prone to rounding error
yO = xO - deltaY;
xO = xO - deltaX;

[xx yy] = ndgrid(xO,yO);
