% MATLAB setup for Frensel apertures. Main sources/tests in top level.
% Barnett 9/4/20
format long g
h = fileparts(mfilename('fullpath'));
addpath(genpath(h))                        % gives access to all subdirs
% don't forget to add finufft/matlab to path
