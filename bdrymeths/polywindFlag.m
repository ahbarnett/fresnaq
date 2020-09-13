%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SISTER: Starshade Imaging Simulation Toolkit for Exoplanet Reconnaissance
%
% Copyright (c) <2019>, California Institute of Technology ("Caltech"). U.S. Government sponsorship acknowledged.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% • Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%
% • Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%
% • Neither the name of Caltech nor its operating division, the Jet Propulsion Laboratory, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Eric Cady. eric.j.cady@jpl.nasa.gov
% Based on function polywind from NR, but vectorized
% One change: input vectors will already by closed, so we just run the full
% loop rather than starting at the "last" vertex.

function wind = polywindFlag(vt, pt, flag)

if nargin < 3 || flag == 0
    pt0 = pt(1);
    pt1 = pt(2);
    
    d0 = vt(2:end, 1);
    d1 = vt(2:end, 2);
    p0 = vt(1:end-1, 1);
    p1 = vt(1:end-1, 2);
    
    windp = (p1 <= pt1) & (d1 >  pt1) & (((p0 - pt0).*(d1-pt1)-(p1-pt1).*(d0-pt0)) > 0);
    windn = (p1 > pt1) & (d1 <= pt1) & (((p0 - pt0).*(d1-pt1)-(p1-pt1).*(d0-pt0)) < 0);
    wind = sum(windp) - sum(windn);
elseif flag == 1 % all points known a priori to be inside occulter.
    wind = 1;
end
