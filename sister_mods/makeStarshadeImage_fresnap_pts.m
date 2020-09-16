% [ u lambdaIn vecPetalArray pupil opt ] = makeStarshadeImage_fresnaq_pts( opt_in, xi, eta)
%
% Hacked version of makeStarshadeImage, with two additional inputs (xi,eta), overriding target list
% Does not save to file, rather, returns u, a Ntarg * n_lmbd complex array of diffraction field at targets.
% Uses fresnaq_pts, after building areal quadrature, under the hood. No rotation is done, because
% targets should be rotated leaving the occulter source the same.
% Only ideal pupil works. The occulter file must have "r" field.
%
% I have not documented this file since there are many aspects I do not understand. However,
% please see the docs for FRESNAQ_PTS and STARSHADEQUAD.

% Barnett 9/15/20.






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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT H/home/alex/physics/starshade/SISTER/softwareOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u lambdaIn vecPetalArray pupil opt ] = makeStarshadeImage_fresnaq_pts( opt_in, xi, eta)

% makeStarshadeImage
% A sample program which creates a locus of edge point from an apodization
% function and propagates the resulting field to a telescope focal plane.
% History:
% 4/25/17: first complete version, Eric J Cady (eric.j.cady@jpl.nasa.gov, JPL/Caltech)
% Modified: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com. Modification of the interface, pupil handling, further interface modifications (lambda range), UTotL computed separately, and several soft and I/O changes to adapt it to sister.m

% Check some options are being passed
  if ~exist( 'opt_in', 'var' )
  disp( '(makeStarshadeImage) Provide some options: opt.x=y ; makeStarshadeImage( opt ) ; returning ...' )
  return
  end

% Get default options (look inside the function for specific definitions)
opt = get_default_options( opt_in ) ;    % ahb: this seems to create a whole bunch of new opt fields :(

% Main parameters
Nx = opt.nx_pupil_pix ;
dlt_lmbd = opt.delta_lambda_nm ;
r_src = opt.r_source_mas ;
psi_src = opt.psi_source_deg ;
ppl_fl = opt.pupil_filename ;

% ahb gutted the prevention of re-saving

%---------------------------
% Step 1: Load up starshade
%---------------------------
units

% Load occulter file. It may be the on-axis or off-axis occulter. It must be set before calling makeStarshadeImage.m
occulterName = [ opt.path_occulter opt.make_occulter_name '.mat' ] ;
load( occulterName ); % Load up the comparison occulter
disp( sprintf( '(makeStarshadeImage) Using: %s', occulterName ) )

% Range of wavelengths to be produced
lambdaIn = nm * ( opt.lambda_1_nm + dlt_lmbd * ( 0 : 1 : ( opt.lambda_2_nm - opt.lambda_1_nm ) / dlt_lmbd ) ) ;
  if lambdaIn( end ) ~= ( nm * opt.lambda_2_nm ) 
  lambdaIn( end + 1 ) = lambdaIn( end ) + dlt_lmbd * nm ;
  end
n_lmbd = numel( lambdaIn ) ;
disp( sprintf( 'Considering %i wavelengths', numel( lambdaIn ) ) ) 
% Build/load
if strcmp( ppl_fl, '0' )
    secondarySize = opt.secondary_size ; % If you're fine with a generic circular pupil with a secondary (but no struts), set this to the fraction of the radius the secondard covers (e.g. 0.2)
    pupil = makePupil( Nx, Nx, 1, secondarySize, 0, 0 ) ;
disp( sprintf( 'Considering a circular pupil with a secondary blocking a fraction of the radius=%f', secondarySize ) ) ;
else
  if strfind( ppl_fl, '.fits' )            % ahb cant' run this since uses Imaging Toolbox...
    pupil = fitsread( ppl_fl );
    disp( [ 'Reading pupil''s file: ' ppl_fl ] ) 
    else
    load( ppl_fl ) ;
    % Test from NG
    pupil = pupil_mask ;
    disp( [ 'Reading pupil''s file: ' ppl_fl ] )
    end
    n_ppl = sqrt( numel( pupil ) ) ;
    % Check
    disp( sprintf( 'Reference pupil size is %ix%i pixels.', n_ppl, n_ppl ) ) 
    % Reducing the size of the pupil grid (fast<1s)
      if Nx ~= n_ppl
      % finding where the pupil starts
      q_ppl_1 = min( find( pupil ~= 0 ) ) ;
      q_ppl_2 = max( find( pupil ~= 0 ) ) ;
      % Column where the pupil data starts and ends
      clmn_1 = floor( q_ppl_1 / n_ppl ) ;
        if clmn_1 ==0, clmn_1 = 1 ; end
      clmn_2 = floor( q_ppl_2 / n_ppl ) ;
      pupil = imresize( pupil( clmn_1 : clmn_2, clmn_1 : clmn_2 ), Nx / ( clmn_2 - clmn_1 + 1 ), 'bilinear' ) ; 
      % Sharping the edges
      pupil( find( pupil < .7 ) ) = 0 ;
      pupil( find( pupil ~= 0 ) ) = 1 ;
      disp( sprintf( 'New pupil size is reduced to %ix%i pixels.', Nx, Nx ) )
      end
end

% Lateral offsets in meters
deltaX = 0;
deltaY = 0;

% -------------------------------
% Step 2: Build edge                 --------------- ahb gutted
% -------------------------------
  vecPetalArray = NaN ;

%--------------------------
% Step 3: Compute field at 
%--------------------------

% ahb inserts setup for areal quadr needed by fresnaq: assumes r,Profile in occulterName file
[~,r0,r1] = eval_sister_apod(occulterName,0);   % get apodization range [r0,r1]
Afunc = @(r) eval_sister_apod(occulterName,r);  % func handle (reads file when called)
Np = opt.n_ptl;
n = 40; m = 300;  % NI2 areal quadr params (should go in occulter file) ... need to check converged (see eg demo_starshades)
verb = 1;
[xq yq wq] = starshadequad(Np,Afunc,r0,r1,n,m,verb);   % ahb: fill areal quadr
% Propagate to telescope aperture plane with fresnaq, each lambda in turn...
tol = 1e-7;   % accuracy of Fourier bit
t0=tic;
u = nan(numel(xi), numel(lambdaIn));
for l=1:numel(lambdaIn)
  fprintf('lambda=%.3g um...\n',lambdaIn(l)*1e6);
  u_aper = fresnaq_pts(xq,yq,wq, lambdaIn(l)*Z, xi,eta, tol);  % diffract to all targs at once
  u(:,l) = 1-u_aper;    % Babinet and stack
end
t=toc(t0);
fprintf('fresnaq_pts (%d sources to %t targets, %d lambdas) took %3.2f seconds, %3.2f per wavelength bin', numel(xq),numel(xi),n_lmbd, t, t / n_lmbd)
% ahb done

% ahb gutted the save bit, moved to sister_basis hacked.

%disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
