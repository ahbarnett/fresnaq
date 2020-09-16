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
function [ UTotL lambdaIn vecPetalArray pupil opt ] = makeStarshadeImage( opt_in )

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

% Settings for saving fields
useSave = opt.save ; % 1 save, 0 don't
savePath = opt.save_path ;
  if ~isdir( savePath ), mkdir( savePath ) ; end
% Skipping the simulation if it is saved and does not need to be re-done
if useSave == 1
  if ( ~opt.redo ) && ( exist( [ savePath '/' opt.save_filename '.mat' ] ) == 2 )
  disp( sprintf( '(makeStarshadeImage) PSF file %s%s exists. Skipping.', savePath, opt.save_filename ) )
  load( [ savePath '/' opt.save_filename '.mat' ] )
  return
  end
  if ( opt.redo ) && ( exist( [ savePath '/' opt.save_filename '.mat' ] ) == 2 )
  disp( sprintf( '(makeStarshadeImage) Simulation %s exists, but re-doing it.', opt.save_filename ) )
  end
end

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
% Step 2: Build edge
% -------------------------------

% Build edge; this function is overkill for an unaberrated edge but will do
% the job
tic
  % Case without perturbations
  vecPetalArray = NaN ;
  % The number of petals can be changed if opt.n_ptl differs from numPetals in the occulter file and the attenuation profile 'r' is available (e.g., NI2.mat, NW2.mat and TV3.mat have these values defined in the mat file)
  if exist( 'r', 'var' )
  % Some checks of consistency since the file with the apodization profile usually contains other basic parameters
    if isfield( opt, 'diameter_telescope_m' ) && exist( 'telescopeDiameter', 'var' )
      if ( opt.diameter_telescope_m ~= telescopeDiameter )
      disp( sprintf( 'WARNING: the value of the diameter of the primary mirror set in the running configuration file (%2.2f m) differs from the one found in the Starshade file (%2.2f m, %s)', opt.diameter_telescope_m, telescopeDiameter, [ opt.path_occulter opt.make_occulter_name '.mat' ] ) )
      end
    end

    if isfield( opt, 'n_ptl' ) && exist( 'numPetals', 'var' )
      if ( opt.n_ptl ~= numPetals )
      disp( sprintf( 'WARNING: the number of petals set in the running configuration file (%i) differs from the one found in the Starshade file (%i, %s)', opt.n_ptl, numPetals, [ opt.path_occulter opt.make_occulter_name '.mat' ] ) )
      end
    end

    if isfield( opt, 'secondary_size' ) && exist( 'secondarySize', 'var' )
      if ( opt.secondary_size ~= secondarySize )
      disp( sprintf( 'WARNING: the size of the secondary set in the running configuration file (%2.2f) differs from the one found in the Starshade file (%2.2f, %s)', opt.secondary_size, secondarySize, [ opt.path_occulter opt.make_occulter_name '.mat' ] ) )
      end
    end

    if isfield( opt, 'distance_starshade_telescope_m' ) && exist( 'Z', 'var' )
      if ( round( opt.distance_starshade_telescope_m ) ~= round( Z ) ) % 1 meter precision for the comparison between them
      disp( sprintf( 'WARNING: the distance of the Starshade to the telescope in the running configuration file (%2.2f) differs from the one found in the Starshade file (%2.2f, %s)', opt.distance_starshade_telescope_m, Z, [ opt.path_occulter opt.make_occulter_name '.mat' ] ) )
      end
    end

    if exist( 'occulterGeoIWA' , 'var' ) && exist( 'Z', 'var' ) && exist( 'occulterDiameter', 'var' )
      if ( occulterGeoIWA * Z ~= occulterDiameter / 2 )
      disp( sprintf( 'The geometric IWA in the Starshade file %2.2f mas is incompatible with the distance between the Starshade and the telescope %2.2f m and the diameter of the Starshade %2.2f m. Check the inout file. Stopping', occulterGeoIWA * 180 / pi * 3600e3, Z, occulterDiameter ) ) 
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end

    if isfield( opt, 'geo_iwa_mas' ) && exist( 'occulterGeoIWA', 'var' )
      if ( round( opt.geo_iwa_mas ) ~= round( occulterGeoIWA * 180 / pi * 3600e3 ) ) % 1 mas precision for the comparison between them
      disp( sprintf( 'WARNING: the geometric IWA in use is %2.2f mas. This is fine if the wavelengths (%04.2f,%04.2f) nm have that value, instead of the one found in the Starshade file (%2.2f mas, %s)', opt.geo_iwa_mas, lambdaIn( 1 ) / nm, lambdaIn( end ) / nm, occulterGeoIWA * 180 / pi * 3600e3, [ opt.path_occulter opt.make_occulter_name '.mat' ] ) )
      end
    end

  % Petals
    vecPetalArray = createErrorProfileV1(r, Profile, occulterDiameter, petalLength, opt.n_ptl, {});
  t = toc ;
    if ( opt.verbose )
    disp( sprintf( 'createErrorProfileV1 took %3.2f seconds', t ) )
    end

  tic
  tmpxVals = [];
  tmpyVals = [];
  tmpzVals = [];
  for j = 1 : opt.n_ptl
      tmpxVals = [tmpxVals vecPetalArray{j}{1}(1, :)]; 
      tmpyVals = [tmpyVals vecPetalArray{j}{1}(2, :)]; 
      tmpzVals = [tmpzVals vecPetalArray{j}{1}(3, :)]; 
  end
  xVals = [tmpxVals tmpxVals(1)];
  yVals = [tmpyVals tmpyVals(1)];
  zVals = [tmpzVals tmpzVals(1)];
  else
  % Given xVals, yVals
  load( [ opt.path_occulter opt.make_occulter_name ] )
    if ~exist( 'zVals', 'var' )
    zVals = 0 * xVals ;
    end
  disp( sprintf( 'Starshade petal shape loaded directly from %s', opt.make_occulter_name ) )
  % Sometimes the x/y/zVals are given as #x1 arrays instead of 1x#
    if size( xVals, 1 ) ~= 1
    xVals = xVals' ;
    end
    if size( yVals, 1 ) ~= 1
    yVals = yVals' ;
    end
    if size( zVals, 1 ) ~= 1
    zVals = zVals' ;
    end
  end

  % Closing the polygon, if it's open
  if (xVals(1) ~= xVals(end)) || (yVals(1) ~= yVals(end))
  xVals = [xVals(:); xVals(1)];
  yVals = [yVals(:); yVals(1)];
  end


% Simple rotation in the XY plane
  if ( opt.starshade_rotation_rad )
  disp( sprintf( '*** Rotating the Starshade an angle %d deg', opt.starshade_rotation_rad * 180 / pi ) ) ;
  alph_rt = opt.starshade_rotation_rad ;
  xVals_nw = xVals * cos( alph_rt ) + yVals * sin( alph_rt ) ;
  yVals_nw = -xVals * sin( alph_rt ) + yVals * cos( alph_rt ) ;
  xVals = xVals_nw ;
  yVals = yVals_nw ;
  clear xVals_nw yVals_nw
  end

% xVals, yVals, zVals give the 3D edge locus
% Under most circumstances zVals will be all zeros
t = toc ;
  if ( opt.verbose )
  disp( sprintf( 'xyzVals took %3.2f seconds', t ) )
  end

%--------------------------
% Step 3: Compute field at 
%--------------------------

% Propagate to telescope aperture plane with line integral
tic
  UTotL = bdwf(xVals, yVals, zVals, opt.dst_strshd_tlscp_m, lambdaIn, opt.dmtr_tlscp_m/Nx, Nx, r_src * mas, psi_src * degree, deltaX, deltaY);
t = toc;
  if ( opt.verbose )
  disp( sprintf('bdwf took %3.2f seconds, %3.2f per wavelength bin', t, t / n_lmbd ) )
  end

  if useSave == 1
  % These fields do not belong to this computation and may vary later on
    if isfield( opt, 'Nx_image_pix' ), opt = rmfield( opt, 'Nx_image_pix' ) ; end
    if isfield( opt, 'Nx_img' ), opt = rmfield( opt, 'Nx_img' ) ; end
    if isfield( opt, 'diam_img_mas' ), opt = rmfield( opt, 'diam_img_mas' ) ; end
  % Avoiding issues with the opt
  opt_make = opt ;
  sv_fl_nm = sprintf( '%s%s.mat', savePath, opt.save_filename ) ;
  save( sv_fl_nm, 'lambdaIn', 'opt_make', 'pupil', 'UTotL' ) 
  disp( sprintf( '(makeStarshadeImage) PSF file %s created.', sv_fl_nm ) )
  end

%disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
