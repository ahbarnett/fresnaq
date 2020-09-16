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

function opt = get_default_options( opt )
% Function to define the default options for Starshade imaging
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% Make options case insensitive.
opt = lower_opt( opt ) ;

% Verbose (a way to control when displaying messages)
  if ~isfield( opt, 'verbose' )
  opt.verbose = 0 ;
  end

% Saving some memory. In most applications, 7 digit precision is enough
  if ~isfield( opt, 'single_precision' )
  opt.sngl_prcsn = 1 ; % By default, store products in single precision
  else
  opt.sngl_prcsn = opt.single_precision ;
  end

% Switch to control some work that needs to be re-done or may be skipped if it already exists
  if ~isfield( opt, 'redo' )
  opt.redo = 0 ;
  end

% Occulter path
  if ~isfield( opt, 'path_occulter' )
  opt.path_occulter = './' ;
  end

% Occulter Name
% Nominal Starshade file
  if ~isfield( opt, 'starshade' )
  opt.starshade.nominal_filename = 'NI2' ;
  end

% Starshade mode: non-spinning/spinning
  if ~isfield( opt, 'starshade' )
  opt.starshade.mode = 'spinning' ;
  else
    if isfield( opt.starshade, 'mode' )
    opt.starshade.mode = lower( opt.starshade.mode ) ;
    else
    opt.starshade.mode = 'spinning' ;
    end
  end

% Check
  if ~strcmp( opt.starshade.mode, 'non-spinning' ) && ~strcmp( opt.starshade.mode, 'spinning' )
  disp( 'Either choose opt.starshade.mode=''spinning'' or opt.starshade.mode=''non-spinning''' )
  return
  end

% Number of petals
  if ~isfield( opt, 'n_ptl' )
  opt.n_ptl = 24 ; % Default number (S. Shaklan email "NI2 24 petals files" 05/03/18)
  end

% Force cosnsistency
  if isfield( opt.starshade, 'number_of_petals' )
  opt.n_ptl = opt.starshade.number_of_petals ;
  end

% PS: make the following lines consistent with sister_imaging.m
% Replace by a file in FITS or Matlab format with an Nx x Nx array if you want a specific pupil
  if ~isfield( opt, 'pupil_filename' )
  opt.pupil_filename = [ opt.path_occulter 'pupil_D1Kpix_256.fits' ] ;
  end

% For now, if it is not NI2, or some perturbed locus of NI2, we will assume an ideal pupil (circularly symmetric)
  if ~strfind( opt.starshade.nominal_filename, 'NI2' )
  opt.pupil_filename = 'ideal' ;
  end

% Also for spinning starshades, we will consider an ideal pupil (otherwise the PSF basis cannot be 1-dimensional, since the telescope's pupil does not rotate whereas the starshade does.
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
  opt.pupil_filename = 'ideal' ;
  end

% One may redefine the pupil as an ideal circularly symmetric pupil at any time
  if isfield( opt, 'pupil_filename' )
    if strcmp( opt.pupil_filename, 'ideal' )
    opt.pupil_filename= '0' ;
    end
  end

% Size of the pupil data in pixels (square)
  if ~isfield( opt, 'nx_pupil_pix' )
  opt.nx_pupil_pix = 64 ;
  end

% Diameter of the primary mirror
  if isfield( opt, 'diameter_telescope_m' )
  opt.dmtr_tlscp_m = opt.diameter_telescope_m ;
  end

% Some well-known cases
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'ni2' ) )
  opt.dmtr_tlscp_m = 2.4 ;
  end
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'nw2' ) ) || numel( strfind( lower( opt.starshade.nominal_filename ), 'tv3' ) )
  opt.dmtr_tlscp_m = 4.0 ;
  end

% Check 
  if ~isfield( opt, 'dmtr_tlscp_m' )
  disp( 'It is necessary to set the telescope''s diameter with opt.diameter_telescope_m. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Size of the secondary if a circular obscuration is fine (no struts). See makeStarshadeImage.m
  if ~isfield( opt, 'secondary_size' )
  opt.secondary_size = 0 ; % (Linear obscuration of the telescope entrance pupil (diameter ratio)
  end

% Some particular cases for the secondary obscuration (it will only have an effect if opt.pupil_filename='0', ideal pupil -see makeStarshadeImage.m). The values below must be the same as in sister_imaging.m
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'ni2' ) ) 
  opt.secondary_size = 0.417257964788079 ; ; % WFIRST value is 0.32 but accounting for the struts, which cover almost an 7.17% of the collecting area, this is the effective value of an equivalent secondary without struts. % https://wfirst.ipac.caltech.edu/sims/Param_db.html#telescope 
  end
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'tv3' ) )
  opt.secondary_size = 0.1125 ; % Table 6.2-2, page 6-4 of https://www.jpl.nasa.gov/habex/pdf/HabEx-Final-Report-Public-Release-LINKED-0924.pdf. No information about struts, but it should be minor compared to the effect of the secondary, which is already a small effect (1-0.1125^2)=0.987.
  end

% Number of pixels across image plane
  if ~isfield( opt, 'nx_image_pix' )
    if ~isfield( opt, 'nx_img' )
    opt.nx_img = 400 ;
    end
  else
  opt.nx_img = opt.nx_image_pix ;
  end

% Diameter of image plane in milliarcseconds
  if ~isfield( opt, 'diam_image_mas' )
  opt.diam_img_mas = 2001 ;
  else
  opt.diam_img_mas = opt.diam_image_mas ;
  end

  if ~isfield( opt, 'lambda_1_nm' ) && isfield( opt, 'lambda_band_nm_min' )
  opt.lambda_1_nm = opt.lambda_band_nm_min ;
  end

% Initial wavelength to consider
  if ~isfield( opt, 'lambda_1_nm' )
  opt.lambda_1_nm = 425 ; 
  end

  if ~isfield( opt, 'lambda_2_nm' ) && isfield( opt, 'lambda_band_nm_max' )
  opt.lambda_2_nm = opt.lambda_band_nm_max ;
  end

% Final wavelength to consider
  if ~isfield( opt, 'lambda_2_nm' )
  opt.lambda_2_nm = 552 ;
  end

% Step of wavelength to consider
  if ~isfield( opt, 'delta_lambda_nm' )
  opt.delta_lambda_nm = 10 ; % nm
  end
% Same variable but renamed without psf for some parts of the imaging simulation software
  if isfield( opt, 'delta_lambda_psf_nm' )
  opt.delta_lambda_nm = opt.delta_lambda_psf_nm ;
  end

% Spinning PSF basis. Number of starshade positions: 2*n_rot_half_petal+1.
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
    if ~isfield( opt, 'n_rot_half_petal' )
    % This default step of the rotations is such that it produces a rotation of the tip of the petal as large as half the pixel size set to build the PSF. That should be enough: the pixel size used to sample the PSF spatially is our limit of spatial resolution.
    opt.n_rot_hlf_ptl = round( ( pi / opt.n_ptl ) / ( opt.px_psf_mas / opt.geo_iwa_mas / 2 ) ) ; % =6, for 24 petals and a rotation step of opt.px_psf_mas(=3)/opt.geo_iwa_mas(=72)/2 (NI2), same for TV3. Total of 13 starshade positions (2*6+1, +/-rotations+0)
    else
    opt.n_rot_hlf_ptl = opt.n_rot_half_petal ;
    end
  end

% Distance Starshade telescope
  if ~isfield( opt, 'distance_starshade_telescope_m' )
  % Some default cases that depend on some well known configurations
  opt.dst_strshd_tlscp_m = set_starshade_distance_and_geometric_iwa( opt, opt.verbose ) ;  % ahb
  else
  opt.dst_strshd_tlscp_m = opt.distance_starshade_telescope_m ;
  end

% X/Y positions on the focal plane (instead of r_source_mas and psi_source_deg)
% Matlab: column major. Therefore in an image, X/Y are to be understood as -y/x in usual Cartesain convention (e.g., X=1, Y=0 means (0,-1) in sual x/y axis)
% Convention ox X/Y with respect to r/Psi as in bdwf core function
% s1 = sin(psi1);
% c1 = cos(psi1);
% s2 = sin(psi2);
% c2 = cos(psi2);

% Checks of consistency
  if isfield( opt, 'r_source_mas' ) && ~isfield( opt, 'psi_source_deg' )
  disp( '(get_default_options) r_source_mas set, but psi_source_deg not. Inconsistent. Returning.' )
  return
  end
  if isfield( opt, 'psi_source_deg' ) && ~isfield( opt, 'r_source_mas' )
  disp( '(get_default_options) psi_source_deg set, but r_source_mas not. Inconsistent. Returning.' )
  return
  end

% Separation of the source from the center of the pointing
  if ~isfield( opt, 'r_source_mas' )
  opt.r_source_mas = 0 ; % mas
  end

% Angle of the source with respect the horizontal axis
  if ~isfield( opt, 'psi_source_deg' )
  opt.psi_source_deg = 0 ; % degrees
  end

% Consistency check
  if isfield( opt, 'y_source_mas' ) && ~isfield( opt, 'x_source_mas' )
  disp( '(get_default_options) If opt.y_source_mas is set, opt.x_source_mas must also be set. Returning.' )
  return
  end
  if isfield( opt, 'x_source_mas' ) && ~isfield( opt, 'y_source_mas' )
  disp( '(get_default_options) If opt.x_source_mas is set, opt.y_source_mas must also be set. Returning.' )
  return
  end

% Sentinel for changes in the filename where results are stored
rplc_fl_nm = 0 ;

  if isfield( opt, 'x_source_mas' )
  % Consistency check
    if ~isfield( opt, 'y_source_mas' )
    disp( '(get_default_options) If opt.x_source_mas is set, opt.y_source_mas must also be set. Returning.' ) 
    return
    end
  % Transform into r and psi
  % Different situations where it is required to check the consistency between r_source_mas, psi_source_deg, x_source_mas and y_source_mas:
    if ( opt.x_source_mas ) || ( opt.y_source_mas ) 
      if ( ( sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ~= opt.r_source_mas ) || ( atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ~= opt.psi_source_deg ) ) 
      r_source_mas_tmp = sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ;
      psi_source_deg_tmp = atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ; % deg
        if ( opt.verbose )
        disp( sprintf( '(get_default_options) Changing the values of r_source_mas and psi_source_deg from %3.3f, %3.3f to %3.3f, %3.3f', opt.r_source_mas, opt.psi_source_deg, r_source_mas_tmp, psi_source_deg_tmp ) )
        end
      opt.r_source_mas = r_source_mas_tmp ;
      opt.psi_source_deg = psi_source_deg_tmp ;
      % Update filename where results are stored
      rplc_fl_nm = 1 ;
      end
    end
  end
  
  %% If they are not fields, create them
    if ~isfield( opt, 'r_source_mas' ) 
    opt.r_source_mas = sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ;
    end
    %%% Recall matlab column major convention and atan2(y,x)
    if ~isfield( opt, 'psi_source_deg' )
      if opt.r_source_mas % if the radius is zero, the angle does not matter.
      opt.psi_source_deg = atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ; % deg
      else
      opt.psi_source_deg = 0 ;
      end
    end

% For ease identification of the pixel on the image
opt.x_source_mas = opt.r_source_mas * cos( opt.psi_source_deg * pi / 180 ) ;
opt.y_source_mas = opt.r_source_mas * sin( opt.psi_source_deg * pi / 180 ) ;

% Hot pixels in the electric fields (by default, do not erase them)
  if ~isfield( opt, 'erase_hot_pixels' )
  opt.erase_hot_pixels = 0 ;
  end

% NB: Name of the filename to store the results at the end of the code
% Saving all the output results and images (0=No, 1=Yes)
  if ~isfield( opt, 'save_all' )
  opt.save_all = 0 ;
  end

% Saving only output results
  if ~isfield( opt, 'save' )
  opt.save = 0 ;
  end
% Saving figures
  if ~isfield( opt, 'save_fig' )
  opt.save_fig = 0 ;
  end

% If all saved, then:
  if opt.save_all
  opt.save = 1 ;
  opt.save_fig = 1 ;
  end

% paths to save the results
  if ~isfield( opt, 'save_path' )
  opt.save_path = './output/' ;
  end

% Saving the figures
  if ~isfield( opt, 'save_path_fig' )
  opt.save_path_fig = './fig' ;
  end

% For the interpolation analysis
  if ~isfield( opt, 'star' )
  opt.star = 0 ;
  end

  if ~isfield( opt, 'planet' )
  opt.planet = 0 ;
  end

  if ~isfield( opt, 'polar' )
  opt.polar = 0 ;
  end

  if ~isfield( opt, 'super_resolution' )
  opt.super_resolution.res = 1 ;
  opt.super_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt.super_resolution, 'res' )
  opt.super_resolution.res = 1 ;
  end

  if ~isfield( opt.super_resolution, 'interp_method' )
  opt.super_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt, 'low_resolution' )
  opt.low_resolution.res = 2 ;
  opt.low_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt.low_resolution, 'res' )
  opt.low_resolution.res = 2 ;
  end

  if ~isfield( opt.low_resolution, 'interp_method' )
  opt.low_resolution.interp_method = 'linear' ;
  end

  % Number of exact simulations to derive interpolation results: for instance, 2*N+1, -N, -N+1, ..., -1, 0, 1, ..., N-1, N
  if ~isfield( opt, 'n_basis_interpolation' )
  opt.n_basis_interpolation = 5 ;
  end

  % Step in mas for a series of simulations
  if ~isfield( opt, 'step_mas' )
  opt.step_mas = 5 ;
  end
 
% Rotating Starshade
  if ~isfield( opt, 'starshade_rotation_rad' )
  opt.starshade_rotation_rad = 0 ;
  end

%%%%%%%%% Filename where to store the results %%%%%%%
  if ~isfield( opt, 'save_filename' ) || ( rplc_fl_nm )
  opt.save_filename = sprintf( 'starshade_UtotL_Nx_%i_pix_%04i_%04i_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', ...
  opt.nx_pupil_pix, opt.lambda_1_nm, opt.lambda_2_nm, opt.delta_lambda_nm, opt.r_source_mas, opt.psi_source_deg ) ;
  % For delta_lambda_nm that is not an integer
    if opt.delta_lambda_nm ~= round( opt.delta_lambda_nm )
    opt.save_filename = sprintf( 'starshade_UtotL_Nx_%i_pix_%04i_%04i_dl_%3.1fnm_dr_%3.1f_mas_psi_%3.1f_deg', ...
    opt.nx_pupil_pix, opt.lambda_1_nm, opt.lambda_2_nm, opt.delta_lambda_nm, opt.r_source_mas, opt.psi_source_deg ) ;
    end
    if strcmp( opt.pupil_filename, '0' ) == 1, opt.save_filename = sprintf( '%s_ideal', opt.save_filename ) ; end
    if opt.diam_img_mas ~= 2001, opt.save_filename = sprintf( '%s_diam_%04i', opt.save_filename, opt.diam_img_mas ) ; end
  % Rotating Starshade
    if ( opt.starshade_rotation_rad ~= 0 )
    str_alph = strrep( sprintf( '%1.4f', opt.starshade_rotation_rad ), '.', 'p' ) ;
    opt.save_filename = sprintf( '%s_alpha_%s', opt.save_filename, str_alph ) ; 
    end
  % Different nominal starshade file than NI2
    if ~strcmp( opt.starshade.nominal_filename, 'NI2' )
    opt.save_filename = sprintf( '%s_%s', opt.save_filename, opt.starshade.nominal_filename ) ;
    end
% Not used for now.
%  if opt.nx_img ~= 400, opt.save_filename = sprintf( '%s_Nx_img_%04i', opt.save_filename, opt.nx_img ) ; end
  end

