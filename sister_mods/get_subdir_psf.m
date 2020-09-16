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

function sbdr_psf = get_subdir_psf( opt )
% Function to get the subdirectory associated with the E fields (non-spinning) of the PSF basis (spinning)
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% Checking for cases of sub-bands that may be part of more than one imaging band
sngl_bnd = 0 ;

% If lambda_band_mn_min or max are not defined, but lambda_1_nm and lambda_2_nm are (creating a basis)

  if ~isfield( opt, 'lambda_band_nm_min' ) && ~isfield( opt, 'lambda_band_nm_max' ) && isfield( opt, 'lambda_1_nm' ) && isfield( opt, 'lambda_2_nm' )
  opt.lambda_band_nm_min = opt.lambda_1_nm ;
  opt.lambda_band_nm_max = opt.lambda_2_nm ;
  end

  n_bnd = numel(  opt.lambda_band_nm_min ) ;

    for i_bnd = 1 : n_bnd
    % The array of wavelengths may be defined when either creating a basis or when imaging a scene.
      if isfield( opt, 'lmbd_arry_scn_nm' )
        if ( size( opt.lmbd_arry_scn_nm, 1 ) == 1 ) || ( size( opt.lmbd_arry_scn_nm, 2 ) == 1 )
        lmbd_lcl_1 = opt.lmbd_arry_scn_nm( 1 ) ;
        lmbd_lcl_2 = opt.lmbd_arry_scn_nm( end ) ;
        else
        lmbd_lcl_1 = opt.lmbd_arry_scn_nm( i_bnd, 1 ) ;
        lmbd_lcl_2 = opt.lmbd_arry_scn_nm( i_bnd, end ) ;
        end
      elseif isfield( opt, 'lambda_imaging_1_nm' )
      lmbd_lcl_1 = opt.lambda_imaging_1_nm ;
      lmbd_lcl_2 = opt.lambda_imaging_2_nm ;
      elseif isfield( opt, 'lambda_1_nm' )
      lmbd_lcl_1 = opt.lambda_1_nm ;
      lmbd_lcl_2 = opt.lambda_2_nm ;
    end
    % + 10 is to accomodate a final wavelength used for interpolation at the upper edge of the wavelength array and avoid extrapolation (this is due to the fact that delta_lambda_nm does not divide the bandwitdh in an integer number of steps, and it is fine).
      if (  opt.lambda_band_nm_min( i_bnd ) - 3 <= lmbd_lcl_1 ) && ( lmbd_lcl_2 <=  opt.lambda_band_nm_max( i_bnd ) + 10 )
      sbdr_psf = sprintf( '%i_%i_nm/',  opt.lambda_band_nm_min( i_bnd ),  opt.lambda_band_nm_max( i_bnd ) ) ;
      sngl_bnd = sngl_bnd + 1 ;
      end
    end
    if ~sngl_bnd
    disp( sprintf( 'Imaging wavelength band (%04.2f,%04.2f) nm incompatible with any imaging band in %s. Stopped', lmbd_lcl_1, lmbd_lcl_2, opt.starshade.nominal_filename ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

    if ( sngl_bnd > 1 )
    disp( sprintf( 'WARNING: %i instrument bands are compatible with the wavelength range selected (%04.2f,%04.2f) nm. Defaulting to the directory %s', sngl_bnd, lmbd_lcl_1, lmbd_lcl_2, sbdr_psf ) )
    end

% Check
  if ~exist( 'sbdr_psf', 'var' )
  disp( sprintf( 'The location for the PSF basis has not been identified. Check the setup and the imaging bands of your instrument (%s). Stopped.', opt.starshade.nominal_filename ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

