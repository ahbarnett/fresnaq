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
function [ dst_m geo_iwa_mas ] = set_starshade_distance_and_geometric_iwa( opt, mssg )
% Setting the distance between the Starshade and the telescope
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% No messages by deafult
  if ~exist( 'mssg', 'var' )
  mssg = 0 ;
  end
  
  % If it is set
  if isfield( opt, 'dst_strshd_tlscp_m' )
  dst_m = opt.dst_strshd_tlscp_m ;
  end
  if isfield( opt, 'geo_iwa_mas' )
  geo_iwa_mas = opt.geo_iwa_mas ;
  end

% Case where both geometric IWA and distance have been set by the user
  if isfield( opt, 'dst_strshd_tlscp_m' ) && isfield( opt, 'geo_iwa_mas' )
  dst_m = opt.dst_strshd_tlscp_m ;
  geo_iwa_mas = opt.geo_iwa_mas ;
    if ( mssg )
    disp( sprintf( 'Distance between the telescope and the Starshade is %f km', dst_m / 1e3 ) ) ;
    disp( sprintf( 'The geometric IWA is %2.2f mas', geo_iwa_mas ) )
    end
  return
  end

% If it belongs to some hard coded instrument
  if isfield( opt, 'lmbd_psf_1_nm' )
  lambda_1_nm = opt.lmbd_psf_1_nm ;
  else
  lambda_1_nm = opt.lambda_1_nm ;
  end
  
  if isfield( opt, 'lmbd_psf_2_nm' )
  lambda_2_nm = opt.lmbd_psf_2_nm ;
  else
  lambda_2_nm = opt.lambda_2_nm ;
  end

% WFIRST. There are 3 bands
  if strfind( lower( opt.starshade.nominal_filename ), 'ni2' )
  % For the default band: 425-552 nm, the distance is 3.724225668350351e+07 m, geometric IWA=72 mas
    if ( lambda_1_nm >= 425 ) && ( lambda_2_nm <= 552 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 3.724225668350351e+07 ;
    geo_iwa_mas = 72 ;
    end
  % Longer wavelength, closer distance the distance is 2.612037672633376e+07 m, greater geometric IWA, 102.6571 mas
    if ( opt.lambda_1_nm >= 606 ) && ( opt.lambda_2_nm <= 787 ) && ( opt.lambda_1_nm <= opt.lambda_2_nm )
    dst_m = 3.724225668350351e+07 * ( 425 + 552 ) / ( 606 + 787 ) ;
    geo_iwa_mas = 102.6571 ; % 72 * ( 606 + 787 ) / ( 425 + 552 )
    end
  % The distance is 2.119142969119565e+07, geometric IWA=126.5343 mas
    if ( lambda_1_nm >= 747 ) && ( lambda_2_nm <= 970 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 3.724225668350351e+07 * ( 425 + 552 ) / ( 747 + 970 ) ;
    geo_iwa_mas = 126.5343 ; % 72 * ( 747 + 970 ) / ( 425 + 552 )
    end
  % Probe study (notice the == to avoid duplication with the 606-787 nm band)
    if ( lambda_1_nm == 615 ) && ( lambda_2_nm == 800 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 3.724225668350351e+07 * ( 425 + 552 ) / ( 615 + 800 ) ;
    geo_iwa_mas = 104 ; % 72 * ( 615 + 800 ) / ( 425 + 552 )
    end
  % Other bands
    if ( lambda_1_nm == 450 ) && ( lambda_2_nm == 600 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 3.724225668350351e+07 * ( 425 + 552 ) / ( 450 + 600 ) ;
    geo_iwa_mas = 77.3797 ; % 72 * ( 450 + 600 ) / ( 425 + 552 )
    end
    if ( lambda_1_nm == 590 ) && ( lambda_2_nm == 770 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 3.724225668350351e+07 * ( 425 + 552 ) / ( 590 + 770 ) ;
    geo_iwa_mas = 100.2252 ; % 72 * ( 590 + 770 ) / ( 425 + 552 )
    end
    if ( lambda_1_nm == 770 ) && ( lambda_2_nm == 1000 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 3.724225668350351e+07 * ( 425 + 552 ) / ( 770 + 1000 ) ;
    geo_iwa_mas = 130.4401 ; % 72 * ( 770 + 1000 ) / ( 425 + 552 )
    end

    if ~exist( 'dst_m', 'var' )
    disp( sprintf( 'Walength is not within one of the three expected bands: 425-552 nm, 606-787 nm, 747-970 nm. The minimum wavelength set is %04f nm, and the maximum wavelength is %04f nm', lambda_1_nm, lambda_2_nm ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  end

% HabEx (2018 version). There is one default band defined on these sims
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'nw2' ) )
  % For the default band: 300-1000 nm, the distance is 1.197666616918624e+08 m, geometric IWA=62 mas
    if ( lambda_1_nm >= 300 ) && ( lambda_2_nm <= 1000 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 1.197666616918624e+08 ; % m
    geo_iwa_mas = 62 ; % mas ( radius Starshade 36m/1.197666616918624e+08m=62 mas)
    end
    if ~exist( 'dst_m', 'var' )
    disp( sprintf( 'Walength is not within the expected band: 300-1000 nm. The minimum wavelength set is %04f nm, and the maximum wavelength is %04f nm', lambda_1_nm, lambda_2_nm ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  end

% HabEx (2019 version). There is one default band defined on these sims
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'tv3' ) )
  % For the default band: 300-1000 nm, the distance is 7.661264232035007e+07 m, geometric IWA=70 mas
    if ( lambda_1_nm >= 300 ) && ( lambda_2_nm <= 1000 ) && ( lambda_1_nm <= lambda_2_nm )
    dst_m = 7.661264232035007e+07 ; % m
    geo_iwa_mas = 70 ; % mas ( radius Starshade 26m/7.661264232035007e+07=70 mas)
    end
    if ~exist( 'dst_m', 'var' )
    disp( sprintf( 'Walength is not within the expected band: 300-1000 nm. The minimum wavelength set is %04f nm, and the maximum wavelength is %04f nm', lambda_1_nm, lambda_2_nm ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  end

% If no distance was set, stop
  if ~exist( 'dst_m', 'var' )
  disp( sprintf( 'The distance between the Starshade and the telescope is not defined. Set it with opt.distance_starshade_telescope_m. The occulter in use is ''%s''. Stopped', opt.starshade.nominal_filename ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

  if ~exist( 'geo_iwa_mas', 'var' )
  disp( sprintf( 'The geometric IWA is not defined. Set it with opt.geo_iwa_mas. The occulter in use is ''%s''. Stopped', opt.starshade.nominal_filename ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  
  if ( mssg )
  disp( sprintf( 'Distance between the telescope and the Starshade is %f km', dst_m / 1e3 ) ) ;
  disp( sprintf( 'The geometric IWA is %2.2f mas', geo_iwa_mas ) )
  end

