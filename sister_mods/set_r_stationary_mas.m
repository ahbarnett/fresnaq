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
function r_mas = set_r_stationary_mas( opt, mssg ) ;
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% The non-stationary PSF study is extended up to some distance from the center of the Starshade for WFIRST (Starshade of 13m radius), HabEx (36m radius), ...
% NI2: 
% 425-552 nm: 150 mas
% 606-787 nm: 220 mas
% 615-800 nm: 220 mas
% 747-970 nm: 270 mas
% NW2:
% 300-1000 nm: 62 mas
% TV3:
% 300-1000 nm: 70 mas

% Generic debug mode
dbstop if error

% If it is already defined, return
  if isfield( opt, 'r_st_mas' )
  r_mas = opt.r_st_mas ;
  return
  end
  if isfield( opt, 'r_stationary_mas' )
  r_mas = opt.r_stationary_mas ;
  return
  end

% By default, no messaging
  if ~exist( 'mssg', 'var' )
  mssg = 0 ;
  end

% Dummy value
r_mas = 0 ; 

  if ~isfield( opt.starshade, 'nominal_filename' )
  disp( 'Occulter unidentified. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

  if ~isfield( opt, 'sbdr_psf' )
  opt.sbdr_psf = get_subdir_psf( opt ) ;
  end

q_tmp = strfind( opt.sbdr_psf, '_' ) ;
lambda_1 = str2num( opt.sbdr_psf( 1 : min( q_tmp ) - 1 ) ) ; 
lambda_2 = str2num( opt.sbdr_psf( min( q_tmp ) + 1 : q_tmp( 2 ) - 1 ) ) ;

% Notice the introduction of some allowance in the wavelength limits of the band to deal with interpolation cases at the edges
  if strfind( lower( opt.starshade.nominal_filename ), 'ni2' )
    if ( lambda_1 == 425 ) && ( lambda_2 == 552 )
    r_mas = 150 ;
    end
    if ( lambda_1 == 606 ) && ( lambda_2 == 787 )
    r_mas = 220 ;
    end
    % WFIRST-Starshade probe study
    if ( lambda_1 == 615 ) && ( lambda_2 == 800 )
    r_mas = 220 ;
    end
    if ( lambda_1 == 747 ) && ( lambda_2 == 970 )
    r_mas = 270 ;
    end
    % Additional bands
    if ( lambda_1 == 450 ) && ( lambda_2 == 600 )
    r_mas = 150 ;
    end
    if ( lambda_1 == 590 ) && ( lambda_2 == 770 )
    r_mas = 220 ;
    end
    if ( lambda_1 == 770 ) && ( lambda_2 == 1000 )
    r_mas = 270 ;
    end
  end % NI2
 
  if numel( strfind( lower( opt.starshade.nominal_filename( 1 : 3 ) ), 'nw2' ) )
  % Preliminary. NW2 has a geometric IWA of 62 mas (Starshade of 18 m radius)
    if ( lambda_1 == 300 ) && ( lambda_2 == 1000 )
    r_mas = 130 ; % 62/72*150 = 130
    end
  end % NW2

  if numel( strfind( lower( opt.starshade.nominal_filename( 1 : 3 ) ), 'tv3' ) )
    if ( lambda_1 == 300 ) && ( lambda_2 == 1000 )
    r_mas = 150 ; 
    end
  end % TV3


% If something went wrong
  if ~( r_mas ) 
  disp( sprintf( 'No radius for the stationary region was assigned for the input wavelengths %04.2f, %04.2f nm for the occulter name ''%s''. Set opt.r_stationary_mas or check set_r_stationary_mas.m if this occulter should have a default value for the stationary region.', lambda_1, lambda_2, opt.starshade.nominal_filename ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  
  if ( mssg )
  disp( sprintf( 'The PSF stationary regions begins at %2.2f mas', r_mas ) )
  end
