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
function sister_basis( opt )
% Function to create the imaging basis to be used with sister.m
% It deals with the case of a non-rotating starshade and with the case os a spinning one.
%
% 1) In the non-spinning case, this function creates the complex electric fields associated with 
% any position of the image plane. This is the ingredient used in the sister code to 
% create an image.
%
% 2) In the spinning case, this function generates a PSF basis for each wavelength of interest 
% for all distances to the center of the starshade (notice that beyond the non-stationary region, 
% the PSF is translationally invariant and a single PSF serves for all r>r_non_stationary. See 
% set_r_stationary_mas for more details)
%
% NOTE: there's no need to run (1) before running (2). In fact, (2) checks whether the electric 
% fields corresponding to the non-spinning case exist, and if they don't, they will be created, 
% stored and used to build the spinning PSF basis. That is running (2), will run (1). However, if (1) 
% has already been run, (2) will use any existing files from (1).
%
% REMARK: remember to have updated the installation path in sister_installation_path.m before running this software.
%
% RECOMMENDATIONS: 
% 1) For creating basis in a personal computer, consider adding the option opt.px_psf_mas=3; to the options below. 
%
% 2) Sample the pupil with less pixels, 16x16, instead of the default 64x64.
%
% The software may run well with 64x64 and opt.px_psf_mas=1; on a personal computer, but it may run out of virtual memory for some 
% cases (red bands, for instance, that require larger PSF and more PSF elements to cover the band). Another option is to 
% modify opt.delta_lambda_nm, from 10 nm (default) to greater values. However, the accuracy will decrease if the band is not 
% properly sampled. Think about the wavelength dependence of your astrophysical or engineering scenario before modifying the 
% default settings.
%
% A) Examples for the non-spinning case:
% A1) All default options
% 	clear opt;opt.starshade.mode='non-spinning';sister_basis( opt )
% A2) Create a PSF basis for the non-spinning case and for the band [ 425, 552 ] nm, with a pupil sampled with 16x16 pixels 
% (rest of parameters get default values, like the delta lambda that is 10 nm)
% 	clear opt;opt.starshade.mode='non-spinning';opt.Nx_pupil_pix=16;opt.lambda_1_nm=425;opt.lambda_2_nm=552;sister_basis( opt )
% A3) Similar to (A2) but for an ideal pupil (implies circular symmetry of the pupil):
%	clear opt;opt.pupil_filename='ideal';opt.starshade.mode='non-spinning';opt.Nx_pupil_pix=16;opt.lambda_1_nm=425;opt.lambda_2_nm=552;sister_basis( opt )
%  
% B) Examples for the spinning case:
% B1) All default options
% 	clear opt;opt.starshade.mode='spinning';sister_basis( opt )
% B2) Create a PSF imaging basis for the spinning case and for the band [ 425, 552 ] nm, with a pupil sampled with 16x16 pixels 
% (rest of parameters get default values, like the delta lambda that is 10 nm)
% 	clear opt;opt.Nx_pupil_pix=16;opt.starshade.mode='spinning';opt.lambda_1_nm=425;opt.lambda_2_nm=552;sister_basis( opt )
% B3) Similar to (B2) but for an ideal pupil (implies circular symmetry of the pupil):
% 	clear opt;opt.pupil_filename='ideal';opt.Nx_pupil_pix=16;opt.starshade.mode='spinning';opt.lambda_1_nm=425;opt.lambda_2_nm=552;sister_basis( opt )
%
% C) Completely new occulter
% Run this example and follow the instructions int he running configuration file:
%  opt = new_occulter_non_kepler_config;sister_basis( opt );
%
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% Generic debug mode
  dbstop if error

% Development option to plot the different PSF that build up a spinning PSF (should be 0, unless interested in getting some videos)
opt.psf_video = 0 ;

  if ~exist( 'opt', 'var' )
  opt = [ ] ;
  end

% * Update the path to your local installation by editing sister_installation_path.m
  if ~isfield( opt, 'installation_path' )
  opt.installation_path = sister_installation_path() ;
  end

  if ~isfield( opt.starshade, 'mode' )
  disp( 'Choose either opt.starshade.mode=''non-spinning'' or opt.starshade.mode=''spinning''. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Get default options related with this code and with the computation of E fields, basis, ...
opt = get_local_default_options( opt ) ;

% Storing the bands to be used in the loop below
lambda_1_nm = opt.lambda_1_nm ;
lambda_2_nm = opt.lambda_2_nm ;

% Minimal check
n_bnd = numel( opt.lambda_1_nm ) ;
n_bnd_2 = numel( opt.lambda_2_nm ) ;
  if ( n_bnd ~= n_bnd_2 )
  disp( sprintf( 'The number of bands is inconsistent: %i or %i? Stopped', n_bnd, n_bnd_2 ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Reminding the setup
disp( sprintf( 'Creating an imaging basis for a %s starshade:', opt.starshade.mode ) )
disp( sprintf( 'Nominal Starshade filename: %s', opt.starshade.nominal_filename ) )
disp( sprintf( '# Petals: %i', opt.n_ptl ) ) 
disp( sprintf( 'Pupil with %ix%i pixels', opt.nx_pupil_pix, opt.nx_pupil_pix ) )
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
  disp( sprintf( 'Number of symmetries used to build the PSF basis: %i', opt.n_symm ) )
  end

wl_str = '' ;
  for i_bnd = 1 : n_bnd
  wl_str = sprintf( '%s%04.2f-%04.2f nm ', wl_str, opt.lambda_1_nm( i_bnd ), opt.lambda_2_nm( i_bnd ) ) ;
  end
wl_str = wl_str( 1 : end - 1 ) ;
  if n_bnd > 1
  disp( [ 'Bands: ' wl_str ] )
  else
  disp( [ 'Band: ' wl_str ] )
  end

% Preparing all jobs, for all bands
t_strt = tic ;
  for i_bnd = 1 : n_bnd
  opt.lambda_1_nm = lambda_1_nm( i_bnd ) ;
  opt.lambda_2_nm = lambda_2_nm( i_bnd ) ;
  % Final radial distance (in general, the one defining the transition from the non-stationary region to the stationary one)
    if ~isfield( opt, 'r_1_mas' )
    opt.r_1_mas = set_r_stationary_mas( opt ) ;
    end
  % Path where the E fileds are stored
    if ~isfield( opt, 'sbdr_psf' )
    opt.sbdr_psf = get_subdir_psf( opt ) ;
    end
  opt.save_path = sprintf( '%ssister_basis/non-spinning/%s_%i_%i/%s/', opt.installation_path, opt.starshade.nominal_filename, opt.n_ptl, opt.nx_pupil_pix, opt.sbdr_psf ) ;
    if ~isdir( opt.save_path )
    mkdir( opt.save_path ) ;
    disp( sprintf( '%s created', opt.save_path ) )
    end

    if strcmp( opt.starshade.mode, 'non-spinning' )
    run_non_spinning( opt )
    else
      error('ahb gutted spinning for simplicity!');
     end % End of spinning case
  end % i_bnd

  if n_bnd > 1
  disp( sprintf( 'Creation of the imaging basis for the bands %s took %3.2f min', wl_str, toc( t_strt ) / 60 ) )
  else
  disp( sprintf( 'Creation of the imaging basis for the band %s took %3.2f min', wl_str, toc( t_strt ) / 60 ) )
  end

%%%%%%%%%%%%%%%%
% Sub-functions
%%%%%%%%%%%%%%%%

% Running the non-spinning case
function run_non_spinning( opt )
opt.save = 1 ;
% For a non-ideal pupil (i.e., with struts) we need all the quadrant
  if ~strcmp( opt.pupil_filename, '0' ) && ~strcmp( opt.pupil_filename, 'ideal' )
    
    error('ahb gutted non-ideal pupil!');  % ----------------------------------------------
    
  else
  % For an ideal pupil, we only need to consider the positions within half a petal, including its axis. The reduction is about sin(alph_ptl_rd)/2, which usually is about 7-8 times less basis elements to derive than in the non-ideal case)
  % Angle corresponding to a petal
  alph_ptl_rd = 2 * pi  / opt.n_ptl ;
  % Maximum y_mas necessary to cover half a petal
  y_mx_mas = ceil( opt.r_1_mas * sin( alph_ptl_rd ) ) ;
  % Array of positions (direct product)
  x_lst_mas = opt.r_0_mas : opt.psf_spacing_mas : opt.r_1_mas ;
  idx_x_1 = 1 ;
    for idx_y = opt.r_0_mas : opt.psf_spacing_mas : y_mx_mas 
    alph_tmp_rd = atan2( idx_y, x_lst_mas ) ;
    q_alph = find( alph_tmp_rd <= alph_ptl_rd ) ;
    n_q_alph = numel( q_alph ) ;
      if n_q_alph ~= 0
      x_ptl_mas( idx_x_1 : idx_x_1 + n_q_alph - 1 ) = x_lst_mas( q_alph ) ;
      y_ptl_mas( idx_x_1 : idx_x_1 + n_q_alph - 1 ) = idx_y ;
      idx_x_1 = idx_x_1 + n_q_alph + 1 ;
      end
    end
    n_px_ptl = numel( x_ptl_mas ) ;
    t0 = tic;
  xi = []; eta = [];           % ahb started targ list
  sv_fl_nm_list = cell(n_px_ptl,1);   % ahb started filename list
    for i_px = 1 : n_px_ptl
    opt.x_source_mas = x_ptl_mas( i_px ) ;
    opt.y_source_mas = y_ptl_mas( i_px ) ;
    if opt.verbose, disp( sprintf( 'Considering the PSF at (%i,%i) mas', opt.x_source_mas, opt.y_source_mas ) ); end
    % but did opt.x_source_mas get used by makeStarshadeImage( opt ) ?? appears not, confusing
    % ah, they are used by get_default_options() called inside - ugh.
    
    % ahb now brings out just code from makeStarshadeImage that decides the set of targets...
    opt = get_default_options( opt ) ;   % fills many opt fields :(
    Nx = opt.nx_pupil_pix ;            % grid params
    dxO = opt.dmtr_tlscp_m/Nx;
    deltaX = 0;   % ahb not sure why these even needed if zero
    deltaY = 0;
    [xi0 eta0] = make_grid_bdwf(dxO, Nx, deltaX, deltaY);  % ahb. a set of Nx^2 targs
    r_src = opt.r_source_mas ;
    mas = pi/180/60^2/1e3;            % one milliarcsecond. no idea where this is defined in SISTER
    degree = pi/180;  % ditto
    psi_src = opt.psi_source_deg ;    % where is this set?
    psi1 = r_src * mas; psi2 = psi_src * degree;
    Z = opt.dst_strshd_tlscp_m;
    xi0 = xi0 - Z*psi1*cos(psi2);                % emulate off-axis wave by targ shift
    eta0 = eta0 - Z*psi1*sin(psi2);              % "
    xi = [xi; xi0(:)]; eta = [eta; eta0(:)];     % append to targ list (kinda slow)
    savePath = opt.save_path ;
    sv_fl_nm_list{i_px} = sprintf( '%s%s.mat', savePath, opt.save_filename ) ;
    end % i_px
    
    Ntargs = numel(xi);
    fprintf('sister_basis collected all %d targets (%d of %dx%d grids) in %.3g s\n',Ntargs, n_px_ptl, Nx,Nx,toc(t0))
    if opt.verbose, figure; plot(xi,eta,'.'); axis equal; title('all targs in sister_basis non-spinning'); drawnow; end
    % wasteful that most targets are essentially repeated many times :( Let's handle them all anyway...
    
    % ahb add build quadr & diffraction here...
    opt.make_occulter_name = opt.starshade.nominal_filename ;
    % makeStarshadeImage( opt ) ;     % ahb this reads all params from opt, writes to file (ugh). Why no output args used??
    t0=tic;
    [u,lambdaIn,~,pupil] = makeStarshadeImage_fresnap_pts( opt, xi, eta) ;   % ahb hacked overrides targets with simple (xi,eta) list
    % output u is size Ntargs*Nlambdas but Nlambdas is not yet defined of course...
    n_lmbd = size(u,2);
    t=toc(t0);
    disp( sprintf('makeStarshadeImage_fresnap_pts took %3.2f seconds, %3.2f per wavelength bin', t, t / n_lmbd ) )

    % make the pupil once (no need to redo each i_px)...
    secondarySize = opt.secondary_size ; % If you're fine with a generic circular pupil with a secondary (but no struts), set this to the fraction of the radius the secondard covers (e.g. 0.2)
    pupil = makePupil( Nx, Nx, 1, secondarySize, 0, 0 ) ;

    % finally, explicitly save each grid chunk of targ results...
    for i_px = 1:n_px_ptl
      targ_inds = (1:Nx^2) + (i_px-1)*Nx^2;
      UTotL = reshape(u(targ_inds,:), [Nx Nx n_lmbd]);             % as if bdwf did it
      % write out all lambdas in a Nx^2 grid to one file, grabbing code from end of makeSt...
      opt_make = opt ;
      save( sv_fl_nm_list{i_px}, 'lambdaIn', 'opt_make', 'pupil', 'UTotL' ) 
      if i_px==1 || i_px==n_px_ptl, disp( sprintf( '(makeStarshadeImage) PSF file %s created.', sv_fl_nm_list{i_px} ) ), end
    end    
    
  end % End of the case of an ideal pupil

  
  
% ahb: spinning case code removed for intelligibility ----------------------------------------
  



% Some global and local options
function opt = get_local_default_options( opt )

% Make options case insensitive
opt = lower_opt( opt ) ;

% Path for the occulter file
  if ~isfield( opt, 'path_occulter' )
  opt.path_occulter = [ opt.installation_path 'input_scenes/locus/in/' ] ;
  end

% PSF basis spacing: it should be 1 mas, regardless of the final pixel pitch of the PSF (Here the basis is created at this pixel pitch based upon analysis on absolute/relative changes of the PSF. Then, if a different pixel pitch is used when imaging scenes, the resampled PSF basis is created and stored for a faster turn around). We've seen that 1 mas keeps relative errors of the PSF below 1e-3 relative, and that may be good to deal with a wide variety of astrophysical scenarios (exo-dust emission, bright planets, broad band integration with several wavelength slices, ...).
  if ~isfield( opt, 'psf_spacing_mas' )
  opt.psf_spacing_mas = 1 ;
  end

% PSF pixel pitch
  if ~isfield( opt, 'px_psf_mas' )
  opt.px_psf_mas = 1 ; % mas
  end

% Wavelengths to consider (should be consistent with set_r_stationary_mas)
  if ~isfield( opt, 'lambda_1_nm' ) && isfield( opt, 'lambda_band_nm_min' )
  opt.lambda_1_nm = opt.lambda_band_nm_min ;
  end  

  if ~isfield( opt, 'lambda_1_nm' )
  opt.lambda_1_nm = [ 425, 606, 747 ] ;
  end

  if ~isfield( opt, 'lambda_2_nm' ) && isfield( opt, 'lambda_band_nm_max' )
  opt.lambda_2_nm = opt.lambda_band_nm_max ;
  end

  if ~isfield( opt, 'lambda_2_nm' )
  opt.lambda_2_nm = [ 552, 787, 970 ] ;
  end

% Starting radial distance (in general, it *should* be 0 mas)
  if ~isfield( opt, 'r_0_mas' )
  opt.r_0_mas = 0 ;
  end

  if ( opt.r_0_mas )
  disp( sprintf( 'The closest position for the PSF basis to the origin does not include the origin, but starts at %2.1f mas. Odd. Make sure this is what is intended.', opt.r_0_mas ) )
  end 

% Pixel resolution of the input PSF
  if ~isfield( opt, 'px_in_mas' )
  opt.px_in_mas = 1 ;
  end

% Pixel resolution on the final image in mas (best results are with 1 mas and then in sister the PSF gets reduced to the pixel resolution of the scene in an optimal way. If the scenes would have the same reolution all the time, it could be done here and save running time and memory, but that does not need to be always the case).
  if ~isfield( opt, 'px_psf_mas' )
  opt.px_psf_mas = 1 ;
  end

% Portion of the PSF stored as times l/D
  if ~isfield( opt, 'n_lambda_over_d' )
  opt.n_lambda_over_d = 7 ;
  end

% Get the general default options to compute the E fields, basis, etc ...
opt = get_default_options( opt ) ;

% For a spinning starshade
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
  % Array of angles to be considered. It depends on the overall symmetry. For a non-ideal starshade, one needs to cover all the circle. If it is ideal and struts are present, it will be necessary to deal with -90 deg to 90 deg in most cases. 
    if ~isfield( opt, 'n_symmetries' )
    opt.n_symm = opt.n_ptl ; % Set to 1 if it is non-ideal
    else
    opt.n_symm = opt.n_symmetries ;
    end
  end

% Conversion from lambda over D to mas (it depends on the telescope)
opt.lD2mas_fct = 1e-9  / opt.dmtr_tlscp_m * 180 / pi * 3600 * 1000 ;

% Number of pixels on the original PSF files 
% 1) Number of mas from the PSF center to derive the PSF (PSF=+/- this value. Derived for the longest wavelength in the band. Pixel pitch of the PSF spacing)
hlf_psf_arry_mas = ceil( ceil( opt.n_lambda_over_d * opt.lambda_2_nm * opt.lD2mas_fct ) / opt.px_psf_mas ) * opt.px_psf_mas ;

% 2) Actual number of pixels to sample the whole PSF: PSF is a square array of size n_px_psfxn_px_psf (an odd number easily sets the center of the PSF on the PSF basis files). 
  if ~isfield( opt, 'n_px_psf' )
  opt.n_px_psf = 2 * ceil( hlf_psf_arry_mas / opt.px_psf_mas ) + 1 ;
  end

% Check (if the modified value is not an odd number)
  if ( 2 * floor( opt.n_px_psf / 2 ) + 1 ) ~= opt.n_px_psf
  disp( 'The number of pixels sampling the PSF should be an odd number. Stopped.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Radius where the stationary region begins
opt.r_st = set_r_stationary_mas( opt ) ;
% Geometric IWA
[ dummy_dst_m opt.geo_iwa_mas ] = set_starshade_distance_and_geometric_iwa( opt ) ;

% Delete or preserve the individual electric fields created when generating the spinning PSF basis (it is only used when the PSF basis is for a spinning starshade). The electric fields are stored in sister_basis/non-spinning/"locus_name"/"band_name"/.
  if ~isfield( opt, 'delete_electric_fields' )
  opt.delete_electric_fields = 1 ; % By default, delete the electric fields after having generated the spinning PSF since they are not used anymore. Only for some technical PSF studies, it should be set to 0.
  end

% Checking whether the files for the spinning imaging exist (the ones that are arranged by wavelength)
function [ wrk_lmbd lmbd_arry n_lmbd fl_nm ] = check_spinning_sister_basis( opt )

% Range of wavelength considered
lmbd_arry = opt.lambda_1_nm : opt.delta_lambda_nm : opt.lambda_2_nm ;
  if lmbd_arry( end ) ~= opt.lambda_2_nm
  lmbd_arry( end + 1 ) = lmbd_arry( end ) + opt.delta_lambda_nm ;
  end
n_lmbd = numel( lmbd_arry ) ;
% Checking if there's any work to do
% Looping over wavelength
  for i_lmbd = 1 : n_lmbd
  wrk_lmbd( i_lmbd ) = 0 ;
  fl_nm_1 = sprintf( '%sstarshade_spinning_psf_Nx%i', [ opt.spinning_path '../' ], opt.nx_pupil_pix ) ;
    if isfield( opt, 'pupil_filename' )
      if ( strcmp( opt.pupil_filename, '0' ) || strcmp( opt.pupil_filename, 'ideal' ) )
      fl_nm_1 = [ fl_nm_1 '_ideal' ] ;
      end
    end
  fl_nm{ i_lmbd } = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%04d_nm.mat', fl_nm_1, opt.geo_iwa_mas, opt.px_psf_mas, lmbd_arry( i_lmbd ) ) ;
    % For delta_lambda_nm that is not an integer
    if lmbd_arry( i_lmbd ) ~= round( lmbd_arry( i_lmbd ) )
    fl_nm{ i_lmbd } = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%3.1f_nm.mat', fl_nm_1, opt.geo_iwa_mas, opt.px_psf_mas, lmbd_arry( i_lmbd ) ) ;
    end
  
    if ~strcmp( opt.starshade.nominal_filename, 'NI2' )
    fl_nm_tmp = fl_nm{ i_lmbd } ;
    fl_nm{ i_lmbd } = [ fl_nm_tmp( 1 : end - 4 ) '_' opt.starshade.nominal_filename '.mat' ] ;
    end

    if exist( fl_nm{ i_lmbd }, 'file' ) ~= 2 || ( opt.redo )
    wrk_lmbd( i_lmbd ) = 1 ;
    end
  end % i_lmbd

% Development functions

function fg_nm_tmp = plot_single_PSF( psf, lmbd, r_psf, alph, opt_plt, opt, i_px )
opt_plt.pbaspect = 1 ;
plt_mnmx( psf, opt_plt ) ;
xlabel( 'PIX (1 mas)' )
ylabel( 'PIX (1 mas)' )
title( sprintf( 'PSF at %i mas from the center. Starshade rotated by %2.1f deg', r_psf, alph ) )
% Storing the image
img = getframe( gcf ) ;
pth_fg_tmp = sprintf( '%soutput/video/spinning_psf/', opt.installation_path ) ;
  if ~isdir( pth_fg_tmp )
  mkdir( pth_fg_tmp ) ;
  end
fg_nm_tmp = sprintf( '%spsf_%i_%03i.png', pth_fg_tmp, r_psf, i_px ) ;
imwrite( img.cdata, fg_nm_tmp ) ;

% Creating a video from the individual images
function create_psf_video( fg_nm )

% Video delay
vd_dly = 40 ;
system( sprintf( 'convert -delay %i -loop 0 %s*.png %s.gif', vd_dly, fg_nm( 1 : end - 7 ), fg_nm( 1 : end - 8 ) ) ) ;
delete( sprintf( '%s*.png', fg_nm( 1 : end - 7 ) ) ) ;
