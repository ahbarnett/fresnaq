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

function opt = lower_opt( opt )
% Function to change to lower case all fields and subfields. The code below works up to 3 levels: opt.field_1.field_2.field_3 which should be enough in general. If necessary, it may be extended to subfields of higher order.
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

  if isstruct( opt )
  flds_opt = fields( opt ) ;
    for i_fld = 1 : numel( flds_opt )
    opt_tmp = opt.( flds_opt{ i_fld } ) ;
      if isstruct( opt_tmp )
      flds_opt_2 = fields( opt_tmp ) ;
        for i_fld_2 = 1 : numel( flds_opt_2 )
        opt_tmp_2 = opt_tmp.( flds_opt_2{ i_fld_2 } ) ;
          if isstruct( opt_tmp_2 )
          flds_opt_3 = fields( opt_tmp_2 ) ;
            for i_fld_3 = 1 : numel( flds_opt_3 )
            opt.( lower( flds_opt{ i_fld } ) ).( lower( flds_opt_2{ i_fld_2 } ) ).( lower( flds_opt_3{ i_fld_3 } ) ) = opt.( flds_opt{ i_fld } ).( flds_opt_2{ i_fld_2 } ).( flds_opt_3{ i_fld_3 } ) ;
            end % i_flds_3
          else
          opt.( lower( flds_opt{ i_fld } ) ).( lower( flds_opt_2{ i_fld_2 } ) ) = opt.( flds_opt{ i_fld } ).( flds_opt_2{ i_fld_2 } ) ;
          end
        end % i_fld_2
      else
      opt.( lower( flds_opt{ i_fld } ) ) = opt.( flds_opt{ i_fld } ) ;
      end
    end % i_fld
  end
