% run all component tests in FRESNAQ.
% Makes lots of figures, text accuracy outputs, takes around 30 secs total.
% Barnett 9/16/20

startup  % in case

% main library
fresnaq_pts
fresnaq_grid

% top level
demo_paramcurve
demo_radial
demo_starshades
shadow_multilambda_starshades

% util
curveareaquad
polarareaquad
starshadequad
perispecdiff
endcorrquad
test_Alpert_quad
globalinterp1d

% bdrymeths
crudecurvequad
test_bdwf
nsli_pts
nsli_emulates_bdwf
test_bdwf_starshades

% occulter
eval_sister_apod
show_sister_profiles

% sister_mods
test_makePupil
test_idealnonspinning          % this needs user to point to SISTER dir

dbclear all     % remove annoying error -> debug mode
