%% example of gapffiling with GapMet
% the example contains a dataset of temperature of 6 automatic weather
% stations (AWS) to be gapfilled with one of GapMet methodologies

%% Gapfilling with date of nearby station
%

dt_coord  = readtable('coordenates.csv','PreserveVariableNames',true); %read file with id and coordanates of the AWS
dt_data   = readtable('original.csv','PreserveVariableNames',true);    %read the temperature dataset

%[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type);

% Gapfilling using the multiple linear regression (pree_meth =
% RLM) using data of at least 3 nearby AWS in a radius of 100 Km of the
% gapfilled AWS (min_est = 3,max_dis = 100,ref_type="nearby" - all deafaut
% values). The estimated data must be between 3 standat deviation of the
% mean temperature to be accept.

[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"RLM");


% Gapfilling using the regional weighting method (pree_meth = MPR) using
% data of at least 2 nearby AWS in a radius of 150 Km of the gapfilled AWS
% (min_est = 4,max_dis = 150,ref_type="nearby" - all deafaut values). The
% estimated data must be between the values of maximmum and minimmum
% temperature of the original dataset to be accept.

[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"MPR",2,150,0,"nearby");

%% Gapfilling with date of external dataset
% external dataset were obtained of the ERA5-Land satelite with pixels in
% the same coordenates of the AWS

dt_extern = readtable('external.csv','PreserveVariableNames',true); %read the temperature external dataset

% Gapfilling using the simple linear regression (pree_meth = RLS) using
% data of external dataset as reference (ref_type="external"). The
% estimated data must be between 4 standat deviation of the mean
% temperature to be accept.

[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"RLS",[],[],4,"external",dt_extern);



