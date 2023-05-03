% The example contains a dataset of air temperature from 6 Automatic Weather Stations (AWS) to be gap-filled with one of the GapMet methodologies.
% The example provides the dataset in the CSV file **original.csv** and the weather station coordinates in **coordinates.csv**.
% Additionally, an external dataset obtained from the ERA5-Land satellite with pixels in the same coordinates of the AWS is provided in **external.csv**.


%% Example 1: Gap-filling using nearby stations as reference time series
% In this example, the dataset will be gap-filled using the multiple linear regression method (**pree_meth = "RLM"**) using data of at least 3 nearby AWS within a radius of 100 km of the AWS to be gap-filled (**min_est = 3**,**max_dis = 100**,**ref_type="nearby"** - all default values). The estimated data must be between 3 standard deviations of the mean temperature to be accepted (**SD_lim = 3**).

% Initially, add the datasets:
dt_coord  = readtable('coordinates.csv','PreserveVariableNames',true); %read file with id and coordinates of the AWS
dt_data   = readtable('original.csv','PreserveVariableNames',true);    %read the temperature dataset

% Set the parameters:
min_est   = [3];          % Minimal number of nearby sations used as reference time series.
max_dis   = [100];        % Maximum distance in Kilometers (km) between stations to be accepted as a reference time serie
SD_lim    = [3];          % Maximum standard deviations (SD). Limits the filling values between “D” standard deviations of the mean observed data.(mean+-D*SD).
pree_meth = ["RLM"];      % Gapfilling methodology used: "RLS" , "RLM" , "MPR" , "MUK" , "IID" , "MAS"
ref_type  = ["nearby"];   % Type of reference time serie: "nearby" ,"external".

% Run the GapMet function:
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type);

% Additionally, since all the parameters were set as default, the function can be run as:
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"RLM");


%% Example 2: Gap-filling using nearby stations as reference time series
% In this example, the dataset will be gap-filled using the regional weighting method (**pre_meth = "MPR"**) using data from at least two nearby AWS stations within a radius of 150 km of the AWS station to be gap-filled (**min_est = 2**, **max_dis = 150**, **ref_type = "nearby"**). The estimated data must fall within the maximum and minimum temperature values of the original dataset to be accepted (**SD_lim = 0**).

% Initially, add the datasets:

dt_coord  = readtable('coordinates.csv','PreserveVariableNames',true); %read file with id and coordinates of the AWS
dt_data   = readtable('original.csv','PreserveVariableNames',true);    %read the temperature dataset

% Set the parameters:
min_est   = [2];          % Minimal number of nearby sations used as reference time series.
max_dis   = [150];        % Maximum distance in Kilometers (Km) between stations to be accepted as a reference time serie
SD_lim    = [0];          % Maximum standard deviations (SD). Limits the filling values between “D” standard deviations of the mean observed data.(mean+-D*SD).
pree_meth = ["MPR"];      % Gapfilling methodology used: "RLS" , "RLM" , "MPR" , "MUK" , "IID" , "MAS"
ref_type  = ["nearby"];   % Type of reference time serie: "nearby" ,"external".

% Run the GapMet function:
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type);

% or 

[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"MPR",2,150,0,"nearby");

%% Example 3: Gap-filling using external stations as reference time series
% In this example, the dataset will be gap-filled using simple linear regression (**pree_meth = "RLS"**) with an external dataset obtained from the ERA5-Land satellite with pixels in the same coordinates as the AWS (**ref_type = "external"**). The estimated data must be within 4 standard deviations of the mean temperature to be accepted (**SD_lim = **).

% Initially, add the datasets:

dt_coord  = readtable('coordinates.csv','PreserveVariableNames',true); %read file with id and coordinates of the AWS.
dt_data   = readtable('original.csv','PreserveVariableNames',true);    %read the temperature dataset.
dt_extern = readtable('external.csv','PreserveVariableNames',true); %read the temperature external dataset.

% Set the parameters:
min_est   = [];          % Minimal number of nearby sations used as reference time series.
max_dis   = [];        % Maximum distance in Kilometers (Km) between stations to be accepted as a reference time serie
SD_lim    = [4];          % Maximum standard deviations (SD). Limits the filling values between “D” standard deviations of the mean observed data.(mean+-D*SD).
pree_meth = ["RLS"];      % Gapfilling methodology used: "RLS" , "RLM" , "MPR" , "MUK" , "IID" , "MAS"
ref_type  = ["external"];   % Type of reference time serie: "nearby" ,"external".

% Run the GapMet function:
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type,dt_extern);

% or 

[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"RLS",[],[],4,"external",dt_extern);