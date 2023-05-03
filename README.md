# GapMET

The script **"GapMet.m"** provides multiple methods for gap-filling meteorological data series.

It was originally tested on 33 automatic weather stations located in the state of Mato Grosso, Brazil, by Sabino and Souza (2022). The data series used in the study include maximum and minimum air temperature, relative humidity, downward solar radiation, and mean wind speed.

The six methodologies implemented in the script are described in detail in the publication by Sabino and Souza (2022) and include:

* RLS: Simple linear regression
* RLM: Multiple linear regression
* MPR: Regional weighting method
* MUK: UK traditional method
* IDM: Inverse distance method
* MAS: Simple arithmetic mean

Users can choose which methodologies to use and adjust the parameters as needed.

## Parameters/Variables

The complete code function to runs "GapMet.m" is:

```javascript
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type,dt_extern)
```

The function takes in **five required inputs**: dt_coord, dt_data, pree_meth, min_est, max_dis, and SD_lim, as well as **two optional inputs**: ref_type and dt_extern.

* Required datasets:

  * **dt_coord**: a matrix containing 3 columns (id, latitude, and longitude) of n weather stations (rows). The latitude and longitude must be in decimal degrees.
  * **dt_data**: a matrix containing the time series. Each column represents one station, and each row represents a time step. *Columns 1 to 3 are reserved* for year, month, and day, respectively.

* Required parameters:

  * **pree_meth**: specifies the gap-filling methodology and can take one of six options: RLS, RLM, MPR, MUK, IDM, and MAS.
  * **min_est**: specifies the minimal number of nearby stations used as a reference time series. Must be a number between 1 and the maximum number of reference stations provided. The RLS and MUK methods will only use more than 1 station as reference if there are commitment gaps between the reference and gap-filled station. Must be a positive integer. *Default = 3*.
  * **max_dis**: specifies the maximum distance between the reference and gap-filled station. Must be a positive integer. *Default = 100*.
  * **SD_lim**: specifies the maximum standard deviation and limits the filling values between "D" standard deviations of the mean observed data. Additionally, the user can set the limit as the maximum and minimum values of the observed data by setting “SD_lim = 0”. *Default = 3*

* Optional parameters:

  * **ref_type**: specifies what dataset will be used as a reference time series for the gap-filling. It can take two options: nearby and external. *Default = 'nearby'*

* Optional datasets:

  * **dt_extern**: a matrix containing external time series and is required if ref_type is set to external. The external reference time series in dt_extern must match the column in dt_data.

The function has **four outputs**: dt_fill, dt_flags, dt_dist, and dt_corr.

  * **dt_fill**: a matrix containing the gap-filled time series.
  * **dt_flags**: a matrix containing the flag of dt_pree data where:
    * 0 = original data;
    * 1 = gap-filled on first interaction. Use reference series original data;
    * 2 ... x = gap-filled on second to x interaction. Use reference series gap-filled data;
    * NaN = unfilled gap.
  * **dt_dist**: a matrix containing the distance in kilometers between the weather stations in dt_data.
  * **dt_corr**: a matrix containing the Pearson correlation between the weather stations in dt_data.

### Flowchart of the GapMET gap-filling meteorological data process

![Flowchart of the GapMET gap-filling meteorological data process](https://user-images.githubusercontent.com/95511913/235925196-8e21f253-4f40-4fd1-83b2-69dbf6fdd298.png)


## Usage/Examples
The example contains a dataset of air temperature from 6 Automatic Weather Stations (AWS) to be gap-filled with one of the GapMet methodologies.

The example provides the dataset in the CSV file **original.csv** and the weather station coordinates in **coordinates.csv**.

Additionally, an external dataset obtained from the ERA5-Land satellite with pixels in the same coordinates of the AWS is provided in **external.csv**.


### Example 1: Gap-filling using nearby stations as reference time series
In this example, the dataset will be gap-filled using the multiple linear regression method (**pree_meth = "RLM"**) using data of at least 3 nearby AWS within a radius of 100 km of the AWS to be gap-filled (**min_est = 3**,**max_dis = 100**,**ref_type="nearby"** - all default values). The estimated data must be between 3 standard deviations of the mean temperature to be accepted (**SD_lim = 3**).

Initially, add the datasets:

```javascript
dt_coord  = readtable('coordinates.csv','PreserveVariableNames',true); %read file with id and coordinates of the AWS
dt_data   = readtable('original.csv','PreserveVariableNames',true);    %read the temperature dataset
```

Set the parameters:
```javascript
min_est   = [3];          % Minimal number of nearby sations used as reference time series.
max_dis   = [100];        % Maximum distance in Kilometers (km) between stations to be accepted as a reference time serie
SD_lim    = [3];          % Maximum standard deviations (SD). Limits the filling values between “D” standard deviations of the mean observed data.(mean+-D*SD).
pree_meth = ["RLM"];      % Gapfilling methodology used: "RLS" , "RLM" , "MPR" , "MUK" , "IID" , "MAS"
ref_type  = ["nearby"];   % Type of reference time serie: "nearby" ,"external".
```

Run the GapMet function:
```javascript
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type);
```

Additionally, since all the parameters were set as default, the function can be run as:

```javascript
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"RLM");
```

### Example 2: Gap-filling using nearby stations as reference time series

In this example, the dataset will be gap-filled using the regional weighting method (**pre_meth = "MPR"**) using data from at least two nearby AWS stations within a radius of 150 km of the AWS station to be gap-filled (**min_est = 2**, **max_dis = 150**, **ref_type = "nearby"**). The estimated data must fall within the maximum and minimum temperature values of the original dataset to be accepted (**SD_lim = 0**).

Initially, add the datasets:

```javascript
dt_coord  = readtable('coordinates.csv','PreserveVariableNames',true); %read file with id and coordinates of the AWS
dt_data   = readtable('original.csv','PreserveVariableNames',true);    %read the temperature dataset
```

Set the parameters:
```javascript
min_est   = [2];          % Minimal number of nearby sations used as reference time series.
max_dis   = [150];        % Maximum distance in Kilometers (Km) between stations to be accepted as a reference time serie
SD_lim    = [0];          % Maximum standard deviations (SD). Limits the filling values between “D” standard deviations of the mean observed data.(mean+-D*SD).
pree_meth = ["MPR"];      % Gapfilling methodology used: "RLS" , "RLM" , "MPR" , "MUK" , "IID" , "MAS"
ref_type  = ["nearby"];   % Type of reference time serie: "nearby" ,"external".
```

Run the GapMet function:
```javascript
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type);
```
or 

```javascript
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"MPR",2,150,0,"nearby");
```

### Example 3: Gap-filling using external stations as reference time series
In this example, the dataset will be gap-filled using simple linear regression (**pree_meth = "RLS"**) with an external dataset obtained from the ERA5-Land satellite with pixels in the same coordinates as the AWS (**ref_type = "external"**). The estimated data must be within 4 standard deviations of the mean temperature to be accepted (**SD_lim = **).

Initially, add the datasets:

```javascript
dt_coord  = readtable('coordinates.csv','PreserveVariableNames',true); %read file with id and coordinates of the AWS.
dt_data   = readtable('original.csv','PreserveVariableNames',true);    %read the temperature dataset.
dt_extern = readtable('external.csv','PreserveVariableNames',true); %read the temperature external dataset.
```

Set the parameters:
```javascript
min_est   = [];          % Minimal number of nearby sations used as reference time series.
max_dis   = [];        % Maximum distance in Kilometers (Km) between stations to be accepted as a reference time serie
SD_lim    = [4];          % Maximum standard deviations (SD). Limits the filling values between “D” standard deviations of the mean observed data.(mean+-D*SD).
pree_meth = ["RLS"];      % Gapfilling methodology used: "RLS" , "RLM" , "MPR" , "MUK" , "IID" , "MAS"
ref_type  = ["external"];   % Type of reference time serie: "nearby" ,"external".
```

Run the GapMet function:
```javascript
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type,dt_extern);
```
or 

```javascript
[dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,"RLS",[],[],4,"external",dt_extern);
```

## Known Issues and Updates

### Unrealistic Gapfilled Values
The quality control provided by the parameter SD_lim is sometimes not enough to prevent the estimation of unrealistic values (such as negative wind speed). Users should double-check the results for better quality control. We are working on an update that will allow users to set their own limit ranges for estimations in future versions of GapMET.

### Limit/Choose Weather Stations to be Gapfilled
In this version of GapMET, all weather stations (columns) in the dataset (dt_data) will be gapfilled automatically. We are working on an update that will allow users to choose which columns to gapfill instead of gapfilling the entire dataset.

## Authors

Please cite the publication by Sabino and Souza if you use this script in your research. 
Sabino, M. and Souza, A.P.D., 2022. Gap-filling meteorological data series using the GapMET software in the state of Mato Grosso, Brazil. Revista Brasileira de Engenharia Agrícola e Ambiental, 27, pp.149-156.http://dx.doi.org/10.1590/1807-1929/agriambi.v27n2p149-156.
