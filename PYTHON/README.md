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

The function takes in **six required inputs**: dt_coord, dt_data, pree_meth, min_est, max_dis, limit_type, and interval, as well as **two optional inputs**: ref_type and dt_extern.

* Required datasets:

  * **dt_coord**: a matrix containing 3 columns (id, latitude, and longitude) of n weather stations (rows). The latitude and longitude must be in decimal degrees.
  * **dt_data**: a matrix containing the time series. Each column represents one station, and each row represents a time step. *Columns 1 to 3 are reserved* for year, month, and day, respectively.

* Required parameters:

  * **pree_meth**: specifies the gap-filling methodology and can take one of six options: RLS, RLM, MPR, MUK, IDM, and MAS.
  * **min_est**: specifies the minimal number of nearby stations used as a reference time series. Must be a number between 1 and the maximum number of reference stations provided. The RLS and MUK methods will only use more than 1 station as reference if there are commitment gaps between the reference and gap-filled station. Must be a positive integer. *Default = 3*.
  * **max_dis**: specifies the maximum distance between the reference and gap-filled station. Must be a positive integer. *Default = 100*.
  * **limit_type**: Limits the filling values between an interval which can be defined as:
    * **range**: Uses either the min and max values of the data in *dt_data* or the min and max values define by the user.(Default).
    * **percentile**: Uses a min and max percentile for the interval.
  * **interval**: Set the min and max interval based on the *limit_type*.
    * *Default = [min,max]* for limit_type = "range".
    * *Default = [5,99]* for limit_type = "percentile".

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

![Flowchart of the GapMET 1.0 gap-filling meteorological data process (need update for GapMET 2.0)](https://user-images.githubusercontent.com/95511913/235925196-8e21f253-4f40-4fd1-83b2-69dbf6fdd298.png)


## Usage/Examples

## Known Issues and Updates

### Limit/Choose Weather Stations to be Gapfilled
In this version of GapMET, all weather stations (columns) in the dataset (dt_data) will be gapfilled automatically. We are working on an update that will allow users to choose which columns to gapfill instead of gapfilling the entire dataset.

## Authors

Please cite the publication by Sabino and Souza if you use this script in your research. 
Sabino, M. and Souza, A.P.D., 2022. Gap-filling meteorological data series using the GapMET software in the state of Mato Grosso, Brazil. Revista Brasileira de Engenharia Agrícola e Ambiental, 27, pp.149-156.http://dx.doi.org/10.1590/1807-1929/agriambi.v27n2p149-156.
