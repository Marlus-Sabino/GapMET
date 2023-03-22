function [dt_fill,dt_flags,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type,dt_extern)

%% GAPMET - MULTIPLE METHODS FOR GAP-FILLING METEOROLOGY STATIONS : 

% The script runs 6 methodologies for gap-filling meteorological data
% series. The script was originally tested in 33 automatic weather stations
% on the state of Mato Grosso, Brazil with data series of the variables
% maximum and minimum air temperature, relative humidity, downward solar
% radiation, and mean wind speed) (Sabino et al., 2023).
% - Sabino, M. and Souza, A.P.D., 2022. Gap-filling meteorological data
% series using the GapMET software in the state of Mato Grosso, Brazil.
% Revista Brasileira de Engenharia Agrícola e Ambiental, 27, pp.149-156.
% DOI: http://dx.doi.org/10.1590/1807-1929/agriambi.v27n2p149-156
%
%[dt_pree,dt_flag,dt_dist,dt_corr] = GapMet(dt_coord,dt_data,pree_meth,min_est,max_dis,SD_lim,ref_type,dt_extern)
%
%  INPUTS:
%  Name:        Description:                                    type:
%  dt_coord  = Matriz containing 3 column with:                 matriz(n,3)
%              (id,latitude,longitude) of 
%              n stations (rows).
%              latitude and longitude must be in decimal degrees
%
%  dt_data   = Matriz containing the timeseries.                matriz(m,n)
%              Each columns  represents one stations.
%              Each row represents a time step 
%              *Columns 1 to 3 are reserved to
%              year, month and day, respectively.
%
%  pree_meth  = Gapfilling methodology:                         string[1,1]
%                - RLS : Simple linear regression
%                - RLM : Multiple linear regression 
%                - MPR : Regionalweighting method
%                - MUK : UK traditionalmethod
%                - IDM : Inverse distance method 
%                - MAS : Simple arithmetic mean
%
%  min_est   =  Minimal number os nearby stations used as       vector[1,1]
%               referrence. Must be a number betwwen 1 and 
%               the maximmum number of reference stations 
%               provided.The RLS and MUK methods will only 
%               use more than 1 station as reference if 
%               there are commitment gaps between the 
%               reference and gapfilled station.
%               Defaut = 3. Must be an possitive integer.
%
%  max_dis   = Maximum distance between the reference and       vector[1,1]
%              gapfilled station.
%              Defaut = 100. Must be an possitive integer.
%
%  SD_lim    = Maximum standard deviations (SD).                vector[1,1]
%              Limits the filling values between “D” 
%              standard deviations of the mean observed data.
%              (mean+-D*SD).
%              Defaut = 3. Must be an possitive integer.
%              Additionally, user can set the limit as the 
%              maximum and minimum values of the observed data 
%              by set “SD_lim = 0”
%
%  ref_type  = What dataset will be used as reference station   string[1,1]
%              for the gapffiling:
%               - nearby   : Use neaby stations provided 
%                            on the same dataset (Use this for 
%                            the following methods: RLM; MPR; 
%                            IID or MAS).(default)
%               - external : Use stations provided on 
%                            an external dataset (Use this for
%                            the following methods: RLS or MUK).
%
%  dt_exter = Matriz containing external timeseries            matriz(n,m)
%              Each columns  represents one stations.
%              Each row represents a time step 
%              *Columns 1 to 3 are reserved to
%              year, month and day, respectively.
%              *External reference stations in dt_satel
%              must match the column number of dt_dados
%              *if "ref_type = external" user must provide dt_satel
%
% OUTPUTS:
% Name:        Description:                                    type:
% dt_fill    = matriz containing the gapfiled timeseries        matriz(n,m)  
% dt_flags   = matriz containing the flag of "dt_pree" data    matriz(n,m)
%              - 0     : original data
%              - 1     : gapfilled on first interaction
%                       (use reference serie original data)
%              - 2...x : gapfilled on secont to x interaction
%                       (use reference serie gapfilled data)
%              - NaN: unfilled gap
% dt_dist    = matriz containing the distance in Kilometers     matriz(n,m) 
%              between the stations in dt_dados.
% dt_corr    = matriz containing the pearson correlation rs     matriz(n,m) 
%              between the stations in dt_dados.
%
%--------------------------------------------------------------------------
%% 1. Check and recover inputs
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1.1 Check Input Parameters
%--------------------------------------------------------------------------
if length(pree_meth)~=1 && ~isstring(pree_meth) && contains(("RLS,RLM,MPR,MUK,IDM,MAS"),pree_meth)==0
    error(['''pree_meth'' must be a string containing one of the gapfilling '...
           'methods: (RLS,RLM,MPR,MUK,IDM, or MAS)'])
end

isaninteger = @(x)isfinite(x) & x==floor(x);

if ~exist('min_est','var')
    min_est = 3;
elseif  isempty(min_est)
    min_est = 3;
elseif isaninteger(min_est)==0 && any(min_est<0)
    error(('"SD_lim" must be a possitive integer'))
end

if ~exist('max_dis','var')
    max_dis = 100;
elseif isempty(max_dis)
    max_dis = 100;
elseif isaninteger(max_dis)==0 && any(max_dis<0)
    error(('"SD_lim" must be a possitive integer'))
end

if ~exist('SD_lim','var')
    SD_lim = 3;
elseif isempty(SD_lim)
    SD_lim = 3;
elseif isaninteger(SD_lim)==0 && any(SD_lim<0)
    error(('"SD_lim" must be a possitive integer'))
end

if ~exist('ref_type','var')
    ref_type = 1;
elseif isempty(ref_type)
    ref_type = 1;
elseif contains(("nearby"),ref_type)
    ref_type = 1;
elseif contains(("external"),ref_type)
    ref_type = 2;
end

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1.2 Check Input data
%--------------------------------------------------------------------------
if any(size(dt_data,2)<(min_est+4))
    error(['The number of columns (n) in ''dt_dados'' must be the equal or higher'...
           'than: n = min_est + 4'])
end

if ~any(size(dt_coord,2)== 3)
    error(['The number of columns (m) in ''dt_coord'' must be the equal to 3 '...
           '(station_id,latitude,longitude)'])
end

if size(dt_coord,1)~=(size(dt_data,2)-3)
    error(['The number of stations (rows) in ''dt_coord'' is diffrent than '...
           'The number of stations (columns) in ''dt_dados'''])
end

if ref_type==2
    if ~any((size(dt_data,2))==(size(dt_extern,2)))
        error(['The number of stations (columns) in ''dt_satel'' is diffrent than'...
               'The number of stations (columns) in ''dt_dados'''])
    end
    if ~any((size(dt_data,1))==(size(dt_extern,1)))
        error(['The timeseries size (time steps) in ''dt_satel'' do not match the'...
               'The timeseries size(time steps) in ''dt_dados'''])
    end
end

dt_ori  = table2array(dt_data(:,4:end));
if ref_type==2;dt_ext = table2array(dt_extern(:,4:end));end  

dt_cod = table2array(dt_coord(:,1));
dt_lat = table2array(dt_coord(:,2));
dt_lon = table2array(dt_coord(:,3));

if any(dt_lat(:,1) >  90) && any(dt_lat(:,1) < -90) &&...
   any(dt_lon(:,1) > 180) && any(dt_lon(:,1) < -180)
    lat_inf = unique(dt_cod(all((dt_lat(:,1)<-90),2),:)).';
    lat_sup = unique(dt_cod(all((dt_lat(:,1)> 90),2),:)).';
    lon_inf = unique(dt_cod(all((dt_lon(:,1)<-180),2),:)).';
    lon_sup = unique(dt_cod(all((dt_lon(:,1)> 180),2),:)).';
    est_err = unique([lat_inf,lat_sup,lon_inf,lon_sup]);
    
     if ~any((size(est_err,2)==0))  
         d = [est_err',[repmat({','},numel(est_err)-1,1);{[]}]]';
         error(['The geographic coordinates in the station',d{:},' must be between'...
               '-90<lat>90 e -180<lon>180']);
     end
end

%--------------------------------------------------------------------------
%% 2 Calculate stations distance
%--------------------------------------------------------------------------
% The the distance between the stations calculatio is done by function
% vdist initially proposed by Vincenty (1975) and adapted by Kleder (2021)
%
% - Vincenty, T. "Direct and Inverse Solutions of Geodesics on the Ellipsoid
%with Application of Nested Equations", Survey Review, vol. 23, no. 176,
%April 1975, p 88-93.
% - Michael Kleder (2021). Geodetic distance on WGS84 earth ellipsoid
%(https://www.mathworks.com/matlabcentral/fileexchange/5379-geodetic-distance-on-wgs84-earth-ellipsoid),
%MATLAB Central File Exchange. Retrieved December 4, 2021.

dt_dist = NaN(length(dt_lat));
for i=1:length(dt_lat)
    for j=1:length(dt_lat)
        %Select two stations.
        lat1 = dt_lat(i);
        lon1 = dt_lon(i);
        lat2 = dt_lat(j);
        lon2 = dt_lon(j);
        %Calculate the distance in Kilometers (Km) between the two stations.
        dist_Km  = vdist(lat1,lon1,lat2,lon2)./1000; 
        %Add results to table of distaces "dt_dist"
        dt_dist(i,j) = dist_Km;
    end
end

%Add station_id to "dt_dist" table.     
dt_dist = array2table(dt_dist,'RowNames',dt_cod(:,1));
dt_dist.Properties.VariableNames = dt_cod(:,1);

%----------------------------------------------------------------------
% 2.1 Check min_est and max_dis parameters (ref_type="nearby")
%----------------------------------------------------------------------  
% Check if there are a minimum number of stations (min_est) inside of
% the maximum radius distance between the stations (max_dis)

if ref_type==1
    dt = table2array(dt_dist);
    dt(dt==0)  = NaN;

    min_dis = round(max(min(dt)));
    if max_dis<min_dis
        error(['The shortest distance to all stations have at least'...
               ' 1 reference station within the "max_dis" is '...
               ,num2str(min_dis),' Km'])
    end

    dt1 = dt;dt1(dt1<=max_dis) = 1;dt1(dt1>max_dis) = 0;
    min_est_in_max_dis = min(sum(dt1,'omitnan'));
    dt1 = sort(dt,2);
    max_dis_in_min_est = round(max(dt1(:,min_est)));

    if max_dis<max_dis_in_min_est && min_est>min_est_in_max_dis
        error(['The number of reference stations within the "max_dis" radius'...
               ' is lower than set in "min_est".'...
               ' decrease the "min_est" parameter to ',num2str(min_est_in_max_dis), ' stations'...
               ' or increase "max_dis" to ',num2str(max_dis_in_min_est),' Km'])
    end
end

%----------------------------------------------------------------------
%% 3 Calculate correlation between stations
%----------------------------------------------------------------------
if ref_type==1
    dt_corr = NaN(length(dt_cod));    
    for i=1:size(dt_ori,2)
        for j=1:size(dt_ori,2)

            pont1 = dt_ori(:,i);
            pont2 = dt_ori(:,j);

            [RHO,PVAL] = corr(pont1,pont2,'Type','Pearson','Rows','pairwise');

            %Check if the correlations are significant at 0.05
            if(PVAL<=0.05)
                dt_corr(i,j) = RHO;
            else
                dt_corr(i,j) = NaN;
            end

        end
    end
    %Creat correlations matriz table     
    dt_corr = array2table(dt_corr,'RowNames',dt_cod);
    dt_corr.Properties.VariableNames = dt_cod;
end

if ref_type==2
    dt_corr = NaN(size(dt_cod));    
    for i=1:size(dt_ori,2)
        pont1 = dt_ori(:,i);
        pont2 = dt_ext(:,i);

        [RHO,PVAL] = corr(pont1,pont2,'Type','Pearson','Rows','pairwise');

        %Check if the correlations are significant at 0.05
        if(PVAL<=0.05)
            dt_corr(i,1) = RHO;
        else
            dt_corr(i,1) = NaN;
        end
    end
    %Creat correlations matriz table     
    dt_corr = array2table(dt_corr,'RowNames',dt_cod);
    %dt_corr.Properties.VariableNames = dt_cod;
end

%----------------------------------------------------------------------
%% 4 Classify reference stations based on the closeness and correlation
%----------------------------------------------------------------------    
if ref_type==1
    dt = table2array(dt_dist);dt(dt==0)  = NaN;
    dt1 = table2array(dt_corr);   
    dt1(dt>max_dis) = 0;dt1(isnan(dt)) = 0;dt1(isnan(dt1)) = 0;dt1=abs(dt1);
    [M,I] = sort(dt1,2,'descend');I(M==0)=NaN;
    index_est = I;
end
%--------------------------------------------------------------------------
%% 5 Calculate the long term monthly mean
%--------------------------------------------------------------------------
% Long term mean is used in the MPR and MUK methods
if contains(("MPR,MUK"),pree_meth)==1
    dt_months  = table2array(dt_data(:,2));
    table_mes(:,1)=dt_months(:,1);table_mes(:,2:size(dt_ori,2)+1)=dt_ori;
    table_mes=array2table(table_mes); 
    dt_mean_month = groupsummary(table_mes,{'table_mes1'}, {'mean'});clear table_mes
    dt_mean_month = table2array(dt_mean_month(:,3:end));
    
    if ref_type == 2
        table_mes(:,1)=dt_months(:,1);table_mes(:,2:size(dt_ext,2)+1)=dt_ext;
        table_mes=array2table(table_mes);
        dt_mean_ext = groupsummary(table_mes,{'table_mes1'}, {'mean'});clear table_mes
        dt_mean_ext = table2array(dt_mean_ext(:,3:end));
    end
end

%--------------------------------------------------------------------------
%% 6 Gapfilling 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 6.1 RLS : Simple linear regression
%--------------------------------------------------------------------------
if pree_meth=="RLS"
    if ref_type==1
        [dt_pree,dt_flag] = RLS(dt_ori,SD_lim,ref_type,index_est,[]);
    else
        [dt_pree,dt_flag] = RLS(dt_ori,SD_lim,ref_type,[],dt_ext);
    end
end

%--------------------------------------------------------------------------
% 6.2 RLM : Multiple linear regression
%--------------------------------------------------------------------------
if pree_meth=="RLM"
    if ref_type==1
        [dt_pree,dt_flag] = RLM(dt_ori,SD_lim,min_est,index_est);
    else
       error('This methodology (RLM) can only be used with nearby reference stations (Ref_type="nearby")')
    end
end

%--------------------------------------------------------------------------
% 6.3 MPR : Regional weighting method
%--------------------------------------------------------------------------
if pree_meth=="MPR"
    if ref_type==1
       [dt_pree,dt_flag] = MPR(dt_ori,dt_mean_month,dt_months,SD_lim,min_est,index_est);
    else
       error('This methodology (MPR) can only be used with nearby reference stations (Ref_type="nearby")')
    end
end

%--------------------------------------------------------------------------
% 6.4 MUK : UK traditionalmethod  
%--------------------------------------------------------------------------
if pree_meth=="MUK"
    if ref_type==1
        [dt_pree,dt_flag] = MUK(dt_ori,dt_mean_month,dt_months,SD_lim,ref_type,index_est,[],[]);
    else
        [dt_pree,dt_flag] = MUK(dt_ori,dt_mean_month,dt_months,SD_lim,ref_type,[],dt_ext,dt_mean_ext);
    end
end

%--------------------------------------------------------------------------
% 6.5 IDM : Inverse distance method
%--------------------------------------------------------------------------
if pree_meth=="IDM"
    if ref_type==1
        [dt_pree,dt_flag] = IDM(dt_ori,dt_dist,SD_lim,min_est,index_est);
    else
       error('This methodology (IDM) can only be used with nearby reference stations (Ref_type="nearby")')
    end
end

%--------------------------------------------------------------------------
% 6.6 MAS : Simple arithmetic mean
%--------------------------------------------------------------------------
if pree_meth=="MAS"
    if ref_type==1
        [dt_pree,dt_flag] = MAS(dt_ori,SD_lim,min_est,index_est);
    else
       error('This methodology (MAS) can only be used with nearby reference stations (Ref_type="nearby")')
    end
end

%--------------------------------------------------------------------------
%% 7. SAVE RESULTS
%--------------------------------------------------------------------------

%check/create a subdirectory to save results
if ~exist('gapfilled', 'dir')
   mkdir('gapfilled')
end

% correlation table
matfile = [pwd '\gapfilled\Correlations.csv'];
writetable(dt_corr,matfile,'WriteRowNames',true)
% distance table
matfile = [pwd '\gapfilled\Distances.csv'];
writetable(dt_dist,matfile,'WriteRowNames',true)
% gapfilled dataset
dt_fill = dt_data; dt_fill(:,4:end) = array2table(dt_pree);
matfile = [pwd '\gapfilled\Gapfilled.csv'];
writetable(dt_fill,matfile)
% gapfilling flags
dt_flags = dt_data; dt_flags(:,4:end) = array2table(dt_flag);
matfile = [pwd '\gapfilled\Flags.csv'];
writetable(dt_flags,matfile)

end


