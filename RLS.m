function [dt_pree,dt_flag] = RLS(dt_ori,SD_lim,ref_type,index_est,dt_ext)
%% RLS : Simple linear regression
% RLS is one of the methods used in the GapMet function to gapfilling
% meteorological time series. RLS computes a simple linear regression and
% estimate the data using the linear and angular coefficients of the
% regression between the gap-filled time series and a reference time
% series. The reference time series can be a nearby station or an external
% dataset.
%
% [dt_pree,dt_flag] = RLS(dt_ori,SD_lim,ref_type,index_est,dt_ext)
%
% INPUTS:
%   Name:        Description:                                    type:
%   dt_ori  =  Matriz containing the timeseries.                matriz(m,n)
%              Each columns  represents one stations.
%              Each row represents a time step
%
%   SD_lim    = Maximum standard deviations (SD).                vector[1,1]
%              Limits the filling values between “D” 
%              standard deviations of the mean observed data.
%              (mean+-D*SD).
%              Defaut = 3. Must be an possitive integer.
%              Additionally, user can set the limit as the 
%              maximum and minimum values of the observed data 
%              by set “SD_lim = 0”
%
%   ref_type = What dataset will be used as reference station   string[1,1]
%              for the gapffiling:
%               - nearby  (1): Use neaby stations provided 
%                            on the same dataset (default)
%               - external(2): Use stations provided on 
%                            an external dataset
%
%  index_est = Matriz containing the pririty order of nearby    matriz(m,n)
%              stations to be used as reference time serie.
%
%  dt_ext    = Matriz containing external timeseries            matriz(n,m)
%              Each columns  represents one stations.
%              Each row represents a time step 
%              *External reference stations in dt_satel
%              must match the column number of dt_dados
%              *if "ref_type = external" user must provide dt_satel
%
% OUTPUTS
%
%  dt_pree   = Matriz containing the gapfilled timeseries       matriz(n,m)
%  dt_flag   = matriz containing the flag of "dt_pree" data     matriz(n,m)
%              - 0     : original data
%              - 1     : gapfilled on first iteration
%                       (use reference serie original data)
%              - 2...x : gapfilled on secont to x iteration
%                       (use reference serie gapfilled data)
%              - NaN: unfilled gap
%
%--------------------------------------------------------------------------
%% 1. Check inputs
%--------------------------------------------------------------------------
if isempty(ref_type)
    ref_type = 1;
elseif ref_type==1
elseif ref_type==2
elseif contains(("nearby"),ref_type)
    ref_type = 1;
elseif contains(("external"),ref_type)
    ref_type = 2;
end

isaninteger = @(x)isfinite(x) & x==floor(x);
if isempty(SD_lim)
    SD_lim = 3;
elseif isaninteger(SD_lim)==0 && any(SD_lim<0)
    error(('"SD_lim" must be a possitive integer'))
end

if ref_type==1
    if ~any((size(dt_ori,2))==(size(index_est,1)))
        error(['The number of stations (columns) in ''index_est'' is diffrent than'...
               'The number of stations (columns) in ''dt_ori'''])
    end   
elseif ref_type==2
    if ~any((size(dt_ori,2))==(size(dt_ext,2)))
        error(['The number of stations (columns) in ''dt_sat'' is diffrent than'...
               'The number of stations (columns) in ''dt_ori'''])
    end
    if ~any((size(dt_ori,1))==(size(dt_ext,1)))
        error(['The timeseries size (time steps) in ''dt_sat'' do not match the'...
               'The timeseries size(time steps) in ''dt_ori'''])
    end
end

dt = dt_ori;
dt_flag = zeros(size(dt_ori,1),size(dt_ori,2)-1);
coef_lin=NaN(size(dt_ori,2)); 
coef_ang=NaN(size(dt_ori,2));
R2=NaN(size(dt_ori,2));
%--------------------------------------------------------------------------
%% 1. RLS: Using nearby stations to gapfilling  (ref_type==1)
%--------------------------------------------------------------------------
if ref_type==1
    
    %----------------------------------------------------------------------
    % 1.1 Obtaining simple regression coefficients
    %----------------------------------------------------------------------
   
    for p=1:size(dt_ori,2)
        kf=index_est(p,:);kf=kf(~isnan(kf));
        for k=1:length(kf)
            gap  =  dt(:,p);                
            ref  =  dt(:,index_est(p,k));   
            mdl = fitlm(ref,gap);             
                coef_lin(p,k) = table2array(mdl.Coefficients(1,1));
                coef_ang(p,k) = table2array(mdl.Coefficients(2,1));
                R2(k,p)       = mdl.Rsquared.Adjusted;
        end
        clearvars kf
    end
    
    %----------------------------------------------------------------------
    % 1.2 Gapfilling 
    %----------------------------------------------------------------------       
    %If after the first iteration there are still gaps in the time series
    %new iterations will be performed, in which the previously filled data
    %from nearby stations can be used as reference. The gapfilling will end
    %either when all the gaps are estimated or when the number of gaps at
    %the end of each interaction becomes constant, meaning that gaps are
    %occurring concomitantly between the station to be filled and the
    %reference stations.
    
    numNans_fim = nnz(isnan(dt));
    iter = 1;                        % Iterations  

    if numNans_fim>0
        numNans_inicio = nnz(isnan(dt)); 
        numNans_fim = 0;                         
        while numNans_fim < numNans_inicio
            numNans_inicio = nnz(isnan(dt));
            if numNans_inicio==0;break;end
            disp(['Running Iteration ' num2str(iter,'%i') '. Number of Gaps ' num2str(numNans_inicio,'%i')]);
            for p=1:size(dt_ori,2)
                disp(['Gapfilling station ' num2str(p,'%i') ' of ' num2str(size(dt_ori,2),'%i')]);
                kf=index_est(p,:);kf=kf(~isnan(kf));
                for k=1:length(kf)
                    if iter==1
                        ref   =  dt_ori(:,index_est(p,k));
                    else 
                        ref   =  dt(:,index_est(p,k));
                    end
                    gap   =  dt(:,p);
                    gap_filled = gap;
                    for i=1:length(gap)
                        if(isnan(gap(i,1)))
                            y = coef_lin(p,k)+(coef_ang(p,k).*ref(i,1));
                            gap_filled(i,1) = y;
                            dt_flag(i,p)  = iter;                                             
                        end
                    end
                    [gap_filled_checked]=LimCheck(dt_ori,gap,gap_filled,SD_lim);
                    dt(:,p) = gap_filled_checked;
                end
            end
            numNans_fim = nnz(isnan(dt));
            iter = iter+1;
        end
    end
end

%--------------------------------------------------------------------------
%% 2. RLS: Using external time series   (ref_type==2)
%--------------------------------------------------------------------------
if ref_type==2   
    %----------------------------------------------------------------------
    % 1.2 Obtaining simple regression coefficients
    %----------------------------------------------------------------------
    dt  = dt_ori;
    dt1 = dt_ext;

    for p=1:size(dt_ori,2)
        gap = dt(:,p);
        ref = dt1(:,p);
        mdl = fitlm(ref,gap);
            coef_lin(p,1) = table2array(mdl.Coefficients(1,1));
            coef_ang(p,1) = table2array(mdl.Coefficients(2,1));
            R2(1,p)       = mdl.Rsquared.Adjusted;  
    end
    
    %----------------------------------------------------------------------
    % 2.2 Gapfilling
    %----------------------------------------------------------------------     
    numNans_fim = nnz(isnan(dt));
    iter = 1; 
    
    if numNans_fim>0
        numNans_inicio = nnz(isnan(dt));
        numNans_fim = 0;
        while numNans_fim < numNans_inicio
            numNans_inicio = nnz(isnan(dt));
            if numNans_inicio==0;break;end
                for p=1:size(dt_ori,2)
                    disp(['Gapfilling station ' num2str(p,'%i') ' of ' num2str(size(dt_ori,2),'%i')]);
                    gap = dt(:,p);
                    ref = dt1(:,p);
                    gap_filled = gap;
                    for i=1:length(gap)
                        if(isnan(gap(i,1)))
                            y = coef_lin(p,1)+(coef_ang(p,1).*ref(i,1));
                            gap_filled(i,1) = y;
                            dt_flag(i,p)  = iter;
                        end
                    end
                    [gap_filled_checked]=LimCheck(dt_ori,gap,gap_filled,SD_lim);
                    dt(:,p) = gap_filled_checked;
                end
            numNans_fim = nnz(isnan(dt));
            iter = iter+1;
        end
    end
end

%----------------------------------------------------------------------
%% 3. Finishing gapfilled dataset
%----------------------------------------------------------------------  
dt_pree = dt;
dt_flag(isnan(dt_pree)) = NaN;

if numNans_fim~=0
     warning([num2str(numNans_fim,'%i') ' Gaps lefts on the dataset']);
end
end 