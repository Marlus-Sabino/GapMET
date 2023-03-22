function [dt_pree,dt_flag] = MUK(dt_ori,dt_mean_month,dt_months,SD_lim,ref_type,index_est,dt_ext,dt_mean_ext)
%% MUK : UK traditionalmethod  
% MUK is one of the methods used in the GapMet function to gapfilling
% meteorological time series. MUK estimatives is made assuming a constant
% difference value between the historical mean on the reference series and
% the gap-filled station time series.
%
%[dt_pree,dt_flag] = MUK(dt_ori,dt_mean_month,dt_months,SD_lim,ref_type,index_est,dt_ext,dt_mean_ext)
%
%INPUTS:
%Name:        Description:                                    type:
%dt_ori        = Matriz containing the timeseries.              matriz(m,n)
%                Each columns  represents one stations.
%                Each row represents a time step
%
%dt_mean_month = Matriz containing the historical monthly       matriz(m,n)
%                mean time series.
%
%dt_months     = Array with the month time serie                vector[m,1]
%
%SD_lim        = Maximum standard deviations (SD).              vector[1,1]
%                Limits the filling values between “D” 
%                standard deviations of the mean observed
%                data.(mean+-D*SD).
%                Defaut = 3. Must be an possitive integer.
%                Additionally, user can set the limit as the 
%                maximum and minimum values of the observed 
%                data by set “SD_lim = 0”
%
%ref_type     = What dataset will be used as reference station   string[1,1]
%               for the gapffiling:
%               - nearby  (1): Use neaby stations provided 
%                            on the same dataset (default)
%               - external(2): Use stations provided on 
%                            an external dataset
%
%index_est     = Matriz containing the pririty order            matriz(m,n)
%                of nearby stations to be used as 
%                reference time serie.
%                requires when "ref_type = nearby"
%
%dt_ext       = Matriz containing the external timeseries.      matriz(m,n)
%               Each columns  represents one stations.
%               Each row represents a time step.
%               requires when "ref_type = external"
%
%dt_mean_ext  = Matriz containing the historical monthly       matriz(m,n)
%               mean time series of the reference dataset.
%               requires when "ref_type = external"
%
%OUTPUTS
%
%dt_pree      = Matriz containing the gapfilled timeseries      matriz(n,m)
%dt_flag      = matriz containing the flag of "dt_pree" data    matriz(n,m)
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

%----------------------------------------------------------------------
%% 2. Gapfilling 
%----------------------------------------------------------------------
numNans_inicio = nnz(isnan(dt));
numNans_fim = 0;
iter = 1;

if ref_type==1 
    while numNans_fim < numNans_inicio
        numNans_inicio = nnz(isnan(dt));
        if numNans_inicio==0;break;end
        disp(['Running Iteration ' num2str(iter,'%i') '. Number of Gaps ' num2str(numNans_inicio,'%i')]);
        for p=1:size(dt_ori,2)
            disp(['Gapfilling station ' num2str(p,'%i') ' of ' num2str(size(dt_ori,2),'%i')]);
            kf=index_est(p,:);kf=kf(~isnan(kf));
            for k=1:length(kf)
                gap        = dt(:,p);
                gap_filled = gap;
                if iter==1
                    ref   =  dt_ori(:,index_est(p,k));
                else
                    ref   =  dt(:,index_est(p,k));
                end
                for i=1:length(gap)
                    if(isnan(gap(i,1)))
                        gap_filled(i,1) = ref(i,1)+(dt_mean_month(dt_months(i),p)-dt_mean_month(dt_months(i),index_est(p,k)));     
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

if ref_type==2
    while numNans_fim < numNans_inicio
        numNans_inicio = nnz(isnan(dt));
        if numNans_inicio==0;break;end
        for p=1:size(dt_ori,2)
            disp(['Gapfilling station ' num2str(p,'%i') ' of ' num2str(size(dt_ori,2),'%i')]);
            gap = dt(:,p);
            gap_filled = gap;
            ref = dt_ext(:,p);
            for i=1:length(gap)
                if(isnan(gap(i,1)))
                    gap_filled(i,1) = ref(i,1)+(dt_mean_month(dt_months(i),p)-dt_mean_ext(dt_months(i),p));
                end
            end
            [gap_filled_checked]=LimCheck(dt_ori,gap,gap_filled,SD_lim);
            dt(:,p) = gap_filled_checked;
        end
        numNans_fim = nnz(isnan(dt));
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