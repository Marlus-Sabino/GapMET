function [dt_pree,dt_flag] = MPR(dt_ori,dt_mean_month,dt_months,limit_type,interval,min_est,index_est)
%% MPR : Regional weighting method
% MPR is one of the methods used in the GapMet function to gapfilling
% meteorological time series. MPR estimate the gap using the mean values of
% the neighbour station weighted by their historical time series mean
%
% [dt_pree,dt_flag] = MPR(dt_ori,dt_mean_month,dt_month,SD_lim,min_est,index_est);
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
%limit_type = Limits the filling values between an interval     string[1,1]
%              which can be defined as:
%              - range      : Uses either the min and max values
%                             of the data in dt_data 
%                             or the min and max values define
%                             by the user.(Default)
%              - percentile : Uses a min and max percentile for 
%                             interval.
%interval  = Set the min and max interval based on the           array[2,1]
%             limit_type.
%             Default = [min,max] for limit_type = "range".
%             Default = [5,99] for limit_type = "percentile".
%
%min_est       = Minimal number os nearby stations used as      vector[1,1]
%                referrence. Must be a number betwwen 1 and 
%                the maximmum number of reference stations 
%                provided.
%                Defaut = 3. Must be an possitive integer.
%
%index_est     = Matriz containing the pririty order            matriz(m,n)
%                of nearby stations to be used as 
%                reference time serie.
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
if ~exist('limit_type','var')
    limit_type = "range";
elseif isempty(limit_type)
    limit_type = "range";
elseif ~contains(("range"),limit_type) 
    if ~contains(("percentile"),limit_type)
        error(('"limit_type" must be "range" or "percentile"'))
    end
end

if contains(("percentile"),limit_type)
    if ~exist('interval','var')
        interval = [5 95];
    elseif isempty(interval)
        interval = [5 95];
    elseif ~isnumeric(interval)
        error(('The "percentile" "interval" must be an array with min and max percentile'))
    elseif interval(1)>=interval(2)
        error(('The "percentile" min interval must be lower than max interval'))
    elseif any(interval<0) || any(interval>100)
        error(('The "percentile" min and max interval must be between 0 and 100'))
    end
elseif contains(("range"),limit_type)
    if ~exist('interval','var')
        interval = [min(dt_ori(:)) max(dt_ori(:))];
    elseif isempty(interval)
        interval = [min(dt_ori(:)) max(dt_ori(:))];
    elseif ~isnumeric(interval)
        error(('The "range" "interval" must be an array with min and max values'))
    elseif interval(1)>=interval(2)
        error(('The min interval must be lower than max interval'))
    end   
end

if any(size(dt_ori,2)<(min_est+1))
    error(['The number of columns (n) in ''dt_dados'' must be the equal or higher'...
           'than: n = min_est + 1'])
end

if any(size(dt_mean_month,1)>(12))
    error('The number of rows (m) in ''dt_mean_month'' must be the equal or lower to 12 (months)')
end

if any(size(dt_mean_month,2)<(min_est+1))
    error(['The number of columns (n) in ''dt_mean_month'' must be the equal or higher'...
           'than: n = min_est + 1'])
end

if any(size(dt_months,1)~=(size(dt_ori,1)))
    error(['The number of rows (m) in ''dt_months'' must be the equal to '...
           'The number of rows (m)(time steps) in ''dt_ori'''])
end

dt = dt_ori;
dt_flag = zeros(size(dt_ori,1),size(dt_ori,2)-1);

%----------------------------------------------------------------------
%% 2. Gapfilling 
%----------------------------------------------------------------------
numNans_inicio = nnz(isnan(dt)); %Inicial number of gaps in the dataset
numNans_fim = 0;                 %Final number of gaps in the dataset
iter = 1;                        %Number of iteration

%If after the first iteration there are still gaps in the time series new
%iterations will be performed, in which the previously filled data from
%nearby stations can be used as reference. The gapfilling will end either
%when all the gaps are estimated or when the number of gaps at the end of
%each interaction becomes constant, meaning that gaps are occurring
%concomitantly between the station to be filled and the reference stations.

while numNans_fim < numNans_inicio
    numNans_inicio = nnz(isnan(dt));
    if numNans_inicio==0;break;end
    disp(['Running Iteration ' num2str(iter,'%i') '. Number of Gaps ' num2str(numNans_inicio,'%i')]);   
    for p=1:size(dt_ori,2)
        disp(['Gapfilling station ' num2str(p,'%i') ' of ' num2str(size(dt_ori,2),'%i')]);
        gap  =  dt(:,p);
        gap_filled = gap;
        for j=1:length(index_est(1,:))
            if(isnan(index_est(p,j)))
                X(:,j) = NaN;
            elseif iter==1
                X(:,j) = dt_ori(:,index_est(p,j));
            else
                X(:,j) = dt(:,index_est(p,j));
            end
        end
        for i=1:length(X)
            if(isnan(dt(i,p)))
                X1 = X(i,:);
                ind = find(~isnan(X1));
                if(length(ind)>=min_est)
                    m_ref = dt_months(i);
                    Xm    = dt_mean_month(m_ref,:);
                    Xm    = Xm(~isnan(X1));
                    X2    = X1(~isnan(X1));
                    gap_filled(i,1) = (1/length(Xm))*(sum(X2./Xm))*(dt_mean_month(m_ref,p));
                    dt_flag(i,p)   = iter;
                    clear X1 X2 Xm
                end
            end
        end
        [gap_filled_checked]=LimCheck(dt_ori,gap,gap_filled,limit_type,interval);
        dt(:,p) = gap_filled_checked;
    end
    numNans_fim = nnz(isnan(dt));
    iter = iter+1;
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