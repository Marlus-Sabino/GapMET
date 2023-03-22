function [dt_pree,dt_flag] = RLM(dt_ori,SD_lim,min_est,index_est)
%% RLM : Multiple linear regression
% RLS is one of the methods used in the GapMet function to gapfilling
% meteorological time series. RLS computes a multiple linear regression and
% estimate the data using the linear and angular coefficients of the
% regression between the gap-filled time series and a the reference time
% series.
%
% [dt_pree,dt_flag] = RLM(dt_ori,SD_lim,min_est,index_est)
%
% INPUTS:
%   Name:        Description:                                    type:
%   dt_ori  =  Matriz containing the timeseries.                matriz(m,n)
%              Each columns  represents one stations.
%              Each row represents a time step
%
%   SD_lim  =  Maximum standard deviations (SD).                vector[1,1]
%              Limits the filling values between “D” 
%              standard deviations of the mean observed data.
%              (mean+-D*SD).
%              Defaut = 3. Must be an possitive integer.
%              Additionally, user can set the limit as the 
%              maximum and minimum values of the observed data 
%              by set “SD_lim = 0”
%
%  min_est =  Minimal number os nearby stations used as       vector[1,1]
%             referrence. Must be a number betwwen 1 and 
%             the maximmum number of reference stations 
%             provided.
%             Defaut = 3. Must be an possitive integer.
%
%  index_est = Matriz containing the pririty order of nearby    matriz(m,n)
%              stations to be used as reference time serie.
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
isaninteger = @(x)isfinite(x) & x==floor(x);
if isempty(SD_lim)
    SD_lim = 3;
elseif isaninteger(SD_lim)==0 && any(SD_lim<0)
    error(('"SD_lim" must be a possitive integer'))
end

if any(size(dt_ori,2)<(min_est+1))
    error(['The number of columns (n) in ''dt_dados'' must be the equal or higher'...
           'than: n = min_est + 1'])
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
        gap  =  dt(:,p); %estacao a ser preenchida
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
                    for j=1:length(ind)
                        X0(:,j+1) = dt(:,index_est(p,ind(j)));
                    end
                    X0(:,1) = 1;
                    b = regress(gap,X0);
                    X2(2:length(b),:) = X1(~isnan(X1)).';
                    X2(1,1) = 1;
                    X3 = X2.*b;
                    y = sum(X3);
                    gap_filled(i,1) = y;
                    dt_flag(i,p)   = iter;
                    clear X0 X1 X2 X3 b y ind stats
                end
            end
        end
        [gap_filled_checked]=LimCheck(dt_ori,gap,gap_filled,SD_lim);
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