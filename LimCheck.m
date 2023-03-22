function [gap_filled_checked] = LimCheck(dt_ori,gap,gap_filled,SD_lim)
%% LimCheck : Validate the estimated values of gapfilled data
%
%[dt_pree,dt_flag] = MAS(dt_ori,SD_lim,min_est,index_est)
%
%INPUTS:
%Name:        Description:                                    type:
%dt_ori        = Matriz containing original timeseries.              matriz(m,n)
%                Each columns  represents one stations.
%                Each row represents a time step
%
%min_est       = Minimal number os nearby stations used as      vector[1,1]
%                referrence. Must be a number betwwen 1 and 
%                the maximmum number of reference stations 
%                provided.
%                Defaut = 3. Must be an possitive integer.
%
%gap           = vector containing the station timeserie        vector(m,1)
%                before being gapfilled.
%
%gap_filled    = vector containing the station timeserie        vector(m,1)
%                after being gapfilled.
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
%OUTPUTS
%
%gap_filled_checked  = Matriz containing the timeseries          matriz(n,m)
%                      ith accepted gapfilled data

gap_filled_checked = gap_filled;

if SD_lim==0
    for line_check=1:size(dt_ori,1)
        if(isnan(gap(line_check,1)))
            if gap_filled(line_check,1)>=(min(dt_ori(:),[],'omitnan')) && gap_filled(line_check,1)<=(max(dt_ori(:),[],'omitnan'))
            else
                gap_filled_checked(line_check,1) = NaN;
            end
        end
    end
else
    for line_check=1:size(dt_ori,1)
        if(isnan(gap(line_check,1)))
            if gap_filled(line_check,1)>=(mean(dt_ori(:),'omitnan')-(SD_lim*(std(dt_ori(:),'omitnan')))) && gap_filled(line_check,1)<=(mean(dt_ori(:),'omitnan')+(SD_lim*(std(dt_ori(:),'omitnan'))))
            else
                gap_filled_checked(line_check,1) = NaN;
            end
        end
    end
end
end