function [gap_filled_checked] = LimCheck(dt_ori,gap,gap_filled,limit_type,interval)
%% LimCheck : Validate the estimated values of gapfilled data
%
%[gap_filled_checked] = LimCheck(dt_ori,gap,gap_filled,SD_lim)
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
%limit_type = Limits the filling values between an interval    string[1,1]
%              which can be defined as:
%              - range      : Uses either the min and max values
%                             of the data in dt_data 
%                             or the min and max values define
%                             by the user.(Default)
%              - percentile : Uses a min and max percentile for 
%                             interval.
%interval  = Set the min and max interval based on the          array[2,1]
%             limit_type.
%             Default = [min,max] for limit_type = "range".
%             Default = [5,99] for limit_type = "percentile".
%
%OUTPUTS
%
%gap_filled_checked  = Matriz containing the timeseries          matriz(n,m)
%                      with accepted gapfilled data

gap_filled_checked = gap_filled;

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
    
    min_lim = prctile(dt_ori(:),interval(1));
    max_lim = prctile(dt_ori(:),interval(2));

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
    
    min_lim = interval(1);
    max_lim = interval(2);
    
end

for line_check=1:size(dt_ori,1)
    if(isnan(gap(line_check,1)))
        if gap_filled(line_check,1)<=min_lim || gap_filled(line_check,1)>=max_lim
            gap_filled_checked(line_check,1) = NaN;            
        end
    end
end
    
end