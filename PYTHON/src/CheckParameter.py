import numpy as np

def validate_ref_type(ref_type):
    if ref_type == "":
        ref_type = 1
    elif "nearby" in ref_type:
        ref_type = 1
    elif "external" in ref_type:
        ref_type = 2
    
    return ref_type

def validate_pree_meth(pree_meth):
    if len(pree_meth) != 3 or pree_meth not in ["RLS", "RLM", "MPR", "MUK", "IDM", "MAS"]:
        raise ValueError("'pree_meth' must be a string containing one of the gap-filling methods: (RLS, RLM, MPR, MUK, IDM, or MAS)")

def validate_min_est(min_est):
    if min_est == "":
        min_est = 3
    elif not isinstance(min_est, int) and any(val < 0 for val in min_est):
        raise ValueError('"min_est" must be a positive integer')

def validate_max_dis(max_dis):
    if max_dis == "":
        max_dis = 100
    elif not isinstance(max_dis, int) and any(val < 0 for val in max_dis):
        raise ValueError('"max_dis" must be a positive integer')

def validate_limit_type(limit_type):
    if limit_type == "":
        limit_type = "range"
    elif "range" not in limit_type:
        if "percentile" not in limit_type:
            raise ValueError('"limit_type" must be "range" or "percentile"')

def validate_interval(interval, limit_type, dt_data):
    if "percentile" in limit_type:
        if interval == "":
            interval = [5, 95]
        elif not isinstance(interval, list) or len(interval) != 2 or not all(isinstance(x, (int, float)) for x in interval):
            raise ValueError('The "percentile" "interval" must be an array with min and max percentile')
        elif interval[0] >= interval[1]:
            raise ValueError('The "percentile" min interval must be lower than max interval')
        elif any(x < 0 or x > 100 for x in interval):
            raise ValueError('The "percentile" min and max interval must be between 0 and 100')
    elif "range" in limit_type:
        if interval == "":
            interval = [np.nanmin(dt_data.iloc[:, 4:].values), np.nanmax(dt_data.iloc[:, 4:].values)]
        if interval[0] == -9999:
             interval[0] = np.nanmin(dt_data.iloc[:, 4:].values)
        if interval[1] == -9999:
             interval[1] = np.nanmax(dt_data.iloc[:, 4:].values)        
        elif not isinstance(interval, list) or len(interval) != 2 or not all(isinstance(x, (int, float)) for x in interval):
            raise ValueError('The "range" "interval" must be an array with min and max values')
        elif interval[0] >= interval[1]:
            raise ValueError('The min interval must be lower than max interval')
