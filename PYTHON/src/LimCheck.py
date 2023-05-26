import numpy as np
import pandas as pd

def LimCheck(dt_ori, gap, gap_filled, limit_type, interval):
    gap_filled_checked = gap_filled.copy()

    if limit_type != "range" and limit_type != "percentile":
        raise ValueError('"limit_type" must be "range" or "percentile"')

    if limit_type == "percentile":
        if interval is None:
            interval = [5, 95]
        elif not isinstance(interval, list) or len(interval) != 2:
            raise ValueError('The "percentile" "interval" must be an array with min and max percentile')
        elif interval[0] >= interval[1]:
            raise ValueError('The "percentile" min interval must be lower than max interval')
        elif any(x < 0 or x > 100 for x in interval):
            raise ValueError('The "percentile" min and max interval must be between 0 and 100')

        min_lim = np.nanpercentile(dt_ori.values.flatten(), interval[0])
        max_lim = np.nanpercentile(dt_ori.values.flatten(), interval[1])

    elif limit_type == "range":
        if interval is None:
            interval = [dt_ori.min().min(), dt_ori.max().max()]
        elif not isinstance(interval, list) or len(interval) != 2:
            raise ValueError('The "range" "interval" must be an array with min and max values')
        elif interval[0] >= interval[1]:
            raise ValueError('The min interval must be lower than max interval')

        min_lim = interval[0]
        max_lim = interval[1]

    for line_check in range(len(dt_ori)):
        if pd.isna(gap[line_check]):
            if gap_filled[line_check] <= min_lim or gap_filled[line_check] >= max_lim:
                gap_filled_checked[line_check] = np.nan

    return gap_filled_checked
