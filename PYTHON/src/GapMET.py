import numpy as np
import pandas as pd

from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
from pymap3d.vincenty import vdist

from CheckParameter import validate_interval, validate_limit_type, validate_max_dis, validate_min_est, validate_pree_meth, validate_ref_type
from IDM import IDM
from MAS import MAS
from MPR import MPR
from MUK import MUK
from RLM import RLM
from RLS import RLS
from libraries import check_libraries

def GapMet(dt_coord, dt_data, pree_meth, min_est, max_dis, limit_type, interval, ref_type, dt_extern):

    check_libraries()
    
    # SECTION 1: CHECK PARAMETERS
    if 'ref_type' not in locals():
        ref_type = 1
    else: 
        ref_type = validate_ref_type(ref_type)
    
    if 'min_est' not in locals():
        min_est = 3
    else:
        validate_min_est(min_est)

    validate_pree_meth(pree_meth)
    
    if 'max_dis' not in locals():
        max_dis = 100
    else:
        validate_max_dis(max_dis)

    if 'limit_type' not in locals():
        limit_type = "range"
    else:
        validate_limit_type(limit_type)

    if 'interval' not in locals():
        if "percentile" in limit_type:
            interval = [5, 95]
        elif "range" in limit_type:
            interval = [np.nanmin(dt_data.iloc[:, 4:].values), np.nanmax(dt_data.iloc[:, 4:].values)]
    else:
        validate_interval(interval, limit_type, dt_data)
    
    # SECTION 2: CHECK DATASETS
    if dt_data.shape[1] < (min_est + 4):
        raise ValueError("The number of columns (n) in 'dt_data' must be equal to or higher than n = min_est + 4")

    if dt_coord.shape[1] != 3:
        raise ValueError("The number of columns (m) in 'dt_coord' must be equal to 3 (station_id, latitude, longitude)")

    if dt_coord.shape[0] != (dt_data.shape[1] - 3):
        raise ValueError("The number of stations (rows) in 'dt_coord' is different than the number of stations (columns) in 'dt_data'")

    if ref_type == 2:
        if dt_data.shape[1] != dt_extern.shape[1]:
            raise ValueError("The number of stations (columns) in 'dt_satel' is different than the number of stations (columns) in 'dt_data'")
        if dt_data.shape[0] != dt_extern.shape[0]:
            raise ValueError("The timeseries size (time steps) in 'dt_satel' do not match the timeseries size (time steps) in 'dt_data'")

    dt_ori = dt_data.iloc[:, 3:]
    if ref_type == 2:
        dt_ext = dt_extern.iloc[:, 3:]

    dt_cod = dt_coord.iloc[:, 0].values
    dt_lat = dt_coord.iloc[:, 1].values
    dt_lon = dt_coord.iloc[:, 2].values

    if np.any(dt_lat > 90) and np.any(dt_lat < -90) and np.any(dt_lon > 180) and np.any(dt_lon < -180):
        lat_inf = np.unique(dt_cod[dt_lat < -90])
        lat_sup = np.unique(dt_cod[dt_lat > 90])
        lon_inf = np.unique(dt_cod[dt_lon < -180])
        lon_sup = np.unique(dt_cod[dt_lon > 180])
        est_err = np.unique(np.concatenate((lat_inf, lat_sup, lon_inf, lon_sup)))

        if est_err.size > 0:
            est_err_str = ",".join(est_err.astype(str))
            raise ValueError(f"The geographic coordinates in the station {est_err_str} must be between -90<lat<90 and -180<lon<180")

    # SECTION 3: CALCULATE THE DISTANCE BETWEEN DATA SERIES
    # Calculate the distances 
    dt_dist = np.zeros((len(dt_lat), len(dt_lat)))
    for i in range(len(dt_lat)):
        for j in range(len(dt_lat)):
            if i != j :
                lat1 = dt_lat[i]
                lon1 = dt_lon[i]
                lat2 = dt_lat[j]
                lon2 = dt_lon[j]
                dist_km, _ = vdist(lat1, lon1, lat2, lon2)   # Calculate distance in kilometers
                dt_dist[i, j] = dist_km/1000
            else:
                dt_dist[i, j] = 0         
    dt_dist = pd.DataFrame(dt_dist, index=dt_cod, columns=dt_cod)

    # Validate the distances with seted max_dis and min_est parameters
    if ref_type == 1:
        dt = dt_dist.copy()
        dt[dt == 0] = np.nan

        min_dis = np.round(np.nanmin(dt))
        if max_dis < min_dis:
            raise ValueError(f"The shortest distance to all stations has at least 1 reference station within the 'max_dis' is {min_dis} km")

        dt1 = dt.copy()
        dt1[dt1 <= max_dis] = 1
        dt1[dt1 > max_dis] = 0
        min_est_in_max_dis = np.nanmin(np.nansum(dt1, axis=1))
        dt1 = np.sort(dt, axis=1)
        max_dis_in_min_est = np.round(np.nanmax(dt1[:, min_est]))

        if max_dis < max_dis_in_min_est and min_est > min_est_in_max_dis:
            raise ValueError("The number of reference stations within the 'max_dis' radius is lower than set in 'min_est'. Decrease the 'min_est' parameter to "
                            f"{min_est_in_max_dis} stations or increase 'max_dis' to {max_dis_in_min_est} km")

    # SECTION 4: CALCULATE THE CORRELATION BETWEEN DATA SERIES
    if ref_type == 1:
        dt_corr = np.full((len(dt_cod), len(dt_cod)), np.nan)
        
        for i in range(dt_ori.shape[1]):
            for j in range(dt_ori.shape[1]):
                if i != j :
                    ponts = pd.DataFrame()
                    ponts[0] = dt_ori.iloc[:, i]
                    ponts[1] = dt_ori.iloc[:, j]
                    ponts = ponts.dropna(subset=[0, 1])
                    rho, pval = pearsonr(ponts[0], ponts[1])
                    # Check if the correlations are significant at 0.05
                    if pval <= 0.05:
                        dt_corr[i, j] = rho
                    else:
                        dt_corr[i, j] = np.nan
                else:
                    dt_corr[i, j] = np.nan        
        # Create correlations matrix table
        dt_corr = pd.DataFrame(dt_corr, index=dt_cod, columns=dt_cod)

    elif ref_type == 2:
        dt_corr = np.full((len(dt_cod), 1), np.nan)

        for i in range(dt_ori.shape[1]):
            ponts = pd.DataFrame()
            ponts[0] = dt_ori.iloc[:, i]
            ponts[1] = dt_ext.iloc[:, j]
            ponts = ponts.dropna(subset=[0, 1])
            rho, pval = pearsonr(ponts[0], ponts[1])
            # Check if the correlation is significant at 0.05
            if pval <= 0.05:
                dt_corr[i, 0] = rho
            else:
                dt_corr[i, 0] = np.nan
        # Create correlations matrix table
        dt_corr = pd.DataFrame(dt_corr, index=dt_cod)

    # SECTION 5: DEFINE REFERENCE SERIES BASED ON THE DISTANCE AND CORRELATION
    if ref_type == 1:
        dt = dt_dist.replace(0, np.nan).to_numpy()
        dt1 = dt_corr.replace(np.nan, 0).to_numpy()
        dt1[dt1 == 1] = 0
        dt1[dt > max_dis] = 0
        dt1 = np.abs(dt1)
        M = np.sort(dt1, axis=1)[:, ::-1]
        I = np.argsort(dt1, axis=1)[:, ::-1]
        I = pd.DataFrame(I, dtype=object)
        M = pd.DataFrame(M, dtype=object)
        I[M == 0] = np.nan
        index_est = I

    if pree_meth in ("MPR", "MUK"):
        dt_months = dt_data.iloc[:, 1].to_numpy()
        table_mes=dt_data.iloc[:, 1:]
        table_mes.drop(columns=['day'], inplace=True)
        dt_mean_month = table_mes.groupby(['month']).mean().to_numpy()
        dt_mean_month = pd.DataFrame(dt_mean_month, index=list(range(1, 13)), columns=dt_cod)
        
        if ref_type == 2:
            table_mes=dt_extern.iloc[:, 1:]
            table_mes.drop(columns=['day'], inplace=True)
            dt_mean_month_ext = table_mes.groupby(['month']).mean().to_numpy()
            dt_mean_ext = pd.DataFrame(dt_mean_month_ext , index=list(range(1, 13)), columns=dt_cod)

    # SECTION 6: GAP-FILLING
    # RLS: Simple linear regression
    if pree_meth == "RLS":
        if ref_type == 1:
            dt_pree, dt_flag = RLS(dt_ori, limit_type, interval, ref_type, index_est, [])
        else:
            dt_pree, dt_flag = RLS(dt_ori, limit_type, interval, ref_type, [], dt_ext)

    # RLM: Multiple linear regression
    if pree_meth == "RLM":
        if ref_type == 1:
            dt_pree, dt_flag = RLM(dt_ori, limit_type, interval, min_est, index_est)
        else:
            raise ValueError('This methodology (RLM) can only be used with nearby reference stations (Ref_type="nearby")')
    
    # MUK: UK traditional method
    if pree_meth == "MUK":
        if ref_type == 1:
            dt_pree, dt_flag = MUK(dt_ori, dt_mean_month, dt_months, limit_type, interval, ref_type, index_est, '','')
        else:
            dt_pree, dt_flag = MUK(dt_ori, dt_mean_month, dt_months, limit_type, interval, ref_type, '', dt_ext, dt_mean_ext)

    # MAS: Simple arithmetic mean   
    if pree_meth == "MAS":
        if ref_type == 1:
            dt_pree, dt_flag = MAS(dt_ori, limit_type, interval, min_est, index_est)
        else:
            raise ValueError("This methodology (MAS) can only be used with nearby reference stations (Ref_type='nearby')")

    #  MPR: Regional weighting method
    if pree_meth == "MPR":
        if ref_type == 1:
            dt_pree, dt_flag = MPR(dt_ori, dt_mean_month, dt_months, limit_type, interval, min_est, index_est)
        else:
            raise ValueError('This methodology (MPR) can only be used with nearby reference stations (Ref_type="nearby")')

    # IDM: Inverse distance method
    if pree_meth == "IDM":
        if ref_type == 1:
            dt_pree, dt_flag = IDM(dt_ori, dt_dist, limit_type, interval, min_est, index_est)
        else:
            raise ValueError('This methodology (IDM) can only be used with nearby reference stations (Ref_type="nearby")')

    return dt_pree, dt_flag, dt_dist, dt_corr
