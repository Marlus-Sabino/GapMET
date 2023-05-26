import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from LimCheck import LimCheck

def RLS(dt_ori, limit_type, interval, ref_type, index_est, dt_ext):
        
    dt = dt_ori.copy()
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))
    #dt_flag = pd.DataFrame(np.zeros((dt_ori.shape[0], dt_ori.shape[1])))
    coef_lin = np.full((dt_ori.shape[1], dt_ori.shape[1]), np.nan)
    coef_ang = np.full((dt_ori.shape[1], dt_ori.shape[1]), np.nan)
    R2 = np.full((dt_ori.shape[1], dt_ori.shape[1]), np.nan)

    # RLS: Using nearby stations to gap-filling (ref_type==1)
    if ref_type == 1:
        for p in range(dt_ori.shape[1]):
            kf = index_est.iloc[p]
            kf = kf.dropna()
            for k in range(len(kf)):
                df = pd.DataFrame()
                df[0] = dt.iloc[:, p]
                df[1] = dt.iloc[:, int(index_est.iloc[p, k])]
                df = df.dropna(subset=[0, 1])
                gap = df[0]
                ref = df[1]
                ref_array = ref.values.reshape(-1, 1)
                gap_array = gap.values.reshape(-1, 1)
                mdl = LinearRegression().fit(ref_array, gap_array)
                coef_lin[p][k] = mdl.intercept_[0]
                coef_ang[p][k] = mdl.coef_[0][0]
                R2[p][k] = mdl.score(ref_array, gap_array)

        # Gap-filling
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter = 1

        if numNans_fim > 0:
            numNans_inicio = np.count_nonzero(np.isnan(dt))
            numNans_fim = 0
            while numNans_fim < numNans_inicio:
                numNans_inicio = np.count_nonzero(np.isnan(dt))
                if numNans_inicio == 0:
                    break
                print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
                for p in range(dt_ori.shape[1]):
                    print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
                    kf = index_est.iloc[p]
                    kf = kf.dropna()
                    for k in range(len(kf)):
                        if iter == 1:
                            ref = dt_ori.iloc[:, int(index_est.iloc[p, k])]
                        else:
                            ref = dt.iloc[:, int(index_est.iloc[p, k])]
                        gap = dt.iloc[:, p]
                        gap_filled = gap.copy()
                        for i in range(len(gap)-1):
                            if np.isnan(gap[i]):
                                y = coef_lin[p, k] + (coef_ang[p, k] * ref[i])
                                gap_filled[i] = y
                                dt_flag[i, p] = iter
                        gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                        dt.iloc[:, p] = gap_filled_checked
                numNans_fim = np.count_nonzero(np.isnan(dt))
                iter += 1

    # RLS: Using external time series (ref_type==2)
    elif ref_type == 2:
        dt = dt_ori.copy()
        dt = pd.concat([dt, dt_ext], axis=1)
        dt_flag = np.zeros((dt.shape[0], dt.shape[1]))
        #dt_flag = pd.DataFrame(np.zeros((dt_ori.shape[0], dt_ori.shape[1])))
        coef_lin = np.full((dt.shape[1], dt.shape[1]), np.nan)
        coef_ang = np.full((dt.shape[1], dt.shape[1]), np.nan)
        R2 = np.full((dt.shape[1], dt.shape[1]), np.nan)

        for p in range(dt_ori.shape[1]):
            for k in range(dt_ext.shape[1]):
                df = pd.DataFrame()
                df[0] = dt.iloc[:, p]
                df[1] = dt.iloc[:, dt_ori.shape[1] + k]
                df = df.dropna(subset=[0, 1])
                gap = df[0]
                ref = df[1]
                ref_array = ref.values.reshape(-1, 1)
                gap_array = gap.values.reshape(-1, 1)
                mdl = LinearRegression().fit(ref_array, gap_array)
                coef_lin[p][k] = mdl.intercept_[0]
                coef_ang[p][k] = mdl.coef_[0][0]
                R2[p][k] = mdl.score(ref_array, gap_array)

        # Gap-filling
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter = 1

        if numNans_fim > 0:
            numNans_inicio = np.count_nonzero(np.isnan(dt))
            numNans_fim = 0
            while numNans_fim < numNans_inicio:
                numNans_inicio = np.count_nonzero(np.isnan(dt))
                if numNans_inicio == 0:
                    break
                print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
                for p in range(dt_ori.shape[1]):
                    print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
                    for k in range(dt_ext.shape[1]):
                        if iter == 1:
                            ref = dt_ext.iloc[:, k]
                        else:
                            ref = dt.iloc[:, dt_ori.shape[1] + k]
                        gap = dt.iloc[:, p]
                        gap_filled = gap.copy()
                        for i in range(len(gap)):
                            if np.isnan(gap[i]):
                                y = coef_lin[p, k] + (coef_ang[p, k] * ref[i])
                                gap_filled[i] = y
                                dt_flag[i, p] = iter
                        gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                        dt.iloc[:, p] = gap_filled_checked
                numNans_fim = np.count_nonzero(np.isnan(dt))
                iter += 1

    # Finishing gap-filled dataset
    dt_pree = dt.copy()
    dt_flag[np.isnan(dt_pree)] = np.nan
    #dt_flag = dt_flag.mask(pd.isnull(dt_pree))

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')

    return dt_pree, dt_flag
