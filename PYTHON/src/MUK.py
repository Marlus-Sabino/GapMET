import numpy as np

from LimCheck import LimCheck

def MUK(dt_ori, dt_mean_month, dt_months, limit_type, interval, ref_type, index_est, dt_ext, dt_mean_ext):

    dt = dt_ori
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))

    numNans_inicio = np.count_nonzero(np.isnan(dt))
    numNans_fim = 0
    iter = 1

    if ref_type == 1:
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
                    gap = dt.iloc[:, p]
                    gap_filled = gap.copy()
                    if iter == 1:
                        ref = dt_ori.iloc[:, kf[k]]
                    else:
                        ref = dt.iloc[:, kf[k]]
                    for i in range(len(gap)):
                        if np.isnan(gap[i]):
                            gap_filled[i] = ref[i] + (dt_mean_month.iloc[dt_months[i]-1, p] - dt_mean_month.iloc[dt_months[i]-1, kf[k]])
                            dt_flag[i, p] = iter
                    gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                    dt.iloc[:, p] = gap_filled_checked
            numNans_fim = np.count_nonzero(np.isnan(dt))
            iter += 1

    if ref_type == 2:
        while numNans_fim < numNans_inicio:
            numNans_inicio = np.count_nonzero(np.isnan(dt))
            if numNans_inicio == 0:
                break
            for p in range(dt_ori.shape[1]):
                print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
                gap = dt.iloc[:, p]
                gap_filled = gap.copy()
                ref = dt_ext.iloc[:, p]
                for i in range(len(gap)):
                    if np.isnan(gap[i]):
                        gap_filled[i] = ref[i] + (dt_mean_month.iloc[dt_months[i]-1, p] - dt_mean_ext.iloc[dt_months[i]-1, p])
                        dt_flag[i, p] = iter
                gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                dt.iloc[:, p] = gap_filled_checked
            numNans_fim = np.count_nonzero(np.isnan(dt))

    dt_pree = dt.copy()
    dt_flag[np.isnan(dt_pree)] = np.nan

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')
    
    return dt_pree, dt_flag
