import numpy as np

from LimCheck import LimCheck

def MAS(dt_ori, limit_type, interval, min_est, index_est):

    dt = dt_ori
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))

    numNans_inicio = np.count_nonzero(np.isnan(dt))
    numNans_fim = 0
    iter = 1

    while numNans_fim < numNans_inicio:
        numNans_inicio = np.count_nonzero(np.isnan(dt))
        if numNans_inicio == 0:
            break
        print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
        for p in range(dt_ori.shape[1]):
            print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
            gap = dt.iloc[:, p]
            gap_filled = gap.copy()
            kf = index_est.iloc[p]
            kf = kf.dropna()
            X = np.full((dt_ori.shape[0],len(kf)), np.nan)
            for j in range(len(kf)):  
                if iter == 1:
                    X[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X[:,j] = dt.iloc[:, kf[j]]
                for i in range(len(gap)):
                    if np.isnan(dt.iloc[i, p]):
                        X1 = X[i]
                        ind = np.where(~np.isnan(X1))[0]
                        if len(ind) >= min_est:
                            if np.isnan(X1).all():
                                gap_filled[i] = np.nan
                            else:
                                gap_filled[i] = np.nanmean(X1)
                            dt_flag[i, p] = iter
                            del X1
                gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                dt.iloc[:, p] = gap_filled_checked
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter += 1

    dt_pree = dt.copy()
    dt_flag[np.isnan(dt_pree)] = np.nan

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')
    return dt_pree, dt_flag
