import numpy as np

from LimCheck import LimCheck

def IDM(dt_ori, dt_dist, limit_type='range', interval=None, min_est=3, index_est=None):

    dist = np.array(dt_dist)
    dt = np.array(dt_ori)
    dt_flag = np.zeros((np.size(dt_ori, 0), np.size(dt_ori, 1)))

    # Gap filling
    num_nans_inicio = np.count_nonzero(np.isnan(dt))
    num_nans_fim = 0
    iter = 1

    while num_nans_fim < num_nans_inicio:
        num_nans_inicio = np.count_nonzero(np.isnan(dt))
        if num_nans_inicio == 0:
            break
        print(f'Running Iteration {iter}. Number of Gaps {num_nans_inicio}')
        for p in range(np.size(dt_ori, 1)):
            print(f'Gapfilling station {p + 1} of {np.size(dt_ori, 1)}')
            gap = dt[:, p]
            gap_filled = gap.copy()
            kf_all = index_est.iloc[p]
            kf = kf_all.dropna()
            X = np.full((dt_ori.shape[0],len(kf)), np.nan)   
            for j in range(len(kf)):
                if iter == 1:
                    X[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X[:,j] = dt[:, kf[j]]                
            for i in range(len(X)):
                if np.isnan(dt[i, p]):
                    X1 = X[i]
                    ind = np.where(~np.isnan(X1))[0]
                    if len(ind) >= min_est:
                        Xd = np.full((len(kf)), np.nan)  
                        for n in range(len(X1)):
                            Xd[n] = dist[p, kf[n]]
                        X2 = X1[~np.isnan(X1)]
                        Xd = Xd[~np.isnan(X1)]
                        gap_filled[i] = np.sum(X2 / Xd) / np.sum(1 / Xd)
                        dt_flag[i, p] = iter
            gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
            dt[:, p] = gap_filled_checked
        num_nans_fim = np.count_nonzero(np.isnan(dt))
        iter += 1
    
    # Finishing gap-filled dataset
    dt_pree = dt
    dt_flag[np.isnan(dt_pree)] = np.nan

    if num_nans_fim != 0:
        print(f'{num_nans_fim} Gaps left on the dataset')

    return dt_pree, dt_flag
