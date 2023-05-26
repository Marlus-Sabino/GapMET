import numpy as np
import pandas as pd

from LimCheck import LimCheck

def MPR(dt_ori, dt_mean_month, dt_months, limit_type, interval, min_est, index_est):
    
    dt = dt_ori
    dt_flag = np.zeros((np.size(dt_ori, 0), np.size(dt_ori, 1)))
    
    # Gap filling
    numNans_inicio = np.count_nonzero(np.isnan(dt))
    numNans_fim = 0
    iter = 1
    
    while numNans_fim < numNans_inicio:
        numNans_inicio = np.count_nonzero(np.isnan(dt))
        if numNans_inicio == 0:
            break
        print('Running Iteration', iter, '. Number of Gaps', numNans_inicio)
        
        for p in range(np.size(dt_ori, 1)):
            print('Gapfilling station', p + 1, 'of', np.size(dt_ori, 1))
            gap = dt.iloc[:, p]
            gap_filled = np.copy(gap)
            kf_all = index_est.iloc[p]
            kf = kf_all.dropna()
            X = np.full((dt_ori.shape[0],len(kf)), np.nan)           
            for j in range(len(kf)):
                if iter == 1:
                    X[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X[:,j] = dt.iloc[:, kf[j]]            
            for i in range(np.size(X, 0)):
                if np.isnan(dt.iloc[i, p]):
                    X1 = X[i, :]
                    ind = np.where(~np.isnan(X1))[0]
                    if len(ind) >= min_est:
                        m_ref = dt_months[i]
                        Xm = dt_mean_month.iloc[m_ref - 1,:].values
                        Xm = Xm[~pd.isna(kf_all)].astype(float)
                        Xm = Xm[~np.isnan(X1)]
                        X2 = X1[~np.isnan(X1)]
                        gap_filled[i] = (1 / len(Xm)) * (np.sum(X2 / Xm)) * dt_mean_month.iloc[m_ref - 1, p]
                        dt_flag[i, p] = iter
                        del X1, X2, Xm            
            gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
            dt.iloc[:, p] = gap_filled_checked
        
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter += 1
    
    dt_pree = np.copy(dt)
    dt_flag[np.isnan(dt_pree)] = np.nan
    
    if numNans_fim != 0:
        print(numNans_fim, 'Gaps left in the dataset')
    
    return dt_pree, dt_flag
