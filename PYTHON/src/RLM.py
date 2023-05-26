import numpy as np
import pandas as pd
import statsmodels.api as sm

from LimCheck import LimCheck

def RLM(dt_ori, limit_type, interval, min_est, index_est):

    dt = np.copy(dt_ori)
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))

    numNans_inicio = np.sum(np.isnan(dt))
    numNans_fim = 0
    iter = 1

    while numNans_fim < numNans_inicio:
        numNans_inicio = np.sum(np.isnan(dt))
        if numNans_inicio == 0:
            break

        print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
        for p in range(dt_ori.shape[1]):
            print(f'Gapfilling station {p + 1} of {dt_ori.shape[1]}')
            gap = dt[:, p]
            gap_filled = np.copy(gap)
            kf_all = index_est.iloc[p]
            kf = kf_all.dropna()
            X = pd.DataFrame(np.full((dt_ori.shape[0], len(kf)), np.nan))
            for j in range(len(kf)):
                if iter == 1:
                    X.iloc[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X.iloc[:,j] = dt[:, kf[j]]
            for i in range(len(X)):
                if np.isnan(dt[i, p]):
                    X1 = X.iloc[i, :]
                    ind = np.where(~np.isnan(X1))[0]
                    if len(ind) >= min_est:
                        X0 = X
                        X0['y'] = gap
                        X0 = X0.dropna()
                        X0 = sm.add_constant(X0)
                        X2 = X0.iloc[:,0:-1]
                        y = X0['y']
                        model = sm.OLS(y, X2)
                        results = model.fit()
                        coefficients = results.params[1:]
                        constant = results.params['const']
                        gap_filled[i] = np.nansum(X1 * coefficients) + constant
                        dt_flag[i, p] = iter

            gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
            dt[:, p] = gap_filled_checked

        numNans_fim = np.sum(np.isnan(dt))
        iter += 1

    dt_pree = np.copy(dt)
    dt_flag[np.isnan(dt_pree)] = np.nan

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')

    return dt_pree, dt_flag
