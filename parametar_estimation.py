#------ インプライドボラティリティの推定 ------
from scipy.optimize import brentq
from scipy.stats import norm
def IV(sigma, p, r, T, K, S0):
    d1 = (np.log(S0/K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return S0 * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2) - p

sigma = brentq(IV, 0.0001, 1, args=(115, 0.001, 14/365, 10500, 10395.18))
print(sigma)
