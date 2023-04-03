



#------アメリカンオプションのLSM法を用いた数値計算法-----
import numpy as np
from scipy.stats import norm

def LSMAput(r, sigma, T, K, S0, M, N):
    # S0 > K を仮定する
    # M : [0, T] の分割数
    # N : パスの数の半分
    Deltat = T/M #一区間の長さ

    # 株価パスの生成
    A = np.sqrt(Deltat) * np.tril(np.ones((M, M)))
    t = np.arange(1, M+1) / M * T
    Z = np.random.normal(size=(M, N))
    Z = np.hstack((Z, -Z)) #対称変数法を使う
    W = A @ Z
    S = S0 * np.exp((r - 0.5 * sigma**2) * t + sigma * W) #株価のパス

    # 満期におけるオプション価値
    Value = np.maximum(0, K - S[-1, :])
    for i in range(M-2, -1, -1):
        # オプション価格を割り引く
        Value = np.exp(-r*Deltat) * Value
        # インザマネーのパスの添字を求める
        IntheMoney = np.where(K > S[i, :])[0]
        # インザマネーの株価を抽出
        XData = S[i, IntheMoney]
        # 説明変数
        XData0 = np.exp(-0.5*XData)
        XData1 = np.exp(-0.5*XData) * (1 - XData0)
        XData2 = np.exp(-0.5*XData) * (1 - 2*XData + 0.5*XData**2)
        XData3 = np.exp(-0.5*XData) * (1 - 3*XData + 1.5*XData**2 - 1/6*XData**3)
        # 被説明変数
        YData = Value[IntheMoney]
        # 回帰分析を実行
        X = np.column_stack((XData0, XData1, XData2, XData3))
        coef = np.linalg.lstsq(X, YData, rcond=None)[0]
        # 継続価値を計算
        ContinuationValue = coef[0] + coef[1]*XData0 + coef[2]*XData1 + coef[3]*XData2 + coef[4]*XData3
        # 本源的価値
        IntrinsicValue = K - XData
        # 本源的価値 > 継続価値となるパスのオプション価値を更新
        Exercise = np.where(IntrinsicValue > ContinuationValue)[0]
        k = IntheMoney[Exercise]
        Value[k] = IntrinsicValue[Exercise]

    return np.mean(np.exp(-r * Deltat) * Value)
