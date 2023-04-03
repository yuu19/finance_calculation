#二項ツリーモデルにおけるヨーロピアンコールオプションの計算

import numpy as np

def CRREcall(r, sigma, T, K, S0, N):
    S = np.zeros((N+1, N+1)) #株価行列の初期化
    C = np.zeros((N+1, N+1)) #オプション価格行列の初期化
    delta = T/N #時間間隔
    u = np.exp(sigma*np.sqrt(delta)) #上昇の大きさ
    d = np.exp(-sigma*np.sqrt(delta)) #下降の大きさ
    p = (np.exp(r*delta) - d)/(u - d) #リスク中立確率
    #-------- 株価ツリーの構築 -------------
    for i in range(N+1):
        for j in range(i+1):
            S[i,j] = S0*u**j*d**(i-j)
  
    
    #-------満期におけるオプションの支払いの計算 ----------
    for j in range(N+1):
        C[N,j] = max(S[N,j] - K, 0) #コールオプション
    #----- 後ろ向きにオプション価格を計算 -------------
    for i in range(N-1,-1,-1):
        for j in range(i+1):
            C[i,j] = np.exp(-r*delta)*(p*C[i+1,j+1] + (1-p)*C[i+1,j])
    return C[0,0] 
# CRREcall(0.1,0.2,5/12,60,62,300)


#二項ツリーモデルにおけるアメリカンプットオプションの計算(コールの場合はヨーロピアンと全く同じコードで実装できる)
def CRRAput(r, sigma, T, K, S0, N):
    S = np.zeros((N+1, N+1)) # 株価行列の初期化
    V = np.zeros((N+1, N+1)) # オプション価格行列の初期化
    delta = T/N # 時間間隔
    u = np.exp(sigma*np.sqrt(delta)) # 上昇の大きさ
    d = np.exp(-sigma*np.sqrt(delta)) # 下降の大きさ
    p = (np.exp(r*delta) - d)/(u - d) # リスク中立確率

    # 株価ツリーの構築
    for i in range(N+1):
        for j in range(i+1):
            S[i,j] = S0*u**j*d**(i-j)

    # 満期におけるオプションの支払いの計算
    for j in range(N+1):
        V[N,j] = max(K - S[N,j], 0) # プットオプション

    # 後ろ向きにオプション価格を計算
    for i in range(N-1, -1, -1):
        for j in range(i+1):
            V[i,j] = max(K - S[i,j], np.exp(-r*delta)*(p*V[i+1,j+1] + (1-p)*V[i+1,j])) # アメリカンオプション

    return V[0,0]
# CRRAput(0.1,0.2,5/12,60,62,100)
