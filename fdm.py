#-------有限差分法による数値計算---------

#---陽的解法---#
def exEcall(r, sigma, T, K, smax, M, N):
    Deltat = T/M #時間の刻み幅
    Deltas = smax/N #株価の刻み幅
    f = np.zeros((M+1, N+1)) #オプション価格行列
    #------ 満期におけるオプション価値 --------
    for j in range(1, N+2):
        f[M, j-1] = max((j-1)*Deltas - K, 0) #コールオプション
    #------ 後ろ向きにオプション価格を計算 -----
    for i in range(M-1, -1, -1):
        f[i, 0] = 0 #境界条件
        f[i, N] = smax - np.exp(-r*(T - i*Deltat))*K
        for j in range(1, N):
            a = (r*j + (sigma*j)**2)*Deltat/2
            b = 1 - ((sigma*j)**2 + r)*Deltat
            c = (-r*j + (sigma*j)**2)*Deltat/2
            f[i,j] = a*f[i+1,j+1] + b*f[i+1,j] + c*f[i+1,j-1]
    return f[0,]

#x = exEcall(0.1,0.2,5/12,60,300,5000,300)
#x[62]


#-----陰的解法--------
import numpy as np

def imEcall(r, sigma, T, K, smax, M, N):
    Deltat = T/M #時間の刻み幅
    Deltas = smax/N #株価の刻み幅
    f = np.zeros((M+1, N+1)) #オプション価格行列

    #-------- 行列の作成 ---------------
    A = np.zeros((N-1, N-1)) #行列の初期化
    a1 = (-r - sigma**2)*Deltat/2 #a_1
    b1 = 1 + (sigma**2 + r)*Deltat #b_1
    c1 = (r - sigma**2)*Deltat/2 #c_1
    A[0,0] = b1
    A[0,1] = a1
    for j in range(1, N-2):
        a = (-r*j - (sigma*j)**2)*Deltat/2
        b = 1 + ((sigma*j)**2 + r)*Deltat
        c = (r*j - (sigma*j)**2)*Deltat/2
        A[j,j+1] = a
        A[j,j] = b
        A[j,j-1] = c
    aNm1 = (-r*(N-1) - (sigma*(N-1))**2)*Deltat/2 #a_N - 1
    bNm1 = 1 + ((sigma*(N-1))**2+r)*Deltat #b_N - 1
    cNm1 = (r*(N-1) - (sigma*(N-1))**2)*Deltat/2 #c_N - 1
    A[N-2,N-2] = bNm1
    A[N-2,N-3] = cNm1

    #------ 満期におけるオプション価値 --------
    for j in range(N+1):
        f[M,j] = max(j*Deltas - K, 0) #コールオプション

    #------ 後ろ向きにオプション価格を計算 -----
    x = np.zeros(N-1) #左辺のベクトル
    y = np.zeros(N-1) #右辺のベクトル
    for i in range(M-1, -1, -1):
        bl = 0 #株価空間の境界条件
        bu = smax - np.exp(-r*(T - (i)*Deltat))*K
        f[i,0] = bl
        f[i,N] = bu
        y[0:N-1] = f[i+1,1:N]
        y[0] = y[0] - c1*bl
        y[N-2] = y[N-2] - aNm1*bu
        x = np.linalg.solve(A,y) #連立方程式を解く
        f[i,1:N] = x
    return f[0,:]
#x = imEcall(0.1,0.2,5/12,60,300,300,300)
#x[62]

#--------クランクニコルソン法
def CrankNicolson(r, sigma, T, K, smax, M, N):
    # 時間、株価の刻み幅を計算
    deltat = T/M
    deltas = smax/N
    
    # オプション価格行列を初期化
    f = np.zeros((M+1, N+1))
    
    # 係数の定義
    def a(j):
        return ((-r*j - (sigma*j)**2)*deltat/4)
    def b(j):
        return (1 + ((sigma*j)**2 + r)*deltat/2)
    def c(j):
        return ((r*j - (sigma*j)**2)*deltat/4)
    def d(j):
        return (b(j)-2)
    
    # 行列Aを初期化
    A = np.zeros((N-1, N-1))
    A[0,0] = b(1)
    A[0,1] = a(1)
    for j in range(1, N-2):
        A[j,j+1] = a(j+1)
        A[j,j] = b(j+1)
        A[j,j-1] = c(j+1)
    A[N-2,N-2] = b(N-1)
    A[N-2,N-3] = c(N-1)
    
    # 満期におけるオプション価値を計算
    for j in range(N+1):
        f[M,j] = max(j*deltas - K, 0)
        
    # 後ろ向きにオプション価格を計算
    x = np.zeros(N-1)
    y = np.zeros(N-1)
    for i in range(M-1, -1, -1):
        bl = 0 # 株価空間の境界条件
        bu = smax - np.exp(-r*(T - i*deltat))*K
        
        # yベクトルを生成
        for j in range(1, N):
            y[j-1] = -a(j)*f[i+1,j+1] - d(j)*f[i+1,j] - c(j)*f[i+1,j-1]
        y[0] -= c(1)*bl
        y[-1] -= a(N-1)*bu
        
        # 連立方程式を解いて、xベクトルを求める
        x = np.linalg.solve(A, y)
        
        # オプション価格を更新
        f[i,0] = bl
        f[i,N] = bu
        f[i,1:N] = x
    
    return f[0,:]
#x = CrankNicolson(0.1,0.2,5/12,60,200,200,200)
#x[62]


#------ アメリカンプットオプションの価格を有限差分法により計算 -------
def imAput(r, sigma, T, K, smax, M, N):
    Deltat = T/M 
    Deltas = smax/N 
    f = np.zeros((M+1, N+1)) #オプション価格行列    
    #---- 行列Aを初期化 -----
    A = np.zeros((N-1, N-1)) 
    a1 = (-r - sigma**2)*Deltat/2 # a_1
    b1 = 1 + (sigma**2 + r)*Deltat # b_1
    c1 = (r - sigma**2)*Deltat/2 # c_1
    A[0,0] = b1
    A[0,1] = a1
    
    for j in range(1, N-2):
        a = (-r*j - (sigma*j)**2)*Deltat/2
        b = 1 + ((sigma*j)**2 + r)*Deltat
        c = (r*j - (sigma*j)**2)*Deltat/2
        A[j,j+1] = a
        A[j,j] = b
        A[j,j-1] = c
    
    aNm1 = (-r*(N-1) - (sigma*(N-1))**2)*Deltat/2 # a_N-1
    bNm1 = 1 + ((sigma*(N-1))**2+r)*Deltat # b_N-1
    cNm1 = (r*(N-1) - (sigma*(N-1))**2)*Deltat/2 # c_N-1
    A[N-2,N-2] = bNm1
    A[N-2,N-3] = cNm1
    
    #---- オプション価格の終端条件 ------
    for j in range(N+1):
        f[M, j] = max(K - j*Deltas, 0) # Put option
    
    #---- 後ろ向きにオプション価格を計算 -----
    x = np.zeros(N-1) # left side vector
    y = np.zeros(N-1) # right side vector
    payoff = np.maximum(K - np.arange(0, N*Deltas+Deltas, Deltas), 0) # ペイオフ
    for i in range(M-1, -1, -1):
        bl = np.exp(-r*(T - i*Deltat))*K #
        bu = 0 
        f[i,0] = bl
        f[i,N] = bu
        y[:N-1] = f[i+1,1:N]
        y[0] = y[0] - c1*bl
        y[N-2] = y[N-2] - aNm1*bu
        x = np.linalg.solve(A,y) #連立方程式を解く
        f[i,1:N] = x
        f[i,:] = np.maximum(payoff, f[i,:]) # ここの部分がヨーロピアンとの違い
    
    return f[0,:]
#x = imAput(0.1,0.2,5/12,60,100,100,100)
#x[62]
