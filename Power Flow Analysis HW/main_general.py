import numpy as np
import pandas as pd
import pf_fnc_general as pf

# データ読み込み
r, x, b, bc, P, Q, n = pf.common_parameter()

# 初期化
theta = np.zeros(n)
V = np.zeros(n)
Y = np.zeros((n,n), dtype=complex)
cnt = 0

p = []
for N in range(n):     
    if P[N] != 0: 
        p.append(P[N])      
    if Q[N] != 0: 
        p.append(Q[N])
p = np.array(p)
print(p)


fnc_v = np.zeros(len(p))

# ノードアドミタンス行列
for jj in range(n):
    for ii in range(n):
        summ = 0.0
        temp = 0.0 # テンポラリ変数初期化
        
        if (ii==jj):
            for jj2 in range(n):
                if (ii!=jj2): # 自ノード除去
                    temp = 1.0/(r[ii,jj2] + 1.0j*x[ii,jj2]) + 1.0j*b[ii,jj2]/2
                    summ = summ + temp
            
            Y[ii,jj] = summ + bc[ii]*1.0j
            
        if (ii!=jj):
            Y[ii,jj] = -1.0/(r[ii,jj] + 1.0j*x[ii,jj])
            
            
G = Y.real
B = Y.imag

## 潮流計算
dfPi_dthetaj = np.zeros((n,n))
dfQi_dthetaj = np.zeros((n,n))
dfPi_dVj = np.zeros((n,n))
dfQi_dVj = np.zeros((n,n))
Jacobian = np.zeros((len(p),len(p)))

# フラットスタート
V = np.array([1.0]*n) # 複素電圧の大きさ[pu]
theta = np.array([0.0]*n) # 位相角[rad]

# 未知数ベクトル
v = []
for N in range(n):
    if P[N] != 0: 
        v.append(theta[N])      
    elif Q[N] != 0: 
        v.append(V[N])
v = np.array(v)
print(v)

# ミスマッチベクトルの無限大ノルムが閾値0.001以下になれば計算終了
while np.linalg.norm(p-fnc_v, np.inf) > 0.001:
    # 初期化
    fP = np.zeros(n) # fPiの配列
    fQ = np.zeros(n) # fQiの配列
    
    # fPiの計算
    for ii in range(n):
        summ = 0.0 # テンポラリ変数の初期化
        
        for jj in range(n):
            summ = summ + V[jj]*(G[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy()) + B[ii,jj] * np.sin(theta[ii].copy()-theta[jj].copy()))
        
        fP[ii] = V[ii].copy() * summ
    
    # fQiの計算
    for ii in range(n):
        summ=0.0
        
        for jj in range(n):
            summ = summ + V[jj].copy()*(G[ii,jj]*np.sin(theta[ii].copy()-theta[jj].copy()) - B[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy()))
        
        fQ[ii] = V[ii].copy() * summ
        
    # ヤコビアン行列の各成分の計算
    Jacobian = {}
    Jacobian = pd.DataFrame(np.zeros((p, p)))
    
    for jj in range(n):
        for ii in range(n):
            if (ii==jj):
                dfPi_dthetaj[ii,jj] = -V[ii].copy()**2 * B[ii,jj] - fQ[ii]
                dfQi_dthetaj[ii,jj] = -V[ii].copy()**2 * G[ii,jj] + fP[ii]
                dfPi_dVj[ii,jj] = V[ii].copy() * G[ii,jj] + fP[ii]/V[ii]
                dfQi_dVj[ii,jj] = -V[ii].copy() * B[ii,jj] + fQ[ii]/V[ii]
            if (ii!=jj):
                dfPi_dthetaj[ii,jj] = V[ii].copy()*V[jj].copy()*(G[ii,jj]*np.sin(theta[ii].copy()-theta[jj].copy())-B[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy()))
                dfQi_dthetaj[ii,jj] = -V[ii].copy()*V[jj].copy()*(G[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy())+B[ii,jj]*np.sin(theta[ii].copy()-theta[jj].copy()))
                dfPi_dVj[ii,jj] = V[ii].copy() * (G[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy())-B[ii,jj]*np.sin(theta[ii].copy()-theta[jj].copy()))
                dfQi_dVj[ii,jj] = V[ii].copy() * (-G[ii,jj]*np.sin(theta[ii].copy()-theta[jj].copy())-B[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy()))
    
    
    for jj in range(fnc_v):
        for ii in range(fnc_v):
            Jacobian.iat[jj,ii] = (dfPi_dthetaj[jj,ii])
                        
    
    fnc_v = []
    for N in range(n):
        if P[N] != 0: 
            fnc_v.append(fP[N])      
        elif Q[N] != 0: 
            fnc_v.append(fP[N])
    fnc_v = np.array(fnc_v)
            
    temp = v + np.dot(np.linalg.inv(Jacobian),(p-fnc_v))
    v = temp.copy()
    
    V = [V[0], v[1], v[3], V[3]]
    theta = [theta[0], v[0], v[2], v[4]]
    cnt = cnt+1
    
# 結果表示
print('Y Matrix Node Admittance')
print(Y)
print('Iteration count: ' + str(cnt))
print('Jacobian Matrix 最終値')
print(Jacobian)
print('未知数ベクトルの最終値')
print(v)
print('ミスマッチベクトルの最終値')
print(p-fnc_v)
print('計算結果')
print('ノード電圧')
print(V)
print('位相差')
print(theta)

# ブランチ潮流の計算
I_dash, Power = pf.node_calc(V, theta, r, x, b, n)

print('I_dash:')
print(I_dash)
print('Power: P[i,j]+jQ[i,j]')
print(Power)