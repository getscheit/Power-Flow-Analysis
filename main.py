import numpy as np
import pf_fnc as pf

# データ読み込み
r, x, b, bc, P, Q, n = pf.common_parameter()

# 初期化
theta = np.zeros(4)
V = np.zeros(4)
Y = np.zeros((4,4), dtype=complex)
cnt = 0
p = np.array([P[1], Q[1], P[2], Q[2], P[3]])
fnc_v = np.zeros(5)

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
Jacobian = np.zeros((5,5))

# フラットスタート
V = np.array([1.0, 1.0, 1.0, 1.0]) # 複素電圧の大きさ[pu]
theta = np.array([0.0, 0.0, 0.0, 0.0]) # 位相角[rad]

# 未知数ベクトル
v = [0.0, 1.0, 0.0, 1.0, 0.0]

# ミスマッチベクトルの無限大ノルムが閾値0.001以下になれば計算終了
while np.linalg.norm(p-fnc_v, np.inf) > 0.001:
    # 初期化
    fP = np.zeros(n) # fPiの配列
    fQ = np.zeros(n) # fQiの配列
    # fPiの計算
    for ii in range(n):
        summ=0.0 # テンポラリ変数の初期化

        for jj in range(n):
            summ = summ + V[jj]*(G[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy()) + B[ii,jj] * np.sin(theta[ii].copy()-theta[jj].copy()))
        fP[ii] = V[ii].copy() * summ

    # fQiの計算
    for ii in range(n):
        summ=0.0
        for jj in range(n):
            summ = summ + V[jj].copy()*(G[ii,jj]*np.sin(theta[ii].copy()- theta[jj].copy()) - B[ii,jj]*np.cos(theta[ii].copy()-theta[jj].copy()))
        fQ[ii] = V[ii].copy() * summ

    # ヤコビアン行列の各成分の計算
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

    Jacobian = np.array([[dfPi_dthetaj[1,1], dfPi_dVj[1,1], dfPi_dthetaj[1,2], dfPi_dVj[1,2], 0],
                        [dfQi_dthetaj[1,1], dfQi_dVj[1,1], dfQi_dthetaj[1,2],dfQi_dVj[1,2], 0],
                        [dfPi_dthetaj[2,1], dfPi_dVj[2,1], dfPi_dthetaj[2,2],dfPi_dVj[2,2], dfPi_dthetaj[2,3]],
                        [dfQi_dthetaj[2,1], dfQi_dVj[2,1], dfQi_dthetaj[2,2],dfQi_dVj[2,2], dfQi_dthetaj[2,3]],
                        [0, 0, dfPi_dthetaj[3,2], dfPi_dVj[3,2],dfPi_dthetaj[3,3]]])
    
    fnc_v = np.array([fP[1], fQ[1], fP[2], fQ[2], fP[3]])
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