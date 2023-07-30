import numpy as np
import pandas as pd
import cmath



def common_parameter():
    # いずれファイル読み込みからの一般化をする！

    #Reading csv file
    line_data = pd.read_csv(input("Enter power line data file name: "))
    bus_data = pd.read_csv(input("Enter bus data file name: "))

    r = np.inf * np.ones((n,n)) # 抵抗
    x = np.zeros((n,n)) # インダクタンス
    b = np.zeros((n,n)) # サセプタンス
    bc = np.zeros(n) # 容量サセプタンス
    P = np.zeros(n)
    Q = np.zeros(n)

    n = len(bus_data) # ノード数
    l = len(line_data) # number of transmission lines


    for L in range(l):
        r[line_data['send'].iat[L], line_data['receive'].iat[L]] = line_data['抵抗'].iat[L]  
        x[line_data['send'].iat[L], line_data['receive'].iat[L]] = line_data['インダクタンス'].iat[L]
        b[line_data['send'].iat[L], line_data['receive'].iat[L]] = line_data['サセプタンス'].iat[L]
        
    for N in range(n):
        P[n] = bus_data['有効電力'].iat[N]
        Q[n] = bus_data['無効電力'].iat[N]
        bc[n] = bus_data['容量サセプタンス'].iat[N]
    
    return r, x, b, bc, P, Q, n

def node_calc(V,theta,r,x,b,n):
    I_dash = np.zeros((n,n), dtype=complex) # ノードiからノードjに向かい流出する電流
    Power = np.zeros((n,n), dtype=complex) # ノードiからノードjに向かい流出する電力潮流
    V_dot = np.zeros(n, dtype=complex)

    #V_dot = V * np.exp(1.0j * theta) # ノード電圧（複素表示）
    for ii in range(n):
        V_dot[ii] = V[ii] * cmath.rect(1.0, theta[ii])

    # ブランチ潮流の計算
    for ii in range(n):
        for jj in range(n):
         I_dash[ii,jj] = -1.0j * b[ii,jj] / 2.0 * V_dot[ii] + (V_dot[ii]- V_dot[jj]) / (r[ii,jj] + 1.0j * x[ii,jj])
   
    for ii in range(n):
        for jj in range(n):
            Power[ii,jj] = V_dot[ii] * I_dash[ii,jj].conjugate()

    return I_dash, Power
