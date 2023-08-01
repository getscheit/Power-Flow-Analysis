import numpy as np
import pandas as pd
import cmath



def common_parameter():
    # いずれファイル読み込みからの一般化をする！

    #Reading csv file
    line_data = pd.read_csv("Line Data.csv", header = 0, index_col = 0)
    bus_data = pd.read_csv("Bus Data.csv",  header = 0, index_col = 0)
    print(line_data)
    print(bus_data)

    n = len(bus_data) # ノード数
    l = len(line_data) # number of transmission lines


    r = np.inf * np.ones((n,n)) # 抵抗
    x = np.zeros((n,n)) # インダクタンス
    b = np.zeros((n,n)) # サセプタンス
    bc = np.zeros(n) # 容量サセプタンス
    P = np.zeros(n)
    Q = np.zeros(n)
    bus_type = np.zeros(n)

    for L in range(l):
        r[(line_data['NF'].iat[L] - 1), (line_data['NT'].iat[L] - 1)] = line_data['R'].iat[L]
        r[(line_data['NT'].iat[L] - 1), (line_data['NF'].iat[L] - 1)] = line_data['R'].iat[L]  
        x[(line_data['NF'].iat[L] - 1), (line_data['NT'].iat[L] - 1)] = line_data['X'].iat[L]
        x[(line_data['NT'].iat[L] - 1), (line_data['NF'].iat[L] - 1)] = line_data['X'].iat[L]
        b[(line_data['NF'].iat[L] - 1), (line_data['NT'].iat[L] - 1)] = line_data['BC'].iat[L]
        b[(line_data['NT'].iat[L] - 1), (line_data['NF'].iat[L] - 1)] = line_data['BC'].iat[L]

    print(r)
    print(x)
    print(b)

    for N in range(n):
        bc[N] = bus_data['SC'].iat[N]
        
        if bus_data['TYPE'].iat[N] == 0: #Slack bus
            P[N] = 0
            Q[N] = 0
        elif bus_data['TYPE'].iat[N] == 1: #Generator Node
            P[N] = bus_data['PG'].iat[N]
            Q[N] = 0
            bus_type[N] = 1
        elif bus_data['TYPE'].iat[N] == 2: #Load Node
            P[N] = bus_data['PL'].iat[N]
            Q[N] = bus_data['QL'].iat[N]
            bus_type[N] = 2

    print(P)
    print(Q)

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