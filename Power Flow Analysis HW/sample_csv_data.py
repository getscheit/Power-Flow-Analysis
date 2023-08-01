import pandas as pd


linedata = {'NF': [1, 2, 3, 4],
            'NT': [2, 3, 4, 5],
            'R': [0.3, 0.5, 0.4, 0.4],
            'X':[0.4, 0.6, 0.2, 0.1],
            'BC':[0.05, 0.01, 0.0, 0.0]}

linedata_df = pd.DataFrame.from_dict(linedata)
print(linedata_df)
linedata_df.to_csv('Line Data.csv')

busdata = {'TYPE': [0, 1, 2, 2, 2],
            'VOLT': [1, 1, 1, 1, 1],
            'DEG': [0.0, 0.0, 0.0, 0.0, 0.0],
            'PG':[0.4, 0.6, 0.2, 0.3, 0.1],
            'QG':[0.3, 0.1, 0.2, 0.2, 0.2],
            'PL':[0.5, 0.6, 0.6, 0.2, 0.1],
            'QL':[0.5, 0.3, 0.2, 0.10, 0.12],
            'SC':[0.0, 0.0, 0.5, 0.6, 0.4]
            }

busdata_df = pd.DataFrame.from_dict(busdata)
print(busdata_df)
busdata_df.to_csv('Bus Data.csv')