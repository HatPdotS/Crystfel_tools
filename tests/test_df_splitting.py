import pandas as pd
from crystfel_tools.handling.data_conversion import split_df

df = pd.read_feather('/home/esrf/hans1507/library/crystfel_tools/test_data/test_df.feather')

dfs = split_df(df)





