from crystfel_tools.handling import io_functions
from crystfel_tools.handling import data_conversion

streamin = '/home/esrf/hans1507/library/crystfel_tools/test_data/stream_cut_6.stream'
cell = [14.77,18.65,18.71,89.34,84.33,67.95]
def test_load_crystals(streamin):
    stream = io_functions.read_entire_stream(streamin)
    df, df_meta = io_functions.build_dataframe(stream)
    # df = data_conversion.apply_pg(df, '-1')
    df = data_conversion.get_resolution(df,cell)
    df = data_conversion.relabel_hkl(df)
    df.to_feather('/home/esrf/hans1507/library/crystfel_tools/test_data/test_df.feather')

test_load_crystals(streamin)