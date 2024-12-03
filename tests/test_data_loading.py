from crystfel_tools.handling import io_functions
from crystfel_tools.handling import data_conversion
pathout = '/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/test_df.feather'
streamin = '/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/distance_optimization_run_3_min_pix_2_int_2_4_6.stream'
cell = [14.77,18.65,18.71,89.34,84.33,67.95]
cell = [61.51,91.01,151.11,90,90,90]
def test_load_crystals(streamin):
    stream = io_functions.read_entire_stream(streamin)
    df, df_meta = io_functions.build_dataframe(stream)
    # df = data_conversion.apply_pg(df, '-1')
    df = data_conversion.get_resolution(df,cell)
    df = data_conversion.relabel_hkl(df)
    df.to_feather(pathout)

test_load_crystals(streamin)