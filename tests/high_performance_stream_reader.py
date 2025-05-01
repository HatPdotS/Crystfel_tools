from crystfel_tools.handling import pystream



test_stream = '/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/distance_optimization_run_3_min_pix_2_int_2_4_6.stream'

instance = pystream.pystream(test_stream)

instance.write_pkl('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/test.pkl')
