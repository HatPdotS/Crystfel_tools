from handling import io_functions
from handling import scaling
from handling import data_conversion
from time import time
import sparse
cell = [61.51,91.1,151.1,90,90,90]

def test_read_stream(file='/das/work/p17/p17490/Rhodopsin_data/processing_dark/python_based_processing/distance_optimization/Distance_optimisation_run_3/Distance_optimization/run_3.stream'):
    """
    Test reading the stream.

    """
    t = time()
    chunks = io_functions.read_entire_stream(file)
    print('read entire stream', time()-t)
    t = time()
    print('read stream')
    chunks = io_functions.filter_for_crystals(chunks)
    print('filter for crystals', time()-t)
    t = time()
    data, df = io_functions.build_dataframe(chunks)
    print('build dataframe', time()-t)
    t = time()
    data = data_conversion.apply_pg(data,pg)
    print('apply pg', time()-t)
    t = time()
    data = data_conversion.get_resolution(data,cell)
    print('get resolution', time()-t)
    t = time()
    data = data.loc[(data.resolution > 1.8) & (data.resolution < 15)]
    print('filter resolution', time()-t)
    t = time()
    data = data_conversion.relabel_hkl(data)
    data.to_feather('/das/work/p17/p17490/Peter/Library/crystfel_tools/tests/data.feather')
    print('relabel hkl', time()-t)
    print(data.idx_hkl.unique().shape)
    t = time()
    data, weights = data_conversion.make_sparse(data)
    print('make sparse', time()-t)
    print(data)
    sparse.save_npz('/das/work/p17/p17490/Peter/Library/crystfel_tools/tests/data.npz',data) 
    sparse.save_npz('/das/work/p17/p17490/Peter/Library/crystfel_tools/tests/weights.npz',weights)


pg = 'mmm'
test_read_stream()
