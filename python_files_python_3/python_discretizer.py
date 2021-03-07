def python_discretizer(data, i, layer_depth_vector, bin_number):
    import numpy as np
    import pdb

    # bins = np.arange(data.min(), data.max()+(data.max()-data.min())/bin_number, (data.max()-data.min())/bin_number)
    # bins = np.arange(data.min(), data.max(), (data.max()-data.min())/bin_number)
    if i == 0:
        min_depth = i
        max_depth = layer_depth_vector[i]
        bins = np.arange(min_depth, max_depth + (max_depth - min_depth) / (2*bin_number), (max_depth - min_depth)
                         / bin_number)
    else:
        min_depth = layer_depth_vector[i-1]
        max_depth = layer_depth_vector[i]
        bins = np.arange(min_depth, max_depth + (max_depth - min_depth) / (2*bin_number), (max_depth - min_depth)
                         / bin_number)
    indices = np.digitize(data, bins)
    return indices, bins

