def python_discretizer(data, i, layer_depth_vector, bin_number):
    # Create depth bins to divide data
    import numpy as np
    import pdb

    # Two cases, based on whether the first layer is processed or not. Intervals are defined in order to always include
    # all interval extremes
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

