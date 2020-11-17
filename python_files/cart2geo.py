def cart2geo(cartesian_coordinates):
    import numpy as np

    X = cartesian_coordinates[:, 0]
    Y = cartesian_coordinates[:, 1]
    Z = cartesian_coordinates[:, 2]
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    lat = np.arcsin(Z/R) * 180 / np.pi
    lon = np.arctan2(Y, X) * 180 / np.pi

    return lat, lon
