import numpy as np


def cart2geo(cartesian_coordinates):
    # Convert cartesian coordinates to lat, lon for every element or node
    # Define cartesian coordinate vectors
    X = cartesian_coordinates[:, 0]
    Y = cartesian_coordinates[:, 1]
    Z = cartesian_coordinates[:, 2]

    # Conversion to geographical coordinates
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    lat = np.arcsin(Z/R) * 180 / np.pi
    lon = np.arctan2(Y, X) * 180 / np.pi

    return lat, lon
