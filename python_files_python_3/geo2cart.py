import numpy as np


def geo2cart(lon, lat, r):

    x = r * np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon))
    y = r * np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon))
    z = r * np.sin(np.deg2rad(lat))

    return x, y, z
