import numpy as np


def geo2cart(lon, lat, r, pandas_length):

    if len(lon) != pandas_length:
        original_rows = len(lon)
        original_columns = len(lon[0])
        lon = lon.reshape(-1, 1)
        lat = lat.reshape(-1, 1)
        r = r.reshape(-1, 1)
    x = r * np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon))
    y = r * np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon))
    z = r * np.sin(np.deg2rad(lat))

    if len(lon) != pandas_length:
        x = x.reshape(original_rows, original_columns)
        y = y.reshape(original_rows, original_columns)
        z = z.reshape(original_rows, original_columns)

    return x, y, z
