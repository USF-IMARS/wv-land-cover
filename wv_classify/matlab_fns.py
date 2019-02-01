# matlab functions ported to python

import math
import numpy


def d2r(deg):
    return deg * math.pi / 180.0


def tand(x):
    return math.tan(d2r(x))


def sind(x):
    return math.sin(d2r(x))


def asind(x):
    return math.asin(d2r(x))


def cosd(x):
    return math.cos(d2r(x))


def acosd(x):
    return math.acos(d2r(x))


def mldivide(X, y):
    """Linear regression X\y AKA mldivide(X,y) AKA `inv(X'*X)*X'*y;`
    as defined by MATLAB docs at:
        https://www.mathworks.com/help/matlab/ref/mldivide.html

    inv(X.T @ X) @ X.T @ y

    NOTE: this is not _strictly_ identical to matlab's mldivide for under-
        determined systems. See the following S.O q/a for more info:
        https://stackoverflow.com/a/38228156/1483986

    """
    return numpy.linalg.lstsq(X, y)


def rdivide(A, B):
    # https://www.mathworks.com/help/matlab/ref/rdivide.html
    # NOTE: maybe A & B should be cast to numpy arrays?
    return A / B


def geotiffread(X):
    # https://www.mathworks.com/help/map/ref/geotiffread.html
    # TODO: read geotiff w/ gdal (or other)
    pass
    # return data_grid, spatial_ref


def geotiffwrite():
    # TODO
    # https://www.mathworks.com/help/map/ref/geotiffwrite.html
    pass
