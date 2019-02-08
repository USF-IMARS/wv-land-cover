# matlab functions ported to python

import math
import numpy

# try:
from osgeo import gdal
# except: import gdal


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


def geotiffread(filename):
    """
    Reads geotiff w/ gdal.
    https://www.mathworks.com/help/map/ref/geotiffread.html

    returns:
    --------
    A : array
        3D array of all raster bands
    R : gdal data object
        NOTE: not a SpatialReference like matlab uses, but
            is passed to geotiffwrite in the same way.
    """
    ds = gdal.Open(filename)
    data_grid = numpy.array([
        ds.GetRasterBand(band+1).ReadAsArray()
        for band in range(ds.RasterCount)
    ])
    # # if srcband is None:
    # #     continue
    # # stats = srcband.GetStatistics(True, True)
    # # if stats is None:
    # #     continue
    # data_grid.append(srcband.ReadAsArray())

    # prj = ds.GetProjection()
    # spatial_ref = osr.SpatialReference(wkt=prj)
    # spatial_ref = ds.SpatialReference()

    return data_grid, ds


def geotiffwrite(outFileName, arr_out, ds, CoordRefSysCode):
    """
    https://www.mathworks.com/help/map/ref/geotiffwrite.html

    parameters:
    ----------
    coor_sys :
        code for projection. eg 4326
    ds :
        gdal data object used only to get the GeoTransform & projection info
    """
    driver = gdal.GetDriverByName("GTiff")
    cols, rows = arr_out.shape
    outdata = driver.Create(outFileName, rows, cols, 1, gdal.GDT_UInt16)
    # set same geotransform as input
    outdata.SetGeoTransform(ds.GetGeoTransform())

    # TODO: set projection using CoordRefSysCode instead of:
    # same projection as input
    outdata.SetProjection(ds.GetProjection())
    outdata.GetRasterBand(1).WriteArray(arr_out)

    # if you want these values transparent
    outdata.GetRasterBand(1).SetNoDataValue(numpy.nan)
    outdata.FlushCache()  # saves to disk!!
