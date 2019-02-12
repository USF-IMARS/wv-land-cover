# matlab functions ported to python

import math
import numpy

try:
    from osgeo import gdal
except ImportError:
    import gdal


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
        3D array of all raster bands. Usage A[band, row, col]
    R : gdal data object
        NOTE: not a SpatialReference like matlab uses, but
            is passed to geotiffwrite in the same way.
    """
    print("reading geotiff '{}'".format(filename))
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
    n_bands, n_rows, n_cols = data_grid.shape
    # usage check...
    assert len(data_grid[:, 0, 0]) == n_bands
    assert len(data_grid[0, :, 0]) == n_rows
    assert len(data_grid[0, 0, :]) == n_cols
    print("read {} bands at resolution {}x{}".format(n_bands, n_rows, n_cols))

    # TODO: use numpy.swapaxes to make indices match geotiffwrite?

    return data_grid, ds


def geotiffwrite(
    outFileName, arr_out, ds, CoordRefSysCode=4326,
    row_index=0,
    col_index=1,
    band_index=2,
):
    """
    https://www.mathworks.com/help/map/ref/geotiffwrite.html

    parameters:
    ----------
    arr_out : 3d numpy.array
        values to write to geotiff
    coor_sys :
        code for projection. eg 4326
    ds :
        gdal data object used only to get the GeoTransform & projection info
    """
    DTYPE_MAP = {  # mappings of array cell data types to gdal data types
        numpy.float32: gdal.GDT_Float32,
        numpy.float64: gdal.GDT_Float64,
        numpy.uint16: gdal.GDT_UInt16,
    }
    driver = gdal.GetDriverByName("GTiff")
    n_bands = arr_out.shape[band_index]
    n_rows = arr_out.shape[row_index]
    n_cols = arr_out.shape[col_index]
    print("writing {}x{}, {}-band geotiff to '{}'".format(
        n_rows, n_cols, n_bands, outFileName
    ))
    # === set up datatype for output file:
    cell_dtype = type(arr_out[0, 0, 0])
    if cell_dtype not in DTYPE_MAP.keys():
        raise ValueError(
            "Unable to map array of type {} to gdal type.".format(cell_dtype) +
            " Available mappings are: \n{}".format(DTYPE_MAP)
        )
    # implied else
    outdata = driver.Create(
        outFileName, n_cols, n_rows, n_bands, DTYPE_MAP[cell_dtype]
    )
    assert outdata is not None
    
    # set same geotransform as input
    outdata.SetGeoTransform(ds.GetGeoTransform())

    # TODO: set projection using CoordRefSysCode instead of:
    # same projection as input
    outdata.SetProjection(ds.GetProjection())
    for band in range(n_bands):
        print('\twriting band #{}...'.format(band+1))
        band_arr = arr_out[band, :, :]
        outdata.GetRasterBand(band+1).WriteArray(band_arr)

        # if you want these values transparent
        outdata.GetRasterBand(band+1).SetNoDataValue(numpy.nan)

        # === required dereference?
        # https://trac.osgeo.org/gdal/wiki/PythonGotchas#Savingandclosingdatasetsdatasources
        band_arr = None
        # outdata.FlushCache()  # saves to disk
    outdata.FlushCache()  # saves to disk
    # === required dereference?
    # https://trac.osgeo.org/gdal/wiki/PythonGotchas#Savingandclosingdatasetsdatasources
    outdata = None
