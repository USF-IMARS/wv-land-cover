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
    return math.degrees(math.acos(x))


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


def geotiffread(filename, numpy_dtype=None):
    """
    Reads geotiff w/ gdal.
    https://www.mathworks.com/help/map/ref/geotiffread.html

    returns:
    --------
    A : array
        3D array of all raster bands. Usage A[band, row, col]
    spatial_ref : gdal data object's SpatialReference like matlab uses.
        actually just the output of `ds.GetGeoTransform()` and
        ` ds.GetProjection()` in an array.
    """
    print("reading geotiff '{}'".format(filename))
    ds = gdal.Open(filename)
    data_grid = numpy.array([
        ds.GetRasterBand(band+1).ReadAsArray()
        for band in range(ds.RasterCount)
    ], dtype=numpy_dtype)
    # # if srcband is None:
    # #     continue
    # # stats = srcband.GetStatistics(True, True)
    # # if stats is None:
    # #     continue
    # data_grid.append(srcband.ReadAsArray())

    # prj = ds.GetProjection()
    # spatial_ref = osr.SpatialReference(wkt=prj)
    # spatial_ref = ds.SpatialReference()

    # swap indices around to match MATLAB geotiff[read/write]
    data_grid = numpy.moveaxis(data_grid, 0, 2)  # bands to index 2

    n_rows, n_cols, n_bands = data_grid.shape
    # usage check...
    assert len(data_grid[:, 0, 0]) == n_rows
    assert len(data_grid[0, :, 0]) == n_cols
    assert len(data_grid[0, 0, :]) == n_bands
    print("read {} bands at resolution {}x{}".format(n_bands, n_rows, n_cols))

    spatial_ref = [ds.GetGeoTransform(), ds.GetProjection()]
    del ds  # close dataset
    return data_grid, spatial_ref


def geotiffwrite(
    outFileName, arr_out, spatial_ref, CoordRefSysCode=4326,
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
    spatial_ref :
        gdal data object used only to get the GeoTransform & projection info
    """
    DTYPE_MAP = {  # mappings of array cell data types to gdal data types
        # numpy.uint8: gdal.GDT_UInt8,  # GDT_UInt8 doesn't exist :(
        numpy.uint16: gdal.GDT_UInt16,
        numpy.float32: gdal.GDT_Float32,
        numpy.float64: gdal.GDT_Float64,  # this might be too big
    }
    driver = gdal.GetDriverByName("GTiff")
    if len(arr_out.shape) == 2:
        # single-band geotiff
        n_bands = 1
        cell_dtype = type(arr_out[0, 0])
    else:
        # must have 3D output array
        assert len(arr_out.shape) == 3
        n_bands = arr_out.shape[band_index]
        cell_dtype = type(arr_out[0, 0, 0])

    n_rows = arr_out.shape[row_index]
    n_cols = arr_out.shape[col_index]
    # === set up datatype for output file:
    if cell_dtype not in DTYPE_MAP.keys():
        raise ValueError(
            "Unable to map array of type {} to gdal type.".format(cell_dtype) +
            " Available mappings are: \n{}".format(DTYPE_MAP)
        )
    print("writing {}x{} '{}', {}-band geotiff to '{}'".format(
        n_rows, n_cols, cell_dtype, n_bands, outFileName
    ))

    # === downsize if output is too big (is it?)
    # this is a workaround for:
    # ERROR 6: A 8810 pixels x 7516 lines x 8 bands Float64 image would be
    # larger than 4GB but this is the largest size a TIFF can be, and BigTIFF
    # is unavailable. Creation failed.
    if cell_dtype == numpy.float64:
        print("WARNING: casting float64 to float32 to reduce image size")
        arr_out = numpy.float32(arr_out)
        cell_dtype = numpy.float32

    outdata = driver.Create(
        # utf8_path, xsize,  ysize,  bands=1, eType=GDT_Byte, char options=None
        outFileName, n_cols, n_rows, n_bands, DTYPE_MAP[cell_dtype]
    )

    if outdata is None:
        raise ValueError("gdal driver failed!")

    # set same geotransform as input
    outdata.SetGeoTransform(spatial_ref[0])

    # TODO: set projection using CoordRefSysCode instead of:
    # same projection as input
    outdata.SetProjection(spatial_ref[1])
    # slicer is weirdness to slice out arbitrary band_index.
    #   This is just arr_out[:,:,band], but with "band" in the position
    #   specified by band_index.
    #   ref: https://stackoverflow.com/a/24434707/1483986
    slicer = [slice(None)] * 3
    for band in range(n_bands):
        if n_bands > 1:
            slicer[band_index] = band
            band_arr = arr_out[tuple(slicer)]
        elif len(arr_out.shape) == 2:
            # can't have multiple bands without 3D output array
            assert n_bands == 1
            band_arr = arr_out
        else:
            raise AssertionError("< 1 bands?")
        print('\twriting band #{}...'.format(band+1))
        outdata.GetRasterBand(band+1).WriteArray(band_arr)

        # if you want these values transparent
        outdata.GetRasterBand(band+1).SetNoDataValue(numpy.nan)

        # === required dereference?
        # https://trac.osgeo.org/gdal/wiki/PythonGotchas#Savingandclosingdatasetsdatasources
        del band_arr
        outdata.FlushCache()  # saves to disk
    outdata.FlushCache()  # saves to disk
    # === required dereference?
    # https://trac.osgeo.org/gdal/wiki/PythonGotchas#Savingandclosingdatasetsdatasources
    del outdata
    print(outFileName + " written.")
