#!/bin/bash
#SBATCH --partition=circe
#SBATCH --job-name ="wv_mosaic_gdal"
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=60480
#SBATCH --time=40:00:00
#SBATCH --array=0

module purge
module load apps/python/3.7.3
module load apps/gdal/3.0.1_el7_gcc

#gdalbuildvrt -a_srs EPSG:4326 my_overview_file.vrt $WORK/output/Rrs/NSF_SWTX/.*tif


## TODO: load & edit my_overview_file.vrt xml:
##       ref: https://gist.github.com/7yl4r/d03f9617212db5efded1f8a0d34550d3
## =================================================================
##       replace extant <VRTRasterBand> opening element tag with:
#<VRTRasterBand dataType="Byte" band="1" subClass="VRTDerivedRasterBand">
##       and add the following elements within that block:
#    <PixelFunctionType>maximu</PixelFunctionType>
#    <PixelFunctionLanguage>Python</PixelFunctionLanguage>
#    <PixelFunctionCode><![CDATA[
#    import numpy as np
#    def maximu(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize,
#                 raster_ysize, buf_radius, gt, **kwargs):
#         np.round_(np.clip(np.max(in_ar, axis = 0),0,255),
#            out = out_ar)
#    ]]></PixelFunctionCode>
    
export GDAL_VRT_ENABLE_PYTHON=YES
export BIGTIFF=YES

gdaladdo --config BIGTIFF_OVERIEW YES -ro my_overview_file.vrt 1
