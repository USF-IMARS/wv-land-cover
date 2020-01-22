# wv2-processing
Processing scripts for decision-tree land use classification on WorldView-2 images.

## Software Dependencies
* gdal
* pygdal
* TODO: more here

## Installation
### PSC Bridges
```
$ git clone https://github.com/iceberg-project/Seals.git
# TODO: more here
```

-----------------------------------------------------------------------------------------------------------------

## General Installation
NOTE: if you are running this code on IMaRS's servers (eg userproc or seashell) jump directly to the IMaRS user quickstart document (./docs/imars-local.md).

1. download: `git clone https://github.com/USF-IMARS/wv2-processing.git`
2. install dependencies
    1. OS-level:
        * Ubuntu:
            * for gdal: `sudo apt istall -y python3-gdal`
    2. remaining python packages w/ setup.py `pip3 install -e .`
        * alternatively: `pip3 install -r requirements.txt` or manually install deps listed therein.

## Testing
### test data
Test data is stored internally at IMaRS and mounted at `/srv/imars-objects/homes/common/wv2-processing/test_data/`.
To run tests you should create a symlink from there to a dir named `test_data` in this project root `ln -s /srv/imars-objects/homes/common/wv2-processing/test_data/ test_data`.

Alternatively, you may download a version of these files from google drive [here](https://drive.google.com/file/d/1kWzAIxrhxD_ROwjMSZW1BTJxWGHtoGGd/view?usp=sharing) if you have been granted the appropriate permissions.
These files are restricted to IMaRS and collaborators; please do not share them in any form.
Once the file is downloaded you must extract this file to `wv2-processing/test_data/`.

TODO: add PGC's files here & merge directories.

### running tests
Python tests herein are generally orchestrated by pytest and live alongside the code they are testing with the suffix `_test`.

Note that comparing hashes on output files doesn't work well b/c of variations in the script and floating point errors so the tests are not very robust; they mostly just check things run without throwing exception.
For much of my testing I had to resort to opening the geotiffs with QGIS and confirming that they look right.

# Usage
## Overview & Manual Steps
Processing is broken into a few steps.
Below are examples of how each step might be run.
0. `INPUT_DIR`, `ORTHO_OUTPUT_DIR`, and other variables below must be set (eg `INPUT_DIR=/home/tylar/wv_proc/my_input_files`).
1. create resampled tifs using pgc_ortho:
    * `python ./pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no-pyramids $INPUT_DIR $ORTHO_OUTPUT_DIR`
2. run the wv_classify script on the resampled tifs
    1. python `python ./wv_classify.py $ORTH_FILE $ID $MET $CRD $DT $SGW $FILT $STAT $LOC $ID_N $RRS_OUT $CLASS_OUT`

3. use gdal or similar tools to mosaic multiple outputs together
    * see [this gist](https://gist.github.com/7yl4r/d03f9617212db5efded1f8a0d34550d3)

## Script Parameter Reference

```
-p                 = projection (4326 is the EPSG code for WGS geographic projection)
-c                 = stretch (ns means "no stretch")
-t                 = output bit depth
-f                 = file format
--no-pyramids      =  prevents the code from creating pyramids in the output GeoTIFF
$INPUT_DIR         = directory for NITF staging/input
$ORTHO_OUTPUT_DIR' = directory for output of pgc_ortho code as GeoTIFF
'$ORTH_FILE        = [this variable is outdated and should be deleted]
$ID                = image file name
$MET               = metadata file
$CRD               = coordinate system (e.g. EPSG 4326)
$DT                = input variable for whether to run the decision tree (DT = 2) or to just run Rrs conversion (DT = 0)
$SGW               = [this variable is outdated and should be deleted]
$FILT              = input variable indicating the size of the moving window filter (1 = 3x2, 2 = 5x5, etc.)
$STAT              = [this variable is outdated and should be deleted]
$LOC               = string identifier (e.g. "NSF_Texas")
$ID_N              = identifier based on the file number being run
$RRS_OUT           = directory for output of Rrs GeoTIFFs
$CLASS_OUT         = directory for output of mapped GeoTIFFs
```

More usage details in the `./docs` directory:
* IMaRS local systems use: docs/imars-local.md
* SLURM : docs/slurm.md
* Airflow: docs/airflow.md
