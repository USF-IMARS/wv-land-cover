# wv2-processing
Processing scripts for decision-tree land use classification on WorldView images.
This project funded by NSF South Big Data Hub and then by the RB & JB NERRs.

## Habitat Cover Classes
### 3D wetlands
Below are the habitat classes created by the 3D wetlands version of this code. 
`wv_classification_colormap.txt` Provides a colormap with similar values and associated colors for use in QGIS (and others).

```
BA = bare soil
WA = water
DG = degraded mangrove
MA = marsh
SC = scrub/shrub
FW = forested wetland (this is mangrove forest in southwest Florida)
FU = forested upland
UG = upland grass
dev = developed
```

## NERRS+IMaRS MCC
For the NERRS Mangrove Coast Collective project and related publications see the [MCC mapping class detials gsheet](https://docs.google.com/spreadsheets/d/1ay7N4hZMNwbxTpRnwHpxpMUNkTmbJZGTEDaQGiepDiU/edit?usp=sharing).

---------------------------------------------------------------------------------------------------------------

## Software Dependencies
* gdal
* pygdal

## Installation
### basic installation
This will install all needed scripts and the easy-to-install dependencies.

```
git clone https://github.com/USF-IMARS/wv-land-cover.git
cd wv-land-cover
git submodule update --init --recursive --remote
```

### detailed dependencies setup
If you are getting errors after performing the basic installation, then your system may need more advanced configuration.
For detailed dependency setup you will need to work with your system administrator.

#### MATLAB setup
Installation instructions for MATLAB are elsewhere.
No special configuration is needed

#### setup for PGC/imagery_utils on linux
```
# gdal
sudo apt install libgdal-dev
sudo apt install gdal-bin
sudo apt install -y python3-gdal


# remaining python packages w/ setup.py `pip3 install -e .`
pip3 install -r requirements.txt

# proj
sudo conda install -y -c conda-forge proj
sudo conda install -y -c conda-forge proj-data
```

NOTE: python bindings for gdal need to be setup manually. See the relevant section in requirements.txt.

#### SLURM setup
SLURM dependency setup is managed via `module add [...]` commands. These will be included in the job submission scripts.

### PSC Bridges
```
$ git clone https://github.com/iceberg-project/Seals.git
# TODO: more here
```

-----------------------------------------------------------------------------------------------------------------

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
Below is an overview of how each step might be run manually.
For specific, detailed examples of running the code see the `./docs/examples/` folder.
0. `INPUT_DIR`, `ORTHO_OUTPUT_DIR`, and other variables below must be set (eg `INPUT_DIR=/home/tylar/wv_proc/my_input_files`).
1. create resampled tifs using pgc_ortho:
    * `python ./pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no-pyramids $INPUT_DIR $ORTHO_OUTPUT_DIR`
2. run the wv_classify script on the resampled tifs
    1. the (now outdated) python version: `python ./wv_classify.py $ORTH_FILE $ID $MET $CRD $DT $SGW $FILT $STAT $LOC $ID_N $RRS_OUT $CLASS_OUT`
    2. the MATLAB version: `matlab -nodisplay -nodesktop -r "wv_classify('$image2','$input_img_basename','$met','$crd_sys','$dt','$filt','$loc','$SLURM_ARRAY_TASK_ID','$rrs_out','$class_out')"`
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
