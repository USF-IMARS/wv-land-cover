# wv2-processing
Processing scripts for decision-tree land use classification on WorldView-2 images.

The submit_py.sh file is what I use in Circe to call the pgc_ortho.py script, which has a number of sub-scripts called.
The submit_py.sh also contains the Matlab script call, so you'll want to comment out those lines before testing it.

## Software Dependencies
* pygdal
* TODO: more here

## Installation
### PSC Bridges
```
$ git clone https://github.com/iceberg-project/Seals.git
# TODO: more here
```

-----------------------------------------------------------------------------------------------------------------

# Installation
NOTE: if you are running this code on IMaRS's servers (eg userproc or seashell) jump directly to the [IMaRS user quickstart section below](README.md#IMaRS-User-Quickstart).

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

# github basics

## download repo to local machine
`git clone https://github.com/USF-IMARS/wv2-processing`

## basic git/github workflow
1. `git pull origin master` - this updates your local to match the remote
2. make your file edits
3. `git status` to review the changes you have made
4. (optional) `git diff` to review even more closely
5. `git add my-new-file.py` to add new files to the "staging area"
6. `git commit -a -m "my new commit"` submits a commit with all changes and your commit message "my new commit"
7. `git push origin master` this uploads your commits to github


# Usage
## Overview & Manual Steps
Processing is broken into a few steps.
Below are examples of how each step might be run.
1. create resampled tifs using pgc_ortho:
    * `python ./pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no-pyramids $INPUT_DIR $ORTHO_OUTPUT_DIR`
2. run one of the wv_classify scripts on the resampled tifs
    1. python `python ./wv_classify.py $ORTH_FILE $ID $MET $CRD $DT $SGW $FILT $STAT $LOC $ID_N $RRS_OUT $CLASS_OUT`
    2. matlab:
        ```
        matlab -nodisplay -nodesktop -r "\
            cd('/opt/wv2_processing');\
            wv2_processing(\
                '$ORTH_FILE',\
                '{{params.id}}',\
                '$MET',\
                '{{params.crd_sys}}',\
                '{{params.dt}}',\
                '{{params.sgw}}',\
                '{{params.filt}}',\
                '{{params.stat}}',\
                '{{params.loc}}',\
                '{{params.id_number}}',\
                '$RRS_OUT',\
                '$CLASS_OUT'\
            );\
            exit\
        "
        ```
3. use gdal or similar tools to mosaic multiple outputs together

## Script Parameter Ref

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

## SLURM
These processing tasks have been executed on USF Research Computing's research cluster, [CIRCE](https://wiki.rc.usf.edu/index.php/CIRCE) using the SLURM task scheduler.
The slurm submission bash script is in the root of this repo at [./submit_py.sh](https://github.com/USF-IMARS/wv2-processing/blob/master/submit_py.sh)

## Apache Airflow
This task has been run as an airflow DAG on IMaRS's airflow cluster.
The airflow dag definition file can be viewed at [USF-IMARS/imars_dags//dags/processing/wv2_classification/wv_classification.py](https://github.com/USF-IMARS/imars_dags/blob/master/dags/processing/wv2_classification/wv_classification.py).
Accompanying wrapper scripts are in the [scripts subdir of the same location](https://github.com/USF-IMARS/imars_dags/tree/master/dags/processing/wv2_classification/scripts).
Of particular note is [USF-IMARS/imars_dags//dags/processing/wv2_classification/scripts/ntf_to_rrs.sh](https://github.com/USF-IMARS/imars_dags/blob/master/dags/processing/wv2_classification/scripts/ntf_to_rrs.sh).

USF-IMaRS/imars_dags expects a certain configuration & software suite is expected to exist on each node.
This airflow cluster's software and configuration management is managed via puppet ([IMaRS-private puppet repo ln](https://github.com/usf-imars/imars_puppet)); related documentation can be found there and can be provided on request.
One of the more important dependencies is [imars-etl](https://github.com/USF-IMARS/imars-etl), which wraps IMaRS's underlying object & metadata storage systems.


# IMaRS User Quickstart
This section is for IMaRS researchers running this code on IMaRS's servers (eg userproc or seashell).
Setup and usage is simplified because a lot has already been set up for you.


## Setup
1. ssh to one of the processing servers
2. make sure you are in your home directory and navigate to wherever you want the code (useful commands: `pwd`, `cd`, `ls`, `mkdir`)
3. download a local copy of the code: `git clone https://github.com/USF-IMARS/wv2-processing`
4. create a link to the test data: `ln -s /srv/imars-objects/homes/common/wv2-processing/test_data/ test_data`
5. run the code tests: `python36 -m pytest`
    * all should pass; open an issue in this repo if they do not

Now that you are set up you can start working with real files.

## Extracting Files to Work With
The best way to get the files you want to work with is to use the imars-etl tool.
This tool allows you to copy a file matching a metadata selection to your current working directory.
Common metadata you might select for: the product_id, the area_id, and the date_time.

This processing currently uses product id 11 and 14 (.ntf & .xml files, respectively), so these will be the product_id values you want.
The area_id may vary; in the example below area_id=9 is used, which corresponds to the FCMaP "monroe" region of Florida, which includes the Florida Keys.

For a full list of `area_id` and `product_id` values and descriptions, see [imars_puppet/.../product_metadata_rows.sql](https://github.com/USF-IMARS/imars_puppet/blob/production/modules/role/files/sql/product_metadata_rows.sql).
For more detailed information and more example SQL queries use [IMaRS's Blazer server](http://imars-physalis.marine.usf.edu:3000/).

Example to download xml & ntf files:
```
# let's start by creating an empty input_data directory and moving into it
mkdir input_data
cd input_data

# find files which match an SQL query for a specific product, area, and time range:
imars-etl select 'WHERE product_id=14 AND area_id=9 AND date_time LIKE "2017-01-03%" ORDER BY date_time'

# same as above but print out only date_time and provenance columns
imars-etl select -c date_time,provenance 'WHERE product_id=11 AND area_id=9 AND date_time LIKE "2017-01-03%" ORDER BY date_time'

# copy a specific file to my current directory using the date_time, area, and product
imars-etl extract 'date_time="2017-01-03T15:57:53.549250" AND area_id=9 AND product_id=11'

# also copy the xml file (product_id 14) for this granule
imars-etl extract 'date_time="2017-01-03T15:57:53.549250" AND area_id=9 AND product_id=14'

# once we are done extracting files we can move up one directory back to the root of this project
cd ..
```

You will now have the files `WV02_20170103155753_0000000000000000_17Jan03155753-M1BS-058526494010_01_P005.ntf` and `WV02_20170103155753_0000000000000000_17Jan03155753-M1BS-058526494010_01_P005.xml` in your current working directory and can use them.

See more docs on how to use imars-etl in [USF-IMARS/imars-etl](https://github.com/USF-IMARS/imars-etl).


## Running the Code on Your Files
Assuming we are using the `WV02_20170103155753_0000000000000000_17Jan03155753-M1BS-058526494010_01_P005` files extracted from the previous section we can now run the python scripts on these input files.

```
# Starting with our working directory in the root of the project.
# Our input files are in `./input_data` and we should see them if we do `ls ./input_data`
# Let's also create some directories for our output data
mkdir ortho_data
mkdir output_data

# 1st step is to run pgy ortho
python pgc_ortho.py \
    -p 4326 \
    -c ns \
    -t UInt16 \
    -f GTiff \
    ./input_data \
    ./ortho_data

# if this was successful we should now have files in ortho_data
ls -lh ./ortho_data

# We can now run the python classifier script on the pgc_ortho output

python36 -m wv_classify.wv_classify \
    ./ortho_data/WV02_20170103155753_0000000000000000_17Jan03155753-M1BS-058526494010_01_P005_u16ns4326.tif \
    ./input_data/WV02_20170103155753_0000000000000000_17Jan03155753-M1BS-058526494010_01_P005.xml \
    ./output_data MONROE "EPSG:4326" 2 1

# if this was successful we should now have rrs, Rrs, and classification_map files in ./output_data/
```
