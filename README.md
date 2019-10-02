# wv2-processing
Processing scripts for decision-tree land use classification on WorldView-2 images.

The submit_py.sh file is what I use in Circe to call the pgc_ortho.py script, which has a number of sub-scripts called.
The submit_py.sh also contains the Matlab script call, so you'll want to comment out those lines before testing it.

# Installation
1. download: `git clone git@github.com:USF-IMARS/wv2-processing.git`
2. install dependencies
    1. OS-level:
        * Ubuntu:
            * for gdal: `sudo apt istall -y python3-gdal`
    2. remaining python packages w/ setup.py `pip3 install -e .`
        * alternatively: `pip3 install -r requirements.txt` or manually install deps listed therein.

## github basics

#### download repo to local machine
`git clone https://github.com/USF-IMARS/wv2-processing`

#### basic git/github workflow
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
