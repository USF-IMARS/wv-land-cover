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
    2. remaining python packages w/ setup.py `pip install -e .`

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
