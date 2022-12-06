⚠️ This file is *very* out of date. Please defer to the README.md.

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
