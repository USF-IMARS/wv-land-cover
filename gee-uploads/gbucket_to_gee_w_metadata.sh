#!/usr/bin/env bash
# Usage:
# ./gbucket_to_gee_w_metadata.sh src_bucket xml_filespath dest_asset
#
# Moves .tif files from a GCloud bucket into GEarthEngine including relevant metadata.
# Some metadata is hard-coded in the file below and some metadata is pulled from the 
# `.xml` file that corresponds to each `.tif`.
# The `.xml` file must be on the local machine and have the same filename as the
# `.tif` file in the GCloud bucket.
#
# example GCloud filename:
# 	20200929T162717_03_P008_WV03_ClassificMap_fullClass_Rookery.tif
# example .NTF file from which the .tif was generated:
#   20SEP29162717-M1BS-504649660010_03_P008.NTF
# corresponding .xml filename:
#   20SEP29162717-M1BS-504649660010_03_P008.XML
#
# Examples: 
# ./gbucket_to_gee_w_metadata.sh \
#		seagrass_mosaics/original_mosaics \
# 		./seagrass_mosiacs/xml_files/ \
#		users/lizcanosandoval/Seagrass/Sentinel/01_OriginalMosaics
#
# ./gee-upload/gbucket_to_gee_w_metadata.sh \
# 		rookery-wv-classmaps \
#		/srv/imars-objects/rookery/Processed/wv_classMaps_rgb \       
#               users/tylarmurray/nerrs/rookery
#
# Modified from: https://www.tucson.ars.ag.gov/notebooks/uploading_data_2_gee.html

# hardcoded metadata
country="USA"
satellite="WV0(2|3)"
generator="Tylar Murray & Digna Rueda"
classifier="NERRS-mangroves-decision-tree"
echo_if_test="echo "  # set this to "echo " to test the script, else set to ""
xml_reader_cmd="python3 ./wv_classify/read_wv_xml.py "
filepanther_cmd="python3 -m filepanther "

echo creating the collection "$3"...
result=`${echo_if_test} earthengine create collection $3`
if `test -z "$result"`; then  # exit if creation failed
    echo $result
    exit 1
fi

# In the following loop we get the entire path to all the geotifs using the specified Gcloud bucket. 
# Each file will have a format like this: `gs://my_gee_bucket/FILE_January2000.tif`.
# Each call to earthengine will launch a task that you can monitor in the JS Code editor "tasks" tab.
for geotiff in `gsutil ls gs://$1/*.tif`; do  
	#filename=${geotiff%.*}
	filename=${geotiff##*/} 
	asset_id="${filename%.*}" 
	echo "*** Transfering file " $asset_id "***"

	year=${asset_id:6:4}
	date="${year}-01-01T12:00:00"
	tile=${asset_id:0:6}
	tile_id="${tile%_*}" 
	code=${asset_id:11:4}
	# python3 filepanther -q parse /srv/imars-objects/rookery/Processed/wv_classMaps_rgb/20180501T160614_01_P003_WV02_ClassificMap_fullClass_Rookery.tif --pattern /srv/imars-objects/rookery/Processed/wv_classMaps_rgb/%Y%m%dT%H%M%S_{number}_P{pass_n}_WV{sat_n}_ClassificMap_fullClass_Rookery.tif > metadata.json
	${echo_if_test} $filepanther_cmd -q parse $filename --pattern %Y%m%dT%H%M%S_{number}_P{pass_n}_WV{sat_n}_ClassificMap_fullClass_Rookery.tif > filepath_metadata.json
	xml_filename=`${echo_if_test} $filepanther_cmd format --type=wv_xml --json_file=filepath_metadata.json`
	xml_vars=`${xml_reader_cmd} $2/${filename}.xml`
	echo "xml_vars : ${xml_vars}"

	${echo_if_test} earthengine upload image gs://$1/$filename \
		-f --asset_id=$3/$asset_id \
		--nodata_value=0 --crs="EPSG:4326" -ts=$date \
		-p="year=${year}" \
		-p="name_code=${code}" \
		-p="tile_id=${tile_id}" \
		-p="country=${country}" \
		-p="satellite=${satellite}" \
		-p="generator=${generator}" \
		-p="classifier=${classifier}"
done
