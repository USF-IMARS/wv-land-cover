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
#               
#		users/lizcanosandoval/Seagrass/Sentinel/01_OriginalMosaics
#
# ./gee-uploads/gbucket_to_gee_w_metadata.sh \
# 		rookery-wv-classmaps \
#		/srv/imars-objects/rookery/Processed/wv_classMaps_rgb \       
#               users/tylarmurray/nerrs_rookery
#
# bash gee-uploads/gbucket_to_gee_w_metadata_jobos.sh \
#    jobos-wv-classmaps \
#    /srv/imars-objects/jobos/Processed/wv_ortho_xml \
#    users/tylarmurray/nerrs_jobos \
#    | tee jobos_upload-2022_10.log


# Modified from: https://www.tucson.ars.ag.gov/notebooks/uploading_data_2_gee.html

# hardcoded metadata
country="USA"
generator="Tylar_Murray+Digna_Rueda"

echo_if_test=""  # set this to "echo " to test the script, else set to ""

xml_reader_cmd="python3 ./wv_classify/read_wv_xml.py "
filepanther_cmd="python3 -m filepanther "

echo checking if the collection "$3" exists...
result=`${echo_if_test} earthengine create collection $3`
if `test -z "$result"`; then  # exit if creation failed
    echo collection created.
fi
echo $result

# In the following loop we get the entire path to all the geotifs using the specified Gcloud bucket. 
# Each file will have a format like this: `gs://my_gee_bucket/FILE_January2000.tif`.
# Each call to earthengine will launch a task that you can monitor in the JS Code editor "tasks" tab.
for geotiff in `gsutil ls gs://$1/*.tif`; do  
	#filename=${geotiff%.*}
	filename=${geotiff##*/} 
	asset_id="${filename%.*}" 
	echo ""
	echo "*** Transfering file " $asset_id "***"
	echo "*** parsing metadata..."
	# python3 filepanther -q parse /srv/imars-objects/rookery/Processed/wv_classMaps_rgb/20180501T160614_01_P003_WV02_ClassificMap_fullClass_Rookery.tif --pattern /srv/imars-objects/rookery/Processed/wv_classMaps_rgb/%Y%m%dT%H%M%S_{number}_P{pass_n}_WV{sat_n}_ClassificMap_fullClass_Rookery.tif > metadata.json
	$filepanther_cmd -q parse $filename \
	    --pattern %Y%m%dT%H%M%S_{number}_P{pass_n}_WV{sat_n}_ClassMap_v{algorithm_version}_Jobos_wDEM.tif \
	    --pickle_fpath metadata.pickle

	echo "*** estimating xml filename..."
	# to get the XML filename we need to do a few weird things:
	# * the xml filename contains 12 numbers that we don't know
	#    a * glob is used to capture these unknown digits (\d{12}).
	# * the filename is all upper-case, so %b is not an exact match.
	#    `tr` is used to convert the output to uppercase
	#
	# python3 -m filepanther -vvv format --pattern '%y%b%d%H%M%S-M1BS-504649660010_{number}_P{pass_n}.XML' --pickle_file metadata.pickle | tr '[:lower:]' '[:upper:]' | sed 's/\\D{12/\\d{12/' 
	xml_fileglob=`$filepanther_cmd -q format --pattern '%y%b%d%H%M%S-M1BS-*_{number}_P{pass_n}.XML' --pickle_file metadata.pickle | tr '[:lower:]' '[:upper:]'`
	echo "xml fname is like: ${xml_fileglob}"

	echo "*** searching for xml file..."
	xml_fpath=`find ${2} -name ${xml_fileglob}`
	if [ -z "${xml_fpath}" ]; then
	        echo find ${2} -name ${xml_fileglob}
		echo "xml file not found!"
		# append file to list of failed files & continue
		echo "missing_xml_file, $filename, find ${2} -name ${xml_fileglob}" >> missing_xml_files.log
		exit 1
	else
		echo "found file: ${xml_fpath}"
	fi

	echo "*** extracting properties from .xml..."
	xml_vars=`${xml_reader_cmd} ${xml_fpath}`
	echo "${xml_vars}"

	echo "*** formatting filename-extracted params for gee..."
	datetime=`$filepanther_cmd -q format --pattern '%Y-%m-%dT%H:%M:%S' --pickle_file metadata.pickle`
	classifier="NERRS-mangroves-decision-tree-v{algorithm_version}"
	echo "$datetime"

	echo "*** transferring image and metadata..."
	${echo_if_test} earthengine upload image gs://$1/$filename \
		-f --asset_id=$3/$asset_id \
		--nodata_value=0 \
		--crs="EPSG:4326" \
		--pyramiding_policy=mode \
		-ts=$datetime \
		${xml_vars} \
		-p country=${country} \
		-p generator=${generator} \
		-p classifier=${classifier}
	echo "done!"
	echo ""
done
