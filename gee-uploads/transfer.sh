# The latest version of gbucket_to_gee script that uploads without xml files.
# Developed for https://github.com/USF-IMARS/wv-land-cover/issues/51
#
# usage:
#  transfer.sh wv-rookery-classifications Rookery projects/nerrs-mcc/assets/nerrs_rookery_v04

# hardcoded metadata
country="USA"
generator="Tylar_Murray+Digna_Rueda"
classifier="NERRS-mangroves-decision-tree"

echo_if_test=""  # set this to "echo " to test the script, else set to ""

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
	# ...*ClassificMap_fullClass_ Rookery-wDEM_v3_DEM.tif
        # wv 2025-03 : 20210109T145803_01_P001_WV03_ClassMap_vB1_Jobos_wDEM.tif
        $filepanther_cmd -q parse $filename \
          --pattern "%Y%m%dT%H%M%S_{number}_P{pass_n}_WV{sat_n}_ClassMap_{adjustments_version}_${2}_wDEM_adjustedLatLon.tif" \
          --pickle_fpath metadata.pickle
	    
#	echo "*** estimating xml filename..."
#	# to get the XML filename we need to do a few weird things:
#	# * the xml filename contains 12 numbers that we don't know
#	#    a * glob is used to capture these unknown digits (\d{12}).
#	# * the filename is all upper-case, so %b is not an exact match.
#	#    `tr` is used to convert the output to uppercase
#	#
	# python3 -m filepanther -vvv format --pattern '%y%b%d%H%M%S-M1BS-504649660010_{number}_P{pass_n}.XML' --pickle_file metadata.pickle | tr '[:lower:]' '[:upper:]' | sed 's/\\D{12/\\d{12/' 
#	xml_fileglob=`$filepanther_cmd -q format --pattern '%y%b%d%H%M%S-M1BS-*_{number}_P{pass_n}.XML' --pickle_file metadata.pickle | tr '[:lower:]' '[:upper:]'`
#	echo "xml fname is like: ${xml_fileglob}"
#
#	echo "*** searching for xml file..."
#	xml_fpath=`find ${2} -name ${xml_fileglob}`
#	if [ -z "${xml_fpath}" ]; then
#		echo "xml file not found!"
#		# append file to list of failed files & continue
#		echo "missing_xml_file, $filename, find ${2} -name ${xml_fileglob}" >> missing_xml_files.log
#		exit 1
#	else
#		echo "found file: ${xml_fpath}"
#	fi
#
#	echo "*** extracting properties from .xml..."
#	xml_vars=`${xml_reader_cmd} ${xml_fpath}`
	echo "${xml_vars}"
        # NOTE: other vars extracted from the filename patter will also be included in xml_vars

	echo "*** formatting ts for gee..."
	datetime=`$filepanther_cmd -q format --pattern '%Y-%m-%dT%H:%M:%S' --pickle_file metadata.pickle`
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
