#!/usr/bin/env bash
# Usage:
# ./transfer_seagrass.sh src_bucket dest_asset
# Example: 
# ./transfer_seagrass.sh seagrass_mosaics/original_mosaics users/lizcanosandoval/Seagrass/Sentinel/01_OriginalMosaics
# ./transfer_seagrass.sh seagrass_mosaics/edited_mosaics users/lizcanosandoval/Seagrass/Sentinel/02_EditedMosaics
# Modified from: https://www.tucson.ars.ag.gov/notebooks/uploading_data_2_gee.html#

##Some metadata
country="USA"
satellite="Sentinel-2"
generator="Lizcano-Sandoval"
classifier="SVM"

result=`earthengine create collection $2`
if `test -z "$result"`; then
    echo $result
    exit 1
fi
# In the following loop we get the entire path to all the geotifs using the specified 
# Gcloud bucket. Each file will have a format like this: gs://my_gee_bucket/FILE_January2000.tif
# Each call to earthengine it will launch a task that you can monitor in the JS Code editor 
# at the "tasks" tab.
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
    earthengine upload image gs://$1/$filename -f --asset_id=$2/$asset_id --nodata_value=0 --crs="EPSG:4326" -ts=$date -p="year=${year}" -p="name_code=${code}" -p="tile_id=${tile_id}" -p="country=${country}" -p="satellite=${satellite}" -p="generator=${generator}" -p="classifier=${classifier}"
done
