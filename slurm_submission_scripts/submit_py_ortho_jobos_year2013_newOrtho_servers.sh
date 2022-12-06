#!/bin/bash
# bash file to run wv processing on a single IMaRS server

IMARS_HOMES_DIR=/srv/imars-objects/homes/
IMAGERY_UTILS_PATH=${IMARS_HOMES_DIR}/scratch/WV_ortho/wv-land-cover/pgc_duplication/
IMAGES_ROOT_PATH=${IMARS_HOMES_DIR}/scratch/digna/WV_varios/img/Jobos_2013b
LOGS_OUTPUT_DIR=${IMARS_HOMES_DIR}/scratch/WV_ortho/wv-land-cover/slurm_submission_scripts/SlurmOutput/
ortho_out=${IMARS_HOMES_DIR}/scratch/digna/WV_varios/Ortho/Jobos_2013b/
rrs_out=${IMARS_HOMES_DIR}/scratch/digna/WV_varios/Rrs/Jobos_2013b/

# Python code to check processing time:
#    starttime = datetime.today()
#    LogMsg('Image: %s' %(info.srcfn))

## Setup input arguments & file locations
## images1=`ls /work/d/druedaro/img/Jobos_perYear/Jobos_2013b/*.[nN][tT][fF]`
## met=`ls /work/d/druedaro/img/Jobos_perYear/Jobos_2013b/*.[xX][mM][lL]`
## ortho_out=/work/d/druedaro/output/Ortho/Jobos_perYear/Jobos_2013b/
## rrs_out=/work/d/druedaro/output/Rrs/Jobos_perYear/Jobos_2013b/

images1=`ls ${IMAGES_ROOT_PATH}/*.[nN][tT][fF]`
met=`ls IMAGES_ROOT_PATH/*.[xX][mM][lL]`


# Setup Matlab arguments
crd_sys=EPSG:4326
areaName='Jobos'

## Run Python code
images1a=($images1)

# for each image
for image_path in ${images1}; do
	base_name=$(basename ${image_path})
	echo "orthorectifying $base_name to $ortho_out..."
	python ${IMAGERY_UTILS_PATH}/pgc_ortho.py \
		-p 4326 -c ns -t UInt16 -f GTiff --no-pyramids $image_path $ortho_out \
	    > ${LOGS_OUTPUT_DIR}/${base_name}_log.txt


	input_img_basename=$(basename "${image%.[nN][tT][fF]}")
	echo "$input_img_basename"
	image2="$ortho_out/${input_img_basename}_u16ns4326.tif"
	# echo $image2

	other_ortho_fpath="$ortho_out/${input_img_basename}_u16ns4326.prj"
	# rm $image2 
	rm $other_ortho_fpath


		#### Calculate Total Time
	 #   endtime = datetime.today()
	 #   td = (endtime-starttime)
	 #   LogMsg("Total Processing Time: %s\n" %(td))
 
 done
 

