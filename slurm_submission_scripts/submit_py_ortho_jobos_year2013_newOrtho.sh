#!/bin/bash
#SBATCH --partition=circe
#SBATCH --job-name ="Jobos-per-year"
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=52240
#SBATCH --time=40:00:00
#SBATCH --array=0-2
#SBATCH --output=/work/d/druedaro/wv2_scripts/slurm_submission_scripts/SlurmOutput/output.%j.txt
## Can submit up to 10,000 jobs at once, but only 512 will run concurrently
## SBATCH --mail-type=ALL
## SBATCH --mail-user=druedaro@usf.edu


module purge
module add apps/python/2.7.5
# module add apps/python/3.7.3  # 2.7.5
# module add apps/gdal/3.0.1

# module add apps/proj/4.9.3
# module add apps/proj/6.2.0
# module add apps/proj/6.2.0_el6
module add apps/proj/6.2.0_el7_gcc
# module add apps/proj/backup

module add apps/gdal/3.0.1_el7_gcc

# Python code to check processing time:
#    starttime = datetime.today()
#    LogMsg('Image: %s' %(info.srcfn))

## Setup input arguments & file locations
images1=`ls /work/d/druedaro/img/Jobos_perYear/Jobos_2013b/*.[nN][tT][fF]`
met=`ls /work/d/druedaro/img/Jobos_perYear/Jobos_2013b/*.[xX][mM][lL]`
ortho_out=/work/d/druedaro/output/Ortho/Jobos_perYear/Jobos_2013b/
rrs_out=/work/d/druedaro/output/Rrs/Jobos_perYear/Jobos_2013b/


# Setup Matlab arguments
crd_sys=EPSG:4326
areaName='Jobos'

## Run Python code
images1a=($images1)
image=${images1a[$SLURM_ARRAY_TASK_ID]}
echo "orthorectifying $image to $ortho_out"

python /work/d/druedaro/wv2_scripts/pgc_imagery_utils2/pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no-pyramids $image $ortho_out


## Run Matlab code
module add apps/matlab/r2018a

input_img_basename=$(basename "${image%.[nN][tT][fF]}")
echo $input_img_basename
image2="$ortho_out${input_img_basename}_u16ns4326.tif"
echo $image2
met=($met)
met=${met[$SLURM_ARRAY_TASK_ID]}

# matlab running
matlab -nodisplay -nodesktop -r "WV_Clasific_cleaner_v5_Jobos_v2019.m('$image2','$met','crd_sys','$areaName','$rrs_out','$SLURM_ARRAY_TASK_ID')"


# other_ortho_fpath="$ortho_out${input_img_basename}_u16ns4326.prj"
# rm $image2 
# rm $other_ortho_fpath


    #### Calculate Total Time
 #   endtime = datetime.today()
 #   td = (endtime-starttime)
 #   LogMsg("Total Processing Time: %s\n" %(td))
 
 
 

