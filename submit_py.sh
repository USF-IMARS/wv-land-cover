#!/bin/bash
#SBATCH --job-name ="wv2_classification_py"
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20480
#SBATCH --time=03:00:00
#SBATCH --array=0


## Setup input arguments & file locations
images1=`ls $WORK/tmp/test/sunglint/test2/*.ntf`
met=`ls $WORK/tmp/test/sunglint/test2/*.xml`
output_dir1=/work/m/mjm8/tmp/test/ortho/
rrs_out=/work/m/mjm8/tmp/test/output/
class_out=/work/m/mjm8/tmp/test/output/

crd_sys=EPSG:4326
dt=0
sgw=5
filt=0
stat=3
loc='testnew'

# Run Python code
images1a=($images1)
image=${images1a[$SLURM_ARRAY_TASK_ID]}

python /work/m/mjm8/progs/pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no_pyramids $image $output_dir1


## Run Matlab code
module add apps/matlab/r2013b
input_img_basename=`basename -s .ntf $image`
image2="$output_dir1${input_img_basename}_uint164326.tif"
met=($met)
met=${met[$SLURM_ARRAY_TASK_ID]}

matlab -nodisplay -nodesktop -r "WV2_Processing('$image2','$met','$crd_sys','$dt','$sgw','$filt','$stat','$loc','$SLURM_ARRAY_TASK_ID','$rrs_out','$class_out')"

# Python code to check processing time:
#    starttime = datetime.today()
#    LogMsg('Image: %s' %(info.srcfn))

#    #### Calculate Total Time
#    endtime = datetime.today()
#    td = (endtime-starttime)
#    LogMsg("Total Processing Time: %s\n" %(td))
