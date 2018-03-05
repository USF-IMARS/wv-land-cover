#!/bin/bash
#SBATCH --job-name ="wv2_classification_py"
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20480
#SBATCH --time=03:00:00
#SBATCH --array=0


## Setup input arguments & file locations
WORK_DIR="$WORK/tmp/test/sunglint/test2"
OUT_DIR="/work/m/mjm8/tmp/test"

# uncomment these for testing
SLURM_ARRAY_TASK_ID=1
WORK_DIR="/home/tylar/wv2-processing/test/work"
OUT_DIR="/home/tylar/wv2-processing/test/out"

images1=`ls $WORK_DIR/*.ntf`
met=`ls $WORK_DIR/*.xml`

output_dir1="$OUT_DIR/ortho/"
rrs_out="$OUT_DIR/output/"
class_out="$OUT_DIR/output/"

crd_sys=EPSG:4326
dt=0
sgw=5
filt=0
stat=3
loc='testnew'

# Run Python code
images1a=($images1)
image=${images1a[$SLURM_ARRAY_TASK_ID]}
echo "processing work file $image..."
input_img_basename=`basename -s .ntf $image`
echo "basename is: $input_img_basename"

# python /work/m/mjm8/progs/pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no_pyramids $image $output_dir1


## Run Matlab code
# module add apps/matlab/r2013b
image2="$output_dir1${input_img_basename}_uint164326.tif"
met=($met)
met=${met[$SLURM_ARRAY_TASK_ID]}
echo "proccessing ortho file $image2..."

# matlab -nodisplay -nodesktop -r "WV2_Processing('$image2','$met','$crd_sys','$dt','$sgw','$filt','$stat','$loc','$SLURM_ARRAY_TASK_ID','$rrs_out','$class_out')"

# Python code to check processing time:
#    starttime = datetime.today()
#    LogMsg('Image: %s' %(info.srcfn))

#    #### Calculate Total Time
#    endtime = datetime.today()
#    td = (endtime-starttime)
#    LogMsg("Total Processing Time: %s\n" %(td))
