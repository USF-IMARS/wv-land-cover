#!/bin/bash
#SBATCH --job-name ="wv2_classification_py"
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20480
#SBATCH --TIME=3:00:00
#SBATCH --array=100-199
##0-611%30
##SBATCH --array=0-611%20
## Can submit up to 10,000 jobs at once, but only 512 will run concurrently

# Python code to check processing time:
#    starttime = datetime.today()
#    LogMsg('Image: %s' %(info.srcfn))


## Setup input arguments & file locations
images1=`ls $WORK/tmp/NSF/raw/*.[nN][tT][fF]`
met=`ls $WORK/tmp/NSF/raw/*.[xX][mM][lL]`
output_dir1=/work/m/mjm8/output/Ortho/NSF_SWTX/
rrs_out=/work/m/mjm8/output/Rrs/NSF_SWTX/
class_out=/work/m/mjm8/output/DT/NSF_SWTX/

# Setup Matlab arguments (dt 0 = Rrs, no DT or rrs; dt 1 = Rrs, DT & rrs; dt 2 = DT, Rrs & rrs | filt=moving-window filter)
dt=2
crd_sys=EPSG:4326
filt=2
loc='NSF_SWTX'

## Run Python code
images1a=($images1)
image=${images1a[$SLURM_ARRAY_TASK_ID]}

python /work/m/mjm8/progs/pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no_pyramids $image $output_dir1


## Run Matlab code
#module add apps/matlab/r2013b
module add apps/matlab/r2017a

input_img_basename=$(basename "${image%.[nN][tT][fF]}")
echo $input_img_basename
image2="$output_dir1${input_img_basename}_u16ns4326.tif"
echo $image2
met=($met)
met=${met[$SLURM_ARRAY_TASK_ID]}

matlab -nodisplay -nodesktop -r "WV_Processing('$image2','$input_img_basename','$met','$crd_sys','$dt','$filt','$loc','$SLURM_ARRAY_TASK_ID','$rrs_out','$class_out')"

other_ortho_fpath="$output_dir1${input_img_basename}_u16ns4326.prj"
rm $image2 
rm $other_ortho_fpath

    #### Calculate Total Time
 #   endtime = datetime.today()
 #   td = (endtime-starttime)
 #   LogMsg("Total Processing Time: %s\n" %(td))