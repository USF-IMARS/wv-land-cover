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
ortho_out=/work/m/mjm8/output/Ortho/NSF_SWTX/  # ortho_out
rrs_out=/work/m/mjm8/output/Rrs/NSF_SWTX/
class_out=/work/m/mjm8/output/DT/NSF_SWTX/

# Setup Matlab arguments (dt 0 = Rrs, no DT or rrs; dt 1 = Rrs, DT & rrs; dt 2 = DT, Rrs & rrs | filt=moving-window filter)
dt=2
crd_sys=EPSG:4326
filt=2
loc='NSF_SWTX'

## === Run Python code
images1a=($images1)  # cast to array
input_image=${images1a[$SLURM_ARRAY_TASK_ID]}

# figure out output filepaths
input_img_basename=$(basename "${input_image%.[nN][tT][fF]}")
echo $input_img_basename
image2="$ortho_out${input_img_basename}_u16ns4326.tif"
echo $image2
other_ortho_fpath="$ortho_out${input_img_basename}_u16ns4326.prj"


if [ ! -f $image2 ]; then  # if output file DNE
    module add apps/python/2.7.5
    python /work/m/mjm8/progs/pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no_pyramids $input_image $ortho_out
fi

## === Run Matlab code
met=($met)
met=${met[$SLURM_ARRAY_TASK_ID]}

final_output_path="$rrs_out${input_img_basename}_$loc_SOALCHI_filt_$filt.tif"

if [ ! -f $final_output_path ]; then  # if output file DNE
    module add apps/matlab/r2017a
    matlab -nodisplay -nodesktop -r "wv_classify('$image2','$input_img_basename','$met','$crd_sys','$dt','$filt','$loc','$SLURM_ARRAY_TASK_ID','$rrs_out','$class_out')"
fi

rm $image2 
rm $other_ortho_fpath

    #### Calculate Total Time
 #   endtime = datetime.today()
 #   td = (endtime-starttime)
 #   LogMsg("Total Processing Time: %s\n" %(td))
