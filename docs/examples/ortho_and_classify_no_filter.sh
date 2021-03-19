# === set basic variables
# file paths
ortho_out=/work/m/mjm8/output/Ortho/NSF_SWTX/  # ortho_out
rrs_out=/work/m/mjm8/output/Rrs/NSF_SWTX/
class_out=/work/m/mjm8/output/DT/NSF_SWTX/

input_image=/work/m/mjm8/input/raw/FILENAME.NTF
input_met=/work/m/mjm8/input/raw/FILENAME.XML

# Matlab arguments
dt=2  # dt 0 = Rrs, no DT, no rrs; dt 1 = Rrs, DT, & rrs; dt 2 = DT, Rrs, & rrs
crd_sys=EPSG:4326
filt=2  # filt=moving-window filter. 2 is 5x5.
loc='rookery'

# other
SLURM_ARRAY_TASK_ID=0

## === Run PGC Orthorectification Python code
# figure out output filepaths pgc_ortho will write to
input_img_basename=$(basename "${input_image%.[nN][tT][fF]}")
echo $input_img_basename
image2="$ortho_out${input_img_basename}_u16ns4326.tif"
echo $image2
other_ortho_fpath="$ortho_out${input_img_basename}_u16ns4326.prj"

# run orthorectification
python /work/m/mjm8/progs/pgc_ortho.py -p 4326 -c ns -t UInt16 -f GTiff --no_pyramids $input_image $ortho_out

## === Run Matlab code decision tree
final_output_path="$rrs_out${input_img_basename}_$loc_SOALCHI_filt_$filt.tif"

matlab -nodisplay -nodesktop -r "wv_classify('$image2','$input_img_basename','$met','$crd_sys','$dt','$filt','$loc','$SLURM_ARRAY_TASK_ID','$rrs_out','$class_out')"

# === clean up intermediate ortho files
rm $image2 
rm $other_ortho_fpath
