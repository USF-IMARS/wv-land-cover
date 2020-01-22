## SLURM
These processing tasks have been executed on USF Research Computing's research cluster, [CIRCE](https://wiki.rc.usf.edu/index.php/CIRCE) using the SLURM task scheduler.
The slurm submission bash script is in the root of this repo at [./submit_py.sh](https://github.com/USF-IMARS/wv2-processing/blob/master/submit_py.sh)


The submit_py.sh file is what I use in Circe to call the pgc_ortho.py script, which has a number of sub-scripts called.
The submit_py.sh also contains the Matlab script call, so you'll want to comment out those lines before testing it.
