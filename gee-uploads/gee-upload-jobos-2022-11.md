The existing (131 or 181) files were deleted from gcloud.

Files were uploaded:

```
(base) tylar@manglilloo:/srv/imars-objects/jobos/Processed/JobosFinal/FinalJobosForGEE$ gsutil cp *fullClass*_Jobos.tif gs://jobos-wv-classmaps | tee ~/wv_classMaps_upload_jobos_2022-09-13.log
Copying file://20100710T151246_01_P004_WV02_ClassificMap_fullClass_vB1_Jobos.tif [Content-Type=image/tiff]...
Copying file://20110328T152054_01_P005_WV02_ClassificMap_fullClass_vB3_Jobos.tif [Content-Type=image/tiff]...
Copying file://20110328T152055_01_P006_WV02_ClassificMap_fullClass_vB3_Jobos.tif [Content-Type=image/tiff]...
Copying file://20110817T152608_01_P006_WV02_ClassificMap_fullClass_vB1_Jobos.tif [Content-Type=image/tiff]...
\ [4 files][ 38.8 MiB/ 38.8 MiB]                                                
==> NOTE: You are performing a sequence of gsutil operations that may
run significantly faster if you instead use gsutil -m cp ... Please
see the -m section under "gsutil help options" for further information
about when gsutil -m can be advantageous.

Copying file://20111107T151653_01_P006_WV02_ClassificMap_fullClass_vB2_Jobos.tif [Content-Type=image/tiff]...
Copying file://20120909T151315_01_P003_WV02_ClassificMap_fullClass_vB2_Jobos.tif [Content-Type=image/tiff]...
[...]
Copying file://20201221T145946_01_P006_WV03_ClassificMap_fullClass_vB3_Jobos.tif [Content-Type=image/tiff]...
Copying file://20210109T145803_01_P001_WV03_ClassificMap_fullClass_vB3_Jobos.tif [Content-Type=image/tiff]...
Copying file://20210220T145743_01_P006_WV02_ClassificMap_fullClass_vB3_Jobos.tif [Content-Type=image/tiff]...
- [173 files][  2.1 GiB/  2.1 GiB]   12.3 MiB/s                                 
==> NOTE: You are performing a sequence of gsutil operations that may
run significantly faster if you instead use gsutil -m cp ... Please
see the -m section under "gsutil help options" for further information
about when gsutil -m can be advantageous.

Operation completed over 173 objects/2.1 GiB.                                    
```

Next up: use `gbucket_to_gee_w_metadata_jobos.sh` to set the metadata and move to gee.
