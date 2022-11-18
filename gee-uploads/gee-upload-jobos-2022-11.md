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

Next up: use `gbucket_to_gee_w_metadata_jobos.sh` to set the metadata and move to gee:

```
(base) tylar@manglilloo:~/wv-land-cover$ bash gee-uploads/gbucket_to_gee_w_metadata_jobos.sh jobos-wv-classmaps /srv/imars-objects/jobos/Processed/wv_ortho_xml/ users/tylarmurray/nerrs_jobos | tee jobos_upload-2022_11_18.log
checking if the collection users/tylarmurray/nerrs_jobos exists...
Asset users/tylarmurray/nerrs_jobos already exists.

*** Transfering file  20100710T151246_01_P004_WV02_ClassificMap_fullClass_vB1_Jobos ***
*** parsing metadata...
{"dt_Y": "2010", "dt_m": "07", "dt_d": "10", "dt_H": "15", "dt_M": "12", "dt_S": "46", "number": "01", "pass_n": "004", "sat_n": "02", "algorithm_version": "B1", "dt_a": "Sat", "dt_A": "Saturday", "dt_w": "6", "dt_dd": "10", "dt_b": "Jul", "dt_B": "July", "dt_mm": "7", "dt_y": "10", "dt_HH": "15", "dt_I": "03", "dt_II": "3", "dt_p": "PM", "dt_MM": "12", "dt_SS": "46", "dt_f": "000000", "dt_z": "", "dt_Z": "", "dt_j": "191", "dt_jj": "191", "dt_U": "27", "dt_W": "27", "dt_c": "Sat Jul 10 15:12:46 2010", "dt_x": "07/10/10", "dt_X": "15:12:46"}
*** estimating xml filename...
xml fname is like: 10JUL10151246-M1BS-*_01_P004.XML
*** searching for xml file...
found file: /srv/imars-objects/jobos/Processed/wv_ortho_xml/10JUL10151246-M1BS-505417665010_01_P004.XML
*** extracting properties from .xml...
 -p IMD_NUMROWS=9216  -p IMD_NUMCOLUMNS=9216  -p ABSCALFACTOR_BAND_C=0.00909474  -p ABSCALFACTOR_BAND_B=0.01773333  -p ABSCALFACTOR_BAND_G=0.01358974  -p ABSCALFACTOR_BAND_Y=0.00671549  -p ABSCALFACTOR_BAND_R=0.01862609  -p ABSCALFACTOR_BAND_RE=0.00598873  -p ABSCALFACTOR_BAND_N=0.02064348  -p ABSCALFACTOR_BAND_N2=0.00888421  -p FIRSTLINETIME=2010-07-10_15:12:46.879750  -p MEANSUNEL=71.4  -p MEANSUNAZ=73.1  -p MEANSATEL=82.4  -p MEANSATAZ=164.8  -p MEANOFFNADIRVIEWANGLE=6.8  -p CLOUDCOVER=0.035  -p MEANINTRACKVIEWANGLE=-6.2  -p MEANCROSSTRACKVIEWANGLE=2.8  -p SATID=WV02  -p MODE=FullSwath  -p SCANDIRECTION=Reverse  -p FILENAME=10JUL10151246-M1BS-505417665010_01_P004.NTF 
*** formatting filename-extracted params for gee...
2010-07-10T15:12:46
*** transferring image and metadata...
Started upload task with ID: B2EVHJXOK35JFPOIMZZECMTM
done!


*** Transfering file  20110328T152054_01_P005_WV02_ClassificMap_fullClass_vB3_Jobos ***
[...]
```
