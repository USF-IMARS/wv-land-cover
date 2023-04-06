## existing images deleted from gcloud bucket
173  images were deleted via the console.cloud.google.com interface.

## upload to gbucket

```
(base) tylar@manglilloo:/srv/imars-objects/jobos/Processed/JobosFinal/FinalJobosForGEE$ gsutil cp *ClassMap_vB3_Jobos_wDEM.tif gs://jobos-wv-classmaps

Copying file://20100201T150638_01_P006_WV02_ClassMap_vB3_Jobos_wDEM.tif [Content-Type=image/tiff]...
Copying file://20100710T151245_01_P003_WV02_ClassMap_vB3_Jobos_wDEM.tif [Content-Type=image/tiff]...
Copying file://20100710T151246_01_P004_WV02_ClassMap_vB3_Jobos_wDEM.tif [Content-Type=image/tiff]...
Copying file://20100922T152144_01_P002_WV02_ClassMap_vB3_Jobos_wDEM.tif [Content-Type=image/tiff]...
- [4 files][ 42.6 MiB/ 42.6 MiB]                                                
==> NOTE: You are performing a sequence of gsutil operations that may
run significantly faster if you instead use gsutil -m cp ... Please
see the -m section under "gsutil help options" for further information
about when gsutil -m can be advantageous.

Copying file://20100922T152205_01_P001_WV02_ClassMap_vB3_Jobos_wDEM.tif [Content-Type=image/tiff]...
[...]
Copying file://20210220T145743_01_P006_WV02_ClassMap_vB3_Jobos_wDEM.tif [Content-Type=image/tiff]...
- [179 files][  2.1 GiB/  2.1 GiB]   24.9 MiB/s                                 
==> NOTE: You are performing a sequence of gsutil operations that may
run significantly faster if you instead use gsutil -m cp ... Please
see the -m section under "gsutil help options" for further information
about when gsutil -m can be advantageous.

Operation completed over 179 objects/2.1 GiB.                                    
```

## transfer files gbucket --> GEE
```
(base) tylar@manglilloo:~/wv-land-cover$ bash gee-uploads/gbucket_to_gee_w_metadata_jobos.sh jobos-wv-classmaps /srv/imars-objects/jobos/Processed/wv_ortho_xml/ users/tylarmurray/nerrs_jobos_v03 | tee jobos_upload-2023_04_06.lo

checking if the collection users/tylarmurray/nerrs_jobos_v03 exists...
collection created.

*** Transfering file  20100201T150638_01_P006_WV02_ClassMap_vB3_Jobos_wDEM ***
*** parsing metadata...
{"dt_Y": "2010", "dt_m": "02", "dt_d": "01", "dt_H": "15", "dt_M": "06", "dt_S": "38", "number": "01", "pass_n": "006", "sat_n": "02", "algorithm_version": "B3", "dt_a": "Mon", "dt_A": "Monday", "dt_w": "1", "dt_dd": "1", "dt_b": "Feb", "dt_B": "February", "dt_mm": "2", "dt_y": "10", "dt_HH": "15", "dt_I": "03", "dt_II": "3", "dt_p": "PM", "dt_MM": "6", "dt_SS": "38", "dt_f": "000000", "dt_z": "", "dt_Z": "", "dt_j": "032", "dt_jj": "32", "dt_U": "05", "dt_W": "05", "dt_c": "Mon Feb  1 15:06:38 2010", "dt_x": "02/01/10", "dt_X": "15:06:38"}
*** estimating xml filename...
xml fname is like: 10FEB01150638-M1BS-*_01_P006.XML
*** searching for xml file...
found file: /srv/imars-objects/jobos/Processed/wv_ortho_xml/10FEB01150638-M1BS-505417676070_01_P006.XML
*** extracting properties from .xml...
 -p IMD_NUMROWS=8192  -p IMD_NUMCOLUMNS=9216  -p ABSCALFACTOR_BAND_C=0.00909474  -p ABSCALFACTOR_BAND_B=0.01257455  -p ABSCALFACTOR_BAND_G=0.00963636  -p ABSCALFACTOR_BAND_Y=0.00501895  -p ABSCALFACTOR_BAND_R=0.01862609  -p ABSCALFACTOR_BAND_RE=0.00447579  -p ABSCALFACTOR_BAND_N=0.02064348  -p ABSCALFACTOR_BAND_N2=0.00888421  -p FIRSTLINETIME=2010-02-01_15:06:38.679250  -p MEANSUNEL=48.4  -p MEANSUNAZ=145.7  -p MEANSATEL=70.6  -p MEANSATAZ=171.2  -p MEANOFFNADIRVIEWANGLE=17.1  -p CLOUDCOVER=0.154  -p MEANINTRACKVIEWANGLE=-16.3  -p MEANCROSSTRACKVIEWANGLE=5.3  -p SATID=WV02  -p MODE=FullSwath  -p SCANDIRECTION=Forward  -p FILENAME=10FEB01150638-M1BS-505417676070_01_P006.NTF 
*** formatting filename-extracted params for gee...
2010-02-01T15:06:38
*** transferring image and metadata...
Started upload task with ID: V3YQS37UW4KIK4X64CTCMWNZ
done!

[...]


*** Transfering file  20210220T145743_01_P006_WV02_ClassMap_vB3_Jobos_wDEM ***
*** parsing metadata...
{"dt_Y": "2021", "dt_m": "02", "dt_d": "20", "dt_H": "14", "dt_M": "57", "dt_S": "43", "number": "01", "pass_n": "006", "sat_n": "02", "algorithm_version": "B3", "dt_a": "Sat", "dt_A": "Saturday", "dt_w": "6", "dt_dd": "20", "dt_b": "Feb", "dt_B": "February", "dt_mm": "2", "dt_y": "21", "dt_HH": "14", "dt_I": "02", "dt_II": "2", "dt_p": "PM", "dt_MM": "57", "dt_SS": "43", "dt_f": "000000", "dt_z": "", "dt_Z": "", "dt_j": "051", "dt_jj": "51", "dt_U": "07", "dt_W": "07", "dt_c": "Sat Feb 20 14:57:43 2021", "dt_x": "02/20/21", "dt_X": "14:57:43"}
*** estimating xml filename...
xml fname is like: 21FEB20145743-M1BS-*_01_P006.XML
*** searching for xml file...
found file: /srv/imars-objects/jobos/Processed/wv_ortho_xml/21FEB20145743-M1BS-505417666010_01_P006.XML
*** extracting properties from .xml...
 -p IMD_NUMROWS=7168  -p IMD_NUMCOLUMNS=9216  -p ABSCALFACTOR_BAND_C=0.00909474  -p ABSCALFACTOR_BAND_B=0.01257455  -p ABSCALFACTOR_BAND_G=0.00963636  -p ABSCALFACTOR_BAND_Y=0.00501895  -p ABSCALFACTOR_BAND_R=0.01098462  -p ABSCALFACTOR_BAND_RE=0.00447579  -p ABSCALFACTOR_BAND_N=0.01217436  -p ABSCALFACTOR_BAND_N2=0.00888421  -p FIRSTLINETIME=2021-02-20_14:57:43.821850  -p MEANSUNEL=52.1  -p MEANSUNAZ=136.9  -p MEANSATEL=57.8  -p MEANSATAZ=167.3  -p MEANOFFNADIRVIEWANGLE=28.3  -p CLOUDCOVER=0.264  -p MEANINTRACKVIEWANGLE=-26.2  -p MEANCROSSTRACKVIEWANGLE=11.0  -p SATID=WV02  -p MODE=FullSwath  -p SCANDIRECTION=Forward  -p FILENAME=21FEB20145743-M1BS-505417666010_01_P006.NTF 
*** formatting filename-extracted params for gee...
2021-02-20T14:57:43
*** transferring image and metadata...
Started upload task with ID: TR7RI5PGXE2ICLIVX7C6ACXH
done!

```

## updated GEE module & assets
* /users/tylarmurray/nerrs_jobos_v03 set to be publicly visible
* variables updated in `/users/tylarmurray/nerrs_jobos/classmap_helpers`
  * `habitat_jobos` imported asset set to use `/users/tylarmurray/nerrs_jobos_v03`
  * `data.collection_id` variable for jobos set to `"users/tylarmurray/nerrs_jobos_v03"`

