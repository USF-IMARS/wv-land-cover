* renamed files using rename_files.sh script created via `find . -name '*.tif' -exec echo mv {} {} \;` and then manually edited
* uploaded files using `gsutil cp *.tif gs://jobos-wv-classmaps
    * `Operation completed over 131 objects/1.6 GiB.`
* transfer gbucket -> GEE started:
```
bash gee-uploads/gbucket_to_gee_w_metadata_jobos.sh \
    jobos-wv-classmaps \
    /srv/imars-objects/jobos/Processed/wv_ortho_xml \
    users/tylarmurray/nerrs_jobos \
    | tee jobos_upload-2022_10.log
```
