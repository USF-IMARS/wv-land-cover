# general gist of how to update the dataset
1. create a bucket on google cloud.
2. upload files from server using `gsutil cp *.tif gs://{{bucket_name}}` (use bucket name from (1))
3. transfer gbucket files into GEE using something like:
    ```
    bash gee-uploads/gbucket_to_gee_w_metadata_jobos.sh \
        jobos-wv-classmaps \
        /srv/imars-objects/jobos/Processed/wv_ortho_xml \
        users/tylarmurray/nerrs_jobos \
        | tee jobos_upload-2022_10.log
    ```
