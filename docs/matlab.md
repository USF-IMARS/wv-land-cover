

2. matlab:
    ```
    matlab -nodisplay -nodesktop -r "\
        cd('/opt/wv2_processing');\
        wv2_processing(\
            '$ORTH_FILE',\
            '{{params.id}}',\
            '$MET',\
            '{{params.crd_sys}}',\
            '{{params.dt}}',\
            '{{params.sgw}}',\
            '{{params.filt}}',\
            '{{params.stat}}',\
            '{{params.loc}}',\
            '{{params.id_number}}',\
            '$RRS_OUT',\
            '$CLASS_OUT'\
        );\
        exit\
    "
    ```
