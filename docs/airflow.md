## Apache Airflow
This task has been run as an airflow DAG on IMaRS's airflow cluster.
The airflow dag definition file can be viewed at [USF-IMARS/imars_dags//dags/processing/wv2_classification/wv_classification.py](https://github.com/USF-IMARS/imars_dags/blob/master/dags/processing/wv2_classification/wv_classification.py).
Accompanying wrapper scripts are in the [scripts subdir of the same location](https://github.com/USF-IMARS/imars_dags/tree/master/dags/processing/wv2_classification/scripts).
Of particular note is [USF-IMARS/imars_dags//dags/processing/wv2_classification/scripts/ntf_to_rrs.sh](https://github.com/USF-IMARS/imars_dags/blob/master/dags/processing/wv2_classification/scripts/ntf_to_rrs.sh).

USF-IMaRS/imars_dags expects a certain configuration & software suite is expected to exist on each node.
This airflow cluster's software and configuration management is managed via puppet ([IMaRS-private puppet repo ln](https://github.com/usf-imars/imars_puppet)); related documentation can be found there and can be provided on request.
One of the more important dependencies is [imars-etl](https://github.com/USF-IMARS/imars-etl), which wraps IMaRS's underlying object & metadata storage systems.
