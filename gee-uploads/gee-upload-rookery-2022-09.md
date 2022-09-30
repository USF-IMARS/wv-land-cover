Uploading classmaps to google earth engine using the instructions from [my GEE cheatsheet](https://github.com/7yl4r/cheatsheets/blob/master/googleEarthEngine.md#imagecollection-upload).

notes from Luis:
> In the first step, you can either use gsutil or drop files manually to a gcloud bucket.
> 
> For the second step, it is important to transfer files to gee including metadata/properties for each file. 
> It might be year, location, product, satellite, etc. 
> See attached an example of a script doing it. 
> The images will be stored in an imageCollection folder you may create before to transfer files. 
> You can save the files in the project’s assets, e.g. ‘projects/imars-3d-wetlands/*new_folder*/*new_imageCollection*/’ 
> 
> Once the images are in a collection use gee for doing mosaics or any other analysis you require.

The mentioned script is [here](https://github.com/USF-IMARS/wv-land-cover/blob/a752e1e45b16357ed2a97db89f0f345db46abfc9/gee-uploads/gbucket_to_gee_w_metadata.sh).

The parameters being passed with the `-p` flag could be augmented using additional metadata from the `.xml` file corresponding to each `.tif`.
[This script](https://github.com/USF-IMARS/wv-land-cover/blob/master/wv_classify/read_wv_xml.py) could be very helpful with that...
