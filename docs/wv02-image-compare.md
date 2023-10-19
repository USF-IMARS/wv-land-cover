Looking at a time-series extraction of some land-cover points, some temporal features can be easily seen:

![image](https://github.com/USF-IMARS/wv-land-cover/assets/1051390/8f9b1fce-2f87-49d1-98a5-c122e38e4ae2)

Looking specifically at the "mud flats" (lower time-series plot), we can inspect the two images nearest to the sudden decrease in multiple bands seen around the 2018 mark:

![image](https://github.com/USF-IMARS/wv-land-cover/assets/1051390/7cbc3689-03cd-4414-827c-b3d31d41d1af)

![image](https://github.com/USF-IMARS/wv-land-cover/assets/1051390/009347b4-fc16-4e81-862f-06b718904642)

In these screenshots: 
* right, top, older (2017-06-25)
* left, bottom, newer (2018-01-06)

The difference between these true-color views is huge.
The max shown on each histogram x-axis is approximately the same (around .06).
I cannot explain the difference.

## Scaling?
Associated with each image is a set of `ABSCALFACTOR_BAND_*` properties for each band.
These images should have been scaled prior to the upload by processing done on CIRCE, however, we can try applying the scale factor in GEE to check.
The following is a mapping between band layer names in the images and ABSCALFACTOR names:

```js
  var bandNameMap = {
    b1: "ABSCALFACTOR_BAND_C",  // CA
    b2: "ABSCALFACTOR_BAND_B",  // Blue
    b3: "ABSCALFACTOR_BAND_G",  // Green
    b4: "ABSCALFACTOR_BAND_Y",  // Yellow
    b5: "ABSCALFACTOR_BAND_R",  // Red
    b6: "ABSCALFACTOR_BAND_RE", // Red Edge
    b7: "ABSCALFACTOR_BAND_N",  // NIR
    b8: "ABSCALFACTOR_BAND_N2", // NIR2
  };
  var bandNameMap = {
    b8: "ABSCALFACTOR_BAND_C",  // CA
    b7: "ABSCALFACTOR_BAND_B",  // Blue
    b6: "ABSCALFACTOR_BAND_G",  // Green
    b5: "ABSCALFACTOR_BAND_Y",  // Yellow
    b4: "ABSCALFACTOR_BAND_R",  // Red
    b3: "ABSCALFACTOR_BAND_RE", // Red Edge
    b2: "ABSCALFACTOR_BAND_N",  // NIR
    b1: "ABSCALFACTOR_BAND_N2", // NIR2
  };


```

After trying either scaling map above on the images major differences between the images remain, and the color views look terrible.
This implies that the scaling has already been done as expected.
