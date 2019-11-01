# std modules:
from unittest import TestCase
import filecmp
import os.path
import warnings

from wv_classify.matlab_fns import geotiffread
from wv_classify.matlab_fns import geotiffwrite


class Test_geotiff_io(TestCase):
    def test_geotiff_cp(self):
        """read then write geotiff; output matches input."""
        INFILEPATH = "test_data/Rrs/16FEB12162517-M1BS-_RB_Rrs.tif"
        OUTFILEPATH = "/tmp/cp_test.tif"
        if os.path.isfile(INFILEPATH):
            data_array, data_obj = geotiffread(INFILEPATH)
            geotiffwrite(
                OUTFILEPATH, data_array, data_obj, 4326
            )
            print("checking if copied file is EXACTLY the same...")
            self.assertTrue(filecmp.cmp(INFILEPATH, OUTFILEPATH))
            # TODO: cleanup rm OUTFILEPATH
        else:
            warnings.warn("Test data not found; skipping test.")
