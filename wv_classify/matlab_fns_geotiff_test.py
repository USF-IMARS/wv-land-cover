# std modules:
from unittest import TestCase
import filecmp

from matlab_fns import geotiffread
from matlab_fns import geotiffwrite


class Test_geotiff_io(TestCase):
    def test_geotiff_cp(self):
        """read then write geotiff; output matches input."""
        INFILEPATH = "test/WV_Testing/Rrs/16FEB12162517-M1BS-_RB_Rrs.tif"
        OUTFILEPATH = "test/cp_test.tif"
        data_array, data_obj = geotiffread(INFILEPATH)
        geotiffwrite(
            OUTFILEPATH, data_array, data_obj, 4326,
            row_index=1,
            col_index=2,
            band_index=0
        )
        self.assertTrue(filecmp.cmp(INFILEPATH, OUTFILEPATH))
        # TODO: cleanup rm OUTFILEPATH
