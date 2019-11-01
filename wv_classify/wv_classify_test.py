import os
import warnings

from wv_classify.wv_classify import process_file


def test_process_file():
    """
    Compare vs Matt's test images 2019-02.

    I've got 2 test images in the following directory:
        lab_share/mmccarthy/Matt/WV_Testing/

    The files are organized as follows:
    ## INPUTS
    - /Raw = Raw NTF & XML files

    ## OUTPUTS:
    ### for pgc_ortho:
    - /Ortho = resampled tifs output from the pgc_ortho script
    ### for wv_classify:
    - /Rrs = calibrated Rrs tifs -THESE ARE THE IMPORTANT ONES TO COMPARE
    - /rrs_subsurface = output rrs tifs - also important to test, but less so
    - /Mapped = filtered and unfiltered classification maps
    """
    os.chdir("test_data/")
    INFILEPATH = "Ortho/16FEB12162517-M1BS-057380245010_01_P001_u16ns4326.tif"
    if os.path.isfile(INFILEPATH):
        OUTPUT_DIR = "/tmp/wv_test_output"
        try:
            os.mkdir(OUTPUT_DIR)
        except FileExistsError:
            pass
        process_file(
            INFILEPATH,
            "Raw/16FEB12162517-M1BS-057380245010_01_P001.xml",
            OUTPUT_DIR,
            "RB",
            4326,
            d_t=1,
            Rrs_write=1
        )
        print(
            "test output written to {}; it must be verified manually.".format(
                OUTPUT_DIR
            )
        )
        # TODO: assert something?
        raise NotImplementedError(
            "Automated checking of output not yet implemented"
        )
        # TODO: cleanup files.
    else:
        warnings.warn("Test data not found; skipping test.")
