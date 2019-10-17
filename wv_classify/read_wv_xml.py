from xml.etree import ElementTree
from datetime import datetime


def read_wv_xml(filename):
    # ==================================================================
    # === read values from xml file
    # ==================================================================
    # Extract calibration factors & acquisition time from
    # metadata for each band
    tree = ElementTree.parse(filename)
    root = tree.getroot()  # tag == 'isd' or 'IMD'
    if root.tag == 'isd':  # xml from digital globe
        pass
    elif root.tag == 'IMD':  # xml output from pgc_ortho
        root = root.find('SOURCE_IMD')
    else:
        raise ValueError(
            "unrecognized .xml format - root element expected to be one of "
            " (isd, IMD)"
        )
    imd = root.find('IMD')  # assumes only one element w/ 'IMD' tag
    szB = [
        int(imd.find('NUMROWS').text),
        int(imd.find('NUMCOLUMNS').text),
        0
    ]
    kf = [
        float(imd.find(band).find('ABSCALFACTOR').text) for band in [
            'BAND_C', 'BAND_B', 'BAND_G', 'BAND_Y', 'BAND_R', 'BAND_RE',
            'BAND_N', 'BAND_N2'
        ]
    ]
    # Extract Acquisition Time from metadata
    aq_dt = datetime.strptime(
        imd.find('IMAGE').find('FIRSTLINETIME').text,
        # "2017-12-22T16:48:10.923850Z"
        "%Y-%m-%dT%H:%M:%S.%fZ"
    )
    aqyear = aq_dt.year
    aqmonth = aq_dt.month
    aqday = aq_dt.day
    aqhour = aq_dt.hour
    aqminute = aq_dt.minute
    aqsecond = aq_dt.second
    # Extract Mean Sun Elevation angle from metadata.Text(18:26))
    sunel = float(imd.find('IMAGE').find('MEANSUNEL').text)
    # Extract Mean Off Nadir View angle from metadata
    satview = float(imd.find('IMAGE').find('MEANOFFNADIRVIEWANGLE').text)
    sunaz = float(imd.find('IMAGE').find('MEANSUNAZ').text)
    sensaz = float(imd.find('IMAGE').find('MEANSATAZ').text)
    satel = float(imd.find('IMAGE').find('MEANSATEL').text)
    cl_cov = float(imd.find('IMAGE').find('CLOUDCOVER').text)
    # ==================================================================
    return (
        szB, aqmonth, aqyear, aqhour, aqminute, aqsecond, sunaz, sunel,
        satel, sensaz, aqday, satview, kf, cl_cov
    )
