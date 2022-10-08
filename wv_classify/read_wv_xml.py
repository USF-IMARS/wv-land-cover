from xml.etree import ElementTree
from datetime import datetime


def read_wv_xml(filename, output_format="list"):
    """ read metadata values from xml file
    params
    ------
    filename : filepath
        The .xml file to read
    output_format : str
        The format to return output.
        Note that "list" output may not support all metadata.
        valid values:
            "list" - [param1Value, param2Value]
            "dict" - {param1:value, param2:value]
    """
    metadata = {}
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
    metadata["IMD_NUMROWS"]    = int(imd.find('NUMROWS').text)
    metadata["IMD_NUMCOLUMNS"] = int(imd.find('NUMCOLUMNS').text)
    BANDNAMES = [
        'BAND_C', 'BAND_B', 'BAND_G', 'BAND_Y', 'BAND_R', 'BAND_RE',
        'BAND_N', 'BAND_N2'
    ]
    for band in BANDNAMES:
        metadata[f"ABSCALFACTOR_{band}"] = float(
            imd.find(band).find('ABSCALFACTOR').text
        )
    # Extract Acquisition Time from metadata
    metadata["FIRSTLINETIME"] = datetime.strptime(
        imd.find('IMAGE').find('FIRSTLINETIME').text,
        # "2017-12-22T16:48:10.923850Z"
        "%Y-%m-%dT%H:%M:%S.%fZ"
    )
    # Extract Mean Sun Elevation angle from metadata.Text(18:26))
    # Extract Mean Off Nadir View angle from metadata
    for param in [
        "MEANSUNEL", "MEANSUNAZ",
        "MEANSATEL", "MEANSATAZ",
        "MEANOFFNADIRVIEWANGLE", "CLOUDCOVER",
        "MEANINTRACKVIEWANGLE", "MEANCROSSTRACKVIEWANGLE", "MEANOFFNADIRVIEWANGLE"
    ]:
        metadata[param] = float(imd.find("IMAGE").find(param).text)
    for param in [
        "SATID", "MODE", "SCANDIRECTION",
    ]:
        metadata[param] = imd.find("IMAGE").find(param).text

    for param in [
        "FILENAME"
    ]:
        metadata[param] = root.find("TIL").find("TILE").find(param).text
        
    if output_format == "list":
        szB = [
            metadata["IMD_NUMROWS"],
            metadata["IMD_NUMCOLUMNS"],
            0
        ]
        kf = [
            metadata[f"ABSCALFACTOR_{band}"] for band in BANDNAMES
        ]
        aq_dt = metadata["FIRSTLINETIME"]
        aqyear = aq_dt.year
        aqmonth = aq_dt.month
        aqday = aq_dt.day
        aqhour = aq_dt.hour
        aqminute = aq_dt.minute
        aqsecond = aq_dt.second
        sunel = metadata['MEANSUNEL']
        satview = metadata['MEANOFFNADIRVIEWANGLE']
        sunaz = metadata['MEANSUNAZ']
        sensaz = metadata['MEANSATAZ']
        satel = metadata['MEANSATEL']
        cl_cov = metadata['CLOUDCOVER']
        return (
            szB, aqmonth, aqyear, aqhour, aqminute, aqsecond, sunaz, sunel,
            satel, sensaz, aqday, satview, kf, cl_cov
        )
    elif output_format == "dict":
        return metadata
    elif output_format == "gee_props":
        res = ""
        for key in metadata:
            try:
                val = str(metadata[key]).replace(" ", "_")
            except:
                val = metadata[key]
            res += f" -p {key}={val} "
        return res
    else:
        raise ValueError(
            f"user requested unknown output_format '{output_format}'"
        )


if __name__ == "__main__":
    #import pprint
    import sys
    #pp = pprint.PrettyPrinter(indent=2)
    fpath = sys.argv[1]
    #pp.pprint(
    #    read_wv_xml(fpath, output_format="dict")
    #)
    print(read_wv_xml(fpath, output_format="gee_props"))
