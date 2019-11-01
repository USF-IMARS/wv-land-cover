# std modules:
from unittest import TestCase

from wv_classify.read_wv_xml import read_wv_xml


class Test_read_wv_xml(TestCase):
    def test_dig_globe_xml_read(self):
        """test reading xml from digital globe."""
        (
            szB, aqmonth, aqyear, aqhour, aqminute, aqsecond, sunaz, sunel,
            satel, sensaz, aqday, satview, kf, cl_cov
        ) = read_wv_xml("test_data/xml/from_digital_globe.xml")
        self.assertEqual(cl_cov, 1.500000000000000e-02)

    def test_pgc_ortho_xml_read(self):
        """test reading xml post pgc ortho."""
        (
            szB, aqmonth, aqyear, aqhour, aqminute, aqsecond, sunaz, sunel,
            satel, sensaz, aqday, satview, kf, cl_cov
        ) = read_wv_xml("test_data/xml/post_pgc_ortho.xml")
        self.assertEqual(cl_cov, 1.500000000000000e-02)
