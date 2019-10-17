from math import log


def stumpf_relative_depth(band_1, band_2):
    """
    Calculate relative depth
    (Stumpf 2003 ratio transform scaled to 1-10)

    returns:
        Relative depth estimate based on ratio transform
        0 if:
            band_1 or band_2 are <= 0
            resulting depth is < 0 or > 2
    """
    try:
        dp = (log(1000*band_1)/log(1000*band_2))
        assert dp > 0 and dp < 2
        return dp
    except (AssertionError, ValueError):
        return 0
