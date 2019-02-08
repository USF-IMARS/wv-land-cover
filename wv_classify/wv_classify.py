# WV2 Processing
#
# Author: Matt McCarthy
# ported from MATLAB by Tylar Murray
#
# Loads TIFF WorldView-2 image files preprocessed through Polar Geospatial
# Laboratory python code, which orthorectifies and projects .NTF files and
# outputs as
# TIFF files
# Radiometrically calibrates digital count data
# Atmospherically corrects images by subtracting Rayleigh Path Radiance
# Converts image to surface reflectance by accounting for Earth-Sun
# distance, solar zenith angle, and average spectral irradiance
# Tests and optionally corrects for sunglint
# Corrects for water column attenuation
# Runs Decision Tree classification on each image
# Optionally smooths results through moving-window filter
# Outputs images as GEOTIFF files with geospatial information.

# built-in imports:
import sys
from os import path
from glob import glob
from math import pi
from math import exp
from math import log
import numpy
from numpy import zeros
from numpy import mean
from numpy import isnan
from numpy import std
from xml.etree import ElementTree
from datetime import datetime

# dep packages:
from skimage.morphology import square as square_strel
from skimage.morphology import white_tophat as imtophat
from skimage.filters import threshold_otsu as imbinarize

# local imports:
from DT_Filter import DT_Filter
from matlab_fns import geotiffread
from matlab_fns import geotiffwrite
from matlab_fns import cosd
from matlab_fns import sind
from matlab_fns import tand
from matlab_fns import acosd
from matlab_fns import asind
from matlab_fns import mldivide
from matlab_fns import rdivide


# TODO: + printout timing of run

DATA_DIR = '/home1/mmccarthy/Matt/USF/Other/NERRS_Mapping/Processing'
OUTPUT_NaN = numpy.nan
# dst_ds.GetRasterBand(1).SetNoDataValue(OUTPUT_NaN)
# === Assign constants for all images
# Effective Bandwidth per band
# (nm converted to um units; from IMD metadata files)
ebw = [
    e * 0.001 * e for e in [47.3, 54.3, 63.0, 37.4, 57.4, 39.3, 98.9, 99.6]
]
# Band-averaged Solar Spectral Irradiance (W/m2/um units)
irr = [
    1758.2229, 1974.2416, 1856.4104, 1738.4791, 1559.4555, 1342.0695,
    1069.7302, 861.2866
]
# Center wavelength
# (used for Rayleigh correction; from Radiometric Use of WorldView-2
# Imagery)
cw = [0.4273, 0.4779, 0.5462, 0.6078, 0.6588, 0.7237, 0.8313, 0.9080]
# Factor used in Rayleigh Phase Function equation (Bucholtz 1995)
gamma = [
    g * 0.01 for g in
    [1.499, 1.471, 1.442, 1.413, 1.413, 1.413, 1.384, 1.384]
]


def read_xml(filename):
    # ==================================================================
    # === read values from xml file
    # ==================================================================
    # Extract calibration factors & acquisition time from
    # metadata for each band
    tree = ElementTree.parse(filename)
    root = tree.getroot()  # assumes tag == 'isd'
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
    # TODO: why this if/else?
    # if isfield(s, 'IMD') == 1:
    #     c = struct2cell(s.Children(2).Children(:))
    # else
    # # end
    # ==================================================================
    return (
        szB, aqmonth, aqyear, aqhour, aqminute, aqsecond, sunaz, sunel,
        satel, sensaz, aqday, satview, kf, cl_cov
    )


def process_file(
    X,  # MS Tiff input image path
    Z,  # XML met input file path
    loc_out,  # output directory
    loc,  # RoI identifier string
    coor_sys=4326,  # coordinate system code
    d_t=1,  # 0=End after Rrs conversion; 1=rrs, bathy, DT; 2 = rrs, bathy & DT
    Rrs_write=1,  # 1=write Rrs geotiff; 0=do not write
):
    """
    process a single set of files
    """
    fname = path.basename(X)
    id = fname[1:19]

    A, R = geotiffread(X)
    szA = A.shape

    (
        szB, aqmonth, aqyear, aqhour, aqminute, aqsecond, sunaz, sunel,
        satel, sensaz, aqday, satview, kf, cl_cov
    ) = read_xml(Z)

    szB[3] = 8

    # ==================================================================
    # === Calculate Earth-Sun distance and relevant geometry
    # ==================================================================
    if aqmonth == 1 or aqmonth == 2:
        year = aqyear - 1
        month = aqmonth + 12
    else:
        year = aqyear
        month = aqmonth
    # end
    # Convert time to UT
    UT = aqhour + (aqminute/60.0) + (aqsecond/3600.0)
    B1 = int(year/100)
    B2 = 2-B1+int(B1/4)
    # Julian date
    JD = (
        int(365.25*(year+4716)) + int(30.6001*(month+1)) + aqday +
        UT/24.0 + B2 - 1524.5
    )
    D = JD - 2451545.0
    degs = float(357.529 + 0.98560028*D)  # Degrees
    # Earth-Sun distance at given date (should be between 0.983 and 1.017)
    ESd = 1.00014 - 0.01671*cosd(degs) - 0.00014*cosd(2*degs)

    inc_ang = 90.0 - sunel
    # Atmospheric spectral transmittance in solar path with solar
    # zenith angle
    TZ = cosd(inc_ang)
    # Atmospheric spectral transmittance in view path with satellite
    # view angle
    TV = cosd(satview)
    # ==================================================================

    # ==================================================================
    # === Calculate Rayleigh Path Radiance
    # ==================================================================
    # (Dash et al. 2012 and references therein)
    # For the following equations, azimuths should be
    # between -180 and +180 degrees
    if sunaz > 180:
        sunaz = sunaz - 360
    # end
    if sensaz > 180:
        sensaz = sensaz - 360
    # end

    az = abs(sensaz - 180 - sunaz)  # Relative azimuth angle
    # Scattering angles
    thetaplus = acosd(
        cosd(90-sunel)*cosd(90-satel) -
        sind(90-sunel)*sind(90-satel)*cosd(az)
    )
    Pr = [0]*8
    for d in range(8):
        # Rayleigh scattering phase function (described in Bucholtz 1995)
        Pr[d] = (
            (3/(4*(1+2*gamma(d)))) *
            ((1+3*gamma(d))+(1-gamma(d))*cosd(thetaplus)**2)
        )
    # end

    tau = [0]*8
    for d in range(8):
        # Rayleigh optical thickness (assume std pressure of 1013.25 mb)
        tau[d] = (
            0.008569*(cw(d)**-4) *
            (1 + 0.0113*(cw(d)**-2) + 0.00013*cw(d)**-4)
        )
    # end

    # Rayleigh calculation (Dash et al., 2012)
    ray_rad = [0]*8
    for d in range(8):
        # Assume standard pressure (1013.25 mb)
        # TODO: cell array -> list?
        ray_rad.append(  # Rayleigh Radience
            ((irr(1, d)/ESd)*1*tau(d)*Pr(d))/(4*pi*cosd(90-satel))
        )
        # ray_rad{1, 1}(d) = (  # Rayleigh Radience
        #     ((irr(1, d)/ESd)*1*tau(d)*Pr(d))/(4*pi*cosd(90-satel))
        # )
    # end

    # rrs constant calculation (Kerr et al. 2018 and Mobley 1994)
    G = float(1.56)  # constant (Kerr eq. 3)
    na = 1.00029  # Refractive index of air
    nw = 1.34  # Refractive index seawater
    # Incident angle for water-air from Snell's Law
    inc_ang2 = (asind(sind(90-satel)*nw/na))
    # Transmission angle for air-water incident light from Snell's Law
    trans_aw = (asind(sind(inc_ang)*na/nw))
    # Transmission angle for water-air incident light from Snell's Law
    trans_wa = 90-satel
    # Fresnel reflectance for air-water incident light (Mobley 1994)
    pf1 = (0.5*(
        (sind(inc_ang - trans_aw)/(sind(inc_ang + trans_aw)))**2 +
        (tand(inc_ang - trans_aw)/(tand(inc_ang + trans_aw)))**2
    ))
    pf2 = (0.5*(
        (sind(inc_ang2 - trans_wa)/(sind(inc_ang2 + trans_wa)))**2 +
        (tand(inc_ang2 - trans_wa)/(tand(inc_ang2 + trans_wa)))**2
    ))
    # rrs constant (~0.52) from Mobley 1994
    zeta = (float((1-pf1)*(1-pf2)/(nw**2)))
    # ==================================================================

    # Adjust file size: Input file (A) warped may contain more or fewer
    # columns/rows than original NITF file, and some may be corrupt.
    sz = [0]*3
    sz[1] = min(szA(1), szB(1))
    sz[2] = min(szA(2), szB(2))
    sz[3] = 8

    # === Assign NaN to no-data pixels and radiometrically calibrate and
    # convert to Rrs
    # Create empty matrix for Rrs output
    # TODO: use numpy.array([interven[-N:]])
    # numpy.array([szA(1), szA(2), 8])
    Rrs = float(zeros(szA(1), szA(2), 8))  # 8 bands x input size
    for j in range(sz(1)):
        # Assign NaN to pixels of no data
        # If a pixel contains data values other than "zero" or
        # "two thousand and forty seven" in any band, it is calibrated;
        # otherwise, it is considered "no-data" - this avoids a
        # problem created during the orthorectification process
        # wherein reprojecting the image may resample data
        for k in range(sz(2)):
            if (
                (A(j, k, 1)) != 0 and
                (A(j, k, 1)) != 2047 or (A(j, k, 2)) != 0 and
                (A(j, k, 2)) != 2047 or (A(j, k, 3)) != 0 and
                (A(j, k, 3)) != 2047 or (A(j, k, 4)) != 0 and
                (A(j, k, 4)) != 2047 or (A(j, k, 5)) != 0 and
                (A(j, k, 5)) != 2047 or (A(j, k, 6)) != 0 and
                (A(j, k, 6)) != 2047 or (A(j, k, 7)) != 0 and
                (A(j, k, 7)) != 2047 or (A(j, k, 8)) != 0 and
                (A(j, k, 8)) != 2047
            ):
                for d in range(8):
                    # Radiometrically calibrate and convert to Rrs
                    # (adapted from Radiometric Use of
                    # WorldView-2 Imagery(
                    Rrs[j, k, d] = float(
                        (
                            pi*(
                                (float(A(j, k, d))*kf(d, 1)/ebw(1, d)) -
                                ray_rad[d]
                                # ray_rad{1, 1}(1, d)
                            )*ESd**2
                        ) / (irr(1, d)*TZ*TV)
                    )
                # end
            else:
                Rrs[j, k, :] = OUTPUT_NaN
            # end
        # end
    # end

    # clear A

    # === Output reflectance image
    if Rrs_write == 1:
        Z = [loc_out, id, '_', loc, '_Rrs']
        geotiffwrite(Z, Rrs, R, CoordRefSysCode=coor_sys)
    # end

    if d_t > 0:
        # Run DT and/or rrs conversion; otherwise end

        # Calculate Kd (water column attenuation coefficient) from
        # Chuanmin Hu's Rrs_Kd_Model.xlsx sheet
        sunzen = 90.0 - sunel
        # c1-4 hard-coded, but v1 and v2 change with modified values of
        # aph(440), adg(440), bbp(440), Sdg, Y
        c1 = 0.005
        c2 = 4.18
        c3 = 0.52
        c4 = 10.8
        v2 = [  # bb (backscatter)
            0.023277868, 0.020883561, 0.018975346, 0.017605058,
            0.016771098, 0.015875283, 0.072734281, 0.068046578
        ]
        v1 = [  # at (total absorption) DEF CHL (0.1 0.1 0.01 0.015 0.5)
            0.22024, 0.142972, 0.099157, 0.286342, 0.443809, 1.491289,
            2.276262, 2.223947
        ]
        v1 = [  # Belize test
            0.069177389, 0.041529229, 0.061338823, 0.268134177,
            0.411563873, 1.489684614, 2.275464564, 2.22325881
        ]
        v2 = [
            0.023277868, 0.020883561, 0.018975346, 0.017605058,
            0.016771098, 0.015875283, 0.049505144, 0.046268752
        ]
        v1 = [  # at (total absorption) LOWER CHL
            0.006921, 0.013933, 0.051193, 0.264353, 0.409819, 1.489006,
            2.275336, 2.223238
        ]
        v2 = [  # bb (backscatter)
            0.023277868, 0.020883561, 0.018975346, 0.017605058,
            0.016771098, 0.015875283, 0.000869138, 0.00067143
        ]

        Kd = [0]*8
        for b in range(8):
            Kd[b] = float(
                (1+c1*sunzen)*v1(b)+c2*(1-c3*exp(-c4*v1(b)))*v2(b)
            )
        # end

#         Kd = [0.036 0.037 0.075 0.32 0.484 1.416]

        # === Setup for Deglint, Bathymetry, and Decision Tree
        # TODO: these counters not needed now that list.append is used(?)
        b = 1  # developed land counter?
        t = 1  # veg counter?
        u = 1  # water counter?
        y = 0
        v = 0
        sum_SD = []  # sand & developed
        num_pix = 0
        sum_veg = [0]
        sum_veg2 = []
        dead_veg = [0]
        sum_water_rrs = []
        sz_ar = sz(1)*sz(2)
        water = zeros(sz_ar, 9)
        c_val = []
        for j in range(sz(1)):
            for k in range(sz(2)):
                if isnan(Rrs[j, k, 1]) is False:
                    num_pix = num_pix + 1  # Count number of non-NaN pixels
                    # Record coastal band value for cloud mask prediction
                    c_val.append(Rrs[j, k, 1])
                    if (
                        (
                            (Rrs[j, k, 7] - Rrs[j, k, 2]) /
                            (Rrs[j, k, 7] + Rrs[j, k, 2])
                        ) < 0.65 and
                        Rrs[j, k, 5] > Rrs[j, k, 4] and
                        Rrs[j, k, 4] > Rrs[j, k, 3]
                    ):  # Sand & Developed
                        sum_SD.append(sum(Rrs[j, k, 6:8]))
                        b = b+1
                    # Identify vegetation (excluding grass)
                    elif (
                        (
                            (Rrs[j, k, 8] - Rrs[j, k, 5]) /
                            (Rrs[j, k, 8] + Rrs[j, k, 5])
                        ) > 0.6 and
                        Rrs[j, k, 7] > Rrs[j, k, 3]
                    ):
                        if (  # Shadow filter
                            (
                                (Rrs[j, k, 7] - Rrs[j, k, 2]) /
                                (Rrs[j, k, 7] + Rrs[j, k, 2])
                            ) > 0.20
                        ):
                            # Sum bands 3-5 for selected veg to distinguish
                            # wetland from upland
                            sum_veg.append(sum(Rrs[j, k, 3:5]))
                            sum_veg2.append(sum(Rrs[j, k, 7:8]))
                            # Compute difference of predicted B5 value from
                            # actual valute
                            dead_veg.append(
                                (
                                    ((Rrs[j, k, 7] - Rrs[j, k, 4])/3) +
                                    Rrs[j, k, 4]
                                ) - Rrs[j, k, 5]
                            )
                            t = t+1
                        # end
                    elif (  # Identify glint-free water
                        Rrs[j, k, 8] < 0.11 and
                        Rrs[j, k, 1] > 0 and
                        Rrs[j, k, 2] > 0 and
                        Rrs[j, k, 3] > 0 and
                        Rrs[j, k, 4] > 0 and
                        Rrs[j, k, 5] > 0 and
                        Rrs[j, k, 6] > 0 and
                        Rrs[j, k, 7] > 0 and
                        Rrs[j, k, 8] > 0
                    ):
                        water[u, 1:8] = float(Rrs[j, k, :])
                        water_rrs = rdivide(
                            Rrs[j, k, 1:6],
                            (zeta + G*Rrs[j, k, 1:6])
                        )
                        if (
                            water_rrs[4] > water_rrs[2] and
                            water_rrs[4] < 0.12 and
                            water_rrs[5] < water_rrs[3]
                        ):
                            sum_water_rrs.append(sum(water_rrs[3:5]))
                        # end
                        # WARN: u increments regardless sum_water_rrs
                        #       append? Is this intentional and what does
                        #       it mean?
                        u = u+1
                        # NDGI to identify glinted water pixels
                        # (some confusion w/ clouds)
                        if (
                            Rrs[j, k, 8] < Rrs[j, k, 7] and
                            Rrs[j, k, 6] < Rrs[j, k, 7] and
                            Rrs[j, k, 6] < Rrs[j, k, 5] and
                            Rrs[j, k, 4] < Rrs[j, k, 5] and
                            Rrs[j, k, 4] < Rrs[j, k, 3]
                        ):
                            v = v+1
                            # Mark array2<array1 glinted pixls
                            water[u, 9] = 2
                        elif(
                            Rrs[j, k, 8] > Rrs[j, k, 7] and
                            Rrs[j, k, 6] > Rrs[j, k, 7] and
                            Rrs[j, k, 6] > Rrs[j, k, 5] and
                            Rrs[j, k, 4] > Rrs[j, k, 5] and
                            Rrs[j, k, 4] > Rrs[j, k, 3]
                        ):
                            v = v+1
                            # Mark array2>array1 glinted pixls
                            water[u, 9] = 3
                        else:
                            # Mark records of glint-free water
                            water[u, 9] = 1
                        # end
                    elif(
                        Rrs[j, k, 8] < Rrs[j, k, 7] and
                        Rrs[j, k, 6] < Rrs[j, k, 7] and
                        Rrs[j, k, 6] < Rrs[j, k, 5] and
                        Rrs[j, k, 4] < Rrs[j, k, 5] and
                        Rrs[j, k, 4] < Rrs[j, k, 3]
                    ):
                        water[u, 1:8] = float(Rrs[j, k, :])
                        # Mark array2<array1 glinted pixels
                        water[u, 9] = 2
                        u = u+1
                        v = v+1
                    elif (
                        Rrs[j, k, 8] > Rrs[j, k, 7] and
                        Rrs[j, k, 6] > Rrs[j, k, 7] and
                        Rrs[j, k, 6] > Rrs[j, k, 5] and
                        Rrs[j, k, 4] > Rrs[j, k, 5] and
                        Rrs[j, k, 4] > Rrs[j, k, 3]
                    ):
                        # Mark array2>array1 glinted pixels
                        water[u, 9] = 3
                        water[u, 1:8] = float(Rrs[j, k, :])
                        u = u + 1
                        v = v + 1
                    # elif (
                    #     (Rrs(j,k,4)-Rrs(j,k,8)) /
                    #     (Rrs(j,k,4)+Rrs(j,k,8)) < 0.55
                    #     and Rrs(j,k,8) < 0.2
                    #     and (Rrs(j,k,7)-Rrs(j,k,2)) /
                    #       (Rrs(j,k,7)+Rrs(j,k,2)) < 0.1
                    #     and (Rrs(j,k,8)-Rrs(j,k,5)) /
                    #       (Rrs(j,k,8)+Rrs(j,k,5)) < 0.3
                    #     and Rrs(j,k,1) > 0
                    #     and Rrs(j,k,2) > 0
                    #     and Rrs(j,k,3) > 0
                    #     and Rrs(j,k,4) > 0
                    #     and Rrs(j,k,5) > 0
                    #     and Rrs(j,k,6) > 0
                    #     and Rrs(j,k,7) > 0
                    #     and Rrs(j,k,8) > 0
                    # ):
                    #
                    #     water(u, 1:8) = float(Rrs(j, k, :))
                    #     u = u + 1
                    #     v = v + 1
                    # end
                # end
            # end
        # end
        # Number of water pixels used to derive E_glint relationships
        n_water = u
        n_glinted = v  # Number of glinted water pixels

        water[water[:, 1] == 0] = numpy.nan
        water7 = water[:, 7]
        water8 = water[:, 8]
        # Positive minimum Band 7 value used for deglinting
        mnNIR1 = min(water7(water7 > 0))
        # Positive minimum Band 8 value used for deglinting
        mnNIR2 = min(water8(water8 > 0))

        # idx_gf = find(water[:, 9] == 1)  # Glint-free water

        if v > 0.25 * u:
            Update = 'Deglinting'
            id2 = 'deglinted'
            # idx_w1 = find(water(:, 9)==2) # Glinted water array1>array2
            # idx_w2 = find(water(:, 9)==3) # Glinted water array2>array1
            # water1 = [water(idx_gf, 1:8);water(idx_w1, 1:8)];
            # water2 = [water(idx_gf, 1:8);water(idx_w2, 1:8)];
            # Calculate linear fitting of all MS bands vs NIR1 & NIR2
            # for deglinting in DT (Hedley et al. 2005)
            E_glint = [0]*6
            for b in range(6):
                if b == 1 or b == 4 or b == 6:
                    # slope1 = water(:, b)\water(:, 8)
                    slope1 = mldivide(
                        water[:, b],
                        water[:, 8]
                    )
                else:
                    # slope1 = water(:, b)\water(:, 7)
                    slope1 = mldivide(
                        water[:, b],
                        water[:, 7]
                    )
                # end
            E_glint[b] = float(slope1)
            # end
            # E_glint  # = [0.8075 0.7356 0.8697 0.7236 0.9482 0.7902]
        else:
            Update = 'Glint-free'
            id2 = 'glintfree'
        # end

        # === Edge Detection
        img_sub = Rrs[:, :, 5]
        # TODO: align imtophat usage w/ docs here:
        # http://scikit-image.org/docs/dev/auto_examples/xx_applications/plot_morphology.html#white-tophat
        # and here:
        # http://scikit-image.org/docs/dev/auto_examples/xx_applications/plot_thresholding.html
        # IE:
        # img_sub = data.camera()
        # BWbin = img_as_ubyte(io.imread(png_path),as_gray=True))
        BWbin = imbinarize(img_sub)
        BW = imtophat(BWbin, square_strel(10))
#        BW1 = edge(BWtop, 'canny')
#        seDil = strel('square', 1)
#        BWdil = imdilate(BW1, seDil)
#        BW = imfill(BWdil, 'holes')
#
#        seDer = strel('', [5 5])
#        BWer = imerode(BW, seDer)

#         # === Depth scaling
#         water10(:, 1:2) = water(idx_gf, 2:3)
#         water10(:, 1:2) = rdivide(
#             water10(:, 1:2),
#             (zeta + G*water10(:, 1:2))
#         )
#         waterdp = rdivide(
#             (log(1000*(water10(:, 1))),
#             log(1000*(water10(:, 2))))
#         )
#         water_dp = waterdp(waterdp>0 & waterdp<2)
#         [N, X] = hist(water_dp)
#         med_dp = median(water_dp)
#         low = X(2) #avg_dp - 5*std(water_dp) #min(water_dp)
#         scale_dp = scale/(med_dp-low)
#
#         clear water10
        #         std_dp = std(water_dp)
#         low = avg_dp - 2*std_dp # Assumed represents 0 depth or min depth
#         high = avg_dp + std_dp

        # === Determine Rrs-infinite from glint-free water pixels
#         water_gf = water(idx_gf, 1:8)
# Sort all values in water by NIR2 column
# (assumes deepest water is darkest is NIR2)
#         dp_max_sort = sortrows(water_gf, 8, 'ascend')
#         # Use "deepest" 0.1# pixels
#         idx_dp = round(size(dp_max_sort, 1)*0.001)
#         dp_pct = dp_max_sort(1:idx_dp, :)
#         # Convert to subsurface rrs
#         dp_rrs = rdivide(
#             dp_pct(:, 1:8),
#             (zeta + G*dp_pct(:, 1:8))
#         )
#         # Mean and Median values too high
#         #median(dp_rrs(:, 1:8)) - 2*std(dp_rrs(:, 1:8))
#         rrs_inf = min(dp_rrs(:, 1:8))
#           # Derived from Rrs_Kd_Model.xlsx for Default values
# #         rrs_inf = [0.00512 0.00686 0.008898 0.002553 0.001506 0.000403]
# #         plot(rrs_inf)
        # === Calculate target class metrics
        avg_SD_sum = mean(sum_SD)
        stdev_SD_sum = std(sum_SD)
        avg_veg_sum = mean(sum_veg)
        avg_dead_veg = mean(dead_veg)
        avg_mang_sum = mean(sum_veg2)

        # exclude sum_water_rrs == 0 in avg calculations
        sum_water_rrs[sum_water_rrs == 0] = numpy.nan
        avg_water_sum = mean(sum_water_rrs)

        if cl_cov > 0:
            # Number of cloud pixels (rounded down to nearest integer)
            # based on metadata-reported percent cloud cover
            num_cld_pix = round(num_pix*cl_cov*0.01)
            # Sort all pixel blue-values in descending order. Cloud mask
            # threshold will be num_cld_pix'th highest value
            srt_c = list(c_val).sort(reverse=True)
            cld_mask = srt_c(num_cld_pix)  # Set cloud mask threshold
        else:
            cld_mask = max(c_val)+1
        # end

        Bathy = float(zeros(szA(1), szA(2)))  # Preallocate for Bathymetry
        Rrs_deglint = float(zeros(5, 1))  # Preallocate for deglinted Rrs
        # Preallocate water-column corrected Rrs
        Rrs_0 = float(zeros(5, 1))
        # Create empty matrix for classification output
        map = zeros(szA(1), szA(2), 'uint8')

    if d_t == 1:  # Execute Deglinting rrs and Bathymetry
        if v > u*0.25:
            # Deglint equation
            Rrs_deglint[1, 1] = (
                Rrs[j, k, 1] - (E_glint[1]*(Rrs[j, k, 8] - mnNIR2))
            )
            Rrs_deglint[2, 1] = (
                Rrs[j, k, 2] - (E_glint[2]*(Rrs[j, k, 7] - mnNIR1))
            )
            Rrs_deglint[3, 1] = (
                Rrs[j, k, 3] - (E_glint[3]*(Rrs[j, k, 7] - mnNIR1))
            )
            Rrs_deglint[4, 1] = (
                Rrs[j, k, 4] - (E_glint[4]*(Rrs[j, k, 8] - mnNIR2))
            )
            Rrs_deglint[5, 1] = (
                Rrs[j, k, 5] - (E_glint[5]*(Rrs[j, k, 7] - mnNIR1))
            )
            Rrs_deglint[6, 1] = (
                Rrs[j, k, 6] - (E_glint[6]*(Rrs[j, k, 8] - mnNIR2))
            )

            # Convert above-surface Rrs to below-surface rrs
            # (Kerr et al. 2018)
            # Was Rrs_0=
            Rrs[j, k, 1:6] = rdivide(
                Rrs_deglint[1:6],
                (zeta + G*Rrs_deglint[1:6])
            )

            # Relative depth estimate
            # Calculate relative depth
            # (Stumpf 2003 ratio transform scaled to 1-10)
            dp = (log(1000*Rrs_0(2))/log(1000*Rrs_0(3)))
            if dp > 0 and dp < 2:
                Bathy[j, k] = dp
            else:
                dp = 0
            # end
            # for d = 1:5
            #     # Calculate water-column corrected benthic reflectance
            #     # (Traganos 2017 & Maritorena 1994)
            #     Rrs(j, k, d) = (
            #        ((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1, d)*dp_sc)) +
            #        rrs_inf(d)
            #     )
            # end

        else:  # For glint-free/low-glint images
            # Convert above-surface Rrs to subsurface rrs
            # (Kerr et al. 2018, Lee et al. 1998)
            Rrs[j, k, 1:6] = rdivide(
                Rrs[j, k, 1:6],
                (zeta + G*Rrs[j, k, 1:6])
            )
            # Calculate relative depth (Stumpf 2003 ratio transform)
            dp = (log(1000*Rrs_0(2))/log(1000*Rrs_0(3)))
            if dp > 0 and dp < 2:
                Bathy[j, k] = dp
            else:
                dp = 0
            # end
        # end

    # Execute Deglinting rrs, Bathymetery, and Decision Tree
    elif d_t == 2:
        update = 'Running DT'
        for j in range(1, szA(1)):
            for k in range(1, szA(2)):
                if isnan(Rrs[j, k, 1]) == 0:
                    # === Mud, Developed and Sand
                    if (
                        (Rrs[j, k, 7] - Rrs[j, k, 2]) /
                        (Rrs[j, k, 7] + Rrs[j, k, 2]) < 0.60 and
                        Rrs[j, k, 5] > Rrs[j, k, 4] and
                        Rrs[j, k, 4] > Rrs[j, k, 3]
                    ):
                        if (
                            Rrs[j, k, 7] < Rrs[j, k, 2] and
                            Rrs[j, k, 8] > Rrs[j, k, 5]
                        ):
                            map[j, k] = 0  # Shadow
                        elif (  # Buildings & bright sand
                            (Rrs[j, k, 8] - Rrs[j, k, 5]) /
                            (Rrs[j, k, 8] + Rrs[j, k, 5]) < 0.01 and
                            Rrs[j, k, 8] > 0.05
                        ):
                            if BW(j, k) == 1:
                                map[j, k] = 11  # Developed
                            elif sum(Rrs[j, k, 6:8]) < avg_SD_sum:
                                map[j, k] = 22  # Mud (intertidal?)
                            else:
                                map[j, k] = 21  # Beach/sand/soil
                            # end
                        elif (
                            Rrs[j, k, 5] >
                            (
                                Rrs[j, k, 2] +
                                ((Rrs[j, k, 7]-Rrs[j, k, 2])/5)*2
                            )
                        ):
                            map[j, k] = 21  # Beach/sand/soil
                        elif (
                            Rrs[j, k, 5] < (
                                ((Rrs[j, k, 7] - Rrs[j, k, 2])/5)*3 +
                                Rrs[j, k, 2]
                            )*0.60 and Rrs[j, k, 7] > 0.2
                        ):
                            map[j, k] = 31  # Marsh grass
                        else:
                            map[j, k] = 22  # Mud
                        # end
                    elif (
                        Rrs[j, k, 2] > Rrs[j, k, 3] and
                        Rrs[j, k, 7] > Rrs[j, k, 3] and
                        Rrs[j, k, 2] < 0.1 and
                        (Rrs[j, k, 8] - Rrs[j, k, 5]) /
                        (Rrs[j, k, 8] + Rrs[j, k, 5]) < 0.20 or
                        Rrs[j, k, 8] > 0.05 and
                        Rrs[j, k, 7] > Rrs[j, k, 2] and
                        (Rrs[j, k, 8] - Rrs[j, k, 5]) /
                        (Rrs[j, k, 8] + Rrs[j, k, 5]) < 0.1
                    ):
                        if BW(j, k) == 1:
                            map[j, k] = 11  # Shadow/Developed
                        else:
                            map[j, k] = 22  # Mud
                        # end
                    # === Vegetation
                    elif (  # Vegetation pixels (NDVI)
                        (Rrs[j, k, 8] - Rrs[j, k, 5]) /
                        (Rrs[j, k, 8] + Rrs[j, k, 5]) > 0.20 and
                        Rrs[j, k, 7] > Rrs[j, k, 3]
                    ):
                        # Shadowed-vegetation filter
                        # (B7/B8 ratio excludes marsh, which tends
                        # to have very similar values here)
                        if (
                            Rrs[j, k, 7] > Rrs[j, k, 2] and
                            (
                                (Rrs[j, k, 7] - Rrs[j, k, 2]) /
                                (Rrs[j, k, 7] + Rrs[j, k, 2])
                            ) < 0.20 and
                            (Rrs[j, k, 7] - Rrs[j, k, 8]) /
                            (Rrs[j, k, 7] + Rrs[j, k, 8]) > 0.01
                        ):
                            map[j, k] = 0  # Shadow
                        elif sum(Rrs[j, k, 3:5]) < avg_veg_sum:
                            # Agriculture filter based on elevated Blue
                            # band values
                            if (
                                (Rrs[j, k, 2] - Rrs[j, k, 5]) /
                                (Rrs[j, k, 2] + Rrs[j, k, 5]) < 0.4
                            ):
                                if (
                                    Rrs[j, k, 7] > 0.12 and
                                    sum(Rrs[j, k, 7:8]) /
                                    sum(Rrs[j, k, 3:5]) > 2
                                ):
                                    map[j, k] = 33  # Forested Wetland
                                # Dead vegetation or Marsh
                                else:
                                    map[j, k] = 31
                                # end
                            else:
                                # Forested Upland
                                # (most likely agriculture)
                                map[j, k] = 32
                            # end
                        elif sum(Rrs[j, k, 7:8]) < avg_mang_sum:
                            # Agriculture filter based on elevated
                            # blue band values
                            if (
                                (
                                    (Rrs[j, k, 2] - Rrs[j, k, 5]) /
                                    (Rrs[j, k, 2] + Rrs[j, k, 5])
                                ) < 0.4
                            ):
                                if (
                                    Rrs[j, k, 7] > 0.12 and
                                    sum(Rrs[j, k, 7:8]) /
                                    sum(Rrs[j, k, 3:5]) > 2
                                ):
                                    map[j, k] = 33  # Forested Wetland
                                else:  # Marsh or Dead Vegetation
                                    map[j, k] = 31
                                # end
                            else:
                                # Forested Upland
                                # (most likely agriculture)
                                map[j, k] = 32
                            # end
                        elif (  # NDVI for high upland values
                            (Rrs[j, k, 8] - Rrs[j, k, 5]) /
                            (Rrs[j, k, 8] + Rrs[j, k, 5]) > 0.65
                        ):
                            map[j, k] = 32  # Upland Forest/Grass
                        elif (

                            Rrs[j, k, 5] > (
                                ((Rrs[j, k, 7] - Rrs[j, k, 2])/5)*3 +
                                Rrs[j, k, 2]
                                )*0.60 and Rrs[j, k, 7] < 0.2
                        ):
                            # Difference of B5 from predicted B5 by
                            # slope of B7:B4 to distinguish marsh
                            # (old: live vs dead trees/grass/marsh)
                            map[j, k] = 31  # Marsh grass
                        elif Rrs[j, k, 7] < 0.12:
                            map[j, k] = 30  # Dead vegetation
                        else:
                            map[j, k] = 32  # Upland Forest/Grass
                        # end
                    # === Water
                    elif (  # Identify all water (glinted & glint-free)
                        Rrs[j, k, 8] < 0.2 and Rrs[j, k, 8] > 0 or
                        Rrs[j, k, 8] < Rrs[j, k, 7] and
                        Rrs[j, k, 6] < Rrs[j, k, 7] and
                        Rrs[j, k, 6] < Rrs[j, k, 5] and
                        Rrs[j, k, 4] < Rrs[j, k, 5] and
                        Rrs[j, k, 4] < Rrs[j, k, 3] and
                        Rrs[j, k, 8] > 0 or
                        Rrs[j, k, 8] > Rrs[j, k, 7] and
                        Rrs[j, k, 6] > Rrs[j, k, 7] and
                        Rrs[j, k, 6] > Rrs[j, k, 5] and
                        Rrs[j, k, 4] > Rrs[j, k, 5] and
                        Rrs[j, k, 4] > Rrs[j, k, 3] and
                        Rrs[j, k, 8] > 0
                    ):
                        # map[j, k] = 5
                        if v > u*0.25:
                            # Deglint equation
                            Rrs_deglint[1, 1] = (
                                Rrs[j, k, 1] -
                                (E_glint[1]*(Rrs[j, k, 8] - mnNIR2))
                            )
                            Rrs_deglint[2, 1] = (
                                Rrs[j, k, 2] -
                                (E_glint[2]*(Rrs[j, k, 7] - mnNIR1))
                            )
                            Rrs_deglint[3, 1] = (
                                Rrs[j, k, 3] -
                                (E_glint[3]*(Rrs[j, k, 7] - mnNIR1))
                            )
                            Rrs_deglint[4, 1] = (
                                Rrs[j, k, 4] -
                                (E_glint[4]*(Rrs[j, k, 8] - mnNIR2))
                            )
                            Rrs_deglint[5, 1] = (
                                Rrs[j, k, 5] -
                                (E_glint[5]*(Rrs[j, k, 7] - mnNIR1))
                            )
                            Rrs_deglint[6, 1] = (
                                Rrs[j, k, 6] -
                                (E_glint[6]*(Rrs[j, k, 8] - mnNIR2))
                            )

                            # Convert above-surface Rrs to
                            # below-surface rrs (Kerr et al. 2018)
                            Rrs[j, k, 1:6] = rdivide(
                                Rrs_deglint[1:6],
                                # Was Rrs_0=
                                (zeta + G*Rrs_deglint[1:6])
                            )

                        # Relative depth estimate
                        # Calculate relative depth
                        # (Stumpf 2003 ratio transform scaled to 1-10)
                        dp = (
                            log(1000*Rrs_0(2))/log(1000*Rrs_0(3))
                        )
                        if dp > 0 and dp < 2:
                            Bathy[j, k] = dp
                        else:
                            dp = 0
                        # end
                        # dp_sc = (dp-low)*scale_dp

                        # for d = 1:5:
                        #     # Calculate water-column corrected
                        #     # benthic reflectance (Traganos 2017 &
                        #     # Maritorena 1994)
                        #     Rrs(j, k, d) = (
                        #         ((Rrs_0(d)-rrs_inf(d)) /
                        #         exp(-2*Kd(1, d)*dp_sc))+rrs_inf(d))
                        # end

                        # === DT
                        if Rrs[j, k, 6] < Rrs[j, k, 7]:
                            map[j, k] = 0  # Shadow
                        elif (
                            (Rrs[j, k, 3] - Rrs[j, k, 4]) /
                            (Rrs[j, k, 3] + Rrs[j, k, 4]) < 0.10
                            # (Rrs[j, k, 2] - Rrs[j, k, 4]) /
                            # (Rrs[j, k, 2]+Rrs[j, k, 4]) < 0
                        ):
                            if (
                                Rrs[j, k, 4] > Rrs[j, k, 3] or
                                Rrs[j, k, 5] > Rrs[j, k, 3]
                            ):
                                map[j, k] = 53  # Soft bottom
                            elif (  # NEW from 0.05
                                sum(Rrs[j, k, 3:5]) > avg_water_sum and
                                (Rrs[j, k, 5] - Rrs[j, k, 2]) /
                                (Rrs[j, k, 5] + Rrs[j, k, 2]) > 0.1
                            ):
                                map[j, k] = 52  # Soft bottom
                            # Separate seagrass from dark water NEW
                            elif (
                                Rrs[j, k, 4] > Rrs[j, k, 2] and
                                (Rrs[j, k, 3] - Rrs[j, k, 6]) /
                                (Rrs[j, k, 3] + Rrs[j, k, 6]) < 0.60
                            ):
                                # Separate seagrass from turbid water
                                # NEW
                                if (
                                    (Rrs[j, k, 3] - Rrs[j, k, 5]) /
                                    (Rrs[j, k, 3] + Rrs[j, k, 5]) > 0.1
                                ):
                                    map[j, k] = 54  # Seagrass
                                else:
                                    map[j, k] = 55  # Turbid water
                                # end
                            else:
                                map[j, k] = 51  # Deep water
                            # end
                        else:
                            map[j, k] = 51  # Deep water
                        # end
                    else:  # For glint-free/low-glint images
                        # Convert above-surface Rrs to subsurface rrs
                        # (Kerr et al. 2018,  Lee et al. 1998)
                        Rrs[j, k, 1:6] = rdivide(
                            Rrs[j, k, 1:6],
                            (zeta + G*Rrs[j, k, 1:6])
                        )
                        # Calculate relative depth
                        # (Stumpf 2003 ratio transform)
                        dp = (
                            log(1000*Rrs_0(2))/log(1000*Rrs_0(3))
                        )
                        if dp > 0 and dp < 2:
                            Bathy[j, k] = dp
                        else:
                            dp = 0
                        # end
                        # dp_sc = (dp-low)*scale_dp
                        # for d = 1:5
                        #     # Calculate water-column corrected
                        #     # benthic reflectance (Traganos 2017 &
                        #     # Maritorena 1994)
                        #     Rrs(j, k, d) = (
                        #         ((Rrs_0(d)-rrs_inf(d)) /
                        #         exp(-2*Kd(1, d)*dp_sc))+rrs_inf(d)
                        #     )
                        # end
                        # === DT
                        if Rrs[j, k, 6] < Rrs[j, k, 7]:
                            map[j, k] = 0  # Shadow
                        elif (
                            (Rrs[j, k, 3] - Rrs[j, k, 4]) /
                            (Rrs[j, k, 3] + Rrs[j, k, 4]) < 0.10
                            # (Rrs[j, k, 2] - Rrs[j, k, 4]) /
                            # (Rrs[j, k, 2]+Rrs[j, k, 4]) < 0
                        ):
                            if (
                                Rrs[j, k, 4] > Rrs[j, k, 3] or
                                Rrs[j, k, 5] > Rrs[j, k, 3]
                            ):
                                map[j, k] = 53  # Soft bottom
                            elif (
                                sum(Rrs[j, k, 3:5]) > avg_water_sum and
                                (Rrs[j, k, 5] - Rrs[j, k, 2]) /
                                (Rrs[j, k, 5] + Rrs[j, k, 2]) > 0.1
                            ):
                                map[j, k] = 52  # Soft bottom
                            elif (  # Separate seagrass from dark water
                                Rrs[j, k, 4] > Rrs[j, k, 2] and
                                (Rrs[j, k, 3] - Rrs[j, k, 6]) /
                                (Rrs[j, k, 3] + Rrs[j, k, 6]) < 0.60
                            ):
                                # Separate seagrass from turbid water
                                if (
                                    (Rrs[j, k, 3] - Rrs[j, k, 5]) /
                                    (Rrs[j, k, 3] + Rrs[j, k, 5]) >
                                    0.10
                                ):
                                    map[j, k] = 54  # Seagrass
                                else:
                                    map[j, k] = 55  # Turbid water
                                # end
                            else:
                                map[j, k] = 51  # Deep water
                            # end
                        else:
                            map[j, k] = 51  # Deep water
                        # end
                    # end  # if v>u
                # end  # If water/land
            # end  # If isnan
        # end  # k
        # if j == szA(1)/4
        #     update = 'DT 25# Complete'
        # end
        # if j == szA(1)/2
        #     update = 'DT 50# Complete'
        # end
        # if j == szA(1)/4*3
        #     update = 'DT 75# Complete'
        # end
    # end  # j

            # === Classes:
            # 1 = Developed
            # 2 = Vegetation
            # 3 = Soil/sand/beach
            # 41 = Deep water
            # 42 = Benthic Sand
            # 43 = Benthic Seagrass
            # 44 = Benthic Coral
            # 45 = Benthic patch coral

            # === DT Filter
            if filter > 0:
                dt_filt = DT_Filter(map, filter, sz(1), sz(2))
                AA = [
                    loc_out, id, '_', loc, '_Map_filt_', str(filter),
                    '_benthicnew'
                ]
                geotiffwrite(
                    AA, dt_filt, R, CoordRefSysCode=coor_sys
                )
            else:
                Z1 = [loc_out, id, '_', loc, '_Map_benthicnew']
                geotiffwrite(Z1, map, R, CoordRefSysCode=coor_sys)
            # end

            # === Output images
            # Z = [loc_out, id, '_', loc, '_Bathy1']
            # geotiffwrite(Z, Bathy, R(1, 1), CoordRefSysCode=coor_sys)
            Z2 = [loc_out, id, '_', loc, '_rrssub']  # last=52
            geotiffwrite(Z2, Rrs, R, CoordRefSysCode=coor_sys)
        # end  # If dt = 1
    # end  # If dt>0
# end


def main(
    input_tiff, input_xml, output_dir, roi_name
):
    crd_sys = "EPSG:4326"
    # === parse arguments:
    if crd_sys == "EPSG:4326":
        coor_sys = 4326  # Change coordinate system code here
    else:
        raise ValueError("unknown coord sys: '{}'".format(crd_sys))

    # sgwid =  num2str(sgw)

    process_file(
        input_tiff, input_xml, output_dir, roi_name, coor_sys, 0, 0
    )


def process_files_in_dir(
    loc_in=DATA_DIR + '/Ortho/',
    _id=0,  # NOTE: unused?
    met_in=DATA_DIR + '/Raw/',
    coor_sys=4326,  # coordinate system code
    d_t=1,  # 0=End after Rrs conversion; 1=rrs, bathy, DT; 2 = rrs, bathy & DT
    sgw=0,  # Sunglint moving-window box = sgw*2 +1 (i.e. 2 = 5x5 box)
    filter=3,  # 0=None, 1=3x3, 3=7x7, 5=11x11
    _stat=3,  # NOTE: unused?
    loc='RB',  # Typically the estuary acronym,
    id_number=0,  # (prev SLURM_ARRAY_TASK_ID) TODO: rm this?
    loc_out=DATA_DIR + '/Output/'
):
    """
    Process a lot of files in directories.

    !!! DEFUNCT
    """
    raise NotImplementedError("This function not yet fully ported to python.")

    # === get list of all product files in directory
    matfiles = glob(path.join(
        'Matt', 'USF', 'Other', 'NERRS_Mapping', 'Processing', 'Ortho', '*.tif'
    ))
    # TODO: Revise this to find both all-caps and all lower-case extensions
    # matfiles2 = glob(path.join(
    #     'Matt', 'USF', 'Other', 'NERRS_Mapping', 'Processing', 'Raw', '*.xml'
    # ))

    # loc_in = ['/home1/mmccarthy/Matt/USF/Other/Seagrass/test/']
    # met_in = ['/home1/mmccarthy/Matt/USF/Other/Seagrass/test/']
    # loc_out = ['/home1/mmccarthy/Matt/USF/Other/Seagrass/test/Rrs/']
    # matfiles = path.join'Matt', 'USF', 'Other', 'Seagrass', 'test', '*.tif'))
    # matfiles2 = path.join'Matt', 'USF', 'Other', 'Seagrass', 'test','*.xml'))

    sz_files = len(matfiles)

    for z in range(sz_files):  # for each file
        process_file()


if __name__ == "__main__":
    main(*sys.argv[1:])
