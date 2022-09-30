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

import numpy
from numpy import isnan
from numpy import zeros
# from numpy import std

# dep packages:
# from skimage.morphology import square as square_strel
# from skimage.morphology import white_tophat as imtophat
# from skimage.filters import threshold_otsu as imbinarize

# local imports:
# from DT_Filter import DT_Filter
from wv_classify.matlab_fns import geotiffread
from wv_classify.matlab_fns import geotiffwrite
from wv_classify.matlab_fns import cosd
from wv_classify.matlab_fns import sind
from wv_classify.matlab_fns import tand
from wv_classify.matlab_fns import acosd
from wv_classify.matlab_fns import asind
from wv_classify.matlab_fns import rdivide
from wv_classify.read_wv_xml import read_wv_xml
from wv_classify.run_rrs import run_rrs
from wv_classify.stumpf_relative_depth import stumpf_relative_depth

OUTPUT_NaN = numpy.nan
BASE_DATATYPE = numpy.float32
# dst_ds.GetRasterBand(1).SetNoDataValue(OUTPUT_NaN)
# === Assign constants for all images
# Effective Bandwidth per band
# (nm converted to um units; from IMD metadata files)
ebw = [0.0473, 0.0543, 0.0630, 0.0374, 0.0574, 0.0393, 0.0989, 0.0996]

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
gamma = [0.0150, 0.0147, 0.0144, 0.0141, 0.0141, 0.0141, 0.0138, 0.0138]


def process_file(
    X,  # MS Tiff input image path
    Z,  # XML met input file path
    loc_out,  # output directory
    loc,  # RoI identifier string
    coor_sys=4326,  # coordinate system code
    d_t=2,  # 0=End after Rrs conversion; 1=rrs, bathy, DT; 2 = rrs, bathy & DT
    Rrs_write=1,  # 1=write Rrs geotiff; 0=do not write
):
    """
    process a single set of files
    """
    if d_t == 1:  # this is here to catch it quickly
        raise NotImplementedError("rrs output only not yet supported")

    if not loc_out.endswith("/"):
        loc_out += "/"

    fname = path.basename(X)
    id = fname[0:18]

    A, R = geotiffread(X, numpy_dtype=BASE_DATATYPE)
    print("\tinput size: {}".format(A.shape))
    szA = [A.shape[0], A.shape[1], A.shape[2]]

    (
        szB, aqmonth, aqyear, aqhour, aqminute, aqsecond, sunaz, sunel,
        satel, sensaz, aqday, satview, kf, cl_cov
    ) = read_wv_xml(Z)

    szB[2] = 8

    print(" === calculating coefficients...")
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
    UT = aqhour + (aqminute/60) + (aqsecond/3600)
    B1 = int(year/100)
    B2 = 2-B1+int(B1/4)
    # Julian date
    JD = (
        int(365.25*(year+4716)) + int(30.6001*(month+1)) + aqday +
        UT/24.0 + B2 - 1524.5
    )
    D = JD - 2451545.0
    degs = float(357.529 + 0.98560028*D)  # Degrees
    # Earth-Sun distance at given date
    ESd = 1.00014 - 0.01671*cosd(degs) - 0.00014*cosd(2*degs)
    # (should be between 0.983 and 1.017)
    assert 0.983 < ESd and ESd < 1.017
    inc_ang = 90 - sunel
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
            (3/(4*(1+2*gamma[d]))) *
            ((1+3*gamma[d])+(1-gamma[d])*cosd(thetaplus)**2)
        )
    # end

    tau = [0]*8
    for d in range(8):
        # Rayleigh optical thickness
        # (Hansen and Travis); Dash et al. 2012 eq 7
        # P_0 = 1013.25
        # rayleigh_optical_thickness = (
        #     (P / P_O) * 0.008569 * wavelength**-4 *
        #     (1 + 0.0113*wavelength**-2 + 0.00013*wavelength**-4)
        # )
        # assuming std pressure of 1013.25 mb (P == P_0)
        # rayleigh_optical_thickness = (
        #     1 * 0.008569 * wavelength**-4 *
        #     (1 + 0.0113*wavelength**-2 + 0.00013*wavelength**-4)
        # )
        tau[d] = (
            1 * 0.008569*(cw[d]**-4) *
            (1 + 0.0113*(cw[d]**-2) + 0.00013*cw[d]**-4)
        )

    # end

    # Rayleigh calculation (aerosol path radiance)
    # (Dash et al., 2012) eq 16
    w_0 = 1  # single_scattering_albedo
    ray_rad = [0]*8
    for d in range(8):
        ray_rad[d] = (
            ((irr[d] / ESd) * w_0 * tau[d] * Pr[d]) /
            (4 * pi * cosd(90-satel))
        )

    # rrs constant calculation (Kerr et al. 2018 and Mobley 1994)
    G = 1.56  # constant (Kerr eq. 3)
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
    sz = [0]*2
    sz[0] = min(szA[0], szB[0])
    sz[1] = min(szA[1], szB[1])
    n_bands = 8

    print("\tszA: {}".format(szA))
    print("\tszB: {}".format(szB))
    print("\tsz : {}".format(sz))

    # TODO: this diagnostic could be made pretting using
    #    https://pypi.org/project/tabulate/
    print("ESd:{}\tTZ:{}\tTV:{}".format(ESd, TZ, TV))
    print("irr\t", irr)
    print("tau\t", tau)
    print("Pr \t", Pr)
    print("rrd\t", ray_rad)
    print("kf \t", kf)
    print("ebw\t", ebw)
    print("gamma\t", gamma)
    print("thetaplus\t", thetaplus)
    # === Radiometrically calibrate and convert to Rrs
    # === optimze calculation by pre-computing coefficients for each band
    # (A * KF / - RAY_RAD) * pi * ESd**2 / ( IRR * tz * tv)
    # (A * KF / - RAY_RAD) * PI_ESD_etc
    # (A * C1       - C2     ) where
    #   C1 = (KF / EBW)*pi*ESd**2 / (IRR*tz*tv)
    #   C2 = RAY_RAD   *pi*ESd**2 / (IRR*tz*tv)
    C1 = numpy.array(
        [
            (pi * ESd**2 * kf[d]) / (irr[d] * TZ * TV * ebw[d])
            for d in range(n_bands)
        ],
        BASE_DATATYPE
    )
    C2 = numpy.array(
        [
            (pi * ray_rad[d] * ESd**2) / (irr[d] * TZ * TV)
            for d in range(n_bands)
        ],
        BASE_DATATYPE
    )
    print("C1\t", C1)
    print("C2\t", C2)

    # === Assign NaN to no-data pixels
    print(" === clearing invalid pixels...")
    # the equation below does invalidity_mask = numpy.where(A in [0, 2047])
    # why abs & 1023.5? See https://stackoverflow.com/a/16343791/1483986
    print("calc mask")
    invalidity_mask = abs(A - 1023.5) == 1023.5  # or should this be <=
    # sum across band dimension, resulting in 2d boolean array of only x,y
    print("reduce ", invalidity_mask.shape)
    invalidity_mask = numpy.add.reduce(invalidity_mask, 2, dtype=bool)
    # get x,y indicies for all pixels who failed the test
    print("index ", invalidity_mask.shape)
    invalid_pixel_indicies = numpy.nonzero(invalidity_mask)
    n_pixels = A.shape[0] * A.shape[1]
    n_invalid = len(invalid_pixel_indicies[0])
    n_valid = n_pixels - n_invalid
    print("{} invalid pixels found at x,y:\n\t{}".format(
        n_invalid, invalid_pixel_indicies
    ))
    print("percent of good pixels in image: {:2.2f}%".format(
        100 * n_valid/n_pixels
    ))
    # set all bands in invalid pixels to NAN
    A[invalid_pixel_indicies] = [OUTPUT_NaN]*8

    print(" === calculating Rrs...")
    # === calculate all at once w/ numpy element-wise broadcasing:
    Rrs = A * C1 - C2

    # TODO: rm less efficient alternatives below:
    # === calculate all at once w/ list comprehension
    # Rrs = [[[
    #     C1[d] * A[j, k, d] - C2[d]
    #     for d in range(8)] for j in range(sz[0])] for k in range(sz[1])
    # ]  # or...
    # === Preallocate & calculate each pixel:
    # Rrs = zeros((sz[0], sz[1], n_bands), dtype=float)  # 8 bands x input size
    # good_pixels = invalid_pixels = 0
    # for j in range(sz[0]):
    #     if j % 50 == 0:  # print every Nth row number to entertain the user
    #         print(j, end='\t', flush=True)
    #     # Assign NaN to pixels of no data
    #     # If a pixel contains data values other than "zero" or
    #     # "two thousand and forty seven" in any band, it is calibrated;
    #     # otherwise, it is considered "no-data" - this avoids a
    #     # problem created during the orthorectification process
    #     # wherein reprojecting the image may resample data
    #     for k in range(sz[1]):
    #         # print(k, end='|')
    #         if any(band_val not in [0, 2047] for band_val in A[j, k, :]):
    #             # Radiometrically calibrate and convert to Rrs
    #             # (adapted from Radiometric Use of
    #             # WorldView-2 Imagery(
    #             Rrs[j, k, :] = [
    #                 A[j, k, d] * C1[d] - C2[d]
    #                 for d in range(n_bands)
    #             ]
    #             good_pixels += 1
    #         else:
    #             Rrs[j, k, :] = OUTPUT_NaN
    #             invalid_pixels += 1
    # print(
    #     "\n\tDone. {} px calculated. {} px skipped.".format(
    #         good_pixels, invalid_pixels
    #     )
    # )
    del A  # clear A
    print("\t  Rrs size: {}".format(Rrs.shape))
    # === Output reflectance image
    if Rrs_write == 1:
        Z = ''.join([loc_out, id, '_', loc, '_Rrs.tif'])
        geotiffwrite(Z, Rrs, R, CoordRefSysCode=coor_sys)
    # end

    if d_t > 0:
        (
            v, u, E_glint_slope, E_glint_y_int, BW,
            avg_SD_sum, avg_veg_sum, avg_mang_sum, avg_water_sum
        ) = run_rrs(sz, Rrs, zeta, G)

        # Preallocate for Bathymetry
        Bathy = numpy.zeros((szA[0], szA[1]), dtype=numpy.float)
        Rrs_deglint = zeros((5, 1))  # Preallocate for deglinted Rrs
        # Preallocate water-column corrected Rrs
        # Rrs_0 = zeros((5, 1))

    if d_t == 1:  # Execute Deglinting rrs and Bathymetry
        raise NotImplementedError("rrs output only not yet supported")
        # print('Executing Deglinting rrs and Bathymetry...')
        # if v > u*0.25:
        #     # Deglint equation
        #     Rrs_deglint[0, 0] = (
        #         Rrs[j, k, 0] - (E_glint[0]*(Rrs[j, k, 7] - mnNIR2))
        #     )
        #     Rrs_deglint[1, 1] = (
        #         Rrs[j, k, 1] - (E_glint[1]*(Rrs[j, k, 6] - mnNIR1))
        #     )
        #     Rrs_deglint[2, 1] = (
        #         Rrs[j, k, 2] - (E_glint[2]*(Rrs[j, k, 6] - mnNIR1))
        #     )
        #     Rrs_deglint[3, 1] = (
        #         Rrs[j, k, 3] - (E_glint[3]*(Rrs[j, k, 7] - mnNIR2))
        #     )
        #     Rrs_deglint[4, 1] = (
        #         Rrs[j, k, 4] - (E_glint[4]*(Rrs[j, k, 6] - mnNIR1))
        #     )
        #     Rrs_deglint[5, 1] = (
        #         Rrs[j, k, 5] - (E_glint[5]*(Rrs[j, k, 7] - mnNIR2))
        #     )
        #
        #     # Convert above-surface Rrs to below-surface rrs
        #     # (Kerr et al. 2018)
        #     # Was Rrs_0=
        #     Rrs[j, k, 0:5] = rdivide(
        #         Rrs_deglint[0:5],
        #         (zeta + G*Rrs_deglint[0:5])
        #     )
        #     dp = stumpf_relative_depth(Rrs_0(1), Rrs_0(2))
        #     Bathy[j, k] = dp
        #     # for d = 1:5
        #     #     # Calculate water-column corrected benthic reflectance
        #     #     # (Traganos 2017 & Maritorena 1994)
        #     #     Rrs(j, k, d) = (
        #     #        ((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1, d)*dp_sc)) +
        #     #        rrs_inf(d)
        #     #     )
        #     # end
        #
        # else:  # For glint-free/low-glint images
        #     # Convert above-surface Rrs to subsurface rrs
        #     # (Kerr et al. 2018, Lee et al. 1998)
        #     Rrs[j, k, 0:5] = rdivide(
        #         Rrs[j, k, 0:5],
        #         (zeta + G*Rrs[j, k, 0:5])
        #     )
        #     # Calculate relative depth (Stumpf 2003 ratio transform)
        #     dp = stumpf_relative_depth(Rrs_0(2), Rrs_0(3))
        #     Bathy[j, k] = dp
        # # end

    elif d_t == 2:
        # Execute Deglinting rrs, Bathymetery, and Decision Tree
        print('Executing Deglinting rrs, Bathymetery, and Decision Tree...')
        # Create empty matrix for classification output
        classif_map = numpy.zeros((szA[0], szA[1]), dtype='uint16')
        # map = zeros(szA[0], szA[1], 'uint8')

        # update = 'Running DT'
        for j in range(1, sz[0]):
            print("\trow {}".format(j), end="\r")
            for k in range(1, sz[1]):
                if isnan(Rrs[j, k, 0]) == 0:
                    # === Mud, Developed and Sand
                    if (
                        (Rrs[j, k, 6] - Rrs[j, k, 1]) /
                        (Rrs[j, k, 6] + Rrs[j, k, 1]) < 0.60 and
                        Rrs[j, k, 4] > Rrs[j, k, 3] and
                        Rrs[j, k, 3] > Rrs[j, k, 2]
                    ):
                        if (
                            Rrs[j, k, 6] < Rrs[j, k, 1] and
                            Rrs[j, k, 7] > Rrs[j, k, 4]
                        ):
                            classif_map[j, k] = 0  # Shadow
                        elif (  # Buildings & bright sand
                            (Rrs[j, k, 7] - Rrs[j, k, 4]) /
                            (Rrs[j, k, 7] + Rrs[j, k, 4]) < 0.01 and
                            Rrs[j, k, 7] > 0.05
                        ):
                            if BW[j, k] == 1:
                                classif_map[j, k] = 11  # Developed
                            elif sum(Rrs[j, k, 5:8]) < avg_SD_sum:
                                classif_map[j, k] = 22  # Mud (intertidal?)
                            else:
                                classif_map[j, k] = 21  # Beach/sand/soil
                            # end
                        elif (
                            Rrs[j, k, 4] >
                            (
                                Rrs[j, k, 1] +
                                ((Rrs[j, k, 6]-Rrs[j, k, 1])/5)*2
                            )
                        ):
                            classif_map[j, k] = 21  # Beach/sand/soil
                        elif (
                            Rrs[j, k, 4] < (
                                ((Rrs[j, k, 6] - Rrs[j, k, 1])/5)*3 +
                                Rrs[j, k, 1]
                            )*0.60 and Rrs[j, k, 6] > 0.2
                        ):
                            classif_map[j, k] = 31  # Marsh grass
                        else:
                            classif_map[j, k] = 22  # Mud
                        # end
                    elif (
                        Rrs[j, k, 1] > Rrs[j, k, 2] and
                        Rrs[j, k, 6] > Rrs[j, k, 2] and
                        Rrs[j, k, 1] < 0.1 and
                        (Rrs[j, k, 7] - Rrs[j, k, 4]) /
                        (Rrs[j, k, 7] + Rrs[j, k, 4]) < 0.20 or
                        Rrs[j, k, 7] > 0.05 and
                        Rrs[j, k, 6] > Rrs[j, k, 1] and
                        (Rrs[j, k, 7] - Rrs[j, k, 4]) /
                        (Rrs[j, k, 7] + Rrs[j, k, 4]) < 0.1
                    ):
                        if BW[j, k] == 1:
                            classif_map[j, k] = 11  # Shadow/Developed
                        else:
                            classif_map[j, k] = 22  # Mud
                        # end
                    # === Vegetation
                    elif (  # Vegetation pixels (NDVI)
                        (Rrs[j, k, 7] - Rrs[j, k, 4]) /
                        (Rrs[j, k, 7] + Rrs[j, k, 4]) > 0.20 and
                        Rrs[j, k, 6] > Rrs[j, k, 2]
                    ):
                        # Shadowed-vegetation filter
                        # (B7/B8 ratio excludes marsh, which tends
                        # to have very similar values here)
                        if (
                            Rrs[j, k, 6] > Rrs[j, k, 1] and
                            (
                                (Rrs[j, k, 6] - Rrs[j, k, 1]) /
                                (Rrs[j, k, 6] + Rrs[j, k, 1])
                            ) < 0.20 and
                            (Rrs[j, k, 6] - Rrs[j, k, 7]) /
                            (Rrs[j, k, 6] + Rrs[j, k, 7]) > 0.01
                        ):
                            classif_map[j, k] = 0  # Shadow
                        elif sum(Rrs[j, k, 2:4]) < avg_veg_sum:
                            # Agriculture filter based on elevated Blue
                            # band values
                            if (
                                (Rrs[j, k, 1] - Rrs[j, k, 4]) /
                                (Rrs[j, k, 1] + Rrs[j, k, 4]) < 0.4
                            ):
                                if (
                                    Rrs[j, k, 6] > 0.12 and
                                    sum(Rrs[j, k, 6:7]) /
                                    sum(Rrs[j, k, 2:4]) > 2
                                ):
                                    classif_map[j, k] = 33  # Forested Wetland
                                # Dead vegetation or Marsh
                                else:
                                    classif_map[j, k] = 31
                                # end
                            else:
                                # Forested Upland
                                # (most likely agriculture)
                                classif_map[j, k] = 32
                            # end
                        elif sum(Rrs[j, k, 6:7]) < avg_mang_sum:
                            # Agriculture filter based on elevated
                            # blue band values
                            if (
                                (
                                    (Rrs[j, k, 1] - Rrs[j, k, 4]) /
                                    (Rrs[j, k, 1] + Rrs[j, k, 4])
                                ) < 0.4
                            ):
                                if (
                                    Rrs[j, k, 6] > 0.12 and
                                    sum(Rrs[j, k, 6:7]) /
                                    sum(Rrs[j, k, 2:4]) > 2
                                ):
                                    classif_map[j, k] = 33  # Forested Wetland
                                else:  # Marsh or Dead Vegetation
                                    classif_map[j, k] = 31
                                # end
                            else:
                                # Forested Upland
                                # (most likely agriculture)
                                classif_map[j, k] = 32
                            # end
                        elif (  # NDVI for high upland values
                            (Rrs[j, k, 7] - Rrs[j, k, 4]) /
                            (Rrs[j, k, 7] + Rrs[j, k, 4]) > 0.65
                        ):
                            classif_map[j, k] = 32  # Upland Forest/Grass
                        elif (

                            Rrs[j, k, 4] > (
                                ((Rrs[j, k, 6] - Rrs[j, k, 1])/5)*3 +
                                Rrs[j, k, 1]
                                )*0.60 and Rrs[j, k, 6] < 0.2
                        ):
                            # Difference of B5 from predicted B5 by
                            # slope of B7:B4 to distinguish marsh
                            # (old: live vs dead trees/grass/marsh)
                            classif_map[j, k] = 31  # Marsh grass
                        elif Rrs[j, k, 6] < 0.12:
                            classif_map[j, k] = 30  # Dead vegetation
                        else:
                            classif_map[j, k] = 32  # Upland Forest/Grass
                        # end
                    # === Water
                    elif (  # Identify all water (glinted & glint-free)
                        Rrs[j, k, 7] < 0.2 and Rrs[j, k, 7] > 0 or
                        Rrs[j, k, 7] < Rrs[j, k, 6] and
                        Rrs[j, k, 5] < Rrs[j, k, 6] and
                        Rrs[j, k, 5] < Rrs[j, k, 4] and
                        Rrs[j, k, 3] < Rrs[j, k, 4] and
                        Rrs[j, k, 3] < Rrs[j, k, 2] and
                        Rrs[j, k, 7] > 0 or
                        Rrs[j, k, 7] > Rrs[j, k, 6] and
                        Rrs[j, k, 5] > Rrs[j, k, 6] and
                        Rrs[j, k, 5] > Rrs[j, k, 4] and
                        Rrs[j, k, 3] > Rrs[j, k, 4] and
                        Rrs[j, k, 3] > Rrs[j, k, 2] and
                        Rrs[j, k, 7] > 0
                    ):
                        # classif_map[j, k] = 5
                        if v > u*0.25:
                            print("deglint:")
                            # Deglint equation
                            print("\t{}x{} - {}".format(
                                float(E_glint_slope[0]),
                                float(Rrs[j, k, 7]),
                                float(E_glint_y_int[0]),
                            ))
                            Rrs_deglint[0, 0] = (
                                Rrs[j, k, 0] -
                                (
                                    float(E_glint_slope[0]) *
                                    float(Rrs[j, k, 7]) -
                                    float(E_glint_y_int[0])
                                )
                            )
                            print("\t{}x{} - {}".format(
                                float(E_glint_slope[1]),
                                float(Rrs[j, k, 6]),
                                float(E_glint_y_int[1]),
                            ))
                            Rrs_deglint[1, 0] = (
                                Rrs[j, k, 1] -
                                (
                                    float(E_glint_slope[1]) *
                                    float(Rrs[j, k, 6]) -
                                    float(E_glint_y_int[1])
                                )
                            )
                            print("\t{}x{} - {}".format(
                                float(E_glint_slope[2]),
                                float(Rrs[j, k, 6]),
                                float(E_glint_y_int[2]),
                            ))
                            Rrs_deglint[2, 0] = (
                                Rrs[j, k, 2] -
                                (
                                    float(E_glint_slope[2]) *
                                    float(Rrs[j, k, 6]) -
                                    float(E_glint_y_int[2])
                                )
                            )
                            print("\t{}x{} - {}".format(
                                float(E_glint_slope[3]),
                                float(Rrs[j, k, 7]),
                                float(E_glint_y_int[3]),
                            ))
                            Rrs_deglint[3, 0] = (
                                Rrs[j, k, 3] -
                                (
                                    float(E_glint_slope[3]) *
                                    float(Rrs[j, k, 7]) -
                                    float(E_glint_y_int[3])
                                )
                            )
                            print("\t{}x{} - {}".format(
                                float(E_glint_slope[4]),
                                float(Rrs[j, k, 6]),
                                float(E_glint_y_int[4]),
                            ))
                            Rrs_deglint[4, 0] = (
                                Rrs[j, k, 4] -
                                (
                                    float(E_glint_slope[4]) *
                                    float(Rrs[j, k, 6]) -
                                    float(E_glint_y_int[4])
                                )
                            )
                            print("\t{}x{} - {}".format(
                                float(E_glint_slope[5]),
                                float(Rrs[j, k, 7]),
                                float(E_glint_y_int[5]),
                            ))
                            Rrs_deglint[5, 0] = (
                                Rrs[j, k, 5] -
                                (
                                    float(E_glint_slope[5]) *
                                    float(Rrs[j, k, 7]) -
                                    float(E_glint_y_int[5])
                                )
                            )

                            # Convert above-surface Rrs to
                            # below-surface rrs (Kerr et al. 2018)
                            Rrs[j, k, 0:5] = rdivide(
                                Rrs_deglint[0:5],
                                # Was Rrs_0=
                                (zeta + G*Rrs_deglint[0:5])
                            )
                        dp = stumpf_relative_depth(Rrs[j, k, 1], Rrs[j, k, 2])
                        Bathy[j, k] = dp
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
                        if Rrs[j, k, 5] < Rrs[j, k, 6]:
                            classif_map[j, k] = 0  # Shadow
                        elif (
                            (Rrs[j, k, 2] - Rrs[j, k, 3]) /
                            (Rrs[j, k, 2] + Rrs[j, k, 3]) < 0.10
                            # (Rrs[j, k, 1] - Rrs[j, k, 3]) /
                            # (Rrs[j, k, 1]+Rrs[j, k, 3]) < 0
                        ):
                            if (
                                Rrs[j, k, 3] > Rrs[j, k, 2] or
                                Rrs[j, k, 4] > Rrs[j, k, 2]
                            ):
                                classif_map[j, k] = 53  # Soft bottom
                            elif (  # NEW from 0.05
                                sum(Rrs[j, k, 2:4]) > avg_water_sum and
                                (Rrs[j, k, 4] - Rrs[j, k, 1]) /
                                (Rrs[j, k, 4] + Rrs[j, k, 1]) > 0.1
                            ):
                                classif_map[j, k] = 52  # Soft bottom
                            # Separate seagrass from dark water NEW
                            elif (
                                Rrs[j, k, 3] > Rrs[j, k, 1] and
                                (Rrs[j, k, 2] - Rrs[j, k, 5]) /
                                (Rrs[j, k, 2] + Rrs[j, k, 5]) < 0.60
                            ):
                                # Separate seagrass from turbid water
                                # NEW
                                if (
                                    (Rrs[j, k, 2] - Rrs[j, k, 4]) /
                                    (Rrs[j, k, 2] + Rrs[j, k, 4]) > 0.1
                                ):
                                    classif_map[j, k] = 54  # Seagrass
                                else:
                                    classif_map[j, k] = 55  # Turbid water
                                # end
                            else:
                                classif_map[j, k] = 51  # Deep water
                            # end
                        else:
                            classif_map[j, k] = 51  # Deep water
                        # end
                    else:  # For glint-free/low-glint images
                        # Convert above-surface Rrs to subsurface rrs
                        # (Kerr et al. 2018,  Lee et al. 1998)
                        Rrs[j, k, 0:5] = rdivide(
                            Rrs[j, k, 0:5],
                            (zeta + G*Rrs[j, k, 0:5])
                        )
                        # Calculate relative depth
                        # (Stumpf 2003 ratio transform)
                        dp = stumpf_relative_depth(Rrs[j, k, 1], Rrs[j, k, 2])
                        Bathy[j, k] = dp
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
                        if Rrs[j, k, 5] < Rrs[j, k, 6]:
                            classif_map[j, k] = 0  # Shadow
                        elif (
                            (Rrs[j, k, 2] - Rrs[j, k, 3]) /
                            (Rrs[j, k, 2] + Rrs[j, k, 3]) < 0.10
                            # (Rrs[j, k, 1] - Rrs[j, k, 3]) /
                            # (Rrs[j, k, 1]+Rrs[j, k, 3]) < 0
                        ):
                            if (
                                Rrs[j, k, 3] > Rrs[j, k, 2] or
                                Rrs[j, k, 4] > Rrs[j, k, 2]
                            ):
                                classif_map[j, k] = 53  # Soft bottom
                            elif (
                                sum(Rrs[j, k, 2:4]) > avg_water_sum and
                                (Rrs[j, k, 4] - Rrs[j, k, 1]) /
                                (Rrs[j, k, 4] + Rrs[j, k, 1]) > 0.1
                            ):
                                classif_map[j, k] = 52  # Soft bottom
                            elif (  # Separate seagrass from dark water
                                Rrs[j, k, 3] > Rrs[j, k, 1] and
                                (Rrs[j, k, 2] - Rrs[j, k, 5]) /
                                (Rrs[j, k, 2] + Rrs[j, k, 5]) < 0.60
                            ):
                                # Separate seagrass from turbid water
                                if (
                                    (Rrs[j, k, 2] - Rrs[j, k, 4]) /
                                    (Rrs[j, k, 2] + Rrs[j, k, 4]) >
                                    0.10
                                ):
                                    classif_map[j, k] = 54  # Seagrass
                                else:
                                    classif_map[j, k] = 55  # Turbid water
                                # end
                            else:
                                classif_map[j, k] = 51  # Deep water
                            # end
                        else:
                            classif_map[j, k] = 51  # Deep water
                        # end
                    # end  # if v>u
                # end  # If water/land
            # end  # If isnan
            # end  # k
            # if j == szA[0]/4
            #     update = 'DT 25# Complete'
            # end
            # if j == szA[0]/2
            #     update = 'DT 50# Complete'
            # end
            # if j == szA[0]/4*3
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
        # if filter > 0:
        #     dt_filt = DT_Filter(map, filter, sz[0], sz[1])
        #     AA = ''.join([
        #         loc_out, id, '_', loc, '_Map_filt_', str(filter),
        #         '_benthicnew.tif'
        #     ])
        #     geotiffwrite(
        #         AA, dt_filt, R, CoordRefSysCode=coor_sys
        #     )
        # else:
        Z1 = ''.join([loc_out, id, '_', loc, '_Map_pytest.tif'])
        geotiffwrite(Z1, classif_map, R, CoordRefSysCode=coor_sys)
        # end

        # === Output images
        # Z = [loc_out, id, '_', loc, '_Bathy1']
        # geotiffwrite(Z, Bathy, R(1, 1), CoordRefSysCode=coor_sys)
        Z2 = ''.join([loc_out, id, '_', loc, '_rrssub.tif'])  # last=52
        geotiffwrite(Z2, Rrs, R, CoordRefSysCode=coor_sys)
    # end  # If dt == 2
# end


def main(
    input_tiff, input_xml, output_dir, roi_name, crd_sys, dt_out, rrs_out
):
    crd_sys = "EPSG:4326"
    # === parse arguments:
    if crd_sys == "EPSG:4326":
        coor_sys = 4326  # Change coordinate system code here
    else:
        raise ValueError("unknown coord sys: '{}'".format(crd_sys))

    # sgwid =  num2str(sgw)

    process_file(
        input_tiff, input_xml, output_dir, roi_name, coor_sys,
        int(dt_out), int(rrs_out)
    )

# TODO: update/rm this:
DATA_DIR = '/home1/mmccarthy/Matt/USF/Other/NERRS_Mapping/Processing'


def process_files_in_dir(
    loc_in=DATA_DIR + '/Ortho/',
    _id=0,  # NOTE: unused?
    met_in=DATA_DIR + '/Raw/',
    coor_sys=4326,  # coordinate system code
    d_t=2,  # 0=End after Rrs conversion; 1=rrs, bathy ; 2 = rrs, bathy & DT
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
