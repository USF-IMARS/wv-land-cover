import numpy
from numpy import zeros
from numpy import mean
from numpy import isnan
# from memory_profiler import profile

from wv_classify.matlab_fns import rdivide


# @profile
def run_rrs(sz, Rrs, zeta, G):
    # Run DT and/or rrs conversion;
    print('Running DT and/or rrs conversion...')

    # Setup for Deglint, Bathymetry, and Decision Tree
    b = 1  # developed land counter?
    t = 1  # veg counter?
    u = 0  # water counter?
    # y = 0
    v = 0
    sum_SD = []  # sand & developed
    num_pix = 0  # count of good pixels
    nan_pix = 0  # count of nan pixels
    sum_veg = [0]
    sum_veg2 = []
    dead_veg = [0]
    sum_water_rrs = []
    sz_ar = sz[0]*sz[1]
    water = zeros((sz_ar, 9))
    # c_val = []
    for j in range(sz[0]):
        print("\trow {}".format(j), end="\r")
        for k in range(sz[1]):
            if isnan(Rrs[j, k, 0]):
                nan_pix += 1
            else:
                num_pix = num_pix + 1  # Count number of non-NaN pixels
                # Record coastal band value for cloud mask prediction
                # c_val.append(Rrs[j, k, 0])
                if (
                    (
                        (Rrs[j, k, 6] - Rrs[j, k, 1]) /
                        (Rrs[j, k, 6] + Rrs[j, k, 1])
                    ) < 0.65 and
                    Rrs[j, k, 4] > Rrs[j, k, 3] and
                    Rrs[j, k, 3] > Rrs[j, k, 2]
                ):  # Sand & Developed
                    sum_SD.append(sum(Rrs[j, k, 5:7]))
                    b = b+1
                # Identify vegetation (excluding grass)
                elif (
                    (
                        (Rrs[j, k, 7] - Rrs[j, k, 4]) /
                        (Rrs[j, k, 7] + Rrs[j, k, 4])
                    ) > 0.6 and
                    Rrs[j, k, 6] > Rrs[j, k, 2]
                ):
                    if (  # Shadow filter
                        (
                            (Rrs[j, k, 6] - Rrs[j, k, 1]) /
                            (Rrs[j, k, 6] + Rrs[j, k, 1])
                        ) > 0.20
                    ):
                        # Sum bands 3-5 for selected veg to distinguish
                        # wetland from upland
                        sum_veg.append(sum(Rrs[j, k, 2:4]))
                        sum_veg2.append(sum(Rrs[j, k, 6:7]))
                        # Compute difference of predicted B5 value from
                        # actual valute
                        dead_veg.append(
                            (
                                ((Rrs[j, k, 6] - Rrs[j, k, 3])/3) +
                                Rrs[j, k, 3]
                            ) - Rrs[j, k, 4]
                        )
                        t = t+1
                    # end
                elif (  # Identify glint-free water
                    Rrs[j, k, 7] < 0.11 and
                    Rrs[j, k, 0] > 0 and
                    Rrs[j, k, 1] > 0 and
                    Rrs[j, k, 2] > 0 and
                    Rrs[j, k, 3] > 0 and
                    Rrs[j, k, 4] > 0 and
                    Rrs[j, k, 5] > 0 and
                    Rrs[j, k, 6] > 0 and
                    Rrs[j, k, 7] > 0
                ):
                    water[u, 0:8] = Rrs[j, k, :]
                    water_rrs = rdivide(
                        Rrs[j, k, 0:5],
                        (zeta + G*Rrs[j, k, 0:5])
                    )
                    if (
                        water_rrs[3] > water_rrs[1] and
                        water_rrs[3] < 0.12 and
                        water_rrs[4] < water_rrs[2]
                    ):
                        sum_water_rrs.append(sum(water_rrs[2:4]))
                    # end
                    # WARN: u increments regardless sum_water_rrs
                    #       append? Is this intentional and what does
                    #       it mean?
                    u = u+1
                    # NDGI to identify glinted water pixels
                    # (some confusion w/ clouds)
                    if (
                        Rrs[j, k, 7] < Rrs[j, k, 6] and
                        Rrs[j, k, 5] < Rrs[j, k, 6] and
                        Rrs[j, k, 5] < Rrs[j, k, 4] and
                        Rrs[j, k, 3] < Rrs[j, k, 4] and
                        Rrs[j, k, 3] < Rrs[j, k, 2]
                    ):
                        v = v+1
                        # Mark array2<array1 glinted pixls
                        water[u, 8] = 2
                    elif(
                        Rrs[j, k, 7] > Rrs[j, k, 6] and
                        Rrs[j, k, 5] > Rrs[j, k, 6] and
                        Rrs[j, k, 5] > Rrs[j, k, 4] and
                        Rrs[j, k, 3] > Rrs[j, k, 4] and
                        Rrs[j, k, 3] > Rrs[j, k, 2]
                    ):
                        v = v+1
                        # Mark array2>array1 glinted pixls
                        water[u, 8] = 3
                    else:
                        # Mark records of glint-free water
                        water[u, 8] = 1
                    # end
                elif(
                    Rrs[j, k, 7] < Rrs[j, k, 6] and
                    Rrs[j, k, 5] < Rrs[j, k, 6] and
                    Rrs[j, k, 5] < Rrs[j, k, 4] and
                    Rrs[j, k, 3] < Rrs[j, k, 4] and
                    Rrs[j, k, 3] < Rrs[j, k, 2]
                ):
                    water[u, 0:8] = Rrs[j, k, :]
                    # Mark array2<array1 glinted pixels
                    water[u, 8] = 2
                    u = u+1
                    v = v+1
                elif (
                    Rrs[j, k, 7] > Rrs[j, k, 6] and
                    Rrs[j, k, 5] > Rrs[j, k, 6] and
                    Rrs[j, k, 5] > Rrs[j, k, 4] and
                    Rrs[j, k, 3] > Rrs[j, k, 4] and
                    Rrs[j, k, 3] > Rrs[j, k, 2]
                ):
                    # Mark array2>array1 glinted pixels
                    water[u, 8] = 3
                    water[u, 0:8] = Rrs[j, k, :]
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
                #     water(u, 1:8) = Rrs(j, k, :)
                #     u = u + 1
                #     v = v + 1
    # Number of water pixels used to derive E_glint relationships
    n_water = u
    n_glinted = v  # Number of glinted water pixels

    print("% good pixels by nan-count: {:05.2}".format(
        nan_pix/(num_pix+nan_pix)
    ))

    print("n_water", n_water)
    print("n_glinted", n_glinted)

    # if band_0 == 0, remove row from water
    # ```matlab
    #   idx = find(water(:,1) == 0);
    #   water(idx,:) = [];
    # ```
    water_len = len(water)
    water = water[water[:, 0] != 0]
    print("{} px removed with band 0 == 0...".format(water_len - len(water)))

    water_len = len(water)
    water = water[water[:, 6] > 0]
    print("{} px removed w/ band 6 < 0".format(water_len - len(water)))

    water_len = len(water)
    water = water[water[:, 7] > 0]
    print("{} px removed w/ band 7 < 0".format(water_len - len(water)))

    # idx_gf = find(water[:, 9] == 1)  # Glint-free water
    water_len = len(water)
    print("{} px remain".format(water_len))
    E_glint_slope = [0]*6
    E_glint_y_int = [0]*6
    if v > 0.25 * u:
        print("Deglinting")
        # idx_w1 = find(water(:, 9)==2) # Glinted water array1>array2
        # idx_w2 = find(water(:, 9)==3) # Glinted water array2>array1
        # water1 = [water(idx_gf, 1:8);water(idx_w1, 1:8)];
        # water2 = [water(idx_gf, 1:8);water(idx_w2, 1:8)];
        # === Calculate linear fitting of all MS bands vs NIR1 & NIR2
        # for deglinting in DT (Hedley et al. 2005)
        for b in range(6):
            if b == 0 or b == 3 or b == 5:
                # slope1 = water(:, b)\water(:, 7)
                correction_ind = 7
            else:
                assert b == 1 or b == 2 or b == 4
                correction_ind = 6
                # slope1 = water(:, b)\water(:, 6)
            # end
            E_glint_slope[b], E_glint_y_int[b] = numpy.linalg.lstsq(
                numpy.vstack([water[:, b], numpy.ones(len(water))]).T,
                water[:, correction_ind]
            )[0]

        # end
        # E_glint  # = [0.8075 0.7356 0.8697 0.7236 0.9482 0.7902]
        print("least-squares glint correction:\n\tslope:{}\n\ty-int:{}".format(
            E_glint_slope, E_glint_y_int
        ))
    else:
        print("Glint-free")
    # end

    # === Edge Detection
    # img_sub = Rrs[:, :, 5]
    # img_sub = img_sub[numpy.logical_not(numpy.isnan(img_sub))]^M
    # TODO: align imtophat usage w/ docs here:
    # http://scikit-image.org/docs/dev/auto_examples/xx_applications/plot_morphology.html#white-tophat
    # and here:
    # http://scikit-image.org/docs/dev/auto_examples/xx_applications/plot_thresholding.html
    # IE:
    # img_sub = data.camera()
    # BWbin = img_as_ubyte(io.imread(png_path),as_gray=True))
    # BWbin = imbinarize(img_sub)
    # BW = imtophat(BWbin, square_strel(10))
    BW = zeros((sz[0], sz[1]))
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
    print("Calculating target class metrics...")
    avg_SD_sum = mean(sum_SD)
    # stdev_SD_sum = std(sum_SD)
    avg_veg_sum = mean(sum_veg)
    # avg_dead_veg = mean(dead_veg)
    avg_mang_sum = mean(sum_veg2)

    # exclude sum_water_rrs == 0 in avg calculations
    sum_water_rrs = list(filter((0).__ne__, sum_water_rrs))
    avg_water_sum = mean(sum_water_rrs)

    if numpy.isnan(avg_water_sum):
        avg_water_sum = [0]
    # if cl_cov > 0:
    #     # Number of cloud pixels (rounded down to nearest integer)
    #     # based on metadata-reported percent cloud cover
    #     num_cld_pix = round(num_pix*cl_cov*0.01)
    #     # Sort all pixel blue-values in descending order. Cloud mask
    #     # threshold will be num_cld_pix'th highest value
    #     srt_c = list(c_val).sort(reverse=True)
    #     cld_mask = srt_c(num_cld_pix)  # Set cloud mask threshold
    # else:
    # cld_mask = max(c_val)+1
    # end

    return (
        v, u, E_glint_slope, E_glint_y_int, BW,
        avg_SD_sum, avg_veg_sum, avg_mang_sum, avg_water_sum
    )
