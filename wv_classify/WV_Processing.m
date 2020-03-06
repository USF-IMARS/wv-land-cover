%% WV2 Processing
% Loads TIFF WorldView-2 image files preprocessed through Polar Geospatial
% Laboratory python code, which orthorectifies and projects .NTF files and outputs as
% TIFF files
% Radiometrically calibrates digital count data
% Atmospherically corrects images by subtracting Rayleigh Path Radiance
% Converts image to surface reflectance by accounting for Earth-Sun
% distance, solar zenith angle, and average spectral irradiance
% Tests and optionally corrects for sunglint
% Corrects for water column attenuation
% Runs Decision Tree classification on each image
% Optionally smooths results through moving-window filter
% Outputs images as GEOTIFF files with geospatial information.

function dt_filt = WV_Processing(images,id,met,crd_sys,dt,filt,loc,idnumber,rrs_out,class_out);

tic
d_t = str2num(dt);
n = num2str(idnumber);
id
met
coor_sys = crd_sys; % Change coordinate system code here
filter = str2num(filt);
loc_out = rrs_out;

% Assign constants for all images
ebw1 = 0.001*[47.3 54.3 63.0 37.4 57.4 39.3 98.9 99.6]; % Effective Bandwidth per WV2 band (nm converted to um units; from IMD metadata files)
ebw2 = 0.001*[40.5 54.0 61.8 38.1 58.5 38.7 100.4 88.9]; % WV3
irr1 = [1758.2229 1974.2416 1856.4104 1738.4791 1559.4555 1342.0695 1069.7302 861.2866]; % Band-averaged Solar Spectral Irradiance (W/m2/um units)
irr2 = [1757.89 2004.61 1830.18 1712.07 1535.33 1348.08 1055.94 858.77]; % WV3 (from Radiometric Use of WorldView-3 Imagery, Thuiller 2003 column Table 3)
cw1 = [.4273 .4779 .5462 .6078 .6588 .7237 .8313 .9080]; % Center wavelength (used for Rayleigh correction; from Radiometric Use of WorldView-2 Imagery)
cw2 = [.4274 .4819 .5471 .6043 .6601 .7227 .8240 .9136]; % WV3
gamma = 0.01*[1.499 1.471 1.442 1.413 1.413 1.413 1.384 1.384]; % Factor used in Rayleigh Phase Function equation (Bucholtz 1995)

    [A, R] = geotiffread(images);
    szA = size(A);
     s = xml2struct(met);
%    save XMLtest.mat s
        % Extract calibration factors and acquisition time from metadata for each band
        if isfield(s,'IMD') == 1
             szB(1) = str2num(s.IMD.SOURCE_IMD.IMD.NUMROWS.Text); %#ok<*ST2NM>
             szB(2) = str2num(s.IMD.SOURCE_IMD.IMD.NUMCOLUMNS.Text);
        	 kf(1,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_C.ABSCALFACTOR.Text);
   	         kf(2,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_B.ABSCALFACTOR.Text);
	         kf(3,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_G.ABSCALFACTOR.Text);
	         kf(4,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_Y.ABSCALFACTOR.Text);
	         kf(5,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_R.ABSCALFACTOR.Text);
	         kf(6,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_RE.ABSCALFACTOR.Text);
	         kf(7,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_N.ABSCALFACTOR.Text);
	         kf(8,1) = str2num(s.IMD.SOURCE_IMD.IMD.BAND_N2.ABSCALFACTOR.Text);
	         aqyear = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.FIRSTLINETIME.Text(12:15)); % Extract Acquisition Time from metadata
	         aqmonth = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.FIRSTLINETIME.Text(17:18)); % Extract Acquisition Time from metadata
    	   	 aqday = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.FIRSTLINETIME.Text(20:21)); % Extract Acquisition Time from metadata
           	 aqhour = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.FIRSTLINETIME.Text(23:24)); % Extract Acquisition Time from metadata
             aqminute = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.FIRSTLINETIME.Text(26:27)); % Extract Acquisition Time from metadata
	    	 aqsecond = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.FIRSTLINETIME.Text(29:37)); % Extract Acquisition Time from metadata
	    	 sunel = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.MEANSUNEL.Text); % Extract Mean Sun Elevation angle from metadata
             satview = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.MEANOFFNADIRVIEWANGLE.Text); % Extract Mean Off Nadir View angle from metadata
	         sunaz = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.MEANSUNAZ.Text);
             sensaz = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.MEANSATAZ.Text);
             satel = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.MEANSATEL.Text);
             cl_cov = str2num(s.IMD.SOURCE_IMD.IMD.IMAGE.CLOUDCOVER.Text);

	elseif isfield(s,'isd') == 1
		 szB(1) = str2num(s.isd.IMD.NUMROWS.Text);
                 szB(2) = str2num(s.isd.IMD.NUMCOLUMNS.Text);
        	 kf(1,1) = str2num(s.isd.IMD.BAND_C.ABSCALFACTOR.Text);
   	         kf(2,1) = str2num(s.isd.IMD.BAND_B.ABSCALFACTOR.Text);
	         kf(3,1) = str2num(s.isd.IMD.BAND_G.ABSCALFACTOR.Text);
	         kf(4,1) = str2num(s.isd.IMD.BAND_Y.ABSCALFACTOR.Text);
	         kf(5,1) = str2num(s.isd.IMD.BAND_R.ABSCALFACTOR.Text);
	         kf(6,1) = str2num(s.isd.IMD.BAND_RE.ABSCALFACTOR.Text);
	         kf(7,1) = str2num(s.isd.IMD.BAND_N.ABSCALFACTOR.Text);
	         kf(8,1) = str2num(s.isd.IMD.BAND_N2.ABSCALFACTOR.Text);
	         aqyear = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(1:4)); % Extract Acquisition Time from metadata
	         aqmonth = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(6:7)); % Extract Acquisition Time from metadata
    	   	 aqday = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(9:10)); % Extract Acquisition Time from metadata
           	 aqhour = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(12:13)); % Extract Acquisition Time from metadata
             aqminute = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(15:16)); % Extract Acquisition Time from metadata
	    	 aqsecond = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(18:26)); % Extract Acquisition Time from metadata
	    	 sunel = str2num(s.isd.IMD.IMAGE.MEANSUNEL.Text); % Extract Mean Sun Elevation angle from metadata
             satview = str2num(s.isd.IMD.IMAGE.MEANOFFNADIRVIEWANGLE.Text); % Extract Mean Off Nadir View angle from metadata
	         sunaz = str2num(s.isd.IMD.IMAGE.MEANSUNAZ.Text);
             sensaz = str2num(s.isd.IMD.IMAGE.MEANSATAZ.Text);
             satel = str2num(s.isd.IMD.IMAGE.MEANSATEL.Text);
             cl_cov = str2num(s.isd.IMD.IMAGE.CLOUDCOVER.Text);
	else
            c = struct2cell(s.Children(2).Children(:));
	    idx{1} = strfind(c(1,:),'NUMROWS');
            idx{2} = strfind(c(1,:),'NUMCOLUMNS');
            idx{3} = strfind(c(1,:),'BAND_C');
            idx{4} = strfind(c(1,:),'BAND_B');
    	    idx{5} = strfind(c(1,:),'BAND_G');
    	    idx{6} = strfind(c(1,:),'BAND_Y');
    	    idx{7} = strfind(c(1,:),'BAND_R');
    	    idx{8} = strfind(c(1,:),'BAND_RE');
    	    idx{9} = strfind(c(1,:),'BAND_N');
    	    idx{10} = strfind(c(1,:),'BAND_N2');
            idx{11} = strfind(c(1,:),'IMAGE');
            for i = 1:11;
                idxb(i,1:2) = find(not(cellfun('isempty',idx{i})));
            end
            szB(1) = str2num(s.Children(2).Children(idxb(1)).Children.Data);
    	    szB(2) = str2num(s.Children(2).Children(idxb(2)).Children.Data);
	    kf(1,1) = str2num(s.Children(2).Children(idxb(3)).Children(26).Children.Data);
	    kf(2,1) = str2num(s.Children(2).Children(idxb(4)).Children(26).Children.Data);
	    kf(3,1) = str2num(s.Children(2).Children(idxb(5)).Children(26).Children.Data);
	    kf(4,1) = str2num(s.Children(2).Children(idxb(6)).Children(26).Children.Data);
	    kf(5,1) = str2num(s.Children(2).Children(idxb(7,1)).Children(26).Children.Data);
	    kf(6,1) = str2num(s.Children(2).Children(idxb(8)).Children(26).Children.Data);
	    kf(7,1) = str2num(s.Children(2).Children(idxb(9,1)).Children(26).Children.Data);
	    kf(8,1) = str2num(s.Children(2).Children(idxb(10)).Children(26).Children.Data);
	    aqyear = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(1:4));
	    aqmonth = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(6:7));
	    aqday = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(9:10));
	    aqhour = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(12:13));
	    aqminute = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(15:16));
	    aqsecond = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(18:26));
	    sunel = str2num(s.Children(2).Children(idxb(11,2)).Children(56).Children.Data);
	    sunaz = str2num(s.Children(2).Children(idxb(11,2)).Children(50).Children.Data);
	    satview = str2num(s.Children(2).Children(idxb(11,2)).Children(86).Children.Data);
	    sensaz = str2num(s.Children(2).Children(idxb(11,2)).Children(62).Children.Data);
	    satel = str2num(s.Children(2).Children(idxb(11,2)).Children(68).Children.Data);
            cl_cov = str2num(s.Children(2).Children(idxb(11,2)).Children(90).Children.Data);
        end
        
        szB(3) = 8;

       % Assign WV2 vs WV3 constant calibration factors
	if id(4) == '3'
       		 ebw = ebw2;
       		 irr = irr2;
	         cw = cw2;
	else ebw = ebw1;
	         irr = irr1;
        	 cw = cw1;
	end

	% Identify growing season vs senesced
	if aqmonth == 11 || aqmonth == 12 || aqmonth == 1 || aqmonth == 2
		season = 0;
	else season = 1;
	end
	    %% Calculate Earth-Sun distance and relevant geometry
	    if aqmonth == 1 || aqmonth == 2;
	        year = aqyear -1;
	        month = aqmonth + 12;
	    else year = aqyear;
	        month = aqmonth;
	    end
	    UT = aqhour + (aqminute/60.0) + (aqsecond/3600.0); % Convert time to UT
	    B1 = int64(year/100);
	    B2 = 2-B1+int64(B1/4);
	    JD = (int64(365.25*(year+4716)) +int64(30.6001*(month+1)) + aqday + UT/24.0 + B2 - 1524.5); % Julian date
	    D = JD - 2451545.0;
	    degs = double(357.529 + 0.98560028*D); % Degrees
	    ESd = 1.00014 - 0.01671*cosd(degs) - 0.00014*cosd(2*degs); % Earth-Sun distance at given date (should be between 0.983 and 1.017)
	
	    inc_ang = 90.0 - sunel;
	    TZ = cosd(inc_ang); % Atmospheric spectral transmittance in solar path with solar zenith angle
	    TV = cosd(satview); % Atmospheric spectral transmittance in view path with satellite view angle

	    %% Calculate Rayleigh Path Radiance (Dash et al. 2012 and references therein)
	    if sunaz > 180 % For the following equations, azimuths should be between -180 and +180 degrees
	        sunaz = sunaz - 360;
	    end
	    if sensaz > 180
	        sensaz = sensaz - 360;
	    end
	    
	    az = abs(sensaz - 180 - sunaz); % Relative azimuth angle
	    thetaplus = acosd(cosd(90-sunel)*cosd(90-satel) - sind(90-sunel)*sind(90-satel)*cosd(az)); % Scattering angles

        for d = 1:8;
            Pr(d) = (3/(4*(1+2*gamma(d))))*((1+3*gamma(d))+(1-gamma(d))*cosd(thetaplus)^2); % Rayleigh scattering phase function (described in Bucholtz 1995)
        end
        
	    for d = 1:8;
	        tau(d) =(0.008569*(cw(d)^-4)*(1 + 0.0113*(cw(d)^-2) + 0.00013*cw(d)^-4)); % Rayleigh optical thickness (assume standard pressure of 1013.25 mb)
	    end
	    
	    % Rayleigh calculation (Dash et al., 2012)
	    for d = 1:8;
	        ray_rad{1,1}(d) = ((irr(1,d)/ESd)*1*tau(d)*Pr(d))/(4*pi*cosd(90-satel)); % Assume standard pressure (1013.25 mb)
	    end

	    % rrs constant calculation (Kerr et al. 2018 and Mobley 1994)
	G = single(1.7); % constant Li et al. 2019
        na = 1.00029; % Refractive index of air
        nw = 1.34; % Refractive index seawater
        inc_ang2 = real(asind(sind(90-satel)*nw/na)); % Incident angle for water-air from Snell's Law
        trans_aw = real(asind(sind(inc_ang)*na/nw)); % Transmission angle for air-water incident light from Snell's Law
	    trans_wa = 90-satel; % Transmission angle for water-air incident light from Snell's Law
	    pf1 = real(0.5*((sind(inc_ang - trans_aw)/(sind(inc_ang + trans_aw)))^2 + (tand(inc_ang - trans_aw)/(tand(inc_ang + trans_aw)))^2)); % Fresnel reflectance for air-water incident light (Mobley 1994)
	    pf2 = real(0.5*((sind(inc_ang2 - trans_wa)/(sind(inc_ang2 + trans_wa)))^2 + (tand(inc_ang2 - trans_wa)/(tand(inc_ang2 + trans_wa)))^2));
	    zeta = real(single((1-pf1)*(1-pf2)/(nw^2))); % rrs constant (~0.52) from Mobley 1994


        % Adjust file size: Input file (A) warped may contain more or fewer columns/rows than original NITF file, and some may be corrupt.
        sz(1) = min(szA(1),szB(1));
        sz(2) = min(szA(2),szB(2));
        sz(3) = 8;

            %% Assign NaN to no-data pixels and radiometrically calibrate and convert to Rrs
            Rrs = single(zeros(szA(1),szA(2),8)); % Create empty matrix for Rrs output
            for j = 1:sz(1); % Assign NaN to pixels of no data
                for k = 1:sz(2); % If a pixel contains data values other than "zero" or "two thousand and forty seven" in any band, it is calibrated; otherwise, it is considered "no-data" - this avoids a problem created during the orthorectification process wherein reprojecting the image may resample data
                    if (A(j,k,1)) ~= 0 && (A(j,k,1)) ~= 2047 || (A(j,k,2)) ~= 0 && (A(j,k,2)) ~= 2047 || (A(j,k,3)) ~= 0 && (A(j,k,3)) ~= 2047 || (A(j,k,4)) ~= 0 && (A(j,k,4)) ~= 2047 || (A(j,k,5)) ~= 0 && (A(j,k,5)) ~= 2047 || (A(j,k,6)) ~= 0 && (A(j,k,6)) ~= 2047 || (A(j,k,7)) ~= 0 && (A(j,k,7)) ~= 2047 || (A(j,k,8)) ~= 0 && (A(j,k,8)) ~= 2047;
                        for d = 1:8;
                            Rrs(j,k,d) = single((pi*((single(A(j,k,d))*kf(d,1)/ebw(1,d)) - ray_rad{1,1}(1,d))*ESd^2)/(irr(1,d)*TZ*TV)); % Radiometrically calibrate and convert to Rrs (adapted from Radiometric Use of WorldView-2 Imagery(
                        end
                    else Rrs(j,k,:) = NaN;
                    end
                end
            end

            clear A

            %% Output reflectance image
%            if Rrs_write == 1;
%		if id(4) == '3'
%			info = geotiffinfo(images);
%			geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
%			tiffTags = struct('TileLength',1024,'TileWidth',1024);
%			Z = [loc_out,id,'_',loc,'_RrsBT']
%			geotiffwrite(Z,Rrs,R(1,1),'GeoKeyDirectoryTag',geoTags,'TiffType','bigtiff','TiffTags',tiffTags);
%		else
%	                Z = [loc_out,id,'_',loc,'_Rrs']
%	                geotiffwrite(Z,Rrs,R(1,1),'CoordRefSysCode',coor_sys);
%		end
%            end

	if d_t > 0; % Run DT and/or rrs conversion; otherwise end

	%% Setup for Deglint, Bathymetry, and Decision Tree
	b = 1;
        t = 1;
	    u = 1;
        y = 0;
	    v = 0;
	    num_pix = 0;
	    sum_SD(b) = 0;
	    sum_veg(t) = 0;
	    sum_veg2(t) = 0;
	    dead_veg(t) = 0;
	    sum_water_rrs(u) = 0;
	    sz_ar = sz(1)*sz(2);
	    water = zeros(sz_ar,9);
	    for j = 1:sz(1);
	        for k = 1:sz(2);
	            if isnan(Rrs(j,k,1)) == 0
                    num_pix = num_pix +1; % Count number of non-NaN pixels
                    c_val(num_pix) = Rrs(j,k,1); % Record coastal band value for use in cloud mask prediction
                    if (Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2)) < 0.65 && Rrs(j,k,5) > Rrs(j,k,4) && Rrs(j,k,4) > Rrs(j,k,3) % Sand & Developed
                        sum_SD(b) = sum(Rrs(j,k,6:8));
                        b = b+1;
                    elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.65 && Rrs(j,k,7) > Rrs(j,k,3); % Identify vegetation (excluding grass)
                        if ((Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2))) > 0.20; % Shadow filter
	                    sum_veg(t) = sum(Rrs(j,k,3:5)); % Sum bands 3-5 for selected veg to distinguish wetland from upland
			    sum_veg2(t) = sum(Rrs(j,k,7:8));
                            dead_veg(t) = (((Rrs(j,k,7) - Rrs(j,k,4))/3) + Rrs(j,k,4)) - Rrs(j,k,5); % Compute difference of predicted B5 value from actual valute
	                        t = t+1;
                        end
 		    elseif Rrs(j,k,8) < 0.11 && Rrs(j,k,1) > 0 && Rrs(j,k,2) > 0 && Rrs(j,k,3) > 0 && Rrs(j,k,4) > 0 && Rrs(j,k,5) > 0 && Rrs(j,k,6) > 0 && Rrs(j,k,7) > 0 && Rrs(j,k,8) > 0; % Identify glint-free water
                        water(u,1:8) = double(Rrs(j,k,:));
			water_rrs(1:6) = Rrs(j,k,1:6)./(zeta + G.*Rrs(j,k,1:6));
			if water_rrs(4) > water_rrs(2) && water_rrs(4) < 0.12 && water_rrs(5) < water_rrs(3)
				sum_water_rrs(u) = sum(water_rrs(3:5));
			end
                        u = u+1;
        		if Rrs(j,k,8)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,3)% NDGI to identify glinted water pixels (some confusion w/ clouds)
                            v = v+1;
                            water(u,9) = 2; % Mark array2<array1 glinted pixels
                        elseif Rrs(j,k,8)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,3)
                            v = v+1;
                            water(u,9) = 3; % Mark array2>array1 glinted pixels
                        else water(u,9) = 1; % Mark records of glint-free water
                        end
                    elseif Rrs(j,k,8)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,3)
                        water(u,1:8) = double(Rrs(j,k,:));
                        water(u,9) = 2; % Mark array2<array1 glinted pixels
                        u = u+1;
                        v = v+1;
                    elseif Rrs(j,k,8)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,3)
                        water(u,9) = 3; % Mark array2>array1 glinted pixels
                        water(u,1:8) = double(Rrs(j,k,:));
                        u = u+1;
            			v = v+1;
%                     elseif (Rrs(j,k,4)-Rrs(j,k,8))/(Rrs(j,k,4)+Rrs(j,k,8)) < 0.55 && Rrs(j,k,8) < 0.2 && (Rrs(j,k,7)-Rrs(j,k,2))/(Rrs(j,k,7)+Rrs(j,k,2)) < 0.1 && (Rrs(j,k,8)-Rrs(j,k,5))/(Rrs(j,k,8)+Rrs(j,k,5)) < 0.3 && Rrs(j,k,1) > 0 && Rrs(j,k,2) > 0 && Rrs(j,k,3) > 0 && Rrs(j,k,4) > 0 && Rrs(j,k,5) > 0 && Rrs(j,k,6) > 0 && Rrs(j,k,7) > 0 && Rrs(j,k,8) > 0; % Identify glinted water
%                         water(u,1:8) = double(Rrs(j,k,:));
%                         u = u+1;
%             			  v = v+1;

                    end
	            end
	        end
	    end
		n_water = u; % Number of water pixels used to derive E_glint relationships
		n_glinted = v; % Number of glinted water pixels

		idx = find(water(:,1) == 0);
		water(idx,:) = [];
        water7 = water(:,7);
        water8 = water(:,8);
        mnNIR1 = min(water7(water7>0)); % Positive minimum Band 7 value used for deglinting
        mnNIR2 = min(water8(water8>0)); % Positive minimum Band 8 value used for deglinting
        
        idx_gf = find(water(:,9)==1); % Glint-free water
	water_gf = water(idx_gf,1:8);

% Identify optically deep water average spectrum
        bn = 7; % Band number
        pctl_l = 5; % Percentile (5th percentile value of glint-free water n-band values chosen based on visual analysis of density slicing of Rrs image)
        pctl_u = 15;
        clear water_gfidx water_odw m0 m1
        water_gfidx = find(water_gf(:,bn) == prctile(water_gf(:,bn),pctl_l) & water_gf(:,bn) <= prctile(water_gf(:,bn),pctl_u));
        water_odw(:,1:8) = (water_gf(water_gfidx(1:end),1:8)); % Li et al. Dove BGR corresponds to WV2 BGY center wavelengths

%         Equations from Li et al. 2019 & Hu et al. 2012
        for h = 1:size(water_odw,1)
%             w1(h) = water_odw(h,3) - (water_odw(h,1) + (546-427)/(659-427)*(water_odw(h,5) - water_odw(h,1))); % Hu et al. 2012
            w2(h) = water_odw(h,3) - 0.46*water_odw(h,4) - 0.54*water_odw(h,1); % Li et al. 2019
        end

	if exist('w2')==1 
	        w = median(w2(w2<0));
	else w = 0;
	end

        if w > -0.0005
            m0 = 0;
            m1 = 0;
            Update = 'Too Turbid for Benthic Mapping'
        else
            chla = 10^(-0.4909 + 191.659*w) % Hu et al. 2012 (Kerr limited chla to 1.0mg/m3; 0.1 mg/m3 WV Cay Sal most accurate value used)
            m0 = 52.083*exp(2.711*chla) % Revised from Li et al. 2019 with exponential scalar derived from Kerr FK WV image field data tuning parameters
            m1 = 50.156*exp(2.711*chla) % TARGET: 64.3 +/- 0.5 & 62.6 +/- 0.5, Predicted: 67.2 & 64.7
        end

	Kd = [0.036 0.037 0.075 0.25 0.415]; %1.416]; %(Based on Kerr 2018 Fig 7a chl-conc 0.1 mg/m3 i.e. lowest RMSE water-depth predictor values)

		if v > 0.25*u
			Update = 'Deglinting'
            		id2 = 'deglinted';
                	for b = 1:6 %% Calculate linear fitting of all MS bands vs NIR1 & NIR2 for deglinting in DT (Hedley et al. 2005)
                        	if b == 1 || b == 4 || b == 6
                                slope1 = water(:,b)\water(:,8);
               		        else slope1 = water(:,b)\water(:,7);
                        	end
                	E_glint(1,b) = single(slope1);
                	end
			E_glint % = [0.8075    0.7356    0.8697    0.7236    0.9482    0.7902]
		else Update = 'Glint-free'
            id2 = 'glintfree';
        end
        
        %% Edge Detection via Morphological Index (improved over Huang & Zhang 2011, Ma et al. 2019)
       waterind = uint16((Rrs(:,:,3)-Rrs(:,:,8))./(Rrs(:,:,3)+Rrs(:,:,8)) > 0.15);
       img_sub2 = Rrs(:,:,2);
       img_sub5 = Rrs(:,:,5);
       img_sub7 = Rrs(:,:,7);

       Rrs_cloud = img_sub2./img_sub7;
       Rrs_cl2 = Rrs_cloud;
       Rrs_cl3 = Rrs_cloud;
       Rrs_cl2(Rrs_cloud >= 0.7) = 1;
       Rrs_cl2(Rrs_cloud < 0.7) = 0;
       Rrs_cl3(Rrs_cloud <= 0.9) = 1;
       Rrs_cl3(Rrs_cloud > 0.9) = 0;
       Rrs_clf = Rrs_cl2 + Rrs_cl3;
       Rrs_clf(Rrs_clf < 2) = 0;
       CLrrs = imbinarize(Rrs_clf);
       CL1 = uint16(imtophat(CLrrs,strel('disk',100))) - waterind;
       CL1(CL1<0) = 0;
       CLe = imerode(CL1,strel('disk',20));
       CLed = imdilate(CLe,strel('disk',150));
       Cloud = imfill(CLed,'holes');
	clear Rrs_cl2 Rrs_cl3 Rrs_clf CLrrs CL1 CLe CLed

	Rrs_sh1 = Rrs_cloud;
	Rrs_sh2 = Rrs_cloud;
	Rrs_sh1(Rrs_cloud >= 1.3) = 1;
	Rrs_sh1(Rrs_cloud < 1.3) = 0;
	Rrs_sh2(Rrs_cloud <= 1.7) = 1;
	Rrs_sh2(Rrs_cloud > 1.7) = 0;
	Rrs_shf = Rrs_sh1 + Rrs_sh2;
	Rrs_shf(Rrs_shf < 2) = 0;
	SHrrs = imbinarize(Rrs_shf);
	Shadow = uint16(imtophat(SHrrs,strel('square',20)));
	clear Rrs_sh1 Rrs_sh2 Rrs_shf SHrrs

        Rrs_map = img_sub5./img_sub7;
	Rrs_map2 = Rrs_map;
	Rrs_map3 = Rrs_map;
	Rrs_map2(Rrs_map >= 0.7) = 1;
	Rrs_map2(Rrs_map < 0.7) = 0;
	Rrs_map3(Rrs_map <= 1.1) = 1;
	Rrs_map3(Rrs_map > 1.1) = 0;
	Rrs_mapf = Rrs_map2 + Rrs_map3;
	Rrs_mapf(Rrs_mapf < 2) = 0;
        BWrrs = imbinarize(Rrs_mapf);

        BW1 = uint16(imtophat(BWrrs,strel('square',30))) - waterind;
	BW1 = imdilate(BW1,strel('square',5)); % Expand developed to include shadows
	BW1(BW1<0) = 0;

        Cloud = Cloud - BW1;
        Cloud(Cloud<0) = 0;
        cld_idx = 0;
        if size(find(Cloud ==1),1) > 0.060*szA(1)*szA(2)
                cld_idx = 1;
        end


%	ns = 2000;
%	BW = uint16(imtophat(BWrrs,strel('square',ns)));
%	CC = bwconncomp(BW);
%	numPixels = cellfun(@numel,CC.PixelIdxList);
%	BW1idx = find(numPixels > 1000);
%	CC.PixelIdxList = CC.PixelIdxList(BW1idx);
%	CC.NumObjects = size(BW1idx,2);
%	BW3 = uint16(labelmatrix(CC));
%	BW3(BW3>0) = 1;
%	BW3e = uint16(imerode(BW3,strel('disk',100)));
%	BW3ed = uint16(imdilate(BW3e,strel('square',200)));
%	BW4 = imfill(BW3ed,'holes');

	BAI = (img_sub2 - img_sub7)./(img_sub2 + img_sub7); % Built Area Index
	BAI = BAI * -1; % Dev & soil negative, soil more negative (water high positive)
	BAI = imbinarize(BAI);
	BAI = imerode(BAI,strel('square',5));

	clear BW3 BW3e BW3ed BW2 BWrrs BWnew BWnewe BW1idx

	Ztest = [loc_out,id,'_',loc,'_BW1']
        geotiffwrite(Ztest,BW1,R(1,1),'CoordRefSysCode',coor_sys);
        Ztest = [loc_out,id,'_',loc,'_BAI']
        geotiffwrite(Ztest,BAI,R(1,1),'CoordRefSysCode',coor_sys);



        %% Determine Rrs-infinite from glint-free water pixels
        rrs_inf = [0.00512 0.00686 0.008898 0.002553 0.001506 0.000403]; % Derived from Rrs_Kd_Model.xlsx for Default values

        %% Calculate target class metrics
        avg_SD_sum = mean(sum_SD(:));
        stdev_SD_sum = std(sum_SD(:));
	avg_veg_sum = mean(sum_veg(:))
	avg_dead_veg = mean(dead_veg(:));
	avg_mang_sum = mean(sum_veg2(:));
	idx_water2 = find(sum_water_rrs==0);
	sum_water_rrs(idx_water2) = [];
	avg_water_sum = mean(sum_water_rrs(:));

	    if cl_cov > 0
		    num_cld_pix = round(num_pix*cl_cov*0.01); % Number of cloud pixels (rounded down to nearest integer) based on metadata-reported percent cloud cover
		    srt_c = sort(c_val,'descend'); % Sort all pixel blue-values in descending order. Cloud mask threshold will be num_cld_pix'th highest value
		    cld_mask = srt_c(num_cld_pix); % Set cloud mask threshold
	    else cld_mask = max(c_val)+1;
	    end


	    Bathy = single(zeros(szA(1),szA(2))); % Preallocate for Bathymetry
	    Rrs_deglint = single(zeros(5,1)); % Preallocate for deglinted Rrs
	    Rrs_0 = single(zeros(5,1)); %Preallocation for water-column corrected Rrs
 	    map = zeros(szA(1),szA(2),'uint8'); % Create empty matrix for classification output

	if d_t == 1; % Execute Deglinting, rrs, Bathymetry
            if v > u*0.25
            for j = 1:szA(1)
                for k = 1:szA(2)
                    if isnan(Rrs(j,k,1)) == 0 && Rrs(j,k,8)<0.2
                        % Deglint equation
                        Rrs_deglint(1,1) = (Rrs(j,k,1) - (E_glint(1)*(Rrs(j,k,8) - mnNIR2)));
                        Rrs_deglint(2,1) = (Rrs(j,k,2) - (E_glint(2)*(Rrs(j,k,7) - mnNIR1)));
                        Rrs_deglint(3,1) = (Rrs(j,k,3) - (E_glint(3)*(Rrs(j,k,7) - mnNIR1)));
                        Rrs_deglint(4,1) = (Rrs(j,k,4) - (E_glint(4)*(Rrs(j,k,8) - mnNIR2)));
                        Rrs_deglint(5,1) = (Rrs(j,k,5) - (E_glint(5)*(Rrs(j,k,7) - mnNIR1)));
                        Rrs_deglint(6,1) = (Rrs(j,k,6) - (E_glint(6)*(Rrs(j,k,8) - mnNIR2)));

                        % Convert above-surface Rrs to below-surface rrs (Kerr et al. 2018)
                        Rrs_0(1:5) = Rrs_deglint(j,k,1:5)./(zeta + G.*Rrs_deglint(j,k,1:5)); % Convert above-surface Rrs to subsurface rrs (Kerr et al. 2018, Lee et al. 1998)
                        b1 = 63.6; % Turning parameters (Kerr 2018)
                        b0 = -60.25;
                        dp = b1*real(log(1000*Rrs_0(2))/log(1000*Rrs_0(3))) + b0; % Calculate depth (Stumpf 2003 ratio transform with Kerr et al. 2018 coefficients)
                        if dp < 15 && dp > 0 % Parameters based on Kerr 2018 RMSE-based recommended constraints (depths greater than 15m inaccurate)
                            Bathy(j,k) = dp;
                        end
                        for d = 1:5
                            Rrs(j,k,d) = real(((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1,d)*dp))+rrs_inf(d)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
                        end
                    end
                end
            end
        else % For glint-free/low-glint images
            for j = 1:szA(1)
                for k = 1:szA(2)
                    if isnan(Rrs(j,k,1)) == 0 && Rrs(j,k,8)<0.2
                        Rrs_0(1:5) = Rrs(j,k,1:5)./(zeta + G.*Rrs(j,k,1:5)); % Convert above-surface Rrs to subsurface rrs (Kerr et al. 2018, Lee et al. 1998)
                        b1 = 63.6; % Turning parameters (Kerr 2018 Table 6 average of 2 forward-modeling WorldView-2 results)
                        b0 = -60.25;
                        dp = b1*real(log(1000*Rrs_0(2))/log(1000*Rrs_0(3))) + b0; % Calculate depth (Stumpf 2003 ratio transform with Kerr et al. 2018 coefficients)
                        if dp < 15 && dp > 0 % Parameters based on Kerr 2018 RMSE-based recommended constraints (depths greater than 15m inaccurate)
                            Bathy(j,k) = dp;
                        else dp = 0;
                        end
                        for d = 1:5
                            Rrs(j,k,d) = real(((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1,d)*dp))+rrs_inf(d)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
                        end
                    end
                end
            end
        end
	elseif d_t == 2; % Only run for Deglinted Rrs and Bathymetry, not Decision Tree
	    update = 'Running DT'
	BS = 2;
	WA = 3;
	DG = 5;
	MA = 6;
	SC = 7;
	FW = 10;
	FU = 9;
	UG = 8;
	dev = 11;
        p = 1;
            for j = 1:szA(1)
               for k = 1:szA(2)
                   if isnan(Rrs(j,k,1)) == 0
                       %% Cloud Cover
		       if Cloud(j,k) == 1 && BW1(j,k) ~= 1
                           map(j,k) = 1; % Cloud
                       %% Vegetation
                       elseif (Rrs(j,k,7) - Rrs(j,k,5))/(Rrs(j,k,7) + Rrs(j,k,5)) > 0.20 && Rrs(j,k,7) > Rrs(j,k,3) % Vegetation pixels (NDVI)
			    if ((Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2))) < 0.20 && (Rrs(j,k,7) - Rrs(j,k,8))/(Rrs(j,k,7) + Rrs(j,k,8)) > 0.01; % Shadowed-vegetation filter (B7/B8 ratio excludes marsh, which tends to have very similar values here)
                            	map(j,k) = 0; % Shadow
                            elseif sum(Rrs(j,k,3:5)) < avg_veg_sum
                                if (Rrs(j,k,3) - Rrs(j,k,8))/(Rrs(j,k,3) + Rrs(j,k,8)) > -0.75 % ML
                                    if (Rrs(j,k,7) - Rrs(j,k,5))/(Rrs(j,k,7) + Rrs(j,k,5)) > 0.75 % M
	                                map(j,k) = FW; % Forested Wetland
                                    elseif sum(Rrs(j,k,3:5)) > 0.12 && sum(Rrs(j,k,7:8)) > 0.45 % ML
                                        map(j,k) = FU; % FORESTED UPLAND
                                    elseif (Rrs(j,k,7) - Rrs(j,k,5))/(Rrs(j,k,7) + Rrs(j,k,5)) > 0.60
                                        map(j,k) = FW; % Forested Wetland
                                    elseif Rrs(j,k,7) < 0.3 && sum(Rrs(j,k,7:8)) > 0.25
					if (Rrs(j,k,5) - Rrs(j,k,3))/(Rrs(j,k,5) + Rrs(j,k,3)) > 0.1
						map(j,k) = DG; % Dead Grass
					elseif Rrs(j,k,7) < 0.27 && sum(Rrs(j,k,7:8)) < 0.5
						map(j,k) = MA; % Marsh
					else map(j,k) = FU; % Forested Upland
					end
                                    end
				elseif (Rrs(j,k,4) - Rrs(j,k,5))/(Rrs(j,k,4) + Rrs(j,k,5)) > 0.08
					map(j,k) = 6; % Marsh (was algal flat)
				else  map(j,k) = FU; % Forested Upland
                                end
			    elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.65
				map(j,k) = FU; % Forested Upland
                            elseif Rrs(j,k,7) < 0.4 % Marsh, Scrub, Grass, Dead Veg
				if (Rrs(j,k,4) - Rrs(j,k,5))/(Rrs(j,k,4) + Rrs(j,k,5)) > 0.08
					map(j,k) = 6; % Marsh (was algal flat)
				elseif (Rrs(j,k,5) - Rrs(j,k,3))/(Rrs(j,k,5) + Rrs(j,k,3)) > 0.05 %&& Rrs(j,k,7) < 0.27 % Agriculture or senesced veg/grass
					map(j,k) = DG; % Dead veg
				else map(j,k) = UG; % Grass
				end
%			    elseif sum(Rrs(j,k,7:8)) < 0.8 && sum(Rrs(j,k,7:8)) > 0.65 % Live grass high, dead grass low
%				map(j,k) = 10; % Upland Forest
			    else map(j,k) = SC; % Scrub/shrub
                            end
                       %% Developed and Soil
                       elseif (Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2)) < 0.60 && Rrs(j,k,5) > Rrs(j,k,4) && waterind(j,k) == 0 %Rrs(j,k,8) > 0.1 % && Rrs(j,k,4) > Rrs(j,k,3)
                           if Rrs(j,k,5)/Rrs(j,k,7) > 0.7 && Rrs(j,k,5)/Rrs(j,k,7) < 1.1
			       if BAI(j,k) == 0 && BW1(j,k) == 1 %BW4(j,k) == 1
					map(j,k) = dev; %Developed. Was: BS; % Soil (fallow field)
			       elseif BAI(j,k) == 1 && BW1(j,k) == 0
					map(j,k) = BS; % Soil
                               elseif BW1(j,k) == 1
                                   if sum(Rrs(j,k,1:2))<0.35
                                       if sum(Rrs(j,k,6:8)) < 0.85%avg_SD_sum
                                           map(j,k) = dev; % Developed
                                       else map(j,k) = BS; % Soil
                                       end
                                   elseif sum(Rrs(j,k,1:2)) > 0.6
                                       map(j,k) = dev;
				else map(j,k) = dev;
                                   end
                               elseif sum(Rrs(j,k,6:8)) < avg_SD_sum
                                   map(j,k) = dev;
                               else map(j,k) = dev; % Developed
                               end
                           else map(j,k) = BS; % Soil
                           end
                       %% Water
                       elseif Rrs(j,k,8)<0.2 && Rrs(j,k,8)>0|| Rrs(j,k,8)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,3) && Rrs(j,k,8)>0 || Rrs(j,k,8)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,3) && Rrs(j,k,8)>0% Identify all water (glinted and glint-free)
                           if v > u*0.25 && u>0.1*num_pix
                                % Deglint equation
                                Rrs_deglint(1,1) = (Rrs(j,k,1) - (E_glint(1)*(Rrs(j,k,8) - mnNIR2)));
                                Rrs_deglint(2,1) = (Rrs(j,k,2) - (E_glint(2)*(Rrs(j,k,7) - mnNIR1)));
                                Rrs_deglint(3,1) = (Rrs(j,k,3) - (E_glint(3)*(Rrs(j,k,7) - mnNIR1)));
                                Rrs_deglint(4,1) = (Rrs(j,k,4) - (E_glint(4)*(Rrs(j,k,8) - mnNIR2)));
                                Rrs_deglint(5,1) = (Rrs(j,k,5) - (E_glint(5)*(Rrs(j,k,7) - mnNIR1)));
                                Rrs_deglint(6,1) = (Rrs(j,k,6) - (E_glint(6)*(Rrs(j,k,8) - mnNIR2)));

                                % Convert above-surface Rrs to below-surface rrs (Kerr et al. 2018)
                                Rrs_0(1:5) = Rrs_deglint(1:5)./(zeta + G.*Rrs_deglint(1:5)); % Was Rrs_0=
                                % Relative depth estimate
                                dp = m0*real(log(1000*Rrs_0(1))/log(1000*Rrs_0(3))) - m1; % Calculate depth (Stumpf 2003 ratio transform with Kerr et al. 2018 coefficients)

                                if dp < 15 && dp > 0 % Parameters based on Kerr 2018 RMSE-based recommended constraints (depths greater than 15m inaccurate)
                                    Bathy(j,k) = dp;
                                else dp = 0;
                                end

%                                    for d = 1:5
%                                        Rrs(j,k,d) = real(((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1,d)*dp))+rrs_inf(d)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
%                                    end

                                    %% DT
                                   if  Shadow(j,k) == 1 && max(Rrs(j,k,:)) == Rrs(j,k,2) % Max band3-6 = turbid/shallow water
                                       map(j,k) = 0; % Shadow
                                   else map(j,k) = WA; % Deep water
                                   end
                            else % For glint-free/low-glint images
                                Rrs_0(1:5) = Rrs(j,k,1:5)./(zeta + G.*Rrs(j,k,1:5)); % Convert above-surface Rrs to subsurface rrs (Kerr et al. 2018, Lee et al. 1998)
                                dp = m0*real(log(1000*Rrs_0(2))/log(1000*Rrs_0(3))) - m1; % Calculate depth (Stumpf 2003 ratio transform with Kerr et al. 2018 coefficients)
                                if dp < 15 && dp > 0 % Parameters based on Kerr 2018 RMSE-based recommended constraints (depths greater than 15m inaccurate)
                                    Bathy(j,k) = dp;
                                else dp = 0;
                                end
                                   %% DT
                                   if  Shadow(j,k) == 1  && max(Rrs(j,k,:)) == Rrs(j,k,2)  % Max band3-6 = turbid/shallow water
                                       map(j,k) = 0; % Shadow/Unclassified
                                   else map(j,k) = WA; % Deep water
%                                    end
                                   end
                             end % if v>u
                       end % If water/land
                   end % If isnan
               end % k

                end % j
    end


%%      DT Filter
         if filter > 0
                update = 'Filtering'
                dt_filt = DT_Filter(map,filter,sz(1),sz(2),dev,FW,FU,UG,WA);
                if cld_idx == 1
                         AA = [loc_out,id,'_',loc,'_SOALCHI_filt_',num2str(filter),'_Cloudy'];
                else AA = [loc_out,id,'_',loc,'_SOALCHI_filt_',num2str(filter)];
                end
                 geotiffwrite(AA,dt_filt,R(1,1),'CoordRefSysCode',coor_sys);
         else
             Z1 = [loc_out,id,'_',loc,'_Map_nofilt'];
             geotiffwrite(Z1,map,R(1,1),'CoordRefSysCode',coor_sys);
         end

%         TP(z,1) = m0;
%         TP(z,2) = m1;
%         TP(z,3) = chla;

        %% Output images
%         Z = [loc_out,id,'_',loc,'_Bathy_MAv1'];
%             geotiffwrite(Z,Bathy,R(1,1),'CoordRefSysCode',coor_sys);

%            Z2 = [Rrs_out,id,'_',loc,'_Rrs']; % last=52
%          geotiffwrite(Z2,Rrs,R(1,1),'CoordRefSysCode',coor_sys);
%
   end % If dt>0


        wtime = toc;
        time_min = wtime/60;
        fprintf(1,'Matlab CPU time (minutes) = %f\n', time_min);

end

