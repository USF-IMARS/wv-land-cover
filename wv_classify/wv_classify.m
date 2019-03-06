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


clear all
tic

%% Assign input and output locations
loc = 'RB'; % Typically the estuary acronym
coor_sys = 4326; % Change coordinate system code here
Rrs_write = 0; % 1=write Rrs geotiff; 0=do not write
d_t = 1; % 0=End after Rrs conversion; 1 = rrs, bathy, DT; 2 = rrs, bathy and DT
filter = 3; % 0=None, 1=3x3, 3=7x7, 5=11x11
% sgw = 0; % Sunglint moving-window box = sgw*2 +1 (i.e. 2 = 5x5 box)
% sgwid = num2str(sgw)

loc_in = '/home1/mmccarthy/Matt/USF/Other/NERRS_Mapping/Processing/Ortho/';
met_in = '/home1/mmccarthy/Matt/USF/Other/NERRS_Mapping/Processing/Raw/';
loc_out = '/home1/mmccarthy/Matt/USF/Other/NERRS_Mapping/Processing/Output/';
matfiles = dir(fullfile('Matt','USF','Other','NERRS_Mapping','Processing','Ortho','*.tif'));
matfiles2 = dir(fullfile('Matt','USF','Other','NERRS_Mapping','Processing','Raw','*.xml')); %% Revise this to find both all-caps and all lower-case extensions

% loc_in = ['/home1/mmccarthy/Matt/USF/Other/Seagrass/test/'];
% met_in = ['/home1/mmccarthy/Matt/USF/Other/Seagrass/test/'];
% loc_out = ['/home1/mmccarthy/Matt/USF/Other/Seagrass/test/Rrs/'];
% matfiles = dir(fullfile('Matt','USF','Other','Seagrass','test','*.tif'));
% matfiles2 = dir(fullfile('Matt','USF','Other','Seagrass','test','*.xml'));
 

sz_files = size(matfiles(:,1),1)

% Assign constants for all images
ebw = 0.001*[47.3 54.3 63.0 37.4 57.4 39.3 98.9 99.6]; % Effective Bandwidth per band (nm converted to um units; from IMD metadata files)
irr = [1758.2229 1974.2416 1856.4104 1738.4791 1559.4555 1342.0695 1069.7302 861.2866]; % Band-averaged Solar Spectral Irradiance (W/m2/um units)
cw = [.4273 .4779 .5462 .6078 .6588 .7237 .8313 .9080]; % Center wavelength (used for Rayleigh correction; from Radiometric Use of WorldView-2 Imagery)
gamma = 0.01*[1.499 1.471 1.442 1.413 1.413 1.413 1.384 1.384]; % Factor used in Rayleigh Phase Function equation (Bucholtz 1995)

for z = 1;%:sz_files;
id = matfiles(z,1).name(1:19)

X = [loc_in, matfiles(z,1).name]; % Change location of MS Tiff images here
Z = [met_in, matfiles2(z,1).name]; % Change location of XML files here

    [A, R] = geotiffread(X);
    szA = size(A);
     s = xml2struct(Z);
%    save XMLtest.mat s
        % Extract calibration factors and acquisition time from metadata for each band
        if isfield(s,'IMD') == 1
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
        else
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
        end
        
        szB(3) = 8;

       
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

        for d = 1:8
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
	    G = single(1.56); % constant (Kerr eq. 3)
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
            if Rrs_write == 1;
                Z = [loc_out,id,'_',loc,'_Rrs']
                geotiffwrite(Z,Rrs,R(1,1),'CoordRefSysCode',coor_sys);
            end

	if d_t > 0; % Run DT and/or rrs conversion; otherwise end
        
	    % Calculate Kd (water column attenuation coefficient) from Chuanmin Hu's Rrs_Kd_Model.xlsx sheet
% 		sunzen = 90.0-sunel;
% 		c1 = 0.005; % c1-4 hard-coded, but v1 and v2 change with modified values of aph(440), adg(440),bbp(440), Sdg, Y
% 		c2 = 4.18;
% 		c3 = 0.52;
% 		c4 = 10.8;
% 	 	v2 = [0.023277868 0.020883561 0.018975346 0.017605058 0.016771098 0.015875283 0.072734281 0.068046578]; % bb (backscatter)
%        v1 = [0.22024 0.142972 0.099157 0.286342 0.443809 1.491289 2.276262 2.223947]; % at (total absorption) DEFAULT CHL (0.1 0.1 0.01 0.015 0.5)
%         v1 = [0.069177389 0.041529229 0.061338823 0.268134177 0.411563873 1.489684614 2.275464564 2.22325881]; % Belize test
%         v2 = [0.023277868 0.020883561 0.018975346 0.017605058 0.016771098 0.015875283 0.049505144 0.046268752];
%        v1 = [0.006921 0.013933 0.051193 0.264353 0.409819 1.489006 2.275336 2.223238]; % at (total absorption) LOWER CHL
%        v2 = [0.023277868 0.020883561 0.018975346 0.017605058 0.016771098 0.015875283 0.000869138 0.00067143]; % bb (backscatter)
% 
% 		for b = 1:8
% 			Kd(b) = single((1+c1*sunzen)*v1(b)+c2*(1-c3*exp(-c4*v1(b)))*v2(b));
% 		end
	
%         Kd = [0.036 0.037 0.075 0.32 0.484 1.416];


	    %% Setup for Deglint, Bathymetry, and Decision Tree
	    b = 1;
        t = 1;
	    u = 1;
        y = 0;
	    v = 0;
	    num_pix = 0;
	    sum_veg(t) = 0;
	    dead_veg(t) = 0;
	    sz_ar = sz(1)*sz(2);
	    water = zeros(sz_ar,9);
	    for j = 1:sz(1)
	        for k = 1:sz(2)
	            if isnan(Rrs(j,k,1)) == 0
                    num_pix = num_pix +1; % Count number of non-NaN pixels
                    c_val(num_pix) = Rrs(j,k,1); % Record coastal band value for use in cloud mask prediction
                    if (Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2)) < 0.65 && Rrs(j,k,5) > Rrs(j,k,4) && Rrs(j,k,4) > Rrs(j,k,3) % Sand & Developed
                        sum_SD(b) = sum(Rrs(j,k,6:8));
                        b = b+1;
                    elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.6 && Rrs(j,k,7) > Rrs(j,k,3); % Identify vegetation (excluding grass)
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

		if v > 0.25*u
			Update = 'Deglinting'
            id2 = 'deglinted';
%             idx_w1 = find(water(:,9)==2); % Glinted water array1>array2
%             idx_w2 = find(water(:,9)==3); % Glinted water array2>array1
%             water1 = [water(idx_gf,1:8);water(idx_w1,1:8)];
%             water2 = [water(idx_gf,1:8);water(idx_w2,1:8)];
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
        
        %% Edge Detection
       img_sub = Rrs(:,:,5);
       BWbin = imbinarize(img_sub);
       BW = imtophat(BWbin,strel('square',10));
%        BW1 = edge(BWtop,'canny');
%        seDil = strel('square',1);
%        BWdil = imdilate(BW1,seDil);
%        BW = imfill(BWdil,'holes');
%        
%        seDer = strel('',[5 5]);
%        BWer = imerode(BW,seDer);

       
%         %% Depth scaling
%         water10(:,1:2) = water(idx_gf,2:3);
%         water10(:,1:2) = water10(:,1:2)./(zeta + G.*water10(:,1:2));
%         waterdp = real(log(1000*(water10(:,1)))./log(1000*(water10(:,2))));
%         water_dp = waterdp(waterdp>0 & waterdp<2);
%         [N,X] = hist(water_dp);
%         med_dp = median(water_dp)
%         low = X(2); %avg_dp - 5*std(water_dp) %min(water_dp)
%         scale_dp = scale/(med_dp-low)
% 
%         clear water10
        %         std_dp = std(water_dp)
%         low = avg_dp - 2*std_dp % Assumed represents 0 depth or min depth
%         high = avg_dp + std_dp;

        %% Determine Rrs-infinite from glint-free water pixels
%         water_gf = water(idx_gf,1:8);
%         dp_max_sort = sortrows(water_gf,8,'ascend'); % Sort all values in water by NIR2 column (assumes deepest water is darkest is NIR2)
%         idx_dp = round(size(dp_max_sort,1)*0.001); % Use "deepest" 0.1% pixels
%         dp_pct = dp_max_sort(1:idx_dp,:);
%         dp_rrs = dp_pct(:,1:8)./(zeta + G.*dp_pct(:,1:8)); % Convert to subsurface rrs
%         rrs_inf = min(dp_rrs(:,1:8)); %median(dp_rrs(:,1:8)) - 2*std(dp_rrs(:,1:8)); % Mean and Median values too high
% %         rrs_inf = [0.00512 0.00686 0.008898 0.002553 0.001506 0.000403]; % Derived from Rrs_Kd_Model.xlsx for Default values
% %         plot(rrs_inf)
        %% Calculate target class metrics
        avg_SD_sum = mean(sum_SD(:));
        stdev_SD_sum = std(sum_SD(:));
	    avg_veg_sum = mean(sum_veg(:));
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

	if d_t == 1; % Execute Deglinting rrs and Bathymetry
        if v > u*0.25
            % Deglint equation
            Rrs_deglint(1,1) = (Rrs(j,k,1) - (E_glint(1)*(Rrs(j,k,8) - mnNIR2)));
            Rrs_deglint(2,1) = (Rrs(j,k,2) - (E_glint(2)*(Rrs(j,k,7) - mnNIR1)));
            Rrs_deglint(3,1) = (Rrs(j,k,3) - (E_glint(3)*(Rrs(j,k,7) - mnNIR1)));
            Rrs_deglint(4,1) = (Rrs(j,k,4) - (E_glint(4)*(Rrs(j,k,8) - mnNIR2)));
            Rrs_deglint(5,1) = (Rrs(j,k,5) - (E_glint(5)*(Rrs(j,k,7) - mnNIR1)));
            Rrs_deglint(6,1) = (Rrs(j,k,6) - (E_glint(6)*(Rrs(j,k,8) - mnNIR2)));
            
            % Convert above-surface Rrs to below-surface rrs (Kerr et al. 2018)
            Rrs(j,k,1:6) = Rrs_deglint(1:6)./(zeta + G.*Rrs_deglint(1:6)); % Was Rrs_0=

            % Relative depth estimate
            dp = real(log(1000*Rrs_0(2))/log(1000*Rrs_0(3))); % Calculate relative depth (Stumpf 2003 ratio transform scaled to 1-10)
            if dp > 0 && dp < 2
                Bathy(j,k) = dp;
            else dp = 0;
            end
%                                 for d = 1:5
%                                     Rrs(j,k,d) = real(((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1,d)*dp_sc))+rrs_inf(d)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
%                                 end
                                
        else % For glint-free/low-glint images
            Rrs(j,k,1:6) = Rrs(j,k,1:6)./(zeta + G.*Rrs(j,k,1:6)); % Convert above-surface Rrs to subsurface rrs (Kerr et al. 2018, Lee et al. 1998)
            dp = real(log(1000*Rrs_0(2))/log(1000*Rrs_0(3))); % Calculate relative depth (Stumpf 2003 ratio transform)
            if dp > 0 && dp < 2
                Bathy(j,k) = dp;
            else dp = 0;
            end
        end


	elseif d_t == 2; % Execute Deglinting rrs, Bathymetery, and Decision Tree
	    update = 'Running DT'
      	    for j = 1:szA(1)
               for k = 1:szA(2)
                   if isnan(Rrs(j,k,1)) == 0
                       %% Mud, Developed and Sand
                       if (Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2)) < 0.60 && Rrs(j,k,5) > Rrs(j,k,4) && Rrs(j,k,4) > Rrs(j,k,3)
                            if Rrs(j,k,7) < Rrs(j,k,2) && Rrs(j,k,8) > Rrs(j,k,5)
                                map(j,k) = 0; % Shadow
                            elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) < 0.01 && Rrs(j,k,8) > 0.05 % Buildings & bright sand
                                if BW(j,k) == 1
                                    map(j,k) = 11; % Developed
                                elseif sum(Rrs(j,k,6:8)) < avg_SD_sum
                                    map(j,k) = 22; % Mud (intertidal?)
                                else map(j,k) = 21; % Beach/sand/soil
                                end
                            elseif Rrs(j,k,5) > (Rrs(j,k,2)+((Rrs(j,k,7)-Rrs(j,k,2))/5)*2)
                                map(j,k) = 21; % Beach/sand/soil
                            elseif Rrs(j,k,5) < (((Rrs(j,k,7) - Rrs(j,k,2))/5)*3+Rrs(j,k,2))*0.60 && Rrs(j,k,7) > 0.2
                                map(j,k) = 61; % Marsh grass
                            else map(j,k) = 22; % Mud
                            end
                       elseif Rrs(j,k,2) > Rrs(j,k,3) && Rrs(j,k,7) > Rrs(j,k,3) && Rrs(j,k,2) < 0.1 && (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) < 0.20|| Rrs(j,k,8) > 0.05 && Rrs(j,k,7) > Rrs(j,k,2) && (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) < 0.1
                           if BW(j,k) == 1
                               map(j,k) = 11; % Shadow/Developed
                           else map(j,k) = 22; % Mud
                           end
                       %% Vegetation
                       elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.20 && Rrs(j,k,7) > Rrs(j,k,3) % Vegetation pixels (NDVI)
                           if Rrs(j,k,7) > Rrs(j,k,2) && ((Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2))) < 0.20 && (Rrs(j,k,7) - Rrs(j,k,8))/(Rrs(j,k,7) + Rrs(j,k,8)) > 0.01; % Shadowed-vegetation filter (B7/B8 ratio excludes marsh, which tends to have very similar values here)
                               map(j,k) = 0; % Shadow
                           elseif sum(Rrs(j,k,3:5)) < avg_veg_sum
                                if ((Rrs(j,k,2) - Rrs(j,k,5))/(Rrs(j,k,2) + Rrs(j,k,5))) < 0.4% Agriculture filter based on elevated Blue band values
                                    if Rrs(j,k,7) > 0.12 && sum(Rrs(j,k,7:8))/sum(Rrs(j,k,3:5)) > 2
                                        map(j,k) = 63; % Forested Wetland
                                    else map(j,k) = 61; % Dead vegetation or Marsh
                                    end
                                else map(j,k) = 62; % Forested Upland (most likely agriculture)
                                end
                           elseif sum(Rrs(j,k,7:8)) < avg_mang_sum
                               if ((Rrs(j,k,2) - Rrs(j,k,5))/(Rrs(j,k,2) + Rrs(j,k,5))) < 0.4% Agriculture filter based on elevated Blue band values
                                   if Rrs(j,k,7) > 0.12 && sum(Rrs(j,k,7:8))/sum(Rrs(j,k,3:5)) > 2
                                       map(j,k) = 63; % Forested Wetland
                                   else map(j,k) = 61; % Marsh or Dead Vegetation
                                   end
                                else map(j,k) = 62; % Forested Upland (most likely agriculture)
                                end
                           elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.65 % NDVI for high upland values
                                map(j,k) = 62; % Upland Forest/Grass;
                           elseif Rrs(j,k,5) > (((Rrs(j,k,7) - Rrs(j,k,2))/5)*3+Rrs(j,k,2))*0.60 && Rrs(j,k,7) < 0.2 % Difference of B5 from predicted B5 by slope of B7:B4 to distinguish marsh (old: live vs dead trees/grass/marsh)
                               map(j,k) = 61; % Marsh grass
                           elseif Rrs(j,k,7) < 0.12
                               map(j,k) = 60; % Dead vegetation
                           else
                               map(j,k) = 62; % Upland Forest/Grass
                           end
                       %% Water
                       elseif Rrs(j,k,8)<0.2 && Rrs(j,k,8)>0|| Rrs(j,k,8)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,7) && Rrs(j,k,6)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,5) && Rrs(j,k,4)<Rrs(j,k,3) && Rrs(j,k,8)>0 || Rrs(j,k,8)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,7) && Rrs(j,k,6)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,5) && Rrs(j,k,4)>Rrs(j,k,3) && Rrs(j,k,8)>0% Identify all water (glinted and glint-free)
%                            map(j,k) = 5;

                           if v > u*0.25
                                % Deglint equation
                                Rrs_deglint(1,1) = (Rrs(j,k,1) - (E_glint(1)*(Rrs(j,k,8) - mnNIR2)));
                                Rrs_deglint(2,1) = (Rrs(j,k,2) - (E_glint(2)*(Rrs(j,k,7) - mnNIR1)));
                                Rrs_deglint(3,1) = (Rrs(j,k,3) - (E_glint(3)*(Rrs(j,k,7) - mnNIR1)));
                                Rrs_deglint(4,1) = (Rrs(j,k,4) - (E_glint(4)*(Rrs(j,k,8) - mnNIR2)));
                                Rrs_deglint(5,1) = (Rrs(j,k,5) - (E_glint(5)*(Rrs(j,k,7) - mnNIR1)));
                                Rrs_deglint(6,1) = (Rrs(j,k,6) - (E_glint(6)*(Rrs(j,k,8) - mnNIR2)));
                                
                                % Convert above-surface Rrs to below-surface rrs (Kerr et al. 2018)
                                Rrs(j,k,1:6) = Rrs_deglint(1:6)./(zeta + G.*Rrs_deglint(1:6)); % Was Rrs_0=

                                % Relative depth estimate
                                dp = real(log(1000*Rrs_0(2))/log(1000*Rrs_0(3))); % Calculate relative depth (Stumpf 2003 ratio transform scaled to 1-10)
                                if dp > 0 && dp < 2
                                    Bathy(j,k) = dp;
                                else dp = 0;
                                end
%                                 dp_sc = (dp-low)*scale_dp;
                                
%                                 for d = 1:5
%                                     Rrs(j,k,d) = real(((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1,d)*dp_sc))+rrs_inf(d)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
%                                 end
                                
                                %% DT
                               if  Rrs(j,k,6) < Rrs(j,k,7) 
                                   map(j,k) = 0; % Shadow
                               elseif (Rrs(j,k,3) - Rrs(j,k,4))/(Rrs(j,k,3) + Rrs(j,k,4)) < 0.10 %(Rrs(j,k,2) - Rrs(j,k,4))/(Rrs(j,k,2)+Rrs(j,k,4)) < 0
                                    if Rrs(j,k,4) > Rrs(j,k,3) || Rrs(j,k,5) > Rrs(j,k,3)
                                        map(j,k) = 53; % Soft bottom
                                    elseif sum(Rrs(j,k,3:5)) > avg_water_sum && (Rrs(j,k,5) - Rrs(j,k,2))/(Rrs(j,k,5) + Rrs(j,k,2)) > 0.1 % NEW from 0.05
                                        map(j,k) = 52; % Soft bottom
                                    elseif Rrs(j,k,4) > Rrs(j,k,2) && (Rrs(j,k,3) - Rrs(j,k,6))/(Rrs(j,k,3) + Rrs(j,k,6)) < 0.60 % Separate seagrass from dark water NEW
                                        if (Rrs(j,k,3) - Rrs(j,k,5))/(Rrs(j,k,3) + Rrs(j,k,5)) > 0.10 % Separate seagrass from turbid water NEW
                                            map(j,k) = 54; % Seagrass
                                        else map(j,k) = 55; % Turbid water
                                        end
                                    else map(j,k) = 51; % Deep water
                                    end
                                else map(j,k) = 51; % Deep water
                                end
                                                                     
                                    
                            else % For glint-free/low-glint images
                                Rrs(j,k,1:6) = Rrs(j,k,1:6)./(zeta + G.*Rrs(j,k,1:6)); % Convert above-surface Rrs to subsurface rrs (Kerr et al. 2018, Lee et al. 1998)
                                dp = real(log(1000*Rrs_0(2))/log(1000*Rrs_0(3))); % Calculate relative depth (Stumpf 2003 ratio transform)
                                if dp > 0 && dp < 2
                                    Bathy(j,k) = dp;
                                else dp = 0;
                                end
%                                 dp_sc = (dp-low)*scale_dp;
                                
%                                 for d = 1:5
%                                     Rrs(j,k,d) = real(((Rrs_0(d)-rrs_inf(d))/exp(-2*Kd(1,d)*dp_sc))+rrs_inf(d)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
%                                 end
                                %% DT
                               if  Rrs(j,k,6) < Rrs(j,k,7) 
                                   map(j,k) = 0; % Shadow
                               elseif (Rrs(j,k,3) - Rrs(j,k,4))/(Rrs(j,k,3) + Rrs(j,k,4)) < 0.10 %(Rrs(j,k,2) - Rrs(j,k,4))/(Rrs(j,k,2)+Rrs(j,k,4)) < 0
                                    if Rrs(j,k,4) > Rrs(j,k,3) || Rrs(j,k,5) > Rrs(j,k,3)
                                        map(j,k) = 53; % Soft bottom
                                    elseif sum(Rrs(j,k,3:5)) > avg_water_sum && (Rrs(j,k,5) - Rrs(j,k,2))/(Rrs(j,k,5) + Rrs(j,k,2)) > 0.1
                                        map(j,k) = 52; % Soft bottom
                                    elseif Rrs(j,k,4) > Rrs(j,k,2) && (Rrs(j,k,3) - Rrs(j,k,6))/(Rrs(j,k,3) + Rrs(j,k,6)) < 0.60 % Separate seagrass from dark water
                                        if (Rrs(j,k,3) - Rrs(j,k,5))/(Rrs(j,k,3) + Rrs(j,k,5)) > 0.10 % Separate seagrass from turbid water
                                            map(j,k) = 54; % Seagrass
                                        else map(j,k) = 55; % Turbid water
                                        end
                                    else map(j,k) = 51; % Deep water
                                    end
                                else map(j,k) = 51; % Deep water
                                end

                             end % if v>u
                       end % If water/land
                   end % If isnan
               end % k
%                 if j == szA(1)/4
%                     update = 'DT 25% Complete'
%                 end
%                 if j == szA(1)/2
%                     update = 'DT 50% Complete'
%                 end
%                 if j == szA(1)/4*3
%                     update = 'DT 75% Complete'
%                 end
	        end % j

%Classes:
% 1 = Developed
% 2 = Vegetation
% 3 = Soil/sand/beach
% 41 = Deep water
% 42 = Benthic Sand
% 43 = Benthic Seagrass
% 44 = Benthic Coral
% 45 = Benthic patch coral

%%      DT Filter
        if filter > 0
                dt_filt = DT_Filter(map,filter,sz(1),sz(2));
                AA = [loc_out,id,'_',loc,'_Map_filt_',num2str(filter),'_benthicnew'];
                geotiffwrite(AA,dt_filt,R(1,1),'CoordRefSysCode',coor_sys);
        else
            Z1 = [loc_out,id,'_',loc,'_Map_benthicnew'];
            geotiffwrite(Z1,map,R(1,1),'CoordRefSysCode',coor_sys);
        end


        %% Output images
%          Z = [loc_out,id,'_',loc,'_Bathy1'];
% 	     geotiffwrite(Z,Bathy,R(1,1),'CoordRefSysCode',coor_sys);
	     Z2 = [loc_out,id,'_',loc,'_rrssub']; % last=52
         geotiffwrite(Z2,Rrs,R(1,1),'CoordRefSysCode',coor_sys);
% 
    end % If dt = 1
   end % If dt>0
end
   

	wtime = toc;
	time_min = wtime/60;
	fprintf(1,'Matlab CPU time (minutes) = %f\n', time_min);


