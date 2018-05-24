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



function dt_filt = WV2_Processing(images,id,met,crd_sys,dt,sgwin,filt,stat,loc,idnumber,rrs_out,class_out)

tic
d_t = str2num(dt)
sgw = str2num(sgwin)
n = num2str(idnumber)
id
met

%% Assign input and output locations

coor_sys = crd_sys % Change coordinate system code here
filter = str2num(filt); % None (0), 3x3 (3), 5x5 (5)
stat = str2num(stat); % Mean, Median, Mode (3)
loc = loc; % Typically the estuary acronym

loc_out = rrs_out;
class_out = class_out;
%loc_in = ['/work/m/mjm8/output/'];
%loc_out = ['/work/m/mjm8/tmp/test/output/'];
%class_out = ['/work/m/mjm8/tmp/test/output/'];

% Assign constants for all images
ebw = 0.001*[47.3 54.3 63.0 37.4 57.4 39.3 98.9 99.6]; % Effective Bandwidth per band (nm converted to um units; from IMD metadata files)
irr = [1758.2229 1974.2416 1856.4104 1738.4791 1559.4555 1342.0695 1069.7302 861.2866]; % Band-averaged Solar Spectral Irradiance (W/m2/um units)
cw = [.4273 .4779 .5462 .6078 .6588 .7237 .8313 .9080]; % Center wavelength (used for Rayleigh correction; from Radiometric Use of WorldView-2 Imagery)
gamma = 0.01*[1.499 1.471 1.442 1.413 1.413 1.413 1.384 1.384]; % Factor used in Rayleigh Phase Function equation (Bucholtz 1995)

    [A, R] = geotiffread(images);
    szA = size(A);
     s = xml2struct(met);
%    save XMLtest.mat s
        % Extract calibration factors and acquisition time from metadata for each band
	    if isfield(s,'IMD') == 0
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
	         kf(8,1) = str2num(s.Children(2).Children(idxb(10)).Children(26).Children.Data)
	        aqyear = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(1:4))
	        aqmonth = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(6:7))
	        aqday = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(9:10))
	        aqhour = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(12:13))
	        aqminute = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(15:16))
	        aqsecond = str2num(s.Children(2).Children(idxb(11,2)).Children(16).Children.Data(18:26))
	        sunel = str2num(s.Children(2).Children(idxb(11,2)).Children(56).Children.Data)
	        sunaz = str2num(s.Children(2).Children(idxb(11,2)).Children(50).Children.Data)
	        satview = str2num(s.Children(2).Children(idxb(11,2)).Children(86).Children.Data)
	        sensaz = str2num(s.Children(2).Children(idxb(11,2)).Children(62).Children.Data)
	        satel = str2num(s.Children(2).Children(idxb(11,2)).Children(68).Children.Data)
		cl_cov = str2num(s.Children(2).Children(idxb(11,2)).Children(90).Children.Data)

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
	         kf(8,1) = str2num(s.isd.IMD.BAND_N2.ABSCALFACTOR.Text)
	         aqyear = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(1:4)) % Extract Acquisition Time from metadata
	         aqmonth = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(6:7)) % Extract Acquisition Time from metadata
    	   	 aqday = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(9:10)) % Extract Acquisition Time from metadata
           	 aqhour = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(12:13)) % Extract Acquisition Time from metadata
             aqminute = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(15:16)) % Extract Acquisition Time from metadata
	    	 aqsecond = str2num(s.isd.IMD.IMAGE.FIRSTLINETIME.Text(18:26)) % Extract Acquisition Time from metadata
	    	 sunel = str2num(s.isd.IMD.IMAGE.MEANSUNEL.Text) % Extract Mean Sun Elevation angle from metadata
            satview = str2num(s.isd.IMD.IMAGE.MEANOFFNADIRVIEWANGLE.Text) % Extract Mean Off Nadir View angle from metadata
	         sunaz = str2num(s.isd.IMD.IMAGE.MEANSUNAZ.Text)
            sensaz = str2num(s.isd.IMD.IMAGE.MEANSATAZ.Text)
            satel = str2num(s.isd.IMD.IMAGE.MEANSATEL.Text)
	    cl_cov = str2num(s.isd.IMD.IMAGE.CLOUDCOVER.Text)
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
	    ESd = 1.00014 - 0.01671*cosd(degs) - 0.00014*cosd(2*degs) % Earth-Sun distance at given date (should be between 0.983 and 1.017)
	
	    TZ = cosd(90.0 - sunel); % Atmospheric spectral transmittance in solar path with solar zenith angle
	    TV = cosd(satview); % Atmospheric spectral transmittance in view path with satellite view angle

	    %% Calculate Rayleigh Path Radiance (Dash et al. 2012 and references therein)
	    if sunaz > 180 % For the following equations, azimuths should be between -180 and +180 degrees
	        sunaz = sunaz - 360;
	    end
	    if sensaz > 180
	        sensaz = sensaz - 360;
	    end
	    
	    az = abs(sensaz - 180 - sunaz); % Relative azimuth angle
	    thetaplus = acos(cosd(90-sunel)*cosd(90-satel) - sind(90-sunel)*sind(90-satel)*cosd(az)); % Scattering angles

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


        % Adjust file size: Input file (A) warped may contain more or fewer columns/rows than original NITF file, and some may be corrupt.
        sz(1) = min(szA(1),szB(1));
        sz(2) = min(szA(2),szB(2));
        sz(3) = 8;

            %% Assign NaN to no-data pixels and radiometrically calibrate and convert to Rrs
%            rad_cal{1,1} = single(zeros(szA(1),szA(2),szA(3))); % Create single-precision empty matrix for radiance image
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
            Z = [loc_out,id,'_',loc,'_Rrs'];
            geotiffwrite(Z,Rrs,R(1,1),'CoordRefSysCode',coor_sys);

	if d_t > 0; % Run DT and/or rrs conversion; otherwise end

	    % Calculate Kd (water column attenuation coefficient) from Chuanmin Hu's Rrs_Kd_Model.xlsx sheet
		sunzen = 90.0-sunel;
		c1 = 0.005;
		c2 = 4.18;
		c3 = 0.52;
		c4 = 10.8;
		v1 = [0.22024 0.142972 0.099157 0.286342 0.443809 1.491289 2.276262 2.223947]; % at (total absorption)
		v2 = [0.023277868 0.020883561 0.018975346 0.017605058 0.016771098 0.015875283 0.014661438 0.013602014]; % bb (backscatter)

		for b = 1:8
			Kd(b) = single((1+c1*sunzen)*v1(b)+c2*(1-c3*exp(-c4*v1(b)))*v2(b));
		end
	
		Kd

%	E_glint = [0.829 0.893 0.941 0.958 0.979 0.984]; % ENDir(band)/ENDir(NIR) from Martin et al. 2016


	    %% Setup for Deglint, Bathymetry, and Decision Tree
	    t = 1;
	    u = 1;
	    v = 0;
	    num_pix = 0;
	    sum_veg(t) = 0;
	    dead_veg(t) = 0;
	    sz_ar = sz(1)*sz(2);
	    water = zeros(sz_ar,8);
	    for j = 1:sz(1);
	        for k = 1:sz(2);
	            if isnan(Rrs(j,k,1)) == 0
			num_pix = num_pix +1; % Count number of non-NaN pixels
			c_val(num_pix) = Rrs(j,k,1); % Record coastal band value for use in cloud mask prediction
	                if (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.6 && Rrs(j,k,7) > Rrs(j,k,3); % Identify vegetation (excluding grass)
	                    if ((Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2))) > 0.20; % Shadow filter
	                        sum_veg(t) = sum(Rrs(j,k,3:5)); % Sum bands 3-5 for selected veg to distinguish wetland from upland
				dead_veg(t) = (((Rrs(j,k,7) - Rrs(j,k,4))/3) + Rrs(j,k,4)) - Rrs(j,k,5); % Compute difference of predicted B5 value from actual valute
	                        t = t+1;
	                    end
			elseif Rrs(j,k,8) < 0.11 && Rrs(j,k,1) > 0 && Rrs(j,k,2) > 0 && Rrs(j,k,3) > 0 && Rrs(j,k,4) > 0 && Rrs(j,k,5) > 0 && Rrs(j,k,6) > 0 && Rrs(j,k,7) > 0 && Rrs(j,k,8) > 0; % Identify water
				water(u,:) = Rrs(j,k,:);
				u = u+1;
				if (Rrs(j,k,4) - Rrs(j,k,8))/(Rrs(j,k,4) + Rrs(j,k,8)) < 0.3 % NDGI to identify glinted water pixels
					v = v+1;
				end
	                end
	            end
	        end
	    end
		n_Eglint = u % Print number of water pixels used to derive E_glint relationships
		n_glinted = v % Print number of glinted water pixels

		idx = find(water(:,1) == 0);
		water(idx,:) = [];
                
		if v > 0.33*u
			Update = 'Deglinting'
                	for b = 1:6 %% Calculate linear fitting of all MS bands vs NIR1 & NIR2 for deglinting in DT
                        	if b == 1 || b ==4 || b == 6
                                	slope = water(:,b)\water(:,8);
               		        else slope = water(:,b)\water(:,7);
                        	end
                	E_glint(b) = slope;
                	end
			E_glint
		else Update = 'Glint-free'
		end

	    avg_veg_sum = mean(sum_veg(:))
	    avg_dead_veg = mean(dead_veg(:))

%	    num_pix = szB(1)*szB(2);
	    num_pix
	    if cl_cov > 0
		    num_cld_pix = round(num_pix*(cl_cov*0.01)); % Number of cloud pixels (rounded down to nearest integer) based on metadata-reported percent cloud cover
		    srt_c = sort(c_val,'descend'); % Sort all pixel blue-values in descending order. Cloud mask threshold will be num_cld_pix'th highest value
		    cld_mask = srt_c(num_cld_pix); % Set cloud mask threshold
	    else cld_mask = max(c_val)+1;
	    end


	    Bathy{1,1} = zeros(szA(1),szA(2),1); % Preallocate for Bathymetry
	    Rrs_deglint = single(zeros(6,1)); % Preallocate for deglinted Rrs
	    Rrs_wc = single(zeros(1,6)); %Preallocation for water-column corrected Rrs
        wn_sz = (sgw*2+1)^2;
	    block_sz = [(sz(1)/10) (sz(2)/10)];

	if d_t == 1; % Execute Deglinting Rrs, Bathymetry, and Decision Tree
	    update = 'Running DT'
	    output{1,1} = zeros(szA(1),szA(2),1); % Create empty matrix for classification output
	   for j = sgw+1:sz(1)-sgw-1; % Start at sunglint moving window +/-1 pixel in for sunglint filter
	        for k = sgw+1:sz(2)-sgw-1;
	            if isnan(Rrs(j,k,1)) == 1
	                if filter > 0
	                    output{1,1}(j,k,1) = NaN; % NaN when using filter
	                else output{1,1}(j,k,1) = 0; % 0 when not using filter
	                end
                    elseif Rrs(j,k,1) >= cld_mask; % Cloud filter
                        output{1,1}(j,k,1) = 1; % Cloud
		    elseif Rrs(j,k,8) < 0.1;
    			output{1,1}(j,k,1) = 3; % Water
	            elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.3 && Rrs(j,k,7) > Rrs(j,k,3); % Identify vegetation (including grass; band 7>3 criterion excludes sunglinted pixels)
	                if ((Rrs(j,k,7) - Rrs(j,k,2))/(Rrs(j,k,7) + Rrs(j,k,2))) < 0.20 && (Rrs(j,k,7) - Rrs(j,k,8))/(Rrs(j,k,7) + Rrs(j,k,8)) > 0.01; % Shadowed-vegetation filter (B7/B8 ratio excludes marsh, which tends to have very similar values here)
	                    output{1,1}(j,k,1) = 0; % Shadow
	                elseif sum(Rrs(j,k,3:5)) < avg_veg_sum;
	                        if ((Rrs(j,k,2) - Rrs(j,k,5))/(Rrs(j,k,2) + Rrs(j,k,5))) < 0.4; % Agriculture filter based on elevated Blue band values
            					output{1,1}(j,k,1) = 8; % Forested Wetland
	                        else output{1,1}(j,k,1) = 4; % Forested Upland (most likely agriculture)
	                        end
			elseif ((Rrs(j,k,7) - Rrs(j,k,4))/3 + Rrs(j,k,4) - Rrs(j,k,5)) < avg_dead_veg
				output{1,1}(j,k,1) = 7; % Dead vegetation
			elseif (Rrs(j,k,8) - Rrs(j,k,5))/(Rrs(j,k,8) + Rrs(j,k,5)) > 0.65; % NDVI for high upland values
				output{1,1}(j,k,1) = 5; % Upland Forest/Grass;
%			elseif (((Rrs(j,k,7) - Rrs(j,k,4))/3)+Rrs(j,k,4)) - Rrs(j,k,5) > 0.035; % Difference of B5 from predicted B5 by slope of B7:B4 to distinguish live vs dead trees/grass/marsh
			else	output{1,1}(j,k,1) = 6; % Live Marsh
	                end


%			f = 1;
%			for p = -sgw:sgw; % Moving-window sunglint extraction (11x11 = 65 minutes)
%				for q = -sgw:sgw;
%					if (Rrs(j+p,k+q,8)) < 0.1; % Record NIR values of adjacent water pixels
%						G(f,1) = Rrs(j+p,k+q,7); % Use NIR1 for array 1, NIR2 for array 2 due to 0.2 second collection gap
%						G(f,2) = Rrs(j+p,k+q,8);
%						f = f+1;
%					end
%				end
%			end
%			mnNIR1 = mean(G(:,1));
%			mnNIR2 = mean(G(:,2));
%			for p = 1:6
%				if p == 1 || p == 4 || p == 6;
%					Rrs_deglint(p) = Rrs(j,k,p) - (E_glint(p).*(Rrs(j,k,8)-mnNIR2)); % Deglint algorithm from Martin et al. 2016 revised to use linear fit of MS bands vs NIR
%				else Rrs_deglint(p) = Rrs(j,k,p) - (E_glint(p).*(Rrs(j,k,7)-mnNIR1)); % Use linear fit of MS bands vs NIR1 or NIR2 
%				end
%			end
%			Z = log(1000*Rrs_deglint(2))/log(1000*Rrs_deglint(3)); % Calculate relative depth from deglinted bands (Stumpf 2003 ratio transform)
%			Bathy{1,1}(j,k,1) = Z;
%			for d = 1:6; % Change Rrs values for water pixels
%				Rrs(j,k,d) = Rrs_deglint(d);
%				Rrs(j,k,d) = Rrs(j,k,d)/exp(-2*Kd(d)*Z); % Calculate water-column corrected benthic reflectance
%				Rrs(j,k,d) = Rrs_deglint(d)/exp(-2*Kd(d)*Z); % Calculate water-column corrected benthic reflectance (Traganos 2017 based on Maritorena 1994)
%			end
%			if (Rrs(j,k,6) - Rrs(j,k,5))/(Rrs(j,k,6) + Rrs(j,k,5)) > 0.68; % Seagrass NDVI algorithm
%				output{1,1}(j,k,1) = 7; % Seagrass
%			else output{1,1}(j,k,1) = 3; % Water
%			end
%                        if (Rrs(j,k,3) - Rrs(j,k,8))/(Rrs(j,k,3) + Rrs(j,k,8)) < 0.3;
%                           output{1,1}(j,k,1) = 0; % Shadow filter #2
%                        else output{1,1}(j,k,1) = 3; % Water
%			 end
	            else output{1,1}(j,k,1) = 2; % Bare/developed
	            end
	        end
	    end
            if filter > 0
                filter
                dt_filt = DT_Filter(output{1,1},filter,stat,sz(1),sz(2));
%               dt_filt_8bit = im2uint8(dt_filt{1,1});
                AA = [class_out,id,'_',loc,'_DT_filt_',n,'_',num2str(filter),'_',num2str(stat)];
%               geotiffwrite(AA,dt_filt_8bit,R(1,1),'CoordRefSysCode',coor_sys);
                geotiffwrite(AA,dt_filt{1,1},R(1,1),'CoordRefSysCode',coor_sys);
            else
                BB = [class_out,id,'_',loc,'_DT_nofilt_',n];
                geotiffwrite(BB,output{1,1},R(1,1),'CoordRefSysCode',coor_sys);
            end
            clear output
	elseif d_t == 2; % Only run for Deglinted Rrs and Bathymetry, not Decision Tree
	    update = 'Not Running DT' % Print to confirm DT is not being used
%	    Rrs1 = single(zeros(szA(1),szA(2),szA(3)));
%	    block = [szA(1)/5 szA(2)/5];
%	    myfunc = @(block_struct) WV2_block_rrs(Rrs,sz,sgw,cld_mask,E_glint,Kd);
%	    Rrs = blockproc(Rrs,[block(1) block(2)],myfunc,'UseParallel',true);
      	    for j = sgw+1:sz(1)-sgw-1; % Start at sunglint moving window +/-1 pixel in for sunglint filter
               for k = sgw+1:sz(2)-sgw-1;
                  if isnan(Rrs(j,k,1)) == 0 &&  Rrs(j,k,1) < cld_mask && Rrs(j,k,8) < 0.11; % && Rrs(j,k,8) < Rrs(j,k,3); % Only run on non-NaN cloud-free pixels containing water
        			if v > u*0.33 && (Rrs(j,k,4) - Rrs(j,k,8))/(Rrs(j,k,4) + Rrs(j,k,8)) < 0.3 % Run deglinting algorithms if more than 1/3 of water pixels contain glint and target pixel is glinted
        				mw1 = Rrs(j-sgw:j+sgw,k-sgw:k+sgw,7);
                		mw2 = Rrs(j-sgw:j+sgw,k-sgw:k+sgw,8);
        				mnNIR1 = single(min(min(mw1(mw1>0))));%
        				mnNIR2 = single(min(min(mw2(mw2>0))));
        				clear mw1 mw2

    %				f = 1;
    %	                	for p = -sgw:sgw; % Moving-window sunglint extraction
    %	                               for q = -sgw:sgw;
    %	                              	    if (Rrs(j+p,k+q,8)) < 0.11 && Rrs(j+p,k+q,7) > 0 && Rrs(j+p,k+q,8) > 0; % Record NIR values of adjacent water pixels
    %		                                G(f,1:2) = Rrs(j+p,k+q,7:8); % Use NIR1 for array 1, NIR2 for array 2 due to 0.2 second collection gap
    %	                	                f = f+1;
    %	                                   end
    %			                end
    %	                	end
    %			        mnNIR1 = single(min(G(:,1)));
    %	                        mnNIR2 = single(min(G(:,2)));
    %				clear G


                        % Deglint equation
                        for d = 1:6
                            dg2(d) = single(Rrs(j,k,d) + (E_glint(d)*(Rrs(j,k,8) + mnNIR2)));
                            dg1(d) = single(Rrs(j,k,d) + (E_glint(d)*(Rrs(j,k,7) + mnNIR1)));
                        end
                        Rrs_deglint = [dg2(1) dg1(2) dg1(3) dg2(4) dg1(5) dg2(6)];
                        clear dg1 dg2

        %				Rrs_0 = Rrs_deglint./(0.52 + 1.7.*Rrs_deglint);
                        dp = single(log(1000*Rrs_deglint(2))/log(1000*Rrs_deglint(3)));
                        for d = 1:6
                            Rrs(j,k,d) = Rrs_deglint(1,d)./(2.7183.^(2*Kd(1,d)*dp)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
                        end
                    else % For glint-free/low-glint images
                        % Subsurface rrs and Depth equations
%         				Rrs_0 = Rrs(j,k,:)./(0.52 + 1.7.*Rrs(j,k,:)); % Convert above-surface Rrs to below-surface rrs (Lee et al. 2002)
        		        dp = single(log(1000*Rrs(j,k,2))/log(1000*Rrs(j,k,3))); % Calculate relative depth from deglinted bands (Stumpf 2003 ratio transform)
           				for d = 1:6
               				Rrs(j,k,d) = Rrs(j,k,d)./(2.7183.^(2*Kd(1,d)*dp)); % Calculate water-column corrected benthic reflectance (Traganos 2017 & Maritorena 1994)
                        end
                    end % If v>u
                    clear dp
                  end % If isnan
    		end
	    end
            %% Output images
%            Z = [loc_out,id,'_',loc,'_Bathy'];
%	     geotiffwrite(Z,Bathy{1,1},R(1,1),'CoordRefSysCode',coor_sys);
	     Z2 = [loc_out,id,'_',loc,'_rrs'];
             geotiffwrite(Z2,Rrs,R(1,1),'CoordRefSysCode',coor_sys);
    end % If dt = 1
   end % If dt>0

	clear Rrs

%	    data_8bit = im2uint8(output{1,1}); % Convert to 8-bit unsigned integers for smaller output filesize

	wtime = toc;
	time_min = wtime/60;
	fprintf(1,'Matlab CPU time (minutes) = %f\n', time_min);

end
 
