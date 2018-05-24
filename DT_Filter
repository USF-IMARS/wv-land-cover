%% DT_Filter.M
%% Written by Matt McCarthy 8/29/2016

function dt_filt = DT_Filter(file,x,y,sz2,sz3)
filt = x; % 3x3 or 5x5 filter
stat = y; % Mean (1), Median (2), Mode (3)
sz_sm(1) = sz2; % Size of unwarped(smaller) file
sz_sm(2) = sz3;

sz1 = size(file);
dt_filt{1,1} = zeros(sz1(1),sz1(2),1);

if filt == 3;    
    for a = 2:sz_sm(1)-1; % Mode filter or median filter: 3x3 or 5x5
        for b = 2:sz_sm(2)-1;
            if isnan(file(a,b)) == 0;
		 C = [file(a,b) file(a-1,b-1) file(a-1,b) file(a-1,b+1) file(a,b-1) file(a,b+1) file(a+1,b) file(a+1,b-1) file(a+1,b+1)];% ...
                idx = find(C == 0); % If any pixels in C are shadows, they are not included in the mode function
                C(idx) = [];
                if stat == 3;
                    mod = mode(C); % Identify most common value (if more than one value, lower value is selected automatically)
                    if isnan(mod) == 1; % Check if mode of box is NaN (redundancy)
                        dt_filt{1,1}(a,b) = 0; % If NaN, assign zero (Arc won't load DT tiffs w/ NaNs)
                    elseif file(a,b) == 6; % If mode of C indicates wetland, check that wetlands comprise at least 2/3 of adjacent vegetation pixels, otherwise assign upland
                        idx2 = C == 6; % This is justified by the homogeneity of wetland vegetation while upland often occurs as individual stands
                        idx3 = C == 4; % Upland
                        if sum(idx2) >= (2/3)*(sum(idx2) + sum(idx3))
                            dt_filt{1,1}(a,b) = 6; % Wetland
                        else dt_filt{1,1}(a,b) = 4; % Upland
                        end
                    else dt_filt{1,1}(a,b) = mod; % If mod is upland, marsh, water or bare/developed, assign it as such
                    end
                elseif stat ==2;
                    med = median(C);
                    if round(med) ~= mod; % Check if median of 3x3 box is a decimal (e.g 2.5)
                        dt_filt{1,1}(a,b) = file(a,b); % If decimal, assign original number (no filter)
                    else
                        dt_filt{1,1}(a,b) = med; % If not decimal, assign median value (filtered)
                    end
                end
            else
                dt_filt{1,1}(a,b) = 0;
            end
        end
    end
elseif filt == 5;
   	for a = 6:sz_sm(1)-6;
      		for b = 6:sz_sm(2) -6;
			if isnan(file(a,b)) == 0;
				f = 1;
				for d = -5:5
			        	for e = -5:5;
						C(f) = file(a+d,b+e);
						f = f+1;
					end
				end
				idx = find(C == 0); % If any pixels in C are shadows, they are not included in the mode function
			        C(idx) = [];
			        if stat == 3;
			            mod = mode(C); % Identify most common value (if more than one value, lower value is selected automatically)
			            if isnan(mod) == 1; % Check if mode of box is NaN (redundancy)
			                dt_filt{1,1}(a,b) = 0; % If NaN, assign zero
			            elseif file(a,b) == 6; % If mode of C indicates wetland, check that wetlands comprise at least 2/3 of adjacent vegetation pixels, otherwise assign upland
			                idx2 = C == 6; % This is justified by the homogeneity of wetland vegetation while upland often occurs as individual stands
	        		        idx3 = C == 4;
					idx4 = C == 2;
	        		        if sum(idx2) >= (2/3)*(sum(idx2) + sum(idx3) + sum(idx4))
	        		            dt_filt{1,1}(a,b) = 6; % Wetland
	        		        else dt_filt{1,1}(a,b) = 4; % Upland
	        		        end
	        		    else dt_filt{1,1}(a,b) = mod;
	        		    end
	        		end
			else dt_filt{1,1}(a,b) = 0;
			end
		end
	end
end
