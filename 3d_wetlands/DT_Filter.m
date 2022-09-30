%% DT_Filter.M
%% Written by Matt McCarthy 8/29/2016

function dt_filt = DT_Filter(file,x,sz2,sz3,dev,FW,FU,UG,WA);
filt = x;
sz_sm(1) = sz2; % Size of unwarped(smaller) file
sz_sm(2) = sz3;

fwfilt = 75
wafilt = 50
sz1 = size(file);
dt_filt =zeros(sz1(1),sz1(2),'uint8');

    for a = filt+1:sz_sm(1)-filt-1; % Mode filter or median filter: 3x3 or 5x5
        for b = filt+1:sz_sm(2)-filt-1;
	    if file(a,b) == dev;
		dt_filt(a,b) = dev;
	    elseif file(a,b) == FW && a > fwfilt && b > fwfilt && a < sz_sm(1)-100 && b < sz_sm(2)-100; % For FW, use larger window to eliminate erroneous urban misclassifications
		C = file(a-fwfilt:a+fwfilt,b-fwfilt:b+fwfilt);
		idx = find(C == 0); % If any pixels in C are shadows, they are not included in the mode function
		C(idx) = [];
		mod = mode(mode(C)); % Mode of window (by definition, lower value is selected if more than one mode value)
		if isnan(mod) == 1;
			dt_filt(a,b) = 0; % If NaN, assign zero because Arc won't load DT tiffs w/ NaNs
		elseif mod == FW; % Check if mode FW is actually urban tree shadow
			idxfor = find(C == dev | C == FU); % Find upland forest, grass, and developed nearby
			if size(idxfor,1)>0.10*size(C,1)*size(C,2); % > 10% is upland, grass, or developed
				dt_filt(a,b) = FU; % Assumed to be non-wetland forest
			else dt_filt(a,b) = FW;
			end
		else dt_filt(a,b) = FU;
		end
            elseif isnan(file(a,b)) == 0;
		C = file(a-filt:a+filt,b-filt:b+filt);
                idx = find(C == 0); % If any pixels in C are shadows, they are not included in the mode function
                C(idx) = [];
                    mod = mode(mode(C)); % Identify most common value (if more than one value, lower value is selected automatically)
                    if isnan(mod) == 1; % Check if mode of box is NaN (redundancy)
                        dt_filt(a,b) = 0; % If NaN, assign zero (Arc won't load DT tiffs w/ NaNs)
		    elseif mod == FW;
			D = dt_filt(a-filt:a,b-filt:b);
			modD = mode(mode(D));
			dt_filt(a,b) = modD;
                    else dt_filt(a,b) = mod; % If mod is upland, marsh, water or bare/developed, assign it as such
                    end
            else
                dt_filt(a,b) = 0;
            end
        end
    end
end
