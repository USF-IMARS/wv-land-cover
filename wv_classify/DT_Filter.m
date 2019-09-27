%% DT_Filter.M
%% Written by Matt McCarthy 8/29/2016

function dt_filt = DT_Filter(file,x,sz2,sz3)
filt = x;
sz_sm(1) = sz2; % Size of unwarped(smaller) file
sz_sm(2) = sz3;

sz1 = size(file);
dt_filt{1,1} = zeros(sz1(1),sz1(2),1);

    for a = 2:sz_sm(1)-1; % Mode filter or median filter: 3x3 or 5x5
        for b = 2:sz_sm(2)-1;
            if isnan(file(a,b)) == 0;
		C = file(a-filt:a+filt,b-filt:b+filt);
                idx = find(C == 0); % If any pixels in C are shadows, they are not included in the mode function
                C(idx) = [];
                    mod = mode(C); % Identify most common value (if more than one value, lower value is selected automatically)
                    if isnan(mod) == 1; % Check if mode of box is NaN (redundancy)
                        dt_filt{1,1}(a,b) = 0; % If NaN, assign zero (Arc won't load DT tiffs w/ NaNs)
%                    elseif file(a,b) == 6; % If mode of C indicates wetland, check that wetlands comprise at least 2/3 of adjacent vegetation pixels, otherwise assign upland
%                       idx2 = C == 6; % This is justified by the homogeneity of wetland vegetation while upland often occurs as individual stands
%                      idx3 = C == 4; % Upland
%                        if sum(idx2) >= (2/3)*(sum(idx2) + sum(idx3))
%                            dt_filt{1,1}(a,b) = 6; % Wetland
%                        else dt_filt{1,1}(a,b) = 4; % Upland
%                        end
                    else dt_filt{1,1}(a,b) = mod; % If mod is upland, marsh, water or bare/developed, assign it as such
                    end
                end
            else
                dt_filt{1,1}(a,b) = 0;
            end
        end
    end
end
