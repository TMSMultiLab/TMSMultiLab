% function filtered=EMG_filter(data, samplehz [,order] [,lowhz] [,highhz])
% data = columns of emg data
% samplehz = sampling frequency of data
% order = of filter; default = 1 (2 x 1st order)
% lowhz = high-pass cut-off; default=20hz
% highhz = low-pass cut-off; default=500hz
% uses: filtfilt

%% version history
% 27th June 2025 (after discussion with Justin Andrushko & Jason DeFreitas)
%	changed defaults to 20-500Hz
% 	added order as an option, default to 1st

function filtered=EMG_filter(data,samplehz,order,lowhz,highhz)
    if nargin==2
        order=1;
        lowhz=20;
        highhz=500;
    end
    if nargin==3
        lowhz=20;
        highhz=500;
    end 
    if nargin==4
        if isempty(order)
            order=1;
        end
        if isempty(lowhz)
            lowhz=20;
        end
        if isempty(highhz)
            highhz=500;
        end
    end
    if lowhz==0
        lowhz=0.01;
    end
    if highhz==0
        highhz=samplehz.*0.499;                                     % butter does not like exactly zero as low cutoff
    end
    if highhz==samplehz./2
        highhz=highhz.*.999;                                        % butter does not like exactly half as high cutoff
    end
    [b,a]=butter(order,[lowhz highhz]./(samplehz./2));
    filtered=filtfilt(b,a,data);
end