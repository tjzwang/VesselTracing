clc; clear all; close all;

%% Inputs

% Object Input
objinput = "C:\Users\exx\Desktop\Zach - Active Analysis\FINAL DATA OUTPUT\Neuron_MEGASHEET_100\neurons_binned_data_MEGASHEET 100.csv";

% Intensity Input

intinput = "C:\Users\exx\Desktop\Zach - Active Analysis\FINAL DATA OUTPUT\Intensity_MEGASHEET_100\tau_intensity_binned_data_MEGASHEET 100.csv";

% bin width away from the blood vessel

BinW = 10;

% max distance from blood vessel

max_dist = 100;

% min distance from blood vessel

min_dist = 0;

% surface data definition

% lower bound

surf_low = 0;

% upper bound

surf_high = 3;

% Output file name:
outputfilename = 'C:\Users\exx\Desktop\Zach - Active Analysis\FINAL DATA OUTPUT\Intensity_MEGASHEET_100\neurons_10micron_bins_to_100microns_3MicronSurface.csv';


%% Create new bins

objdata = readmatrix(objinput);

intdata = readmatrix(intinput);

[C,iInt,iObj] = intersect(intdata(:,1:6),objdata(:,1:6),'rows');

intdata2 = intdata(iInt,:);

objdata = objdata(iObj,:);

final_data(:,1:6) = objdata(:,1:6);

surface_data = intdata2(:,(6+11+surf_low):(6+11+surf_high));

final_data(:,7) = mean(surface_data,2);

num_bins = (max_dist-min_dist)/BinW; 

for i = 1:num_bins

    bindata = objdata(:,(6+11+min_dist+BinW*(i-1)):(6+11+BinW*(i)));

    binmean = mean(bindata,2);
    
    final_data(:,7+i) = binmean;

end


%% Output binned data

csvwrite(outputfilename,final_data);


