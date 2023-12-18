clc; clear all; close all;

%% Inputs

input = "C:\Users\exx\Desktop\Zach - Active Analysis\FINAL DATA OUTPUT\Tangle_MEGASHEET\tau_tangles_binned_data_MEGASHEET 2.csv";

% bin width away from the blood vessel

BinW = 10;

% max distance from blood vessel

max_dist = 30;

% min distance from blood vessel

min_dist = 0;

% surface data definition

% lower bound

surf_low = 0;

% upper bound

surf_high = 3;

% Output file name:
outputfilename = 'C:\Users\exx\Desktop\Zach - Active Analysis\FINAL DATA OUTPUT\Tangle_MEGASHEET\tangle_10micron_bins_to_30microns_3MicronSurface.csv';


%% Create new bins

data = readmatrix(input);

final_data(:,1:6) = data(:,1:6);

surface_data = data(:,(6+11+surf_low):(6+11+surf_high));

final_data(:,7) = mean(surface_data,2);

num_bins = (max_dist-min_dist)/BinW; 

for i = 1:num_bins

    bindata = data(:,(6+11+min_dist+BinW*(i-1)):(6+11+BinW*(i)));

    binmean = mean(bindata,2);
    
    final_data(:,7+i) = binmean;

end

%% Output binned data

csvwrite(outputfilename,final_data);


