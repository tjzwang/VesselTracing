clc; clear all; close all;

%% Inputs

input = '/Users/zacharyhoglund/Downloads/REAL_Piecewise_tau_intensity_binned_data_MEGASHEET_with_headers.xlsx';

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
outputfilename = '/Users/zacharyhoglund/Downloads/Perc_Bins_Piecewise_intensty_10micron_bins_to_100microns_3MicronSurface.csv';

samples = [2392 2302 2290 2131 2267 24632 23642 2477 2470 24652 2417 2399];

% Percentile bin width:

perc_binW = [ 10 ];

%% Create new bins

data = readmatrix(input);

final_data(:,1:6) = data(:,1:6);

surface_data = data(:,(6+11+surf_low):(6+11+surf_high));

final_data(:,7) = mean(surface_data,2);

data_backup = data;

%%

data = final_data;

for i = 1:length(samples)

    idx = find(data(:,1) == samples(i));

    temp_data = data(idx,:);

    edges = 0:(100/perc_binW):100;

   ADpercentiles = invprctile(temp_data(:,7),temp_data(:,7));
   
   ADbin = discretize(ADpercentiles,edges);

   temp_data = [temp_data ADbin];

   if i == 1

   sample_data(1:length(temp_data(:,1)),:) = temp_data;    

   else

       sample_data(end+1:end+length(temp_data(:,1)),:) = temp_data;
   
   end

end

final_data = sample_data;

data = data_backup;

%%

num_bins = (max_dist-min_dist)/BinW; 

for i = 1:num_bins

    bindata = data(:,(6+11+min_dist+BinW*(i-1)):(6+11+BinW*(i)));

    binmean = mean(bindata,2);
    
    final_data(:,8+i) = binmean;

end

%% Output binned data

csvwrite(outputfilename,final_data);


