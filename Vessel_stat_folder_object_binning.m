clc; clear all; close all;

% Exported Data Format: [ Vessel#  Frame#  MinFrameDistance MaxFrameDistance  ObjectDensity  ]

%% Import Data

% .txt file (Recommended for speed):

% Input data folder name (INCLUDE FOLDER PATH AND NAME):
sourcefolder = '/Users/zacharyhoglund/Documents/test data';

% Max Distance From Vessel (in microns) (Recommended Value 30):

ves_dist = [ 30 ];

min_dist = [ -10 ];

% Resolution in microns? (Recommended Value: 1):

res = [ 1 ];

% Frame width (# Microns):

frame_w = [ 10 ];

% Distance to define surface data (in microns)? (Recommended Value: 3):

surf_d = [ 3 ];

% Data Output Name (True, unshuffled data only) (Include Name + Path):
% Vessel number and data type will be appended, so do not include in names.
% NEED FOLDER WITH SAME NAME FOR EACH DATA TYPE: XDATA, YDATA, AND ERRORS
% ex. you need folders stat_data_xdata, stat_data_ydata, and
% stat_data_errors if stat_data is the folder name.

% Folder:
outputfoldername = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2302\Stat_data_imp_fit_ind_Percent';
% File:
outputname = 'Human2302_Tau_obj_test_data';

% % Set number of bins by increment of intensity
% nbins = round((max(Tau_data(:,4)) - min(Tau_data(:,4)))/50);

%% Load Data

% # Bins along vessel
nbins_length = 5;

f = waitbar(0,'Calculating Statistics: Reading Data 0 %');

fcontent = dir(fullfile(sourcefolder, '*.txt')); %fcontent is a column vector of structures

l1 = 1;

% Open files and load data from folder
% If multiple files exist, their data will be stacked

Tau_data = [];

n = length(fcontent');

b = 1;

for b = 1:n
    sizes = [];
        n1 = b;
        n2 = b;
        Tau_data = [];
        bin_edges = [];
        surface_data = [];
        bin_numbers = [];
        clc
        % close all
for i = n1:n2

%i = 7;

   filename = fcontent(i).name;
   Tau_data_single = [];
   Tau_data_single = readmatrix(fullfile(sourcefolder,filename));
   [r, c] = size(Tau_data_single);
   l2 = l1 + r - 1;
   Tau_data(l1:l2,:) = Tau_data_single;
   l1 = 1;

    waitbar(0,f,sprintf('Calculating Statistics: Reading Data %d %%',floor(i/n*100)));

end


waitbar(0.1,f,'Calculating Statistics: Preparing Data');

% Remove data where Tau is within blood vessel and greater than ves_dist
% away
Tau_data_backup = Tau_data;
idx = find( Tau_data(:,6) <= min_dist);
Tau_data(idx,:) = [];
idx = find(Tau_data(:,6) > ves_dist );
Tau_data(idx,:) = [];

%% Bin Data by Surface Frames

% b = vessel ID    g = frame ID
frame_data = zeros(ceil(max(Tau_data(:,5)/frame_w)),3);
g = 1;
for i = frame_w:frame_w:(ceil(max(Tau_data(:,5)/frame_w))*frame_w)
    
    idx = find(Tau_data(:,5) < i & Tau_data(:,5) >= (i - frame_w));
     %  [b g (i - frame_w) i Tau_data(idx,:)]
     if length(idx)>=1
     for d = 1:length(idx)
    binned_data{g}(d,:) = [b g (i - frame_w) i (Tau_data(idx(d),:))];
     end
        else
        binned_data{g}(d,:) = [b g (i - frame_w) i zeros(1,length(Tau_data(1,:)))];
     end
    g = g+1;

end

%%
% Get intensity at distances
if b == 1;
    new_data = [];
    endrow = 0;
end
new_data = [new_data;zeros(length(binned_data),((ves_dist-min_dist)+5))];

for i = 1:length(binned_data)

   single_bin_data = binned_data{i};
  
   if ~isempty(single_bin_data)
   new_data(endrow+i,1:4) = single_bin_data(1,1:4);
   else
   new_data(endrow+i,1:4) = [b    i   new_data(endrow+i-1,3)+frame_w      new_data(endrow+i-1,4)+frame_w];
   end

   for d = 0:(ves_dist-min_dist)
       dist_min = d+min_dist-0.5;
       dist_max = d+min_dist+res-0.5;
       if d == 0
       idx = find(single_bin_data(:,10) >= dist_min & single_bin_data(:,10)<= dist_max);
       mean_int = mean(single_bin_data(idx,8));
       if ~isnan(mean_int) 
       volume = length(idx)/((pi*(dist_max)^2 - pi*(dist_min)^2)*frame_w);
       density = length(idx)/volume;
       new_data(i,d+4) = density;
       end

       else
       idx = find(single_bin_data(:,10) > d+min_dist-0.5 & single_bin_data(:,10)<= dist_max);
       mean_int = mean(single_bin_data(idx,8));
       if ~isnan(mean_int) 
       volume = ((pi*(dist_max)^2 - pi*(dist_min)^2)*frame_w);
       density = length(idx)/volume;
       new_data(i,d+4) = density;
       end

       end
  
   end

end

endrow = length(new_data(:,1));
%binned_data = {};
end

%% Save Data
writematrix(new_data,append(outputfoldername,'/',outputname,'.txt'));
xlswrite(append(outputfoldername,'/',outputname),new_data);
