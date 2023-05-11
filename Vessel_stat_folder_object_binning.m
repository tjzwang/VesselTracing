clc; clear all; close all;

% Exported Data Format: [ Vessel#  Frame#  MinFrameDistance MaxFrameDistance  ObjectDensity  ]

%% Import Data

% .txt file (Recommended for speed):

% Input tau data folder name (INCLUDE FOLDER PATH AND NAME):
Tau_sourcefolder = '/Users/zacharyhoglund/Documents/test data';

% Input distance transform data folder name for volume calculation (INCLUDE FOLDER PATH AND NAME):
Distance_sourcefolder = '/Users/zacharyhoglund/Documents/test data';

% Max Distance From Vessel (in microns) (Recommended Value 30):

ves_dist = [ 30 ];

min_dist = [ -10 ];

% Resolution in microns? (Recommended Value: 1):

res = [ 1 ];

% Frame width (# Microns):

frame_w = [ 10 ];

% Distance to define surface data (in microns)? (Recommended Value: 3):

surf_d = [ 3 ];

% Voxel Size

% Voxel Size:  [   X    Y   Z   ]
voxel =        [ 0.32 0.32 1.75 ];

% Data Output Name (True, unshuffled data only) (Include Name + Path):
% Vessel number and data type will be appended, so do not include in names.
% NEED FOLDER WITH SAME NAME FOR EACH DATA TYPE: XDATA, YDATA, AND ERRORS
% ex. you need folders stat_data_xdata, stat_data_ydata, and
% stat_data_errors if stat_data is the folder name.

% Folder:
outputfoldername = '/Users/zacharyhoglund/Documents/test data output';

% File:
outputname = 'Human2302_Tau_obj_test_data';

% % Set number of bins by increment of intensity
% nbins = round((max(Tau_data(:,4)) - min(Tau_data(:,4)))/50);

%% Load Data

% # Bins along vessel
nbins_length = 5;

f = waitbar(0,'Calculating Statistics: Reading Data 0 %');

fcontent = dir(fullfile(Tau_sourcefolder, '*.txt')); %fcontent is a column vector of structures

l1 = 1;

Tau_data = [];

Distance_data = [];

n = length(fcontent');

%b = 1;

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
   Tau_data_single = readmatrix(fullfile(Tau_sourcefolder,filename));
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


% Read Distance Data

fcontent = dir(fullfile(Distance_sourcefolder, '*.txt')); %fcontent is a column vector of structures

l1 = 1;

for i = n1:n2

%i = 7;

   filename = fcontent(i).name;
   Distance_data_single = [];
   Distance_data_single = readmatrix(fullfile(Distance_sourcefolder,filename));
   [r, c] = size(Distance_data_single);
   l2 = l1 + r - 1;
   Distance_data(l1:l2,:) = Distance_data_single;
   l1 = 1;

    waitbar(0,f,sprintf('Calculating Statistics: Reading Data %d %%',floor(i/n*100)));

end


waitbar(0.1,f,'Calculating Statistics: Preparing Data');

% Remove data where Tau is within blood vessel and greater than ves_dist
% away
Distance_data_backup = Distance_data;
idx = find( Distance_data(:,6) <= min_dist);
Distance_data(idx,:) = [];
idx = find(Distance_data(:,6) > ves_dist );
Distance_data(idx,:) = [];

%% Bin Data by Surface Frames

% b = vessel ID    g = frame ID

g = 1;
for i = frame_w:frame_w:(ceil(max(Tau_data(:,5)/frame_w))*frame_w)
    
    idx = find(Tau_data(:,5) < i & Tau_data(:,5) >= (i - frame_w));
    idx2 = find(Distance_data(:,5) < i & Distance_data(:,5) >= (i - frame_w));
     %  [b g (i - frame_w) i Tau_data(idx,:)]
     if ~isempty(idx)
     for d = 1:length(idx)
    Tau_binned_data{g}(d,:) = [b g (i - frame_w) i (Tau_data(idx(d),:))];
    Distance_binned_data{g}(d,:) = [b g (i - frame_w) i (Distance_data(idx(d),:))];
     end
        else
        Tau_binned_data{g}(d,:) = [b g (i - frame_w) i zeros(1,length(Tau_data(1,:)))];
        Distance_binned_data{g}(d,:) = [b g (i - frame_w) i zeros(1,length(Distance_data(1,:)))];
     end
    g = g+1;
    idx = [];
    d = 1;

end

%%
% Get intensity at distances
if b == 1;
    new_data = [];
    endrow = 0;
end
new_data = [new_data;zeros(length(Tau_binned_data),((ves_dist-min_dist)+5))];

for i = 1:length(Tau_binned_data)

   Tau_single_bin_data = Tau_binned_data{i};
   Distance_single_bin_data = Distance_binned_data{i};
  
   if ~isempty(Tau_single_bin_data)
   new_data(endrow+i,1:4) = Tau_single_bin_data(1,1:4);
   else
   new_data(endrow+i,1:4) = [b    i   new_data(endrow+i-1,3)+frame_w      new_data(endrow+i-1,4)+frame_w];
   end

   for d = 0:res:(ves_dist-min_dist)
       dist_min = d+min_dist-0.5;
       dist_max = d+min_dist+res-0.5;
       if d == 0
       idx = find(Tau_single_bin_data(:,10) >= dist_min & Tau_single_bin_data(:,10) <= dist_max);
       idx2 = find(Distance_single_bin_data(:,10) >= dist_min & Distance_single_bin_data(:,10) <= dist_max);
       if ~isempty(idx) 
       volume = length(idx2)*(voxel(1)*voxel(2)*voxel(3));
       density = length(idx)/volume;
       new_data(i,d+4) = density;
       new_data(endrow+i,d+5) = density;
       end

       else
       idx = find(Tau_single_bin_data(:,10) > d+min_dist-0.5 & Tau_single_bin_data(:,10)<= dist_max);
       idx2 = find(Distance_single_bin_data(:,10) > d+min_dist-0.5 & Distance_single_bin_data(:,10)<= dist_max);
       if ~isempty(idx) 
       volume = length(idx2)*(voxel(1)*voxel(2)*voxel(3));
       density = length(idx)/volume;
        new_data(endrow+i,d+5) = density;
       end

       end
  
   end

end

endrow = length(new_data(:,1));
%binned_data = {};
end

%% Headers
intensity_headers = [min_dist:ves_dist];
intensity_headerstemp = {};
for i = 1:length(intensity_headers)
    intensity_headerstemp(:,i) = {append(sprintf('%g',intensity_headers(1,i)),' um Density (#/um^3)')};
end
intensity_headers = intensity_headerstemp;
headers = cell2table({ 'Vessel#'  'Frame#'  'MinFrameDistance' 'MaxFrameDistance' intensity_headers{:}});

%% Save Data (Add header export in separate file)
writetable(headers,append(outputfoldername,'/',outputname,'_HEADERS.csv'));
csvwrite(append(outputfoldername,'/',outputname,'.csv'),new_data);



