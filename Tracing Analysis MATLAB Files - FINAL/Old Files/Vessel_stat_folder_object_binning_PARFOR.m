clc; clear all; close all;

% Exported Data Format: [ Vessel#  Frame#  MinFrameDistance MaxFrameDistance  ObjectDensity  ]

%% Import Data

% .txt file (Recommended for speed):

% Input tau data folder name (INCLUDE FOLDER PATH AND NAME):
Tau_sourcefolder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\Tangle_distance_data';

% Input distance transform data folder name for volume calculation (INCLUDE FOLDER PATH AND NAME):
Distance_sourcefolder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\tangle_volume_data';

% Max Distance From Vessel (in microns) (Recommended Value 30):

ves_dist = [ 100 ];

min_dist = [ -10 ];

% Resolution in microns? (Recommended Value: 1):

res = [ 1 ];

% Frame width (# Microns):

frame_w = [ 10 ];

% Voxel Size

% Voxel Size:  [   X     Y    Z   ]
voxel =        [ 0.621 0.621 2.50 ];

% Data Output Name (True, unshuffled data only) (Include Name + Path):
% Vessel number and data type will be appended, so do not include in names.
% NEED FOLDER WITH SAME NAME FOR EACH DATA TYPE: XDATA, YDATA, AND ERRORS
% ex. you need folders stat_data_xdata, stat_data_ydata, and
% stat_data_errors if stat_data is the folder name.

% Folder:
outputfoldername = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\Tangle_density_data';

% File:
outputname = 'Human2267_tangle_density_data';

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

fcontent = dir(fullfile(Tau_sourcefolder, '*.txt')); %fcontent is a column vector of structures

%i = 7;

   filename = fcontent(b).name;
   Tau_data_single = [];
   Tau_data_single = readmatrix(fullfile(Tau_sourcefolder,filename));
   [r, c] = size(Tau_data_single);
   l2 = l1 + r - 1;
   Tau_data(l1:l2,:) = Tau_data_single;
   l1 = 1;

    waitbar(0,f,sprintf('Calculating Statistics: Reading Data %d %%',floor(i/n*100)));



waitbar(0.1,f,'Calculating Statistics: Preparing Data');

% vesnum = vessel ID    g = frame ID

fcontent = dir(fullfile(Tau_sourcefolder, '*.txt')); %fcontent is a column vector of structures
filename = fcontent(b).name;
[pathstr,name,ext] = fileparts(filename);

name_parts = regexp(name,'_','split');
vesnum = str2num(name_parts{3});

if b == 1;
    new_data = [];
    endrow = 0;
end

if isempty(Tau_data)

maxlength = (ceil(max(Distance_data(:,5)/frame_w))*frame_w);

for i = 1:maxlength/frame_w

    new_data(endrow+i,:) = [ vesnum i (i*frame_w - frame_w) i*frame_w  zeros(1,((ves_dist-min_dist)/res)+1)];

end
endrow = length(new_data(:,1));
    Tau_binned_data = [];
    Distance_binned_data = [];
    continue
end

% Remove data where Tau is within blood vessel and greater than ves_dist
% away
Tau_data_backup = Tau_data;
idx = find( Tau_data(:,6) <= min_dist);
Tau_data(idx,:) = [];
idx = find(Tau_data(:,6) > ves_dist );
Tau_data(idx,:) = [];

if isempty(Tau_data)
    Tau_data = zeros(1,6);
end

% Read Distance Data

fcontent = dir(fullfile(Distance_sourcefolder, '*.txt')); %fcontent is a column vector of structures

   filename = fcontent(b).name;
   Distance_data_single = [];
   Distance_data_single = readmatrix(fullfile(Distance_sourcefolder,filename));
   Distance_data = [Distance_data_single(:,1:3) zeros(length(Distance_data_single(:,1)),1) Distance_data_single(:,5) Distance_data_single(:,4)];

    waitbar(0,f,sprintf('Calculating Statistics: Reading Data %d %%',floor(i/n*100)));




waitbar(0.1,f,'Calculating Statistics: Preparing Data');


% Remove data where Tau is within blood vessel and greater than ves_dist
% away
Distance_data_backup = Distance_data;
idx = find( Distance_data(:,6) <= min_dist);
Distance_data(idx,:) = [];
idx = find(Distance_data(:,6) > ves_dist );
Distance_data(idx,:) = [];

%% Bin Data by Surface Frames


waitbar(0.3,f,sprintf('Binning Data for Vessel %d/%d: Cropping Along Vessel',b,n));


maxlength = (ceil(max(Distance_data(:,5)/frame_w))*frame_w);

idx = cell(1,maxlength/frame_w);
idx2 = cell(1,maxlength/frame_w);

parfor i = 1:(maxlength/frame_w)

    g = i*frame_w;

    idx{i} = find(Tau_data(:,5) < g & Tau_data(:,5) >= (g - frame_w));

if ~isempty(idx{i})

    idx2{i} = find(Distance_data(:,5) < g & Distance_data(:,5) >= (g - frame_w));

end
end


for i = frame_w:frame_w:maxlength

    g = i/frame_w;

    idx_temp = cell2mat(idx(g));
    idx2_temp = cell2mat(idx2(g));

     %  [b g (i - frame_w) i Tau_data(idx,:)]
     if ~isempty(idx_temp)

     for d = 1:length(idx_temp)
     Tau_binned_data{g}(d,:) = [vesnum g (i - frame_w) i (Tau_data(idx_temp(d),:))];
     end

     if ~isempty(idx2_temp)
     for d = 1:length(idx2_temp)
     Distance_binned_data{g}(d,:) = [vesnum g (i - frame_w) i (Distance_data(idx2_temp(d),:))];
     end
     else
     Distance_binned_data{g} = [vesnum g (i - frame_w) i zeros(1,length(Distance_data(1,:)))];
     end

     
        else
        Tau_binned_data{g} = [vesnum g (i - frame_w) i zeros(1,(length(Tau_data(1,:))-1)) min_dist-res-100];
        Distance_binned_data{g} = [vesnum g (i - frame_w) i zeros(1,length(Distance_data(1,:))-1) min_dist-res-100];
     end

    d = 1;

end

%%
% Get intensity at distances

waitbar(0.6,f,sprintf('Binning Data for Vessel %d/%d: Cropping Away From Vessel',b,n));

new_data = [new_data;zeros(length(Tau_binned_data),((ves_dist-min_dist)+5))];

for i = 1:length(Tau_binned_data)

   Tau_single_bin_data = Tau_binned_data{i};
   Distance_single_bin_data = Distance_binned_data{i};
  
   if ~isempty(Tau_single_bin_data)
   new_data(endrow+i,1:4) = Tau_single_bin_data(1,1:4);

   for d = 0:res:(ves_dist-min_dist)
       dist_min = d+min_dist-0.5;
       dist_max = d+min_dist+res-0.5;
       if d == 0
       idx = find(Tau_single_bin_data(:,10) >= dist_min & Tau_single_bin_data(:,10) <= dist_max);
       if ~isempty(idx)
       idx2 = find(Distance_single_bin_data(:,10) >= dist_min & Distance_single_bin_data(:,10) <= dist_max);
       volume = length(idx2)*(voxel(1)*voxel(2)*voxel(3));
       density = length(idx)/volume;
       new_data(endrow+i,d+5) = density;
       else
           new_data(endrow+i,d+5) = 0;
       end
     
       else
       idx = find(Tau_single_bin_data(:,10) > d+min_dist-0.5 & Tau_single_bin_data(:,10) <= dist_max);
       if ~isempty(idx)
       idx2 = find(Distance_single_bin_data(:,10) > d+min_dist-0.5 & Distance_single_bin_data(:,10) <= dist_max);
       volume = length(idx2)*(voxel(1)*voxel(2)*voxel(3));
       density = length(idx)/volume;
       new_data(endrow+i,d+5) = density;
       else
            new_data(endrow+i,d+5) = 0;
       end

       end
  
   end

   else
   new_data(endrow+i,1:4) = [Tau_single_bin_data(1,1:4) zeros(1,((ves_dist-min_dist)/res)+1)];
   end

end

endrow = length(new_data(:,1));
Tau_binned_data = [];
Distance_binned_data = [];
%binned_data = {};
end

%% Headers
sampleID = name_parts{1};
intensity_headers = [min_dist:ves_dist];
intensity_headerstemp = {};
for i = 1:length(intensity_headers)
    intensity_headerstemp(:,i) = {append(sprintf('%g',intensity_headers(1,i)),' um Density (#/um^3)')};
end
intensity_headers = intensity_headerstemp;
headers = cell2table({ sampleID 'Vessel#'  'Frame#'  'MinFrameDistance' 'MaxFrameDistance' intensity_headers{:}});

%% Save Data (Add header export in separate file)
writetable(headers,append(outputfoldername,'/',outputname,'_HEADERS.csv'));
csvwrite(append(outputfoldername,'/',outputname,'.csv'),new_data);



