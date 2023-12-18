%%
clc; close all; clear all;

%% Nancy included this for dealing with parallel pool shutdown
availableGPUs = gpuDeviceCount("available");

parpool('Processes',availableGPUs);

spmd
    gpuDevice
end

%% Exported Data format: 
% Column:     1  2  3      4             5              6
%           [ X  Y  Z   Distance   CorticalLayer   TauIntensity]

%% IMPORT DATA/SETTINGS

% DISTANCE TRANSFORM Images (must be .tif):

distance_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\2131\2465v2_2 transform tiffs';

% Cortical Layer File:

cortical_file = 'C:\Users\exx\Desktop\Zach - Active Analysis\2465\Human2465_layers.tif';

% Tau Intensity File:

intensity_file = 'C:\Users\exx\Desktop\Zach - Active Analysis\2465\2465v2_tau_intensity.tif';

% Save final data? What Format? Yes = [1]  No = [0]  (.txt only is recommended)

%       [  .txt    .xlsx  ]
yesno = [    1       0    ];

% Folder path for output (ex. '/Users/zacharyhoglund/Desktop'):

folder_dir = 'C:\Users\exx\Desktop\Zach - Active Analysis\2465\FINAL_coordinates';

% Data files will be exported as Human####_Vessel_##_DT_COORD

%%
fcontent = dir(fullfile(distance_folder, '*.tif'));

if length(fcontent') == 0
    fcontent =  dir(fullfile(distance_folder, '*.tif'));
end

k = length(fcontent');

if k == 0
    error('ERROR: Distance transform folder is empty');
end

g = waitbar(0,sprintf('Reading image for vessel 1/%d',k));

% Import Cortical Layer Data

info = imfinfo(cortical_file);
numberOfPages = length(info);

for i = 1:numberOfPages
    cortical_image(:,:,i) = imread(cortical_file,i);
end

% Import Intensity Data

info = imfinfo(intensity_file);
numberOfPages = length(info);

for i = 1:numberOfPages
    intensity_image(:,:,i) = imread(intensity_file,i);
end

for h = 1:k

% Import distance transform data

waitbar(0,g,sprintf('Reading image for vessel %d/%d',h,k));
dist_fcontent = dir(fullfile(distance_folder, '*.tif')); %fcontent is a column vector of structures
dist_sourcefolder = distance_folder;

%% Open distance transform tiff

dist_filename = dist_fcontent(h).name;
%%
info = imfinfo(fullfile(dist_sourcefolder,dist_filename));
numberOfPages = length(info);

%%
for i = 1:numberOfPages
    distance_image(:,:,i) = imread(fullfile(dist_sourcefolder,dist_filename),i);
end

if size(distance_image) ~= size(cortical_image) | size(distance_image) ~= size(intensity_image)
    error('ERROR: Image sizes do not match!');
end

waitbar(0.5,g,sprintf('Getting coodinates for vessel %d/%d',h,k));

Tau_data = [];

[y, x, z] = ind2sub(size(distance_image),find(distance_image < 101));
distance_data = double(zeros(length(x),6));

[r,c,t] = size(distance_image);

parfor i = 1:length(x)
distance_data(i,:) = [double(x(i)) double(y(i)) double(z(i)) double(distance_image(y(i),x(i),z(i))) double(cortical_image(y(i),x(i),z(i))) double(intensity_image(y(i),x(i),z(i)))];
end

waitbar(0.9,g,sprintf('Exporting data for vessel %d/%d',h,k));

l = length(distance_data(:,6));

idx = randperm(l);

distance_data_backup = distance_data;

shuffled_distance_data(idx,6) = distance_data(:,6);

distance_data = shuffled_distance_data;

[pathstr,name,ext] = fileparts(dist_filename);

name_parts = regexp(name,'_','split');
vesnum = str2num(name_parts{3});
new_name = append(name_parts{1},'_',name_parts{2},'_',name_parts{3});


txt_name = append(folder_dir,'\',new_name,'_SHUFFLED_DT_CORTICAL_INTENSITY_COORD','.txt');
xlsx_name = append(folder_dir,'\',new_name,'_SHUFFLED_DT_CORTICAL_INTENSITY_COORD', '.xlsx');

% % txt_name = append(final_name,ves_num,'.txt');
% % xlsx_name = append(final_name,ves_num, '.xlsx');

txt_name = char(txt_name);
xlsx_name = char(xlsx_name);

if yesno(1) == 1
writematrix(distance_data,txt_name);
end

if yesno(2) == 1
writematrix(distance_data,xlsx_name);
end


clear distance_image
%clear distance_data

end

waitbar(1,g,'Finished!');