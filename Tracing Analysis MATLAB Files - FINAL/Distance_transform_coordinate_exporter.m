%%
clc; close all; clear all;

%% Exported Data format:
% Column:     1  2  3     4
%           [ X  Y  Z  Distance ]

%% IMPORT DATA/SETTINGS

% DISTANCE TRANSFORM Images (must be .tiff):

distance_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2131 DT test tiff';

% Save final data? What Format? Yes = [1]  No = [0]  (.txt only is recommended)

%       [  .txt    .xlsx  ]
yesno = [    0       0    ];

% Folder path for output (ex. '/Users/zacharyhoglund/Desktop'):

folder_dir = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2302\Intensity_distance_data';

% Data files will be exported as Human####_Vessel_##_DT_COORD


%%
fcontent = dir(fullfile(distance_folder, '*.tif'));

if length(fcontent') == 0
    fcontent =  dir(fullfile(distance_folder, '*.tiff'));
end

k = length(fcontent');

if k == 0
    error('ERROR: Distance transform folder is empty');
end

f = waitbar(0,sprintf('Reading image for vessel 1/%d',k));

for h = 1:k

% Import vessel data

waitbar(0,f,sprintf('Reading image for vessel %d/%d',h,k));
dist_fcontent = dir(fullfile(distance_folder, '*.tif')); %fcontent is a column vector of structures
dist_sourcefolder = distance_folder;

%% Open distance transform tiff

dist_filename = dist_fcontent(h).name;

info = imfinfo(fullfile(dist_sourcefolder,dist_filename));
numberOfPages = length(info);


for i = 1:numberOfPages
    distance_image(:,:,i) = imread(fullfile(dist_sourcefolder,dist_filename),i);
end

waitbar(0.5,f,sprintf('Getting coodinates for vessel %d/%d',h,k));

Tau_data = [];

[y, x, z] = ind2sub(size(distance_image),find(distance_image < 101));
distance_data = zeros(length(x),4);

[r,c,t] = size(distance_image);

parfor i = 1:length(x)
distance_data(i,:) = [x(i) y(i) z(i) distance_image(y(i),x(i),z(i))];
end

waitbar(0.9,f,sprintf('Exporting data for vessel %d/%d',h,k));


[pathstr,name,ext] = fileparts(dist_filename);

name_parts = regexp(name,'_','split');
vesnum = str2num(name_parts{3});
new_name = append(name_parts{1},'_',name_parts{2},'_',name_parts{3});


txt_name = append(folder_dir,'/',new_name,'_COORD','.txt');
xlsx_name = append(folder_dir,'/',new_name,'_COORD', '.xlsx');

% % txt_name = append(final_name,ves_num,'.txt');
% % xlsx_name = append(final_name,ves_num, '.xlsx');

txt_name = char(txt_name);
xlsx_name = char(xlsx_name);

if yesno(1) == 1
writematrix(distance_data,txt_name);
elseif yesno(2) == 1
writematrix(distance_data,xlsx_name);
end

clear distance_image
clear distance_data

end

waitbar(1,f,'Finished!');
