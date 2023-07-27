clc; close all; clear all;
%% READ ME

% DATA OUTPUT FORMAT: 
%      Column: [  1    2    3         4                   5                      6           ]
%              [  X    Y    Z     Intensity      DistanceAlongVessel    DistanceFromSurface  ]

% IMPORTANT NOTES: 
% DATA FROM IMAGEJ MUST BE OF THE SAME RESOLUTION LEVEL
% DATA IN FOLDERS MUST BE IN SAME ORDER 
% ONLY WORKS WITH DATA FROM SAME SOURCE IMAGE
% INCLUDE FOLDER PATH AND NAME

%% IMPORT DATA/SETTINGS

% TAU INTENSITY Image (must be .tiff):

tau_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\Tau Tiff';

% DISTANCE TRANSFORM Images (must be .tiff):

distance_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\Distance Tiff';

% Input Imaris Image Dimentions (Get this from Edit -> Image Properties)

% X-axis lower bound (in microns):

xlow = [ 0 ];

% X-axis upper bound (in microns):

xup = [ 2454 ];

% Y-axis lower bound (in microns):

ylow = [ 0 ];

% Y-axis upper bound (in microns):

yup = [ 2444 ];

% Z-axis lower bound (in microns):

zlow = [ 430 ];

% Z-axis upper bound (in microns):

zup = [ 853 ];

% Save final data? What Format? Yes = [1]  No = [0]  (.txt only is recommended)

%       [  .txt    .xlsx  ]
yesno = [    1       0    ];

% Folder path for output (ex. '/Users/zacharyhoglund/Desktop'):

folder_dir = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2302\Intensity_distance_data';

% Data files will be exported as Human####_Vessel_##_int_dist_data


%%
fcontent = dir(fullfile(distance_folder, '*.tif')); %fcontent is a column vector of structures
k = length(fcontent');

f = waitbar(0,'Reading Tau Image');

tau_fcontent = dir(fullfile(tau_folder, '*.tif'));
tau_sourcefolder = tau_folder;
tau_filename = tau_fcontent.name;

info = imfinfo(fullfile(tau_sourcefolder,tau_filename));
numberOfPages = length(info);
%%
for i = 1:numberOfPages
    tau_image(:,:,i) = imread(fullfile(tau_sourcefolder,tau_filename),i);
end
%%
for h = 1:k

% Import vessel data

waitbar(0,f,'Reading Distance Transform Image');
dist_fcontent = dir(fullfile(distance_folder, '*.tif')); %fcontent is a column vector of structures
dist_sourcefolder = distance_folder;

%% Open distance transform tiff

dist_filename = dist_fcontent(h).name;

info = imfinfo(fullfile(dist_sourcefolder,dist_filename));
numberOfPages = length(info);

distance_data = [];

for i = 1:numberOfPages
    distance_image(:,:,i) = imread(fullfile(dist_sourcefolder,dist_filename),i);
end

waitbar(0,f,'Getting Distance & Tau Coordinates');
distance_data = [];
Tau_data = [];

[x, y, z] = ind2sub(size(distance_image),find(distance_image < 101));
distance_data = [x y z distance_image(x,y,z)];
Tau_data = [x y z tau_image(x,y,z) zeros(length(x),1) distance_image(x,y,z)];

waitbar(0,f,'Getting Vessel Coordinates');

[x, y, z] = ind2sub(size(distance_image),find(distance_image <= 0));
vessel_data = [x y z];

   waitbar(0,f,'Getting Tau Coordinates');
   
Tau_data_orig = Tau_data;

% k = 1;
% h = 1;

waitbar(0,f,sprintf('Getting Distances for Vessel %d/%d: Preparing Data',h,k));

%%
% 
% tau_xdim = max(Tau_data(:,1));
% tau_ydim = max(Tau_data(:,2));
% tau_znum = max(Tau_data(:,3));

tau_xdim = size(Tau_image,1);
tau_ydim = size(Tau_image,2);
tau_zdim = size(Tau_image,3);

ves_xdim = size(distance_image,1);
ves_ydim = size(distance_image,2);
ves_zdim = size(distance_image,3);

Tau_data(:,1) = ((Tau_data(:,1)./tau_xdim).*(xup - xlow))+xlow;

Tau_data(:,2) = ((Tau_data(:,2)./tau_ydim).*(yup - ylow))+ylow;

Tau_data(:,3) = ((Tau_data(:,3)./tau_znum).*(zup - zlow))+zlow;

vessel_data(:,1) = ((vessel_data(:,1)./ves_xdim).*(xup - xlow))+xlow;

vessel_data(:,2) = ((vessel_data(:,2)./ves_ydim).*(yup - ylow))+ylow;

vessel_data(:,3) = ((vessel_data(:,3)./ves_znum).*(zup - zlow))+zlow;

%% Find Line of best fit equation and dataset 
waitbar(0.1,f,sprintf('Getting Distances for Vessel %d/%d: Finding Line of Best Fit',h,k));

% Find minimum axes values of vessel data for improvement of fit line. 
minx = min(vessel_data(:,1));
miny = min(vessel_data(:,2));
minz = min(vessel_data(:,3));

% Finds values of coefficients for polynomial equation of best fit line.
fitdata = [0 0 0];

% finds best fit line for X-Y plane
fxy = fit((vessel_data(:,1)-minx),(vessel_data(:,2)-miny),'poly1');

% finds best fit line for X-Z plane
fxz = fit((vessel_data(:,1)-minx),(vessel_data(:,3)-minz),'poly1');

% Get coefficients from equation output
fxy_values = coeffvalues(fxy);
fxz_values = coeffvalues(fxz);

%fit_res_1 = (max(vessel_data(:,1)) - min(vessel_data(:,1)))/fit_res;

x_fit = ((min(vessel_data(:,1))-minx-5):(max(vessel_data(:,1))-minx+5))';
y_fit = fxy(x_fit);
z_fit = fxz(x_fit);

x_fit = x_fit + minx;
y_fit = y_fit + miny;
z_fit = z_fit + minz;

fitdata = [x_fit y_fit z_fit];

% Remove rows containing all zeros from best fit line dataset. 
fitdata = fitdata(any(fitdata,2),any(fitdata,1));

%% Improve Best Fit Line

vessel_data2 = vessel_data;
nearfit = knnsearch(fitdata(:,1:3),vessel_data2(:,1:3),Distance='Euclidean');

% Calculate arclength between the end of the best fit line and each nearest
% point on the best fit line to get distance along the vessel. 
for i = 1:length(vessel_data2)
    if length(fitdata(1:nearfit(i,:),1)) == 1
        vessel_data2(i,5) = [0];
    else
    % Calculate arclength
    vesdist(i,:) = arclength(fitdata(1:nearfit(i,:),1),fitdata(1:nearfit(i,:),2),fitdata(1:nearfit(i,:),3));
    % Store distance data in new column of tau dataset. 
    vessel_data2(i,5) = vesdist(i,:);
    end
end

vessel_data2 = sortrows(vessel_data2,5);

mindist = min(vessel_data2(:,5));
maxdist = max(vessel_data2(:,5));

%range = maxdist - mindist;

%fit_resolution = range/fit_res;
fit_data_imp = [];
o = mindist;
p = 1;
fit_data_imp = zeros(1,4);
for i = mindist:1:maxdist
    idx = find(vessel_data2(:,5) >= o & vessel_data2(idx,5) <= i);
    fit_data_imp(p,:) = [mean(vessel_data2(idx,1))  mean(vessel_data2(idx,2)) mean(vessel_data2(idx,3)) 0];
o = i;
p = p+1;

end

%%
p = [];
idx = find(isnan(fit_data_imp(:,1)));
fit_data_imp(idx,:) = [];
fitdata = fit_data_imp(:,1:3);

%% Plot Best Fit Line with Vessel Surface Data

% Plot the best fit line
figure()
scatter3(fitdata(:,1), fitdata(:,2), fitdata(:,3),0.3,'.','r');
hold on

% Plot the Blood vessel on same plot as best fit line.
scatter3(vessel_data(:,1),vessel_data(:,2),vessel_data(:,3),0.3,vessel_data(:,3));
hold on
% % % vis_points(:,1) = randperm(max(Tau_data(:,1)),3);
% % % vis_points(:,2) = randperm(max(Tau_data(:,2)),3);
% % % vis_points(:,3) = randperm(max(Tau_data(:,3)),3);
% % % scatter3(vis_points(:,1), vis_points(:,2), vis_points(:,3));
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal
grid on
colormap(jet)
title('Vessel Only Plot');

%%
figure()
scatter3(Tau_data(:,1),Tau_data(:,2),Tau_data(:,3),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
hold on
scatter3(vessel_data(:,1),vessel_data(:,2),vessel_data(:,3));
title('Tau + Vessel Plot');

%% Find Tau distance along Blood Vessel
%waitbar(0.4,f,sprintf('Getting Distances for Vessel %d/%d: Finding Distance Along Vessel',h,k));
% Create matrix for distance along vessel
dist = zeros(length(Tau_data),1);

% Find the nearest point in the best fit line dataset to each point in the
% tau dataset. 
near = knnsearch(fitdata(:,1:3),Tau_data(:,1:3),Distance="fasteuclidean");

% Calculate arclength between the end of the best fit line and each nearest
% point on the best fit line to get distance along the vessel. 
parfor i = 1:length(Tau_data)
    if length(fitdata(1:near(i,:),1)) == 1
        Tau_data(i,5) = [0];
    else
    % Calculate arclength
    dist(i,:) = arclength(fitdata(1:near(i,:),1),fitdata(1:near(i,:),2),fitdata(1:near(i,:),3));
    % Store distance data in new column of tau dataset. 
    Tau_data(i,5) = dist(i,:);
    end
end

% Sort tau dataset by distance down the vessel for plotting. 
Tau_data = sortrows(Tau_data,5);

%% Export Data
f = waitbar(0.9,sprintf('Getting Distances for Vessel %d/%d: Exporting Data',h,k));

fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures
filename = fcontent(h).name;
[pathstr,name,ext] = fileparts(filename);

name_parts = regexp(name,'_','split');
vesnum = str2num(name_parts{3});
new_name = append(name_parts{1},'_',name_parts{2},'_',name_parts{3});

fcontent3 = dir(fullfile(tau_folder, '*.txt')); %fcontent is a column vector of structures
filename3 = fcontent2(h).name;
[pathstr3,name3,ext3] = fileparts(filename3);

name_parts3 = regexp(name3,'_','split');
vesnum3 = str2num(name_parts3{3});

if length(unique([vesnum vesnum3])) > 1
    waitbar(0,f,sprintf('VESSEL MATCH ERROR: Object, vessel, and/or distance transform data are from different vesels (obj: Vessel %d Volume: Vessel %d Vessel: Vessel %d)',vesnum3,vesnum, vesnum2));
    'ERROR: VESSEL MISMATCH (see error window)'
    break
end

Tau_data = Tau_data(:,1:6);
final_data = Tau_data;

txt_name = append(folder_dir,'/',new_name,'intensity_dist_data','.txt');
xlsx_name = append(folder_dir,'/',new_name,'intensity_dist_data', '.xlsx');

% % txt_name = append(final_name,ves_num,'.txt');
% % xlsx_name = append(final_name,ves_num, '.xlsx');

txt_name = char(txt_name);
xlsx_name = char(xlsx_name);

if yesno(1) == 1
writematrix(final_data,txt_name);
else
end
if yesno(2) == 1
writematrix(final_data,xlsx_name);
else
end
close all;
end

waitbar(1,f,'Finished!');

toc

close(f)

