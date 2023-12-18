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

% VESSEL DISTANCE TRANSFORM, CORTICAL, and INTENSITY COORDINATES:

distance_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\2392\Transform_intensity_cortical_coordinates';

% Input Imaris Image Dimentions (Get this from Edit -> Image Properties)

% X-axis lower bound (in microns):

xlow = [ 0 ];

% X-axis upper bound (in microns):

xup = [ 3663 ];

% Y-axis lower bound (in microns):

ylow = [ 0 ];

% Y-axis upper bound (in microns):

yup = [ 2459 ];

% Z-axis lower bound (in microns):

zlow = [ 3575 ];

% Z-axis upper bound (in microns):

zup = [ 3930 ];

% Fiji Image Dimensions:

% Distance Transform X dimension:

dist_xdim = [ 5894 ];

% Distance Transform Y dimension:

dist_ydim = [ 3958 ];

% Number of Distance Transform  z-stack slices:

dist_znum = [ 71 ];

% Save final data? What Format? Yes = [1]  No = [0]  (.txt only is recommended)

%       [  .txt    .xlsx  ]
yesno = [    1       0    ];

% Folder path for tau intensity output (ex. '/Users/zacharyhoglund/Desktop'):

folder_dir = 'C:\Users\exx\Desktop\Zach - Active Analysis\2392\Intensity_distance_data';


% Folder path for tau vessel cortical layer output (ex. '/Users/zacharyhoglund/Desktop'):

folder_dir2 = 'C:\Users\exx\Desktop\Zach - Active Analysis\2392\vessel_cortical_data';

% Data files will be exported as Human####_Vessel_##_int_dist_data


%%
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures
k = length(fcontent');

f = waitbar(0,'Reading Tau Data');

for h = 1:k

% Import distance transform data

waitbar(0,f,'Reading Distance Data');
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = distance_folder;

Distance_data = {};

   filename = fcontent(h).name;
   distance_data_single = [];
   Distance_data = readmatrix(fullfile(sourcefolder,filename));

Distance_data_orig = Distance_data;
Tau_data = [Distance_data(:,1:3) Distance_data(:,6) zeros(length(Distance_data(:,1)),1) Distance_data(:,4)];
Distance_data = [Distance_data(:,1:5)];
% k = 1;
% h = 1;

waitbar(0,f,sprintf('Getting Distances for Vessel %d/%d: Preparing Data',h,k));


%%

n = length(Tau_data(:,1));

%%
% 
% tau_xdim = max(Tau_data(:,1));
% tau_ydim = max(Tau_data(:,2));
% tau_znum = max(Tau_data(:,3));

Tau_data(:,1) = ((Tau_data(:,1)./dist_xdim).*(xup - xlow))+xlow;
Distance_data(:,1) = ((Distance_data(:,1)./dist_xdim).*(xup - xlow))+xlow;

Tau_data(:,2) = ((Tau_data(:,2)./dist_ydim).*(yup - ylow))+ylow;
Distance_data(:,2) = ((Distance_data(:,2)./dist_ydim).*(yup - ylow))+ylow;

Tau_data(:,3) = ((Tau_data(:,3)./dist_znum).*(zup - zlow))+zlow;
Distance_data(:,3) = ((Distance_data(:,3)./dist_znum).*(zup - zlow))+zlow;

%% create vessel data:

vessel_idx = find(Distance_data(:,4) <= 0);

vessel_data = Distance_data(vessel_idx,:);


%%


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
    idx = find(vessel_data2(:,5) >= o);
    idx2 = find(vessel_data2(idx,5) <= i);
    fit_data_imp(p,:) = [mean(vessel_data2(idx(idx2),1))  mean(vessel_data2(idx(idx2),2)) mean(vessel_data2(idx(idx2),3)) 0];
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

%% Find Cortical Layer distance along Blood Vessel
%waitbar(0.4,f,sprintf('Getting Distances for Vessel %d/%d: Finding Distance Along Vessel',h,k));
% Create matrix for distance along vessel
dist = zeros(length(Tau_data),1);

% Find the nearest point in the best fit line dataset to each point in the
% tau dataset. 
near = knnsearch(fitdata(:,1:3),vessel_data(:,1:3),Distance="fasteuclidean");

% Calculate arclength between the end of the best fit line and each nearest
% point on the best fit line to get distance along the vessel. 
parfor i = 1:length(vessel_data)
    if length(fitdata(1:near(i,:),1)) == 1
        vessel_data(i,6) = [0];
    else
    % Calculate arclength
    dist(i,:) = arclength(fitdata(1:near(i,:),1),fitdata(1:near(i,:),2),fitdata(1:near(i,:),3));
    % Store distance data in new column of tau dataset. 
    vessel_data(i,6) = dist(i,:);
    end
end

% Sort vessel dataset by distance down the vessel for plotting. 
vessel_data = sortrows(vessel_data,6);

%% Export Data
f = waitbar(0.9,sprintf('Getting Distances for Vessel %d/%d: Exporting Data',h,k));

fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures
filename = fcontent(h).name;
[pathstr,name,ext] = fileparts(filename);

name_parts = regexp(name,'_','split');
vesnum = str2num(name_parts{3});
new_name = append(name_parts{1},'_',name_parts{2},'_',name_parts{3});

Tau_data = Tau_data(:,1:6);
final_data = Tau_data;
cortical_data = vessel_data;

%% Export Intensity Data

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

%% Export Vessel Cortical Layer Data

txt_name = append(folder_dir2,'/',new_name,'vessel_cortical_data','.txt');
xlsx_name = append(folder_dir2,'/',new_name,'vessel_cortical_data', '.xlsx');

% % txt_name = append(final_name,ves_num,'.txt');
% % xlsx_name = append(final_name,ves_num, '.xlsx');

txt_name = char(txt_name);
xlsx_name = char(xlsx_name);

if yesno(1) == 1
writematrix(cortical_data,txt_name);
else
end
if yesno(2) == 1
writematrix(cortical_data,xlsx_name);
else
end


close all;
end

waitbar(1,f,'Finished!');

toc

close(f)
