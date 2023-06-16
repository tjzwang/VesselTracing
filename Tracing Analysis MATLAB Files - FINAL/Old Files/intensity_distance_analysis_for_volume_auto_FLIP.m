%% READ ME

% DATA OUTPUT FORMAT: 
%      Column: [  1    2    3         4                   5                      6           ]
%              [  X    Y    Z     Intensity       DistanceAlongVessel    DistanceFromSurface ]

% IMPORTANT NOTES: 
% DATA FROM IMAGEJ MUST BE OF THE SAME RESOLUTION LEVEL
% DATA IN FOLDERS MUST BE IN SAME ORDER 
% ONLY WORKS WITH DATA FROM SAME SOURCE IMAGE
% INCLUDE FOLDER PATH AND NAME

%% IMPORT DATA/SETTINGS

clc; clear all; close all;

tic

% INPUT FOLDER NAMES + PATH (ex. (ex. '/Users/zacharyhoglund/Desktop/Human2131/Vessel_Data')


% VESSEL DISTANCE TRANSFORMS:

distance_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\DT_coord_hi_res';

% VESSEL DATA:

vessel_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\Vessel_coord_hi_res';

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

% Fiji Image Dimensions:

% Vessel X dimension:

ves_xdim = [ 3949 ];

% Distance Transform X dimension:

dist_xdim = [ 3949 ];

% Vessel Y dimension:

ves_ydim = [ 3934 ];

% Distance Transform Y dimension:

dist_ydim = [ 3934 ];

% Number of vessel z-stack slices:

ves_znum = [ 169 ];

% Number of Distance Transform  z-stack slices:

dist_znum = [ 169 ];

% Set Moving Average Size (DOES NOT AFFECT EXPORTED DATA)
Int_moving_ave_size = [  1000  ];

Int_surfplot_moving_ave_size = [  500  ];

% Make test plots? Yes = [1] No = [0] (DOES NOT AFFECT EXPORTED DATA):

testplot = [  1  ];

% Plot Intensity vs Distance Metrics? Yes = [1] No = [0] (DOES NOT AFFECT EXPORTED DATA):

intdistplots = [  1  ];

% Axes equal for 3D plot? Yes = [1] No = [0] (DOES NOT AFFECT EXPORTED DATA):

axes = [   1   ];

% Save final data? What Format? Yes = [1]  No = [0]  (.txt only is recommended)

%       [  .txt    .xlsx  ]
yesno = [    1       0    ];

% Folder path for output (ex. '/Users/zacharyhoglund/Desktop'):

folder_dir = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\tangle_volume_data';

% Adjusted Data Export Name? (DO NOT INCLUDE FILE TYPE OR VESSEL NUMBER!)
% ex. Human2131_int_dist_data
% Vessel number will be appended to file name (ex.
% Human2131_int_dist_data_1)

final_name =  'MATLAB_2470v3_volume_output' ;

%%
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures

k = length(fcontent');

%for h = 1:k

k = 1;
h = 8;

% Import distance transform data

f = waitbar(0,'Reading Distance Data');
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = distance_folder;

n = length(fcontent');

   filename = fcontent(h).name;
   dist_filename = filename;
   distance_data_single = [];
   distance_data_single = readmatrix(fullfile(sourcefolder,filename));
   Distance_data = distance_data_single;

% Import vessel data:

waitbar(0,f,'Reading Vessel Data');
fcontent = dir(fullfile(vessel_folder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = vessel_folder;


n = length(fcontent');



   filename = fcontent(h).name;
   ves_filename = filename;
   vessel_data_single = [];
   vessel_data_single = readmatrix(fullfile(sourcefolder,filename));
   vessel_data = vessel_data_single;



%%

n = length(Distance_data(:,1));

Distance_data = sortrows(Distance_data,1);
Distance_data = sortrows(Distance_data,2);
Distance_data = sortrows(Distance_data,3);

%%
% 
% tau_xdim = max(Tau_data(:,1));
% tau_ydim = max(Tau_data(:,2));
% tau_znum = max(Tau_data(:,3));

%find average centerpoint

ves_xmean = mean(vessel_data(:,1));
ves_ymean = mean(vessel_data(:,2));
ves_zmean = mean(vessel_data(:,3));

dist_xmean = mean(Distance_data(:,1));
dist_ymean = mean(Distance_data(:,2));
dist_zmean = mean(Distance_data(:,3));

ves_mean_sum = ves_xmean + ves_ymean + ves_zmean;
dist_mean_sum = dist_xmean + dist_ymean + dist_zmean;

vessel_data(:,1) = ((vessel_data(:,1)./ves_xdim).*(xup - xlow))+xlow;
Distance_data(:,1) = ((Distance_data(:,1)./dist_xdim).*(xup - xlow))+xlow;

if ves_mean_sum ~= dist_mean_sum

% TO FLIP X AND Y AXES
vessel_data(:,2) = yup - ((vessel_data(:,2)./ves_ydim).*(yup - ylow));
Distance_data(:,2) = ((Distance_data(:,2)./dist_ydim).*(yup - ylow))+ylow;

vessel_data(:,3) = zup - ((vessel_data(:,3)./ves_znum).*(zup - zlow));
Distance_data(:,3) = ((Distance_data(:,3)./dist_znum).*(zup - zlow))+zlow;

else
    vessel_data(:,2) = ((vessel_data(:,2)./ves_ydim).*(yup - ylow)) + ylow;
Distance_data(:,2) = ((Distance_data(:,2)./dist_ydim).*(yup - ylow))+ylow;

vessel_data(:,3) = ((vessel_data(:,3)./ves_znum).*(zup - zlow)) + ylow;
Distance_data(:,3) = ((Distance_data(:,3)./dist_znum).*(zup - zlow))+zlow;

end

%%
ves_xmean = mean(vessel_data(:,1));
ves_ymean = mean(vessel_data(:,2));
ves_zmean = mean(vessel_data(:,3));

dist_xmean = mean(Distance_data(:,1));
dist_ymean = mean(Distance_data(:,2));
dist_zmean = mean(Distance_data(:,3));

ves_mean_sum = ves_xmean + ves_ymean + ves_zmean;
dist_mean_sum = dist_xmean + dist_ymean + dist_zmean;


%%

%% Find Line of best fit equation and dataset 
waitbar(0.1,f,sprintf('Getting Distances for Vessel %d/%d: Finding Line of Best Fit',h,k));

% Find minimum axes values of vessel data for improvement of fit line. 
minx = min(vessel_data(:,1));
miny = min(vessel_data(:,2));
minz = min(vessel_data(:,3));

% Finds values of coefficients for polynomial equation of best fit line.
fitdata = [0 0 0];

% finds best fit line fo X-Y plane
fxy = fit((vessel_data(:,1)-minx),(vessel_data(:,2)-miny),'poly3');

% finds best fit line fo X-Z plane
fxz = fit((vessel_data(:,1)-minx),(vessel_data(:,3)-minz),'poly3');

% Get coefficients from equation output
fxy_values = coeffvalues(fxy);
fxz_values = coeffvalues(fxz);

%fit_res_1 = (max(vessel_data(:,1)) - min(vessel_data(:,1)))/fit_res;

x_fit = ((min(vessel_data(:,1))-minx-5):1:(max(vessel_data(:,1))-minx+5))';
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
nearfit = dsearchn(fitdata(:,1:3),vessel_data2(:,1:3));

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

if testplot == 1
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

else
end

figure()
scatter3(Distance_data(:,1),Distance_data(:,2),Distance_data(:,3),'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha',0.1);
hold on
scatter3(vessel_data(:,1),vessel_data(:,2),vessel_data(:,3));

%% Find Tau distance along Blood Vessel
%waitbar(0.4,f,sprintf('Getting Distances for Vessel %d/%d: Finding Distance Along Vessel',h,k));
% Create matrix for distance along vessel
dist = zeros(length(Distance_data),1);

% Find the nearest point in the best fit line dataset to each point in the
% tau dataset. 
near = knnsearch(fitdata(:,1:3),Distance_data(:,1:3));

% Calculate arclength between the end of the best fit line and each nearest
% point on the best fit line to get distance along the vessel. 
for i = 1:length(Distance_data)
    if length(fitdata(1:near(i,:),1)) == 1
        Distance_data(i,5) = [0];
    else
    % Calculate arclength
    dist(i,:) = arclength(fitdata(1:near(i,:),1),fitdata(1:near(i,:),2),fitdata(1:near(i,:),3));
    % Store distance data in new column of tau dataset. 
    Distance_data(i,5) = dist(i,:);
    end
end

% Sort tau dataset by distance down the vessel for plotting. 
Distance_data = sortrows(Distance_data,5);

%% Export Data
waitbar(0.9,f,sprintf('Getting Distances for Vessel %d/%d: Exporting Data',h,k));
% % norm_int = Tau_data(:,4)./median(Tau_data(:,4));
% % Tau_data = [Tau_data(:,1:3) norm_int  Tau_data(:,5:6)];
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures
filename = fcontent(h).name;
[pathstr,name,ext] = fileparts(filename);

final_data = Distance_data;
ves_num = compose('_%02d', h);
txt_name = append(folder_dir,'/',name,'_volume_data','.txt');
xlsx_name = append(folder_dir,'/',name,'_volume_data', '.xlsx');

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
%end

waitbar(1,f,'Finished!');

toc

close(f)

