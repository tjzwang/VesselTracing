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

% TAU TANGLE DATA (.csv):

tau_folder = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\Tangle_data';

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

% Tau X dimension:

tau_xdim = [ 3949 ];

% Distance Transform X dimension:

dist_xdim = [ 3949 ];

% Vessel Y dimension:

ves_ydim = [ 3934 ];

% Tau Y dimension:

tau_ydim = [ 3934 ];

% Distance Transform Y dimension:

dist_ydim = [ 3934 ];

% Number of vessel z-stack slices:

ves_znum = [ 169 ];

% Number of Tau z-stack slices:

tau_znum = [ 169 ];

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

folder_dir = 'C:\Users\exx\Desktop\Zach - Active Analysis\Human2267\Tangle_distance_data';

% Adjusted Data Export Name? (DO NOT INCLUDE FILE TYPE OR VESSEL NUMBER!)
% ex. Human2131_int_dist_data
% Vessel number will be appended to file name (ex.
% Human2131_int_dist_data_1)

final_name =  '2699_tangle_distance_data' ;







%%
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures

k = length(fcontent');

for h = 1:k

f = waitbar(0,'Reading Tau Data');
fcontent = dir(fullfile(tau_folder, '*.csv')); %fcontent is a column vector of structures

sourcefolder = tau_folder;

n = length(fcontent');

   filename = fcontent(1).name;
   Tau_data_single = [];
   Tau_data_single = readmatrix(fullfile(sourcefolder,filename));
   Tau_data = [round(Tau_data_single(:,29:31)) Tau_data_single(:,1)];

Tau_data_orig = Tau_data;

% Import distance transform data

waitbar(0,f,'Reading Distance Data');
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = distance_folder;

n = length(fcontent');

   filename = fcontent(h).name;
   distance_data_single = [];
   distance_data_single = readmatrix(fullfile(sourcefolder,filename));
   Distance_data = distance_data_single;

% Import vessel data:

waitbar(0,f,'Reading Vessel Data');
fcontent = dir(fullfile(vessel_folder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = vessel_folder;


n = length(fcontent');



   filename = fcontent(h).name;
   vessel_data_single = [];
   vessel_data_single = readmatrix(fullfile(sourcefolder,filename));
   vessel_data = vessel_data_single;



% % k = 1;
% % h = 1;

f = waitbar(0,sprintf('Getting Distances for Vessel %d/%d: Preparing Data',h,k));

new_Tau_data = Tau_data;

C = find(~ismember(Tau_data(:,1:3),Distance_data(:,1:3),'rows'));

new_Tau_data(C,:) = [];

C = find(~ismember(Distance_data(:,1:3),Tau_data(:,1:3),'rows'));

Distance_data(C,:) = [];

%%

n = length(new_Tau_data(:,1));
new_Tau_data = sortrows(new_Tau_data,1);
new_Tau_data = sortrows(new_Tau_data,2);
new_Tau_data = sortrows(new_Tau_data,3);
Distance_data = sortrows(Distance_data,1);
Distance_data = sortrows(Distance_data,2);
Distance_data = sortrows(Distance_data,3);

new_Tau_data(:,6) = Distance_data(:,4);

Tau_data = new_Tau_data;

%%
% 
% tau_xdim = max(Tau_data(:,1));
% tau_ydim = max(Tau_data(:,2));
% tau_znum = max(Tau_data(:,3));

vessel_data(:,1) = ((vessel_data(:,1)./ves_xdim).*(xup - xlow))+xlow;
Tau_data(:,1) = ((Tau_data(:,1)./tau_xdim).*(xup - xlow))+xlow;
Distance_data(:,1) = ((Distance_data(:,1)./dist_xdim).*(xup - xlow))+xlow;

% % vessel_data(:,2) = ((vessel_data(:,2)./ves_ydim).*(yup - ylow))+ylow;
% % Tau_data(:,2) = ((Tau_data(:,2)./tau_ydim).*(yup - ylow))+ylow;
% % Distance_data(:,2) = ((Distance_data(:,2)./dist_ydim).*(yup - ylow))+ylow;
% % 
% % vessel_data(:,3) = ((vessel_data(:,3)./ves_znum).*(zup - zlow))+zlow;
% % Tau_data(:,3) = ((Tau_data(:,3)./tau_znum).*(zup - zlow))+zlow;
% % Distance_data(:,3) = ((Distance_data(:,3)./dist_znum).*(zup - zlow))+zlow;

% TO FLIP X AND Y AXES
vessel_data(:,2) = yup - ((vessel_data(:,2)./ves_ydim).*(yup - ylow));
Tau_data(:,2) = ((Tau_data(:,2)./tau_ydim).*(yup - ylow))+ylow;
Distance_data(:,2) = ((Distance_data(:,2)./dist_ydim).*(yup - ylow))+ylow;

vessel_data(:,3) = zup - ((vessel_data(:,3)./ves_znum).*(zup - zlow));
Tau_data(:,3) = ((Tau_data(:,3)./tau_znum).*(zup - zlow))+zlow;
Distance_data(:,3) = ((Distance_data(:,3)./dist_znum).*(zup - zlow))+zlow;

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

x_fit = ((min(vessel_data(:,1))-minx-5):0.1:(max(vessel_data(:,1))-minx+5))';
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
scatter3(Tau_data(:,1),Tau_data(:,2),Tau_data(:,3),'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha',0.1);
hold on
scatter3(vessel_data(:,1),vessel_data(:,2),vessel_data(:,3));

%% Find Tau distance along Blood Vessel
%waitbar(0.4,f,sprintf('Getting Distances for Vessel %d/%d: Finding Distance Along Vessel',h,k));
% Create matrix for distance along vessel
dist = zeros(length(Tau_data),1);

% Find the nearest point in the best fit line dataset to each point in the
% tau dataset. 
near = knnsearch(fitdata(:,1:3),Tau_data(:,1:3),Distance="fasteuclidean");

% Calculate arclength between the end of the best fit line and each nearest
% point on the best fit line to get distance along the vessel. 
parfor i = 1:length(Tau_data(:,1))
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

% Plot the moving average of tau intesity vs. distance down the blood vessel

if intdistplots == 1
% Calculate the moving average of intesity for tau dataset
means4 = movmean(Tau_data(:,4),Int_moving_ave_size);


%Calculate relative moving average based on the minimum value calculated.
relativemeans4 = means4./(min(means4));

%Plot relative moving average vs. distance along vessel
figure()
plot(Tau_data(:,5),relativemeans4);
xlabel('Distance Along Vessel');
ylabel('Relative Tau Intesity');

%Plot absolute moving average vs. distance along vessel.
figure()
plot(Tau_data(:,5),means4);
xlabel('Distance Along Vessel');
ylabel('Absolute Tau Intesity');

else
end

%% IF NO DISTANCE TRANSFORM AVAILABLE

% % %% Find nearest points on surface of blood vessel to tau data 
% % 
% % % Find the nearest point in the blood vessel XYZ dataset and the tau
% % % dataset to find the sortest distance between tau and blood vessel
% % % surface.
% % 
% % nearsurf = dsearchn(gpuArray(vessel_data(:,1:3)),gpuArray(Tau_data(:,1:3)));
% % 
% % %% Calculate the distance between the tau data and their corresponding nearest point on blood vessel surface
% % Tau_data(:,6) = zeros(length(Tau_data), 1);
% %  for i = 1:length(Tau_data)
% %     % Calculate pythagorean distance between nearest points
% %     distsurf = sqrt(    (((Tau_data(i,1)-vessel_data(nearsurf(i),1)))^2) + (((Tau_data(i,2)-vessel_data(nearsurf(i),2)))^2) +   (((Tau_data(i,3)-vessel_data(nearsurf(i),3)))^2) );
% %     % add distance from surface to the tau dataset
% %     Tau_data(i,6) = distsurf;
% %  end

%% Eliminate tau data over 100 microns from vessel surface

waitbar(0.8,f,sprintf('Getting Distances for Vessel %d/%d: Finalizing Data',h,k));
index = find(Tau_data(:,6) > 100);
Tau_data(index,:) = [];
maxtest = max(Tau_data(:,6));
figure()
scatter3(Tau_data(:,1), Tau_data(:,2), Tau_data(:,3));

%% Plot tau intensity vs. distance along blood vessel vs. distance from blood vessel surface

if intdistplots == 1
% Create dataset used for a 3D surface plot containing distance from the
% surface (X), distance along the vessel (Y), and tau intensity (Z).
% distance from surface was rounded so plots can be generated at each
% distance increment away from vessel surface.
surfplotdata = [ round(Tau_data(:,6)) Tau_data(:,5) Tau_data(:,4)];
surfplotdata = sortrows(surfplotdata, 1);

surfplotdata(surfplotdata(:,1)<0,:)=[];

% Find the indexes in the surface plot data where the distance from the
% blood vessel surface is equal to each integer. 
idx = zeros(length(surfplotdata),max(surfplotdata(:,1)));
for i = min(surfplotdata(:,1)):max(surfplotdata(:,1))
    indexes = [];
    indexes = find(surfplotdata(:,1) == i);
   for n = 1:length(indexes)
   idx(n,i+1) = indexes(n);
   end
end

% remove empty rows from index matrix.
idx = idx(any(idx,2),:);
meansurfplotdata = {};

% Group data into cell array by distance from the surface of the blood
% vessel where all data is a certain integer distance from the blood vessel
% surface is placed in a separate cell.
for i = 1:(max(surfplotdata(:,1)+1))
    idx2 = idx(:,i);
    idx2 = idx2(any(idx2,2),:);
meansurfplotdata(i,:,:) = { [surfplotdata(idx2,1)  surfplotdata(idx2,2) movmean(surfplotdata(idx2,3),Int_surfplot_moving_ave_size)] };
end

% Plot a line for each distance from the blood vessel on the same 3D plot.
figure()
% set colormap for use in plot
colors = flipud(jet(length(meansurfplotdata)));
for i = 1:length(meansurfplotdata)
% Plot 3D scatter plots for intensity vs. distance along vessel for each
% increment of distance from the blood vessel surface
scatter3(meansurfplotdata{i,:}(:,1),meansurfplotdata{i,:}(:,2),meansurfplotdata{i,:}(:,3),0.3,meansurfplotdata{i,:}(:,3));
colormap(jet);
hold on
end

grid on

% Set axes to equal increments
if axes == 1
axis equal
else
end

% Label Axes
lx = xlabel('Distance From Blood vessel Surface');
ly = ylabel('Distance Along Blood Vessel');
lz = zlabel('Tau Intensity');

else
end


%% Export Data
waitbar(0.9,f,sprintf('Getting Distances for Vessel %d/%d: Exporting Data',h,k));
% % norm_int = Tau_data(:,4)./median(Tau_data(:,4));
% % Tau_data = [Tau_data(:,1:3) norm_int  Tau_data(:,5:6)];
fcontent = dir(fullfile(distance_folder, '*.txt')); %fcontent is a column vector of structures
filename = fcontent(h).name;
[pathstr,name,ext] = fileparts(filename);

name_parts = regexp(name,'_',split);
new_name = append(name_parts{1},'_',name_parts{2},'_',name_parts{3});

final_data = Tau_data;

txt_name = append(folder_dir,'/',new_name,'tangle_dist_data','.txt');
xlsx_name = append(folder_dir,'/',new_name,'tangle_dist_data', '.xlsx');

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

%close all;

Tau_data = [];
Distance_data = [];
vessel_data = [];
end

waitbar(1,f,'Finished!');

toc

close(f)

