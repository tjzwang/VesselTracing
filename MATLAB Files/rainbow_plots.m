
clc; clear all; close all;

% Input intensity distance calculation data to get rainbow plots for each
% vessel. 

% Input data folder name (INCLUDE FOLDER PATH AND NAME):
sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Intensity_distance_data_imp_fit';

%%
% Set Moving Average Size
Int_moving_ave_size = [  300  ];

Int_surfplot_moving_ave_size = [  500  ];
%%

f = waitbar(0,'Reading Data');
fcontent = dir(fullfile(sourcefolder, '*.txt')); %fcontent is a column vector of structures

Data = {};

n = length(fcontent');
n_ves = n;

for i = 1:n

   filename = fcontent(i).name;
   Data_single = readmatrix(fullfile(sourcefolder,filename));
   Data{i,1} = Data_single;
   Data_single = [];

   waitbar(i/n,f,sprintf('Reading Tau Data: %d %%',floor(i/n*100)));

end

%%
Tau_data = Data;
for i = 1:length(Tau_data)
Tau_data{i}(Tau_data{i}(:,6)>60,:) = [];
end
clear Data;

%% Plot tau intensity vs. distance along blood vessel vs. distance from blood vessel surface

for ii = 1:n_ves
% Create dataset used for a 3D surface plot containing distance from the
% surface (X), distance along the vessel (Y), and tau intensity (Z).
% distance from surface was rounded so plots can be generated at each
% distance increment away from vessel surface.
surfplotdata = [ round(Tau_data{ii}(:,6)) Tau_data{ii}(:,5) Tau_data{ii}(:,4)];
surfplotdata = sortrows(surfplotdata, 1);

surfplotdata(surfplotdata(:,1)<0,:)=[];

% Find the indexes in the surface plot data where the distance from the
% blood vessel surface is equal to each integer. 
idx = zeros(length(surfplotdata),max(surfplotdata(:,1)));
i = [];
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
i = [];
for i = 1:(max(surfplotdata(:,1)+1))
    idx2 = idx(:,i);
    idx2 = idx2(any(idx2,2),:);
meansurfplotdata(i,:,:) = { [surfplotdata(idx2,1)  surfplotdata(idx2,2) movmean(surfplotdata(idx2,3),Int_surfplot_moving_ave_size)] };
end

% Plot a line for each distance from the blood vessel on the same 3D plot.
figure()
% set colormap for use in plot
colors = flipud(jet(length(meansurfplotdata)));
%%
i = [];
for i = 1:length(meansurfplotdata)
% Plot 3D scatter plots for intensity vs. distance along vessel for each
% increment of distance from the blood vessel surface
scatter3(meansurfplotdata{i,:}(:,1),meansurfplotdata{i,:}(:,2),meansurfplotdata{i,:}(:,3),0.3,meansurfplotdata{i,:}(:,3));
colormap(jet);
hold on
end
%%
grid on

% Set axes to equal increments
%axis equal

% Label Axes
lx = xlabel('Distance From Blood Vessel Surface (\mum)');
ly = ylabel('Distance Along Blood Vessel (\mum)');
lz = zlabel('Tau Intensity (Norm)');
title(sprintf('Vessel %g',ii));
i = [];
end
