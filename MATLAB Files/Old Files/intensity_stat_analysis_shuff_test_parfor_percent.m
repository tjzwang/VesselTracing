clc; clearvars; close all;

%% Import Data


% Input data folder name (INCLUDE FOLDER PATH AND NAME):
sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2302\Intensity_distance_data_imp_fit';

% Set Moving Average Size (Recommended Value: 50)
Int_moving_ave_size = [  50  ];

% Max Distance From Vessel (Recommended Value: 30):

ves_dist = [ 30 ];

% Resolution in Microns (Recommended Value: 1):

res = [ 1 ];

% Frame width (# Points) (Recommended Value: 20):

frame_w = [ 20 ];

% Surface Data Definition (Recommended Value: 3): (Within # Microns from Vessel Surface)

surf_dist = [ 3 ];

% Standard Deviation or Percent Binning? 
% (1 = SD, 2 = %)

bintype = [ 2 ];

% Error type? (0 = standard error, 1 = anova group interval)

errortype = [0];

% Number of percent bins:

numbins = [ 10 ];

% Shuffled Data Output Name (Include Name + Path):

outputname = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2290\Shuffled_data_imp_fit\Human2131_percentage';


%%

% # Bins along vessel
nbins_length = 5;

y = 0;
for h = 1:1

%%

f = waitbar(0,'Calculating Statistics: Reading Data 0 %');

fcontent = dir(fullfile(sourcefolder, '*.txt')); %fcontent is a column vector of structures

l1 = 1;

% Open files and load data from folder
% If multiple files exist, their data will be stacked

n = length(fcontent');

for i = 1:length(fcontent')

   filename = fcontent(i).name;
   Tau_data_single = [];
   Tau_data_single = readmatrix(fullfile(sourcefolder,filename));
   [r, c] = size(Tau_data_single);
   l2 = l1 + r - 1;
   Tau_data(l1:l2,:) = Tau_data_single;
   l1 = r + l1;

    waitbar(0,f,sprintf('Calculating Statistics: Reading Data %d %%',floor(i/n*100)));

end

% % % % Set Moving Average Size
% % % Int_moving_ave_size = [  50  ];
% % % 
% % % % Max Distance From Vessel:
% % % 
% % % ves_dist = [ 30 ];
% % % 
% % % % Resolution in Microns?:
% % % 
% % % res = [ 1 ];
% % % 
% % % % Frame width (# Points):
% % % 
% % % frame_w = [ 10 ];
% % % 
% % % % Surface Data Definition: (Within # Microns from Vessel Surface)
% % % 
% % % surf_dist = [ 2 ];
% % % 
% % % % # Bins along vessel
% % % nbins_length = 5;

% % Set number of bins by increment of intensity
% nbins = round((max(Tau_data(:,4)) - min(Tau_data(:,4)))/50);

%% Sort Data into groups by intensity

 waitbar(0.2,f,sprintf('Calculating Statistics: Sorting Data %d %%',20));

% Remove data where Tau is within blood vessel and greater than ves_dist
% away
idx = find( Tau_data(:,6) < 0);
Tau_data(idx,:) = [];
idx = find(Tau_data(:,6) > ves_dist );
Tau_data(idx,:) = [];

nbins = ves_dist/res;

% % Set number of bins by increment of intensity
% nbins = round((max(Tau_data(:,4)) - min(Tau_data(:,4)))/50);

group_Tau_data = Tau_data;

% sort the ungrouped data by intensity
group_Tau_data = sortrows(group_Tau_data,4);

% get bin edges for grouping the tau data by intensity
[N,edges7,bin] = histcounts(group_Tau_data(:,4),nbins);


group_Tau_data = [group_Tau_data(:,1:6) bin];

%% ANOVA on distance from vessel vs. intensity

% % [p,tbl,stats] = anovan(group_Tau_data(:,6),group_Tau_data(:,7),'display','off');
% % 
% % c = multcompare(stats);

%% Group Data by distance from vessel

dist_Tau_data = sortrows(group_Tau_data,6);

% get bin edges for grouping the tau data by intensity
[N,edges8,bin] = histcounts(dist_Tau_data(:,6),nbins);

group_Tau_data = [dist_Tau_data(:,1:7) bin];

%% ANOVA on Intensity vs. distance from vessel

% % [p,tbl,stats] = anovan(group_Tau_data(:,4),group_Tau_data(:,8),'display','off');
% % figure()
% % c = multcompare(stats);

%% Group Data by distance along vessel

group_Tau_data = sortrows(group_Tau_data,5);

% get bin edges for grouping the tau data by intensity
[N,edges9,bin] = histcounts(group_Tau_data(:,5),nbins_length);


group_Tau_data = [group_Tau_data(:,1:8) bin];

%% ANOVA on Intensity vs. distance along vessel

% % [p,tbl,stats] = anovan(group_Tau_data(:,4),{group_Tau_data(:,8),group_Tau_data(:,9)},'display','off');
% % figure()
% % c = multcompare(stats,"Dimension", [1 2]);

%% Group data by intensity standard deviation

 waitbar(0.4,f,sprintf('Calculating Statistics: Grouping Data %d %%',40));

% filter data to only include tau within surf_dist microns of vessel.
surface_idx = find(group_Tau_data(:,6) < surf_dist & group_Tau_data(:,6) >= 0);
surface_data = group_Tau_data(surface_idx,:);
surface_data_back = group_Tau_data(surface_idx,:);

% % [m,n] = size(surface_data);
% % shuffle_idx = [randperm(m)]';
% % surface_data1(shuffle_idx,1:3) = surface_data(:,1:3);
% % surface_data1(shuffle_idx,5:6) = surface_data(:,5:6);
% % 
% % [m,n] = size(surface_data);
% % shuffle_idx = [randperm(m)]';
% % surface_data2(shuffle_idx,1:3) = surface_data(:,1:3);
% % surface_data2(shuffle_idx,5:6) = surface_data(:,5:6);
% % 
% % [m,n] = size(surface_data);
% % shuffle_idx = [randperm(m)]';
% % surface_data3(shuffle_idx,1:3) = surface_data(:,1:3);
% % surface_data3(shuffle_idx,5:6) = surface_data(:,5:6);
% % 
% % surface_data = (surface_data1 + surface_data2 + surface_data3)/3;

% % [m,n] = size(surface_data);
% % shuffle_idx1 = [randperm(m)]';
% % shuffle_idx2 = [randperm(m)]';
% % shuffle_idx3 = [randperm(m)]';
% % 
% % shuffle_idx = round((shuffle_idx1 + shuffle_idx2 + shuffle_idx3)/3);
% % maxshuff = max(shuffle_idx);
% % 
% % surface_data = surface_data(shuffle_idx,:);

group_Tau_data_back = group_Tau_data;

[m,n] = size(group_Tau_data);
shuffle_idx1 = [randperm(m)]';
shuffle_idx2 = [randperm(m)]';
shuffle_idx3 = [randperm(m)]';
shuffle_idx = round((shuffle_idx1 + shuffle_idx2 + shuffle_idx3)/3);

group_Tau_data = group_Tau_data(shuffle_idx,:);
group_Tau_data(:,4) = group_Tau_data_back(:,4);

%%
% % for i = length(surface_data)
% % x(i,:) = surface_data_back(i,4);
% % surface_data(i,4) = surface_data_back(i,4);
% % end
surface_data(:,4) = surface_data_back(:,4);
%%
for i = length(surface_data)
group_Tau_data(surface_idx(i),1:6) = surface_data(i,1:6);
end

surface_data = sortrows(surface_data,5);
% Calculate standard deviation and median
tauStd = std(group_Tau_data(:,4));
% tauStd = std(surface_data(:,4));
med_int = median(group_Tau_data(:,4));

d = surface_data(:,4);
s = tauStd;
m = med_int;

if bintype == 1

% Find indexes of vaules for each bin based on 1/2 standard deviation
idx = padcat(find( d < m-3.75*s),find( d < m-3.25*s & d>= m-3.75*s),find( d < m-2.75*s & d>= m-3.25*s),find( d < m-2.25*s & d>= m-2.75*s),find( d < m-1.75*s & d>= m-2.25*s),find( d < m-1.25*s & d>= m-1.75*s),find( d < m-0.75*s & d>= m-1.25*s),find( d < m-0.25*s & d>= m-0.75*s),find( d < m+0.25*s & d>= m-0.25*s),find( d < m+0.75*s & d>= m+0.25*s),find( d < m+1.25*s & d>= m+0.75*s),find( d < m+1.75*s & d>= m+1.25*s),find( d < m+2.25*s & d>= m+1.75*s),find( d < m+2.75*s & d>= m+2.25*s),find( d < m+3.25*s & d>= m+2.75*s),find( d < m+3.75*s & d>= m+3.25*s),find( d < m+4.25*s & d>= m+3.75*s),find(d>= m+4.25*s));

% Put bin number for each row of tau data in column 10
[r,c] = size(idx);
for i = 1:c
    idx2 = idx(:,i);
    idx2(isnan(idx2))=[];
    surface_data(idx2,10) = i;
end

% Bin #s
% 1 =                        <   -3.75   SD
% 2 =   -3.75    >=    and   <   -3.25   SD
% 3 =   -3.25    >=    and   <   -2.75   SD
% 4 =   -2.75    >=    and   <   -2.25   SD
% 5 =   -2.25    >=    and   <   -1.75   SD
% 6 =   -1.75    >=    and   <   -1.25   SD
% 7 =   -1.25    >=    and   <   -0.75   SD
% 8 =   -0.75    >=    and   <   -0.25   SD
% 9 =   -0.25    >=    and   <    0.25   SD
% 10 =  +0.25    >=    and   <   +0.75   SD
% 11 =  +0.75    >=    and   <   +1.25   SD
% 12 =  +1.25    >=    and   <   +1.75   SD
% 13 =  +1.75    >=    and   <   +2.25   SD
% 14 =  +2.25    >=    and   <   +2.75   SD
% 15 =  +2.75    >=    and   <   +3.25   SD
% 16 =  +3.25    >=    and   <   +3.75   SD
% 17 =  +3.75    >=    and   <   +4.25   SD
% 18 =  +4.25    >=                      SD

%%
figure()
histogram(surface_data(:,10));

%%

% sort tau data by distance along blood vessel
surface_data = sortrows(surface_data,5);

% Find regions of similar tau density lasting a certain number of points
% (f)
f = frame_w; % width of region with similar tau intensity
n = 1;  
z = [];
sim_bin_data = {};
std_sens = 1;
% Find regions of similar tau intensity
% Similar tau intensity is defined as a region with a standard deviation less than 0.01. 
for i = 1:max(surface_data(:,10))
for z = 1:length(surface_data)
    if z <= length(surface_data) - f
    frame_data = surface_data(z:(z+f),10);
    frame_std = std(frame_data);
    frame_median = median(frame_data);
    frame_mean = mean(frame_data);
    % if frame_mean >= i - 0.1 && frame_mean <= i + 0.1
    if frame_median == i
    % if frame_std <= 0.1 && frame_median == i
    % if frame_mean >= i - 0.1 && frame_mean <= i + 0.1 && frame_std <= 0.1
        idx3(n,:) = [i z z+f];
        bin_edges(n,:) = [i surface_data(idx3(n,2),5) surface_data(idx3(n,3),5)];
        n = n+1;
    else
    end
    else
    end
end
end
% % for i = 1:length(idx3)
% %     if i ~= length(idx3)
% %     if idx3(i+1,2) == (idx3(i,2)+1) && idx3(i,1) == idx3(i+1,1)
% %        idx4(i,:) = [ idx3(i,1:2) idx3(i+1,3)];
% %     else
% %     end
% %     else
% %     end
% % end

else if bintype == 2

% Find regions of similar tau density lasting a certain number of points (f)
f = frame_w; % width of region with similar tau intensity
n = 1;  
z = [];
sim_bin_data = {};
std_sens = 1;
    for z = 1:f:length(surface_data)
    if z+f>length(surface_data)
    f = length(surface_data)-z;
    frame_data = surface_data(z:(z+f),4);
    frame_median_int = median(frame_data);
    idx3(z,:) = [frame_median_int z z+f];
    for k = 0:f
    bin_edges(z+k,1:3) = [frame_median_int surface_data(idx3(z,2),5) surface_data(idx3(z,3),5)];
    end

    else
    %if z <= length(surface_data) - f
    frame_data = surface_data(z:(z+f),4);
    frame_median_int = median(frame_data);
    idx3(z,:) = [frame_median_int z z+f];
    for k = 0:f
    bin_edges(z+k,1:3) = [frame_median_int surface_data(idx3(z,2),5) surface_data(idx3(z,3),5)];
    end
%     else
%         
%          bin_edges(z,:) = bin_edges(z-f,:);
%     end
    end
    end
    %%
    [bin_edges,sort_idx] = sortrows(bin_edges,1);
    [r,c] = size(bin_edges);
    bin_size = round(r/numbins); 
    lower_bound = 0;
    upper_bound = bin_size;
    for i = 1:numbins
        if i ~= 10
        bin_edges((lower_bound+1):upper_bound,1) = i;
        lower_bound = upper_bound;
        upper_bound = bin_size + upper_bound;
        else
        bin_edges((lower_bound+1):end,1) = i;
        end
    
    end
    bin_edges = bin_edges(sort_idx,:);
    surface_data(:,10) = bin_edges(:,1);
end
end

%%
% Append group number for distance from vessel to initial tau data
Tau_data = [group_Tau_data(:,1:6)  group_Tau_data(:,8)];
% sort data by distance from vessel
Tau_data = sortrows(Tau_data,6);

%% 
%% 
% Create cells with data between regions of similar intensity using bin
% edges
% mark each cell with its corresponding intensity bin number
sim_bin_data = {zeros(length(bin_edges),2)};
lb = length(bin_edges);
p = 1;

parfor_progress(lb); % Set the total number of iterations
idx4 = {zeros(lb,1)};
parfor i = 1:lb
    Tau_data2 = Tau_data;
    bin_edges2 = bin_edges;
    idx4{i} = {find(Tau_data2(:,5) >= bin_edges2(i,2) & Tau_data2(:,5) <= bin_edges2(i,3))};
    parfor_progress;
end
parfor_progress(0);
%% 

for i = 1:lb
    sim_bin_data(i,1:2) = [{Tau_data(cell2mat(idx4{i}),:)} bin_edges(i,1)];
end


% get bin numbers from cell array
int_bin = cell2mat(sim_bin_data(:,2));


% combine cells of array with the same intensity bin number
for i = 1:max(surface_data(:,10))
    idx5 = find(int_bin == i);
    sim_bin_data_comb(i,1:2) = [{unique(cat(1,sim_bin_data{idx5}),'rows')} i];
    idx5 = [];
end

% remove empty rows
sim_bin_data_comb(all(cellfun(@isempty, sim_bin_data_comb(:,1)),2),:) = [];

% take first column of cell array for ANOVA analysis
anova_bin_data = sim_bin_data_comb(:,1);

% run ANOVA for each bin of cell array comparing intensity with distance
% from vessel
[r,c] = size(sim_bin_data_comb);
if r ~= 1
for i = 1:r

[p,tbl,stats(i,:)] = anovan(anova_bin_data{i}(:,4),anova_bin_data{i}(:,7),'display','off');

figure()
c = multcompare(stats(i,:));

end

else

[p,tbl,stats(1,:)] = anovan(anova_bin_data{1}(:,4),anova_bin_data{1}(:,7),'display','off');

figure()
c = multcompare(stats(1,:));

end
%%
dist_groups = {{zeros(1,8)} zeros(1,1)};

cell = [];
int_group = cell2mat(sim_bin_data_comb(:,2));

% Create new cell array with cells for each distance from vessel bin
% including grouping variable grouping by intensity at vessel surface
for n = 1:nbins
    z = 1;
    for i = 1:length(anova_bin_data)
        idx6 = find(anova_bin_data{i}(:,7) == n);
        l = length(idx6);
        k = int_group(i)*ones(l,1);
        t = anova_bin_data{i}(idx6,:);
        cell(z:(z+l-1),:) = [anova_bin_data{i}(idx6,:) int_group(i)*ones(l,1)];
        dist_groups(n,:) = [{cell} n];
        z = z+l;
    end
    cell = [];
end

%%
title_edge = round(edges8,1);
%%

% run ANOVA comparing intensity at local intensity for each distance from
% vessel bin for the surface intensity bins. 
%%
if bintype == 1
bin_numbers = cell2mat(sim_bin_data_comb(:,2));
bin_numbers = [(bin_numbers - 9)/2];
else
bin_numbers = cell2mat(sim_bin_data_comb(:,2));
bin_numbers = bin_numbers.*10;
end
%%
figure(99)

if r ~= 1
for i = 1:length(dist_groups)
n = i;
[p,tbl,stats(i,:)] = anovan(dist_groups{i}(:,4),dist_groups{i}(:,8),'display','off');
set(0,'DefaultFigureVisible','off')
[c,m] = multcompare(stats(i,:));
set(0,'DefaultFigureVisible','on')
bin_numbers = unique(dist_groups{i}(:,8));

if bintype == 1
bin_numbers = [(bin_numbers - 9)/2];
else
bin_numbers = bin_numbers.*10;
end

m = [m bin_numbers];
[m_length,c] = size(m);
xdata(1:m_length,i) = [m(:,1)];
ydata(1:m_length,i) = [m(:,3)];

if errortype == 0
errors(1:m_length,i) = [m(:,2)];
else
    errors(1:m_length,i) = std(dist_groups{i}(:,4))/sqrt(length(dist_groups{i}(:,4)));
end
end

if h == 1
xdata1 = xdata;
xdata = [];
ydata1 = ydata;
ydata = [];
errors1 = errors;
errors = [];
else
if h == 2
xdata2 = xdata;
xdata = [];
ydata2 = ydata;
ydata = [];
errors2 = errors;
errors = [];
else
if h == 3
xdata3 = xdata;
xdata = [];
ydata3 = ydata;
ydata = [];
errors3 = errors;
errors = [];
else
if h == 4
xdata4 = xdata;
xdata = [];
ydata4 = ydata;
ydata = [];
errors4 = errors;
errors = [];
else
end
end
end
end
end
%%
%xdata = (xdata1 + xdata2 + xdata3 + xdata4)/4;
%ydata = (ydata1 + ydata2 + ydata3 + ydata4)/4;
%errors = (errors1 + errors2 + errors3 + errors4)/4;

xdata = xdata1;
ydata = ydata1;
errors = errors1;
%%
f=waitbar(0.8,sprintf('Calculating Statistics: Exporting Data %d %%',80));

writematrix(xdata,append(outputname,'_shuffle_xdata'));
writematrix(ydata,append(outputname,'_shuffle_ydata'));
writematrix(errors,append(outputname,'_shuffle_errors'));

[r,c] = size(xdata);

n = 1;
tiledlayout(3,ceil(length(xdata(:,1))/3));
for i = 1:r
nexttile
std_range = 0.25*ones(size(ydata(:,i)));
%subplot(ceil(nbins/6),6,i);
% % e = errorbar(xdata(i,:),edges8(1:(length(edges8)-1)),errors(i,:),'horizontal','.');
e = errorbar(edges8(1:(length(edges8)-1)),xdata(i,:),errors(i,:),'.');

%e = errorbar(xdata(:,i),ydata(:,i),std_range,std_range,errors(:,i),errors(:,i),'.');
e.MarkerSize = 15;
e.LineWidth = 1;
e.MarkerFaceColor = 'r';
e.MarkerEdgeColor = 'r';
e.CapSize = 0;
e.AlignVertexCenters = 'on';

e2.MarkerSize = 15;
e2.LineWidth = 1;
e2.Color = 'm';
e2.MarkerFaceColor = 'g';
e2.MarkerEdgeColor = 'g';
e2.CapSize = 0;
e2.AlignVertexCenters = 'on';

% % xlabel('Local Intensity (Norm)','FontSize',8);
% % ylabel('Distance From Vessel','FontSize',8);
% % set(gca, 'YDir','reverse');

ylabel('Local Intensity (Norm)','FontSize',8);
xlabel('Distance From Vessel','FontSize',8);

ylim([0.65 1.5]);

%xlim([0.5 1.4])
% % hold on;
% % e = errorbar(xdata(:,i),ydata(:,i),std_range,'.');
% % xlabel('Local Intensity (Norm)','FontSize',8);
% % ylabel('Surface Intensity (SD)','FontSize',8);
% % e.MarkerSize = 1;
% % e.LineWidth = 0.5;
% % e.Color = 'r';
% % e.MarkerFaceColor = 'r';
% % e.MarkerEdgeColor = 'r';
% % e.CapSize = 2;
% % e.AlignVertexCenters = 'on';
if bintype == 2
title(sprintf('%g - %g %%',(ydata(i,1)-10),ydata(i,1)));
sgtitle('Surface Intensity %') 
else
title(sprintf('%g - %g (SD)',(ydata(i,1)-0.25),(ydata(i,1)+0.25)));
sgtitle('Surface Intensity SD') 
end
end
end


