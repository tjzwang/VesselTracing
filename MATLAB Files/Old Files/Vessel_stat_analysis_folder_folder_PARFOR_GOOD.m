clc; clear all; close all;

%% Import Data

% If a .txt file (Recommended for speed):

% Input data folder name:
sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Intensity_distance_data';

% Input shuffled data for comparison (INCLUDE FOLDER PATH AND NAME)
shuff_xdata = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data\Human2131_shuffle_xdata');
shuff_ydata = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data\Human2131_shuffle_ydata');
shuff_errors = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data\Human2131_shuffle_errors');

% Set Moving Average Size
Int_moving_ave_size = [  50  ];

% Max Distance From Vessel:

ves_dist = [ 30 ];

% Resolution in Microns?:

res = [ 1 ];

% Frame width (# Points):

frame_w = [ 20 ];

% Distance to define surface data (in microns)?:

surf_d = [ 5 ];

% # Bins along vessel
nbins_length = 5;

% % Set number of bins by increment of intensity
% nbins = round((max(Tau_data(:,4)) - min(Tau_data(:,4)))/50);

%% Load Data

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
%% Sort Data into groups by intensity


waitbar(0.1,f,'Calculating Statistics: Preparing Data');

% Remove data where Tau is within blood vessel and greater than ves_dist
% away
idx = find( Tau_data(:,6) <= 0);
Tau_data(idx,:) = [];
idx = find(Tau_data(:,6) > ves_dist );
Tau_data(idx,:) = [];

% calculate number of bins based on specified resolution
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


waitbar(0.2,f,'Calculating Statistics: Running Initial ANOVA');

[p,tbl,stats] = anovan(group_Tau_data(:,6),group_Tau_data(:,7),'display','off');

c = multcompare(stats);

%% Group Data by distance from vessel

waitbar(0.3,f,'Calculating Statistics: Grouping Data');

dist_Tau_data = sortrows(group_Tau_data,6);

% get bin edges for grouping the tau data by intensity
[N,edges8,bin] = histcounts(dist_Tau_data(:,6),nbins);

group_Tau_data = [dist_Tau_data(:,1:7) bin];

%% ANOVA on Intensity vs. distance from vessel

[p,tbl,stats] = anovan(group_Tau_data(:,4),group_Tau_data(:,8),'display','off');
figure()
c = multcompare(stats);

%% Group Data by distance along vessel

group_Tau_data = sortrows(group_Tau_data,5);

% get bin edges for grouping the tau data by intensity
[N,edges9,bin] = histcounts(group_Tau_data(:,5),nbins_length);


group_Tau_data = [group_Tau_data(:,1:8) bin];

%% ANOVA on Intensity vs. distance along vessel

[p,tbl,stats] = anovan(group_Tau_data(:,4),{group_Tau_data(:,8),group_Tau_data(:,9)},'display','off');
figure()
c = multcompare(stats,"Dimension", [1 2]);

%% Group data by intensity standard deviation


waitbar(0.4,f,'Calculating Statistics: Finding Surface Data & Intensity Standard Deviation');

% Calculate standard deviation and median
tauStd = std(group_Tau_data(:,4));
med_int = median(group_Tau_data(:,4));

% filter data to only include tau within stated microns of vessel.
surface_idx = find(group_Tau_data(:,6) < surf_d);
surface_data = group_Tau_data(surface_idx,:);

d = surface_data(:,4);
s = tauStd;
m = med_int;

% Find indexes of values for each bin based on 1/2 standard deviation
idx = padcat(find( d < m-4*s),find( d < m-3.5*s & d>= m-4*s),find( d < m-3*s & d>= m-3.5*s),find( d < m-2.5*s & d>= m-3*s),find( d < m-2*s & d>= m-2.5*s),find( d < m-1.5*s & d>= m-2*s),find( d < m-1*s & d>= m-1.5*s),find( d < m-0.5*s & d>= m-1*s),find( d < m-0*s & d>= m-0.5*s),find( d < m+0.5*s & d>= m+0*s),find( d < m+1*s & d>= m+0.5*s),find( d < m+1.5*s & d>= m+1*s),find( d < m+2*s & d>= m+1.5*s),find( d < m+2.5*s & d>= m+2*s),find( d < m+3*s & d>= m+2.5*s),find( d < m+3.5*s & d>= m+3*s),find( d < m+4*s & d>= m+3.5*s),find(d>= m+4*s));

% Put bin number for each row of tau data in column 10
[r,c] = size(idx);
for i = 1:c
    idx2 = idx(:,i);
    idx2(isnan(idx2))=[];
    surface_data(idx2,10) = i;
end

% Bins #s
% 1 =                     <   -4     SD
% 2 =   -4    >=    and   <   -3.5   SD
% 3 =   -3.5  >=    and   <   -3     SD
% 4 =   -3    >=    and   <   -2.5   SD
% 5 =   -2.5  >=    and   <   -2     SD
% 6 =   -2    >=    and   <   -1.5   SD
% 7 =   -1.5  >=    and   <   -1     SD
% 8 =   -1    >=    and   <   -0.5   SD
% 9 =   -0.5  >=    and   <    0     SD
% 10 =   0    >=    and   <   +0.5   SD
% 11 =  +0.5  >=    and   <   +1     SD
% 12 =  +1    >=    and   <   +1.5   SD
% 13 =  +1.5  >=    and   <   +2     SD
% 14 =  +2    >=    and   <   +2.5   SD
% 15 =  +2.5  >=    and   <   +3     SD
% 16 =  +3    >=    and   <   +3.5   SD
% 17 =  +3.5  >=    and   <   +4     SD
% 18 =  +4    >=                     SD


%%
figure()
histogram(surface_data(:,10));

%%

waitbar(0.6,f,'Calculating Statistics: Binning Data');

% sort tau data by distance along blood vessel
surface_data = sortrows(surface_data,5);

% Find regions of similar tau density lasting a certain number of points
% (f)
fw = frame_w; % width of region with similar tau intensity
n = 1;  
z = [];
sim_bin_data = {};
std_sens = 1;

% Find regions of similar tau intensity
% Similar tau intensity is defined as a region with a standard deviation less than 0.01. 
% % for i = 1:max(surface_data(:,10))
% % for z = 1:length(surface_data)
% %     if z <= length(surface_data) - f
% %     frame_data = surface_data(z:(z+f),10);
% %     frame_std = std(frame_data);
% %     frame_median = median(frame_data);
% %     frame_mean = mean(frame_data);
% %     % if frame_mean >= i - 0.1 && frame_mean <= i + 0.1
% %     if frame_median == i
% %     % if frame_std <= 0.1 && frame_median == i
% %     % if frame_mean >= i - 0.1 && frame_mean <= i + 0.1 && frame_std <= 0.1
% %         idx3(n,:) = [i z z+f];
% %         bin_edges(n,:) = [i surface_data(idx3(n,2),5) surface_data(idx3(n,3),5)];
% %         n = n+1;
% %     else
% %     end
% %     else
% %     end
% % end
% % end

% Mark surface bins with their median intensity
for i = 1:max(surface_data(:,10))
y = 1;
for z = 1:length(surface_data)
    if y <= length(surface_data) - fw
    frame_data = surface_data(y:(y+fw),10);
    frame_std = std(frame_data);
    frame_median = median(frame_data);
    frame_mean = mean(frame_data);
    % if frame_mean >= i - 0.1 && frame_mean <= i + 0.1
    if frame_median == i
    % if frame_std <= 0.1 && frame_median == i
    % if frame_mean >= i - 0.1 && frame_mean <= i + 0.1 && frame_std <= 0.1
        idx3(n,:) = [i y y+fw];
        bin_edges(n,:) = [i surface_data(idx3(n,2),5) surface_data(idx3(n,3),5)];
        n = n+1;
    else
    end
    y = y+fw;
    else
    end
end
end

figure()
histogram(bin_edges(:,1));


% % for i = 1:length(idx3)
% %     if i ~= length(idx3)
% %     if idx3(i+1,2) == (idx3(i,2)+1) && idx3(i,1) == idx3(i+1,1)
% %        idx4(i,:) = [ idx3(i,1:2) idx3(i+1,3)];
% %     else
% %     end
% %     else
% %     end
% % end

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

%%

waitbar(0.6,f,'Calculating Statistics: Running Final ANOVA');

% remove empty rows
sim_bin_data_comb(all(cellfun(@isempty, sim_bin_data_comb(:,1)),2),:) = [];

% take first column of cell array for ANOVA analysis
anova_bin_data = sim_bin_data_comb(:,1);

% run ANOVA for each bin of cell array comparing intensity with distance
% from vessel
for i = 1:length(sim_bin_data_comb)

[p,tbl,stats(i,:)] = anovan(anova_bin_data{i}(:,4),anova_bin_data{i}(:,7),'display','off');

figure()
c = multcompare(stats(i,:));

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

bin_numbers = cell2mat(sim_bin_data_comb(:,2));
bin_numbers = [(bin_numbers - 9.5)/2];

for i = 1:length(dist_groups)
n = i;
[p,tbl,stats(i,:)] = anovan(dist_groups{i}(:,4),dist_groups{i}(:,8),'display','off');
set(0,'DefaultFigureVisible','off')
[c,m] = multcompare(stats(i,:));
set(0,'DefaultFigureVisible','on')
m = [m bin_numbers];
xdata(:,i) = [m(:,1)];
ydata(:,i) = [m(:,3)];
errors(:,i) = [m(:,2)];
end


% Plot ANOVA data

waitbar(0.8,f,'Calculating Statistics: Plotting Data');

figure(99)
% Horizontal
for i = 1:length(dist_groups)
std_range = 0.25*ones(size(ydata(:,i)));
subplot(ceil(nbins/6),6,i);
e = errorbar(xdata(:,i),ydata(:,i),errors(:,i),'horizontal','.');
%e = errorbar(xdata(:,i),ydata(:,i),std_range,std_range,errors(:,i),errors(:,i),'.');
e.MarkerSize = 15;
e.LineWidth = 1;
e.MarkerFaceColor = 'r';
e.MarkerEdgeColor = 'r';
e.CapSize = 0;
e.AlignVertexCenters = 'on';
xlabel('Local Intensity (Norm)','FontSize',8);
ylabel('Surface Intensity (SD)','FontSize',8);

%xlim([0.5 3.5])

xlim([(min(xdata,[],'all')-0.4) (max(xdata,[],'all')+0.4)]);
ylim([(min(ydata,[],'all')-0.4) (max(ydata,[],'all')+0.4)]);

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
title(sprintf('%g - %g \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
sgtitle('Distance From Vessel') 
end


[r,c] = size(xdata);

figure(100)

% % plot_width = 6;
% % tiledlayout(ceil(nbins/6),plot_width);

edge_length = length(edges8);
n = 1;

bin_num = length(sim_bin_data_comb(:,2));

for i = 1:r
subplot(ceil(bin_num/3),3,i);
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


% % xlabel('Local Intensity (Norm)','FontSize',8);
% % ylabel('Distance From Vessel','FontSize',8);
% % set(gca, 'YDir','reverse');

ylabel('Local Intensity (Norm)','FontSize',8);
xlabel('Distance From Vessel','FontSize',8);

ylim([(min(ydata,[],'all')-0.1) (max(ydata,[],'all')+0.1)]);
xlim([-2 (ves_dist+1)]);


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
title(sprintf('%g (SD)',ydata(i,1)));
sgtitle('Surface Intensity SD') 
% % nexttile
end



% % % % Load shuffled data for comparison
% % % shuff_xdata = readmatrix('/Users/zacharyhoglund/Desktop/Arivis/Human2131/MATLAB Test Shuffle Data/Human2131_test_shuffle_xdata');
% % % shuff_ydata = readmatrix('/Users/zacharyhoglund/Desktop/Arivis/Human2131/MATLAB Test Shuffle Data/Human2131_test_shuffle_ydata');
% % % shuff_errors = readmatrix('/Users/zacharyhoglund/Desktop/Arivis/Human2131/MATLAB Test Shuffle Data/Human2131_test_shuffle_errors');


% Plot shuffled data along with the true data for comparison
figure(99)
plot_width = 6;
%tiledlayout(ceil(nbins/plot_width),plot_width);

% Horizontal
for i = 1:length(dist_groups)

subplot(ceil(nbins/6),plot_width,i);
std_range = 0.25*ones(size(ydata(:,i)));
%subplot(ceil(nbins/6),6,i);
e = errorbar(xdata(:,i),ydata(:,i),errors(:,i),'horizontal','.');
hold on
e2 = errorbar(shuff_xdata(:,i),shuff_ydata(:,i),shuff_errors(:,i),'horizontal','.');

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


xlabel('Local Intensity (Norm)','FontSize',8);
ylabel('Surface Intensity (SD)','FontSize',8);

xlim([(min(xdata,[],'all')-0.4) (max(xdata,[],'all')+0.4)]);
ylim([(min(ydata,[],'all')-0.4) (max(ydata,[],'all')+0.4)]);
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
title(sprintf('%g - %g \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
sgtitle('Distance From Vessel') 

if i == plot_width

lg = legend('Real','Shuffled','location','northeast','orientation','horizontal');
%lg.Location = 'northeastoutside';

else
end

%nexttile

end


[r,c] = size(xdata);
[r2,c2] = size(shuff_ydata);

figure(100)
plot_width = 6;
%tiledlayout(ceil(nbins/6),plot_width);

edge_length = length(edges8);
n = 1;

bin_num = length(sim_bin_data_comb(:,2));

for i = 1:r

subplot(ceil(bin_num/3),3,i);

std_range = 0.25*ones(size(ydata(:,i)));
%subplot(ceil(nbins/6),6,i);
% % e = errorbar(xdata(i,:),edges8(1:(length(edges8)-1)),errors(i,:),'horizontal','.');
e = errorbar(edges8(1:(length(edges8)-1)),xdata(i,:),errors(i,:),'.');

if n <= r2
idx = find(shuff_ydata(:,1) == ydata(i,1));
if length(idx) == 1
    hold on
% % e2 = errorbar(shuff_xdata(n,:),edges8(1:(length(edges8)-1)),shuff_errors(n,:),'horizontal','.');
e2 = errorbar(edges8(1:(length(edges8)-1)),shuff_xdata(idx,:),shuff_errors(idx,:),'.');

e2.MarkerSize = 15;
e2.LineWidth = 1;
e2.Color = 'm';
e2.MarkerFaceColor = 'g';
e2.MarkerEdgeColor = 'g';
e2.CapSize = 0;
e2.AlignVertexCenters = 'on';

hold on
else
end
else
end

%e = errorbar(xdata(:,i),ydata(:,i),std_range,std_range,errors(:,i),errors(:,i),'.');
e.MarkerSize = 15;
e.LineWidth = 1;
e.MarkerFaceColor = 'r';
e.MarkerEdgeColor = 'r';
e.CapSize = 0;
e.AlignVertexCenters = 'on';

% % xlabel('Local Intensity (Norm)','FontSize',8);
% % ylabel('Distance From Vessel','FontSize',8);
% % set(gca, 'YDir','reverse');

ylabel('Local Intensity (Norm)','FontSize',8);
xlabel('Distance From Vessel','FontSize',8);

ylim([(min(xdata,[],'all')-0.1) (max(xdata,[],'all')+0.1)]);
xlim([-2 (ves_dist+1)]);


if i == ceil(bin_num/3)

lg = legend('Real','Shuffled','location','northeast','orientation','horizontal');
%lg.Location = 'northeastoutside';

else
end

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
title(sprintf('%g (SD)',ydata(i,1)));
sgtitle('Surface Intensity SD') 
%nexttile
end

%%

% % figure()
% % 
% % for i = 1:length(sim_bin_data_comb(:,2))
% % 
% %     subplot(1,length(sim_bin_data_comb(:,2)),i)
% % 
% %     data = sortrows(cell2mat(sim_bin_data_comb(i,1)),8);
% %     data = [data(:,8) data(:,4)];
% %     data(:,2) = movmean(data(:,2),1000);
% % 
% %     plot(data(:,1),data(:,2));
% % 
% %     sd = ((cell2mat(sim_bin_data_comb(i,2)) - 9.5)/2);
% % 
% %     title(sprintf('%g (SD)', sd));
% % 
% %     ylim([0.8 1.7]);
% %     xlim([0 max(data(:,1))]);
% % 
% %     
% % 
% % end




% % % % Vertical
% % % 
% % % figure(100)
% % % for i = 1:length(dist_groups)
% % % std_range = 0.25*ones(size(ydata(:,i)));
% % % subplot(ceil(nbins/6),6,i);
% % % e = errorbar(ydata(:,i),xdata(:,i),errors(:,i),'.');
% % % %e = errorbar(xdata(:,i),ydata(:,i),std_range,std_range,errors(:,i),errors(:,i),'.');
% % % e.MarkerSize = 15;
% % % e.LineWidth = 1.;
% % % e.MarkerFaceColor = 'r';
% % % e.MarkerEdgeColor = 'r';
% % % e.CapSize = 0;
% % % e.AlignVertexCenters = 'on';
% % % ylabel('Local Intensity (Norm)','FontSize',8);
% % % xlabel('Surface Intensity (SD)','FontSize',8);
% % % xlim([0.6 1.6])
% % % title(sprintf('%g - %g \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
% % % sgtitle('Distance From Vessel') 
% % % end



% % % run ANOVA comparing intensity at local intensity for each distance from
% % % vessel bin for the surface intensity bins. 
% % close all;
% % for i = 1:length(dist_groups)
% % 
% % figure(i)
% % [p,tbl,stats(i,:)] = anovan(dist_groups{i}(:,4),dist_groups{i}(:,8),'display','off');
% % c = multcompare(stats(i,:));
% % end
% % 
% % for i = 1:length(dist_groups)
% % xlabel('Local Intensity');
% % ylabel('Surface Int');
% % title(sprintf('%g - %g. \x3bcm From Vessel',round(edges8(i),1), round(edges8(i+1),1)));
% % if i<=15
% % figure(nbins + 1)
% % subplot(3,5,i)
% % boxplot(dist_groups{i}(:,4),dist_groups{i}(:,8), 'symbol', '');
% % q = quantile(dist_groups{i}(:,4),[.25 .5 .75]);
% % ylim([(0.5*q(1)) (1.5*q(3))]);
% % title(sprintf('%g - %g. \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
% % ylabel('Local Int');
% % xlabel('Surface Int');
% % hold on
% % sgtitle('Distance From Vessel') 
% % else 
% % figure(nbins + 2)
% % subplot(3,5,i-15)
% % %c = mult_plot_multcompare(stats(i,:));
% % boxplot(dist_groups{i}(:,4),dist_groups{i}(:,8), 'symbol', '');
% % q = quantile(dist_groups{i}(:,4),[.25 .5 .75]);
% % ylim([(0.5*q(1)) (1.5*q(3))]);
% % title(sprintf('%g - %g. \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
% % ylabel('Local Int');
% % xlabel('Surface Int');
% % hold on
% % sgtitle('Distance From Vessel') 
% % end
% % 
% % end

waitbar(1,f,'Finished!');


% % % % Vertical
% % % 
% % % figure(100)
% % % for i = 1:length(dist_groups)
% % % std_range = 0.25*ones(size(ydata(:,i)));
% % % subplot(ceil(nbins/6),6,i);
% % % e = errorbar(ydata(:,i),xdata(:,i),errors(:,i),'.');
% % % %e = errorbar(xdata(:,i),ydata(:,i),std_range,std_range,errors(:,i),errors(:,i),'.');
% % % e.MarkerSize = 15;
% % % e.LineWidth = 1.;
% % % e.MarkerFaceColor = 'r';
% % % e.MarkerEdgeColor = 'r';
% % % e.CapSize = 0;
% % % e.AlignVertexCenters = 'on';
% % % ylabel('Local Intensity (Norm)','FontSize',8);
% % % xlabel('Surface Intensity (SD)','FontSize',8);
% % % xlim([0.6 1.6])
% % % title(sprintf('%g - %g \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
% % % sgtitle('Distance From Vessel') 
% % % end



% % % run ANOVA comparing intensity at local intensity for each distance from
% % % vessel bin for the surface intensity bins. 
% % close all;
% % for i = 1:length(dist_groups)
% % 
% % figure(i)
% % [p,tbl,stats(i,:)] = anovan(dist_groups{i}(:,4),dist_groups{i}(:,8),'display','off');
% % c = multcompare(stats(i,:));
% % end
% % 
% % for i = 1:length(dist_groups)
% % xlabel('Local Intensity');
% % ylabel('Surface Int');
% % title(sprintf('%g - %g. \x3bcm From Vessel',round(edges8(i),1), round(edges8(i+1),1)));
% % if i<=15
% % figure(nbins + 1)
% % subplot(3,5,i)
% % boxplot(dist_groups{i}(:,4),dist_groups{i}(:,8), 'symbol', '');
% % q = quantile(dist_groups{i}(:,4),[.25 .5 .75]);
% % ylim([(0.5*q(1)) (1.5*q(3))]);
% % title(sprintf('%g - %g. \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
% % ylabel('Local Int');
% % xlabel('Surface Int');
% % hold on
% % sgtitle('Distance From Vessel') 
% % else 
% % figure(nbins + 2)
% % subplot(3,5,i-15)
% % %c = mult_plot_multcompare(stats(i,:));
% % boxplot(dist_groups{i}(:,4),dist_groups{i}(:,8), 'symbol', '');
% % q = quantile(dist_groups{i}(:,4),[.25 .5 .75]);
% % ylim([(0.5*q(1)) (1.5*q(3))]);
% % title(sprintf('%g - %g. \x3bcm',round(edges8(i),1), round(edges8(i+1),1)));
% % ylabel('Local Int');
% % xlabel('Surface Int');
% % hold on
% % sgtitle('Distance From Vessel') 
% % end
% % 
% % end


