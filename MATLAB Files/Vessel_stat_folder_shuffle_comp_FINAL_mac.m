clc; clear all; close all;

%% Import Data

% .txt file (Recommended for speed):

% Input data folder name (INCLUDE FOLDER PATH AND NAME):
sourcefolder = '/Users/zacharyhoglund/Desktop/Intensity_distance_data_abs';

% Input shuffled data for comparison (INCLUDE FILE PATH AND NAME)
shuff_xdata = readmatrix('/Users/zacharyhoglund/Desktop/Shuffled_data_abs/Human2470_percentage_shuffle_xdata');
shuff_ydata = readmatrix('/Users/zacharyhoglund/Desktop/Shuffled_data_abs/Human2470_percentage_shuffle_ydata');
shuff_errors = readmatrix('/Users/zacharyhoglund/Desktop/Shuffled_data_abs/Human2470_percentage_shuffle_errors');

% Set Moving Average Size (Recommended Value: 50):
Int_moving_ave_size = [  50  ];

% Max Distance From Vessel (in microns) (Recommended Value 30):

ves_dist = [ 30 ];

% Resolution in microns? (Recommended Value: 1):

res = [ 1 ];

% Frame width (# Points):

frame_w = [ 20 ];

% Distance to define surface data (in microns)? (Recommended Value: 3):

surf_d = [ 3 ];

% # Bins along vessel
nbins_length = 5;

% Cortex? (0 = no, 1 = yes):

cortex_yesno = [ 0 ];

% Standard Deviation/Percentile or Percent Binning? 
% (1 = SD, 2 = %)

bintype = [ 2 ];

% Error type? (0 = standard error, 1 = anova group interval)

errortype = [ 0 ];

% Number of percent bins:

numbins = [ 10 ];

% All data or single vessel binning? (0 = all, 1 = single)

vesnum = [ 1 ];

% Data Output Name (True, unshuffled data only) (Include Name + Path):
% Vessel number and data type will be appended, so do not include in names.
% NEED FOLDER WITH SAME NAME FOR EACH DATA TYPE: XDATA, YDATA, AND ERRORS
% ex. you need folders stat_data_xdata, stat_data_ydata, and
% stat_data_errors if stat_data is the folder name.

% Folder:
outputfoldername = '/Users/zacharyhoglund/Desktop/Stat_data';
% File:
outputname = 'Human2470_percentage';

% % Set number of bins by increment of intensity
% nbins = round((max(Tau_data(:,4)) - min(Tau_data(:,4)))/50);

%% Load Data

% # Bins along vessel
nbins_length = 5;

f = waitbar(0,'Calculating Statistics: Reading Data 0 %');

fcontent = dir(fullfile(sourcefolder, '*.txt')); %fcontent is a column vector of structures

l1 = 1;

% Open files and load data from folder
% If multiple files exist, their data will be stacked

Tau_data = [];

n = length(fcontent');

for b = 1:n
    sizes = [];
    if vesnum == 1
        n1 = b;
        n2 = b;
        Tau_data = [];
        bin_edges = [];
        surface_data = [];
        bin_numbers = [];
        clc
        % close all
    else
        n1 = 1;
        n2 = n;
    end
for i = n1:n2

%i = 7;

   filename = fcontent(i).name;
   Tau_data_single = [];
   Tau_data_single = readmatrix(fullfile(sourcefolder,filename));
   if length(Tau_data) >= 1 && vesnum == 0
   Tau_data_single(:,5) = Tau_data_single(:,5)+max(Tau_data(:,5));
   end
   [r, c] = size(Tau_data_single);
   l2 = l1 + r - 1;
   Tau_data(l1:l2,:) = Tau_data_single;
   l1 = r + l1;
   if vesnum ==1
       l1 = 1;
   end

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

 if cortex_yesno == 0
       Tau_data(:,7) = zeros(length(Tau_data(:,1)),1);
 end

Tau_data = Tau_data(:,1:7);

% calculate number of bins based on specified resolution
nbins = ves_dist/res;

% % Set number of bins by increment of intensity
% nbins = round((max(Tau_data(:,4)) - min(Tau_data(:,4)))/50);

group_Tau_data = Tau_data;

% sort the ungrouped data by intensity
group_Tau_data = sortrows(group_Tau_data,4);

% get bin edges for grouping the tau data by intensity
[N,edges7,bin] = histcounts(group_Tau_data(:,4),nbins);

cortex_data = Tau_data(:,7);

group_Tau_data = [group_Tau_data(:,1:6) bin zeros(length(bin),1) zeros(length(bin),1) zeros(length(bin),1) cortex_data];

%% ANOVA on distance from vessel vs. intensity


waitbar(0.2,f,'Calculating Statistics: Running Initial ANOVA');

[p,tbl,stats] = anovan(group_Tau_data(:,6),group_Tau_data(:,7),'display','off');

c = multcompare(stats);

%% Group Data by distance from vessel

waitbar(0.3,f,'Calculating Statistics: Grouping Data');

dist_Tau_data = sortrows(group_Tau_data,6);

% get bin edges for grouping the tau data by intensity
[N,edges8,bin] = histcounts(dist_Tau_data(:,6),nbins);

group_Tau_data = [dist_Tau_data(:,1:7) bin zeros(length(bin),1) zeros(length(bin),1) dist_Tau_data(:,11)];

%% ANOVA on Intensity vs. distance from vessel

[p,tbl,stats] = anovan(group_Tau_data(:,4),group_Tau_data(:,8),'display','off');
figure()
c = multcompare(stats);

%% Group Data by distance along vessel

group_Tau_data = sortrows(group_Tau_data,5);

% get bin edges for grouping the tau data by intensity
[N,edges9,bin] = histcounts(group_Tau_data(:,5),nbins_length);

group_Tau_data = [group_Tau_data(:,1:8) bin zeros(length(bin),1) dist_Tau_data(:,11)];


%% ANOVA on Intensity vs. distance along vessel

[p,tbl,stats] = anovan(group_Tau_data(:,4),{group_Tau_data(:,8),group_Tau_data(:,9)},'display','off');
figure()
c = multcompare(stats,"Dimension", [1 2]);

%% Group data by intensity standard deviation

waitbar(0.4,f,'Calculating Statistics: Finding Surface Data & Intensity Standard Deviation');

% filter data to only include tau within surf_dist microns of vessel.
surface_idx = find(group_Tau_data(:,6) < surf_d & group_Tau_data(:,6) >= 0);
surface_data = group_Tau_data(surface_idx,:);
surface_data_back = group_Tau_data(surface_idx,:);


% Calculate standard deviation and median
tauStd = std(group_Tau_data(:,4));
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
Tau_data = [group_Tau_data(:,1:6)  group_Tau_data(:,8) group_Tau_data(:,11)];

group_Tau_data = [];

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

if vesnum == 0 
parfor_progress(lb); % Set the total number of iterations
end
idx4 = {zeros(lb,1)};
parfor i = 1:lb
   
    Tau_data2 = Tau_data;
    bin_edges2 = bin_edges;
    idx4{i} = {find(Tau_data2(:,5) >= bin_edges2(i,2) & Tau_data2(:,5) <= bin_edges2(i,3))};
    if vesnum == 0
        parfor_progress;
    end
end
if vesnum == 0
parfor_progress(0);
end
Tau_data2 = [];
%% 

for i = 1:lb
    sim_bin_data(i,1:2) = [{Tau_data(cell2mat(idx4{i}),:)} bin_edges(i,1)];
end

Tau_data = [];

% get bin numbers from cell array
int_bin = cell2mat(sim_bin_data(:,2));


% combine cells of array with the same intensity bin number
for i = 1:max(surface_data(:,10))
    idx5 = find(int_bin == i);
    sim_bin_data_comb(i,1:2) = [{unique(cat(1,sim_bin_data{idx5}),'rows')} i];
    idx5 = [];
end

%%

f = waitbar(0.6,'Calculating Statistics: Running Final ANOVA');

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
        cell(z:(z+l-1),:) = [anova_bin_data{i}(idx6,1:7) int_group(i)*ones(l,1)];
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

if bintype == 1
bin_numbers = cell2mat(sim_bin_data_comb(:,2));
bin_numbers = [(bin_numbers - 9)/2];
else
bin_numbers = cell2mat(sim_bin_data_comb(:,2));
bin_numbers = bin_numbers.*10;
end

for i = 1:length(dist_groups)
n = i;
[p,tbl,stats(i,:)] = anovan(dist_groups{i}(:,4),dist_groups{i}(:,8),'display','off');
set(0,'DefaultFigureVisible','off')
[c,m] = multcompare(stats(i,:));
set(0,'DefaultFigureVisible','on')

sizes(i) = length(dist_groups{i});

% % if bintype == 1
% % bin_numbers = [(bin_numbers - 9)/2];
% % else
% % bin_numbers = bin_numbers.*10;
% % end

m = [m bin_numbers];
[m_length,c] = size(m);
xdata(:,i) = [m(:,1)];
ydata(:,i) = [m(:,3)];

% % bin_numbers = [];

if errortype == 0
errors(1:m_length,i) = [m(:,2)];
else
errors(1:m_length,i) = std(dist_groups{i}(:,4))/sqrt(length(dist_groups{i}(:,4)));
end

end



% % % % Load shuffled data for comparison
% % % shuff_xdata = readmatrix('/Users/zacharyhoglund/Desktop/Arivis/Human2131/MATLAB Test Shuffle Data/Human2131_test_shuffle_xdata');
% % % shuff_ydata = readmatrix('/Users/zacharyhoglund/Desktop/Arivis/Human2131/MATLAB Test Shuffle Data/Human2131_test_shuffle_ydata');
% % % shuff_errors = readmatrix('/Users/zacharyhoglund/Desktop/Arivis/Human2131/MATLAB Test Shuffle Data/Human2131_test_shuffle_errors');

if bintype == 1
    bin_num = length(sim_bin_data_comb(:,2));
end

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

f=waitbar(0.8,sprintf('Calculating Statistics: Exporting Data %d %%',80));

outputname2 = append(outputname,sprintf('%g',b));


writematrix(xdata,append(outputfoldername,'_xdata','/',outputname2,'_stat_xdata'));
writematrix(ydata,append(outputfoldername,'_ydata','/',outputname2,'_stat_ydata'));
writematrix(errors,append(outputfoldername,'_errors','/',outputname2,'_stat_errors'));
%writematrix(sizes,append(outputfoldername,'_sizes','\',outputname,'_stat_sizes'));

[r,c] = size(xdata);
[r2,c2] = size(shuff_ydata);

figure(100)

if bintype == 1
plot_width = 6;
tiledlayout(ceil(nbins/6),plot_width);
end

edge_length = length(edges8);
n = 1;

if bintype ==2
bin_num = length(sim_bin_data_comb(:,2));
tiledlayout(3,ceil(length(xdata(:,1))/3));
end

for i = 1:r
nexttile
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
if bintype == 2
title(sprintf('%g - %g %%',(ydata(i,1)-10),ydata(i,1)));
sgtitle('Surface Intensity %') 
else
title(sprintf('%g - %g (SD)',(ydata(i,1)-0.25),(ydata(i,1)+0.25)));
sgtitle('Surface Intensity SD') 
end
%nexttile
end

%%

% figure()

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
if vesnum == 0
    break
end
close all;
end

