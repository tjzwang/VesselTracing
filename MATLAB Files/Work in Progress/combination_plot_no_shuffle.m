clc; close all; clear all;
%% make combination plot

% Input data folder name (INCLUDE FOLDER PATH AND NAME):
xdata_sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Stat_data_in_vessel_percent_xdata';
ydata_sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Stat_data_in_vessel_percent_ydata';
errors_sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Stat_data_in_vessel_percent_errors';
%sizes_sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Stat_data_imp_fit_sizes';

% Input shuffled data for comparison (INCLUDE FOLDER PATH AND NAME)
% shuff_xdata = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data_imp_fit\Human2131_percentage_shuffle_xdata');
% shuff_ydata = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data_imp_fit\Human2131_percentage_shuffle_ydata');
% shuff_errors = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data_imp_fit\Human2131_percentage_shuffle_errors');

% Max Distance From Vessel (in microns) (Recommended Value 30):

ves_dist = [ 30 ];

% Resolution in microns? (Recommended Value: 1):

res = [ 1 ];

% Standard Deviation/Percentile or Percent Binning? 
% (1 = SD, 2 = %)

bintype = [ 2 ];

% Single Vessel Plots?:
% (1 = no, 2 = yes)

singleplots = [ 2 ];

% Combined plots?:
% (1 = no, 2 = yes)

combineplots = [ 1 ];

%%

% read real xdata
f = waitbar(0,'Reading Xdata Data');
fcontent = dir(fullfile(xdata_sourcefolder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = xdata_sourcefolder;

n = length(fcontent');

for i = 1:n

   filename = fcontent(i).name;
   xdata_single = [];
   xdata_single = readmatrix(fullfile(sourcefolder,filename));
   xdata{i,1} = xdata_single;

   waitbar(i/n,f,sprintf('Reading Xdata: %d %%',floor(i/n*100)));

end

% Read real ydata

f = waitbar(0,'Reading Ydata Data');
fcontent = dir(fullfile(ydata_sourcefolder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = ydata_sourcefolder;

n = length(fcontent');

for i = 1:n

   filename = fcontent(i).name;
   ydata_single = [];
   ydata_single = readmatrix(fullfile(sourcefolder,filename));
   ydata{i,1} = ydata_single;

   waitbar(i/n,f,sprintf('Reading ydata: %d %%',floor(i/n*100)));

end

% Read real errors

f = waitbar(0,'Reading errors Data');
fcontent = dir(fullfile(errors_sourcefolder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = errors_sourcefolder;

n = length(fcontent');

for i = 1:n

   filename = fcontent(i).name;
   errors_single = [];
   errors_single = readmatrix(fullfile(sourcefolder,filename));
   errors{i,1} = errors_single;

   waitbar(i/n,f,sprintf('Reading errors: %d %%',floor(i/n*100)));

end
%%
% % % Read real sizes
% % 
% % f = waitbar(0,'Reading errors Data');
% % fcontent = dir(fullfile(sizes_sourcefolder, '*.txt')); %fcontent is a column vector of structures
% % 
% % sourcefolder = sizes_sourcefolder;
% % 
% % n = length(fcontent');
% % 
% % for i = 1:n
% % 
% %    filename = fcontent(i).name;
% %    sizes_single = [];
% %    sizes_single = readmatrix(fullfile(sourcefolder,filename));
% %    sizes{i,1} = sizes_single;
% % 
% %    waitbar(i/n,f,sprintf('Reading sizes: %d %%',floor(i/n*100)));
% % 
% % end
%%

bin_num = length(ydata{1}(:,1));

if singleplots == 2

    [r,c] = size(xdata{1});
[r2,c2] = size(xdata);

ydata_ind = ydata{1};


for i = 1:r2

    dist = xdata{i};

for ii = 1:res:ves_dist
    dist(:,ii) = ii;
end

ydata_comb = zeros(r,c);

vessel_figures{i} = figure();
plot_width = 6;
tiledlayout(ceil(r/6),plot_width);

n = 1;

%bin_num = length(sim_bin_data_comb(:,2));
dim = ceil(length(xdata{i}(:,1))/3);
tiledlayout(3,dim);

    for p = 1:r
        nexttile;
        e = errorbar(dist,xdata{i}(p,:),errors{i}(p,:),'.');
        hold on;
        %e3 = errorbar(dist,shuff_xdata(p,:),shuff_errors(p,:),'.');

        if bintype == 2
        title(sprintf('%g - %g %%',(ydata_ind(p,1)-10),ydata_ind(p,1)));
        sgtitle(sprintf('Surface Intensity %% (Vessel %g)',i)) 
        else
        title(sprintf('%g - %g (SD)',(ydata_ind(p,1)-0.25),(ydata_ind(p,1)+0.25)));
        sgtitle('Surface Intensity SD') 
        end

        ylabel('Local Intensity (Norm)','FontSize',8);
        xlabel('Distance From Vessel (\mum)','FontSize',8);
        pos = get(gcf, 'Position');
        set(gcf, 'Position',pos+[-35 -20 35 20]);

ylim([(min(xdata{i},[],'all')-0.1) (max(xdata{i},[],'all')+0.1)]);
xlim([-2 (ves_dist+1)]);


if i == ceil(bin_num/3)

%lg = legend('Shuffled','Real','location','northeast','orientation','horizontal');
%lg.Location = 'northeastoutside';

else
end

end

    end
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
%%

if combineplots == 2

% Combine individual vessel data
[r,c] = size(xdata{1});
[r2,c2] = size(xdata);
xdata_comb = zeros(r,c);
for i = 1:numel(xdata{1})
    for p = 1:r2
        if i > numel(xdata{p}(i))
            break
        end
        xdata_temp(p) = xdata{p}(i);
    end
xdata_comb(i) = mean(xdata_temp);
xdata_temp = [];
end

[r,c] = size(ydata{1});
[r2,c2] = size(ydata);
ydata_comb = zeros(r,c);
for i = 1:numel(ydata{1})
    for p = 1:r2
        if i > numel(ydata{p}(i))
            break
        end
        ydata_temp(p) = ydata{p}(i);
    end
ydata_comb(i) = mean(ydata_temp);
ydata_temp = [];
end

[r,c] = size(errors{1});
[r2,c2] = size(errors);
errors_comb = zeros(r,c);
for i = 1:numel(errors{1})
    for p = 1:r2
        if i > numel(errors{p}(i))
            break
        end
        errors_temp(p) = errors{p}(i);
    end
errors_comb(i) = sqrt((1/r2)*(sum(errors_temp.^2)));
errors_temp = [];
end

xdata = xdata_comb;
ydata = ydata_comb;
errors = errors_comb;

dist = xdata;

for i = 1:res:ves_dist
    dist(:,i) = i;
end


%%

[r,c] = size(xdata);
% [r2,c2] = size(shuff_ydata);
bin_num = r;
figure(100)
plot_width = 6;
tiledlayout(ceil(r/6),plot_width);

n = 1;

%bin_num = length(sim_bin_data_comb(:,2));
tiledlayout(3,ceil(length(xdata(:,1))/3));

for i = 1:r
nexttile
std_range = 0.25*ones(size(ydata(:,i)));
%subplot(ceil(nbins/6),6,i);
% % e = errorbar(xdata(i,:),edges8(1:(length(edges8)-1)),errors(i,:),'horizontal','.');
e = errorbar(dist,xdata(i,:),errors(i,:),'.');

% if n <= r2
% % idx = find(shuff_ydata(:,1) == ydata(i,1));
% if length(idx) == 1
%     hold on
% % % e2 = errorbar(shuff_xdata(n,:),edges8(1:(length(edges8)-1)),shuff_errors(n,:),'horizontal','.');
% e3 = errorbar(dist,shuff_xdata(idx,:),shuff_errors(idx,:),'.');

% e2.MarkerSize = 15;
% e2.LineWidth = 1;
% e2.Color = 'm';
% e2.MarkerFaceColor = 'g';
% e2.MarkerEdgeColor = 'g';
% e2.CapSize = 0;
% e2.AlignVertexCenters = 'on';

hold on


%e = errorbar(xdata(:,i),ydata(:,i),std_range,std_range,errors(:,i),errors(:,i),'.');
% e.MarkerSize = 15;
% e.LineWidth = 1;
% e.MarkerFaceColor = 'r';
% e.MarkerEdgeColor = 'r';
% e.CapSize = 0;
% e.AlignVertexCenters = 'on';

% % xlabel('Local Intensity (Norm)','FontSize',8);
% % ylabel('Distance From Vessel','FontSize',8);
% % set(gca, 'YDir','reverse');

ylabel('Local Intensity (Norm)','FontSize',8);
xlabel('Distance From Vessel (\mum)','FontSize',8);

ylim([(min(xdata,[],'all')-0.1) (max(xdata,[],'all')+0.1)]);
xlim([-2 (ves_dist+1)]);


if i == ceil(bin_num/3)

lg = legend('Shuffled','Real','location','northeast','orientation','horizontal');
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
sgtitle('Surface Intensity % (All Vessels)') 
else
title(sprintf('%g - %g (SD)',(ydata(i,1)-0.25),(ydata(i,1)+0.25)));
sgtitle('Surface Intensity SD') 
end
%nexttile
end

end
%%

% if singleplots == 2
% for i = 1:length(vessel_figures)
%     exportgraphics(vessel_figures{i},append('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Figures\',sprintf('vessel_%g_ind_percent_plot.jpg',i)));
% end
% end