clc; clear all; close all;

%% Import Data

% .txt file (Recommended for speed):

% Input data folder name (INCLUDE FOLDER PATH AND NAME):

% % input_file = '/Users/zacharyhoglund/Downloads/FINAL DATA OUTPUT/Tangle_MEGASHEET/tangle_10micron_bins_to_100microns_3MicronSurface 1.xlsx';

input_file = '/Users/zacharyhoglund/Downloads/FINAL DATA OUTPUT/Neuron_MEGASHEET/neurons_10micron_bins_to_100microns_3MicronSurface.xlsx';

% Object Type:

obj_type = 'Neuron';


% Max Distance From Vessel (in microns) (Recommended Value 30):

ves_dist = [ 100 ];

% Resolution in microns? (Recommended Value: 1):

res = [ 1 ];

% Distance to define surface data (in microns)? (Recommended Value: 3):

surf_d = [ 3 ];

% Cortex? (0 = no, 1 = yes):

cortex_yesno = [ 0 ];

% Standard Deviation/Percentile or Percent Binning? 
% (1 = SD, 2 = %)

bintype = [ 2 ];

% Error type? (0 = standard error, 1 = anova group interval)

errortype = [ 1 ];

% Percentile bin width:

perc_binW = [ 10 ];

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

data = readmatrix(input_file);

data2 = data;

% % samples = [2392 2302 2290 2131 2267 24632 23642 2477 2470 24652 2417 2399];

samples = [ 2302 2290 2131 2267 23642 2399];

for q = 1:length(samples)

    idx = find(data(:,1) == samples(q));

    if length(idx) == 0
        continue
    end

    ADdata = data(idx,:);

ADdata2 = ADdata;

% % % for i = 1:length(ADdata)
% % % 
% % %     ADdata2(i,:) = [ADdata(i,1:6) mean(ADdata(i,17:(17+surf_d))) ADdata(i,17:(17+ves_dist)) ];
% % % 
% % % end

%% Group data by surface intensity standard deviation

%waitbar(0.4,f,'Calculating Statistics: Finding Surface Data & Intensity Standard Deviation');


% Calculate standard deviation and median
ADstd = std(ADdata2(:,7));
ADmed = median(ADdata2(:,7));


ADd = ADdata(:,7);
ADs = ADstd;
ADm = ADmed;

if bintype == 1

% Find indexes of vaules for each bin based on 1/2 standard deviation
ADidx = padcat(find( ADd < ADm-3.75*ADs),find( ADd < ADm-3.25*ADs & ADd>= ADm-3.75*ADs),find( ADd < ADm-2.75*ADs & ADd>= ADm-3.25*ADs),find( ADd < ADm-2.25*ADs & ADd>= ADm-2.75*ADs),find( ADd < ADm-1.75*ADs & ADd>= ADm-2.25*ADs),find( ADd < ADm-1.25*ADs & ADd >= ADm-1.75*ADs),find( ADd < ADm-0.75*ADs & ADd>= ADm-1.25*ADs),find( ADd < ADm-0.25*ADs & ADd>= ADm-0.75*ADs),find( ADd < ADm+0.25*ADs & ADd>= ADm-0.25*ADs),find( ADd < ADm+0.75*ADs & ADd>= ADm+0.25*ADs),find( ADd < ADm+1.25*ADs & ADd>= ADm+0.75*ADs),find( ADd < ADm+1.75*ADs & ADd>= ADm+1.25*ADs),find( ADd < ADm+2.25*ADs & ADd>= ADm+1.75*ADs),find( ADd < ADm+2.75*ADs & ADd >= ADm+2.25*ADs),find( ADd < ADm+3.25*ADs & ADd>= ADm+2.75*ADs),find( ADd < ADm+3.75*ADs & ADd>= ADm+3.25*ADs),find( ADd < ADm+4.25*ADs & ADd>= ADm+3.75*ADs),find(ADd>= ADm+4.25*ADs));

% Put bin number for each row of tau data in column 10
[r,c] = size(ADidx);
for i = 1:c
    ADidx2 = ADidx(:,i);
    ADidx2(isnan(ADidx2))=[];
    if ~isempty(ADidx2)
    bin = ones(length(ADidx2),1)*i;
    ADdata3(ADidx2,:) = [bin ADdata2(ADidx2,:)];
    end
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


else if bintype == 2

edges = 0:(100/perc_binW):100;

   % % % for i = 1:length(ADdata2)
   % % % 
   % % %    ADpercentiles(i,:) = comp_percentile(ADdata2(:,7),ADdata2(i,7));
   % % % 
   % % % end

   ADpercentiles = invprctile(ADdata2(:,7),ADdata2(:,7));
   
   ADbin = discretize(ADpercentiles,edges);

   ADdata3 = [ADbin ADdata2];

   % % % for i = 1:length(CONTdata2)
   % % %    CONTpercentiles(i,:) = comp_percentile(CONTdata2(:,7),CONTdata2(i,7));
   % % % 
   % % % end


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

%% Bin data into percentile/standard dev groups

for i = 1:max(ADbin)

   
    idx = find(ADdata3(:,1) == i);

    ADplotdatatemp = ADdata3(idx,9:end);

    [ADRcoord,ADCcoord] = ind2sub(size(ADplotdatatemp),find(ADplotdatatemp > -10000));
    
    if length(ADRcoord(:,1))==0
        continue
    end

    for ii = 1:length(ADRcoord)
        ADplotdata2(ii,:) = [ADCcoord(ii) ADplotdatatemp(ADRcoord(ii),ADCcoord(ii))];
    end
    
    ADplotdata{i,1} = ADplotdata2;
    ADplotdata{i,2} = i;

end


ADplotdata(all(cellfun('isempty',ADplotdata),2),:) = [];


%%

% f = waitbar(0.6,'Calculating Statistics: Running Final ANOVA');

% run ANOVA for each bin of cell array comparing intensity with distance
% from vessel
for i = 1:length(ADplotdata(:,1))

[ADp,ADtbl,ADstats(i,:)] = anovan(ADplotdata{i}(:,2),ADplotdata{i}(:,1),'display','off');

ADc = multcompare(ADstats(i,:),'display','off');


end

if bintype == 1

ADbin_numbers = cell2mat(ADplotdata(:,2));
ADbin_numbers = [(ADbin_numbers - 9)/2];

else

ADbin_numbers = cell2mat(ADplotdata(:,2));
ADbin_numbers = ADbin_numbers.*10;

ADbin_numbers = 0:max(ADplotdata{1,1}(:,1))-1;
ADbin_numbers = ADbin_numbers';

end
%%

for i = 1:length(ADplotdata(:,1))

n = i;
[ADp,ADtbl,ADstats(i,:)] = anovan(ADplotdata{i,1}(:,2),ADplotdata{i,1}(:,1),'display','off');
%set(0,'DefaultFigureVisible','off')
[ADcMULT,ADm] = multcompare(ADstats(i,:),'display','off');

% % if bintype == 1
% % bin_numbers = [(bin_numbers - 9)/2];
% % else
% % bin_numbers = bin_numbers.*10;
% % end

ADm = [ADm ADbin_numbers];
[ADm_length, ADcMULT] = size(ADm);
ADxdata(:,i) = [ADm(:,1)];
ADydata(:,i) = [ADm(:,3)];

% % bin_numbers = [];

ADerrors(1:ADm_length,i) = [ADm(:,2)];


end

%%
figure()

ADydata = ADydata.*10 + 5;

if bintype == 1
plot_width = 5;
tiledlayout(ceil(length(ADplotdata(:,1))/plot_width),plot_width);
end

% edge_length = length(edges8);
n = 1;

if bintype ==2
bin_num = length(ADplotdata(:,2));
tiledlayout(2,5);
end

[ADr,ADc] = size(ADxdata);

for i = 1:ADc
nexttile
std_range = 0.25*ones(size(ADydata(:,i)));
%subplot(ceil(nbins/6),6,i);
% % e = errorbar(xdata(i,:),edges8(1:(length(edges8)-1)),errors(i,:),'horizontal','.');
e = errorbar(ADydata(:,i),ADxdata(:,i),ADerrors(:,i),'.');


% % e2.MarkerSize = 15;
% % e2.LineWidth = 1;
% % e2.Color = 'm';
% % e2.MarkerFaceColor = 'g';
% % e2.MarkerEdgeColor = 'g';
% % e2.CapSize = 0;
% % e2.AlignVertexCenters = 'on';

hold on


%e = errorbar(xdata(:,i),ydata(:,i),std_range,std_range,errors(:,i),errors(:,i),'.');
e.MarkerSize = 10;
e.LineWidth = 1;
e.MarkerFaceColor = '#EDB120';
e.MarkerEdgeColor = '#EDB120';
e.Color = '#EDB120';
e.CapSize = 10;
e.AlignVertexCenters = 'on';

% % xlabel('Local Intensity (Norm)','FontSize',8);
% % ylabel('Distance From Vessel','FontSize',8);
% % set(gca, 'YDir','reverse');

ylabel(append(sprintf('%s Density',obj_type),' (\mum^{-3})'),'FontSize',8);
xlabel('Distance From Vessel (\mum)','FontSize',8);

%ylim([(min(ADxdata,[],'all')-abs(min(ADxdata,[],'all'))*1) (max(ADxdata,[],'all')+abs(max(ADxdata,[],'all'))*0.4)]);

ylim([(min(ADxdata,[],'all')-abs(max(ADxdata,[],'all'))*0.2) (max(ADxdata,[],'all')+abs(max(ADxdata,[],'all'))*0.4)]);

xlim([-10 (ves_dist+10)]);


if i == ceil(bin_num/3)

%lg = legend('Real','Shuffled','location','northeast','orientation','horizontal');
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

ay=gca; 
ay.YAxis.Exponent = -6;


if bintype == 2
title(sprintf('%g - %g %%',i*perc_binW-perc_binW,i*perc_binW));
sgtitle(sprintf('Human %g: Surface Intensity %%', samples(q))); 
else
title(sprintf('%g - %g (SD)',(ADydata(i,1)-0.25),(ADydata(i,1)+0.25)));
sgtitle('Surface Intensity SD') 
end
%nexttile



    if q == 1 && i == 1

    exportdata(1:length(ADxdata(:,1)),:) = [ samples(q).*ones(length(ADxdata(:,1)),1) ADxdata(:,i) ADydata(:,i) (i*perc_binW).*ones(length(ADxdata(:,1)),1) ];
    
    else

    exportdata(end+1:end+length(ADxdata(:,1)),:) = [ samples(q).*ones(length(ADxdata(:,1)),1) ADxdata(:,i) ADydata(:,i) (i*perc_binW).*ones(length(ADxdata(:,1)),1) ];

    end



end


clear ADbin ADbin_numbers ADc ADCcoord ADcMULT ADd ADdata ADdata2 ADdata3 ADerrors ADm ADm_length ADmed ADp ADpercentiles ADplotdata ADplotdata2 ADplotdatatemp ADr ADRcoord ADs ADstats ADstd ADtbl ADydata ADxdata 

end

for i = 1:100/perc_binW

    idx = find(exportdata(:,4) == i*perc_binW);

    perc_bin_data = exportdata(idx,:);

    if i == 1

    exportdata2(1,:) = [ 0 0 unique(exportdata(:,3))'];

    end

     endrow = length(exportdata2(:,1));

    for q = 1:length(samples)

    idx = find(perc_bin_data(:,1) == samples(q));

    temp_data = perc_bin_data(idx,:);

    temp_data = sortrows(temp_data,3);

    exportdata2(endrow+q,:) = [ i*perc_binW samples(q) temp_data(:,2)' ];

    end

   exportdata2(end+1,:) = [0];

end

%% Export all data

writematrix(exportdata, append('/Users/zacharyhoglund/Downloads/AD_',obj_type,'_data_plot_coordinates.xlsx'));

writematrix(exportdata2, append('/Users/zacharyhoglund/Downloads/AD_',obj_type,'_data_plot_coordinates_binned.xlsx'));

