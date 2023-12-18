clc; clear all; close all;

%% Import Data

% .txt file (Recommended for speed):

% Input data folder name (INCLUDE FOLDER PATH AND NAME):
input_file = '/Users/zacharyhoglund/Downloads/TEST_Piecewise_tau_intensity_binned_data_MEGASHEET_with_headers.xlsx';

% Max Distance From Vessel (in microns) (Recommended Value 30):

ves_dist = [ 100 ];

% Resolution in microns? (Recommended Value: 1):

res = [ 1 ];

% Distance to define surface data (in microns)? (Recommended Value: 3):

surf_d = [ 3 ];

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


samples = [ 2302 2290 2131 2267  23642 2399];

for q = 1:length(samples)

% for q = 1:1

    idx = find(data(:,1) == samples(q));

    if length(idx) == 0
        continue
    end

    ADdata = data(idx,:);


for i = 1:length(ADdata)

    ADdata2(i,:) = [ADdata(i,1:6) mean(ADdata(i,17:(17+surf_d))) ADdata(i,17:(17+ves_dist)) ];

end


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

    ADplotdatatemp = ADdata3(idx,8:end);

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

if bintype == 1
plot_width = 5;
tiledlayout(ceil(length(ADplotdata(:,1))/plot_width),plot_width);
end

% edge_length = length(edges8);
n = 1;

if bintype ==2
bin_num = length(ADplotdata(:,2));
plot_width = 5;
tiledlayout(2,plot_width);
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

e.MarkerSize = 5;
e.LineWidth = 1;
e.CapSize = 0.1;
e.AlignVertexCenters = 'on';

if samples(q) == 2131 | samples(q) == 2302 | samples(q) == 2290 | samples(q) == 2267 | samples(q) == 23642 | samples(q) == 2399
e.MarkerFaceColor = '#EDB120';
e.MarkerEdgeColor = '#EDB120';
e.Color = '#EDB120';
else
   e.MarkerFaceColor = "#0072BD";
e.MarkerEdgeColor = "#0072BD"; 
e.Color = "#0072BD";
end

e.AlignVertexCenters = 'on';

% % xlabel('Local Intensity (Norm)','FontSize',8);
% % ylabel('Distance From Vessel','FontSize',8);
% % set(gca, 'YDir','reverse');

ylabel('Local Intensity (Norm)','FontSize',8);
xlabel('Distance From Vessel (\mum)','FontSize',8);

% ylim([(min(ADxdata,[],'all')-abs(min(ADxdata,[],'all'))*0.1) (max(ADxdata,[],'all')+abs(max(ADxdata,[],'all'))*0.1)]);
xlim([-2 (ves_dist+1)]);

ylim([0 1500]);


if i == ceil(bin_num/plot_width)

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
title(sprintf('%g - %g %%',i*10-10,i*10));
sgtitle(sprintf('SAMPLE %g Surface Intensity %', samples(q))) 
else
title(sprintf('%g - %g (SD)',(ADydata(i,1)-0.25),(ADydata(i,1)+0.25)));
sgtitle(sprintf('SAMPLE %g Surface Intensity SD', samples(q))); 
end
%nexttile


 if q == 1 && i == 1 
        
        
        if isempty(ADxdata) == 0

    exportdata(1:length(ADxdata(:,1)),:) = [ samples(q).*ones(length(ADxdata(:,1)),1) ADxdata(:,i) ADydata(:,i) (i*perc_binW).*ones(length(ADxdata(:,1)),1) ];

        else

    exportdata(1:length(CONTxdata(:,1)),:) = [ samples(q).*ones(length(CONTxdata(:,1)),1) CONTxdata(:,i) CONTydata(:,i) (i*perc_binW).*ones(length(CONTxdata(:,1)),1) ];
        
        end
    
    else

        if isempty(ADxdata) == 0

    exportdata(end+1:end+length(ADxdata(:,1)),:) = [ samples(q).*ones(length(ADxdata(:,1)),1) ADxdata(:,i) ADydata(:,i) (i*perc_binW).*ones(length(ADxdata(:,1)),1) ];

        else

    exportdata(end+1:end+length(CONTxdata(:,1)),:) = [ samples(q).*ones(length(CONTxdata(:,1)),1) CONTxdata(:,i) CONTydata(:,i) (i*perc_binW).*ones(length(CONTxdata(:,1)),1) ];

        end

    end





end
%%
%saveas(gcf,append('/Users/zacharyhoglund/Downloads/plot_',sprintf('%g',q),'_replot_axis1500.jpg'));

%%

   

    

clear ADbin ADbin_numbers ADc ADCcoord ADcMULT ADd ADdata ADdata2 ADdata3 ADerrors ADm ADm_length ADmed ADp ADpercentiles ADplotdata ADplotdata2 ADplotdatatemp ADr ADRcoord ADs ADstats ADstd ADtbl ADydata ADxdata 
clear CONTbin CONTbin_numbers CONTc CONTCcoord CONTcMULT CONTd CONTdata CONTdata2 CONTdata3 CONTerrors CONTm CONTm_length CONTmed CONTp CONTpercentiles CONTplotdata CONTplotdata2 CONTplotdatatemp CONTr CONTRcoord CONTs CONTstats CONTstd CONTtbl CONTydata CONTxdata 



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

     if isempty(temp_data) == 1
         continue
     end

    temp_data = sortrows(temp_data,3);

    exportdata2(endrow+q,:) = [ i*perc_binW samples(q) temp_data(:,2)' ];

    end

   exportdata2(end+1,:) = [0];

end

for i = 1:length(exportdata2(2:end,1))

% % exportdata2(i,3:end) = exportdata2(i,3:end)./mean(exportdata2(i,101:104));

exportdata2(i,3:end) = exportdata2(i,3:end)./exportdata2(i,104);

exportdata2(i,3:end) = (exportdata2(i,3:end) - 1).*100;

end



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

%% Export all data

writematrix(exportdata, append('/Users/zacharyhoglund/Downloads/AD_','intensity', '_adj3_data_plot_coordinates.xlsx'));

writematrix(exportdata2, append('/Users/zacharyhoglund/Downloads/AD_','intensity','_adj3_data_plot_coordinates_binned.xlsx'));



