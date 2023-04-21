clc; close all; clear all;
%% make combination plot

% Input data folder name (INCLUDE FOLDER PATH AND NAME):
xdata_sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Stat_data_imp_fit_xdata';
ydata_sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Stat_data_imp_fit_ydata';
errors_sourcefolder = 'C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Stat_data_imp_fit_errors';

% Input shuffled data for comparison (INCLUDE FOLDER PATH AND NAME)
shuff_xdata = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data_imp_fit\Human2131_percentage_shuffle_xdata');
shuff_ydata = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data_imp_fit\Human2131_percentage_shuffle_ydata');
shuff_errors = readmatrix('C:\Users\zh624\Desktop\AD Blood Vessel Tracing\Human2131\Shuffled_data_imp_fit\Human2131_percentage_shuffle_errors');

% read real xdata
f = waitbar(0,'Reading Xdata Data');
fcontent = dir(fullfile(tau_folder, '*.txt')); %fcontent is a column vector of structures

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
fcontent = dir(fullfile(tau_folder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = ydata_sourcefolder;

n = length(fcontent');

for i = 1:n

   filename = fcontent(i).name;
   ydata_single = [];
   ydata_single = readmatrix(fullfile(sourcefolder,filename));
   ydata{i,1} = ydata_single;

   waitbar(i/n,f,sprintf('Reading Xdata: %d %%',floor(i/n*100)));

end

% Read real errors

f = waitbar(0,'Reading errors Data');
fcontent = dir(fullfile(tau_folder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = errors_sourcefolder;

n = length(fcontent');

for i = 1:n

   filename = fcontent(i).name;
   errors_single = [];
   errors_single = readmatrix(fullfile(sourcefolder,filename));
   errors{i,1} = errors_single;

   waitbar(i/n,f,sprintf('Reading errors: %d %%',floor(i/n*100)));

end

