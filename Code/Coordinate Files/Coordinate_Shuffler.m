clc; clear all; close all; 

inputCoordFolder = 'C:\Users\exx\Desktop\Zach - Active Analysis\2131\FINAL_coordinates';


fcontent = dir(fullfile(inputCoordFolder, '*.txt')); %fcontent is a column vector of structures

sourcefolder = inputCoordFolder;

for h = 1:length(fcontent)

filename = fcontent(h).name;
%%

hREAD = h

realCoord = readmatrix(fullfile(sourcefolder,filename));

hShuffle = h

l = length(realCoord(:,6));

idx = randperm(l);

shuffleCoord = realCoord;

shuffleCoord(idx,6) = realCoord(:,6);

hWrite = h

%%


writematrix(shuffleCoord,append('C:\Users\exx\Desktop\Zach - Active Analysis\2131\Shuffled_Coord','\shuffled_',filename));

clear shuffleCoord


end