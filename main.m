clear;close;clc;close all;
cd('C:\Users\thmorin\Documents\Methods\Li7810ChamberProcessing\');%Directory where your codes are
dataDir='C:\Users\thmorin\Documents\Projects\2021_LaurenOWC_Microsite\Data\'; %Directory where your data are

userName='Tim Morin'; %Update to your name when you process

dat_num='06022022'; %Date vector of when the measurements were taken - User is expected to fill this out accurately

%% Read the file.
metadata=ImportMetaDataSheet([dataDir 'Trip' dat_num '\metadata_' dat_num '.csv']); %It must be in a folder called "Data" which must then have a folder called "Trip06292022" or something like that. 

%% Makes the calculated flux folder if it does not already exist
fluxDir=[dataDir '\Trip' dat_num '\Flux' dat_num '\'];
f=dir(fluxDir);
if isempty(f)
   mkdir(fluxDir);
end
%% Set up data structure for output file
fname=['FluxOutput' dat_num '.csv'];%fid= fopen(fname,'w'); % w= write access
filepath = fluxDir;
file= fullfile(filepath, fname);
fid = fopen(file, 'wt');
fprintf(fid,['Written by ' userName '\n']);
fprintf(fid,[datestr(now()) '\n']);
fprintf(fid,'Date,Site,CO2_Flux,CO2_NRMSE,CH4_Flux,CH4_NRMSE,Bubble_Flux\n');
fprintf(fid,'Date,Site,umol m2 s-1,%%,umol m2 s-1,%%,umol m2 s-1\n');
fclose(fid);

for i=1:length(metadata.Date)
    %% Chamber analysis for i_th chamber
    [Flux_HM_CO2, NRMSE_CO2,Flux_HM_CH4, NRMSE_CH4,Bubble_CH4]=ProcessLI7810Chamber(... 
        metadata.Date(i),metadata.start_time(i),metadata.end_time(i),...
        [dataDir 'Trip' dat_num '\' metadata.data_sheet{i}],metadata.Start_offset(i),metadata.End_offset(i),... 
        metadata.ChamberVolume(i),metadata.ChamberArea(i),1,metadata.Name{i},metadata.Start_Temp(i),metadata.End_Temp(i),...
        metadata.Pressure(i),metadata.Shift_back(i));
	%% Append data to output file
    fid = fopen(fname,'a+'); % a+ = Open or create new file for reading and writing. Append data to the end of the file.
    fprintf(fid,'%s,%s,%f,%f,%f,%f,%f\n',dat_num,metadata.Name{i},Flux_HM_CO2,NRMSE_CO2,Flux_HM_CH4,NRMSE_CH4,Bubble_CH4);
    fclose(fid);
end
