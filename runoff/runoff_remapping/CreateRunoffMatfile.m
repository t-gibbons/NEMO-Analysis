% Get Runoff data from csv file ready to put into ANHA4 grid
% this script creates the mat file used in remapping the runoff (RemapHBCRunoff.m).

close all
clear all
clc
HYPEfilename='monthlyQ-20170410_MR3';
path='/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HYPEcal_HBC/';

% Load in runoff data file (numbers at end say ignore first two columns and first two rows)
csv_discharge=csvread([path HYPEfilename '.csv'],1,2);

% Get haroIDs
csv_haroID=importdata([path HYPEfilename '.csv']);
csv_haroID=csv_haroID.textdata(1,3:end); % HaroIDs (in cell format)

% Days in Month
csv_DaysInMonth=importdata([path HYPEfilename '.csv']);
csv_DaysInMonth=csv_DaysInMonth.data(:,1);

% load file with HaroIDs and lat lon info
load /mnt/storage0/xhu/PROGRAM/RUNOFF/runoff_HBC/HBSdischarge_geoinfo.csv
HBCdischarge_geoinfo=HBSdischarge_geoinfo;
clear HBSdischarge_geoinfo

% Create matrix for discharge and geo locations
szdischarge=size(csv_discharge);
HYPEdischarge(:,:,1)=csv_discharge;
HYPEdischarge(:,:,2)=zeros(szdischarge);
HYPEdischarge(:,:,3)=zeros(szdischarge);

% Calculate HYPE annual means (used in RemapHBCRunoff.m)
months=1:12;
for ny=1:(szdischarge(1)/12)
    temp=mean(HYPEdischarge(months,:,1));
    temp=temp*3600*24*sum(csv_DaysInMonth(months,1)); % m3/yr
    HYPEdischarge_AnnualMean(ny,:)=temp;
    months=months+12;
end
HYPEdischarge_AnnualMean=mean(HYPEdischarge_AnnualMean);


for i=1:length(HBCdischarge_geoinfo)
    
    temphid_char=csv_haroID{1,i}; % get cell haroid to string
    temphid_num=str2double(temphid_char(3:end-1)); % convert str haroid to number
    ind = find(HBCdischarge_geoinfo(:,1)==temphid_num);
    
    HBCdischarge_geoinfo(ind,5)= HYPEdischarge_AnnualMean(1,i);
    
    lat=HBCdischarge_geoinfo(ind,3); % latitude
    lon=HBCdischarge_geoinfo(ind,4); % longitude
    HYPEdischarge(:,i,2)=lat; % put latitude in each row in the column
    HYPEdischarge(:,i,3)=lon; % put longitude in each row in the column
end

save([path HYPEfilename '.mat'],'HYPEdischarge') % Units of m3/s


HBCdischarge_geoinfo(:,5)= HBCdischarge_geoinfo(:,5)/1000000000; % change from m3 to km3
% save HBSdischarge_geoinfo which has the haroid, km2 drainage area,
% lat, lon, discharge km3/year
save([path HYPEfilename '_geoinfo.mat'],'HBCdischarge_geoinfo','-ascii')

disp('Check:')
disp('Max river discharge:')
max(HBCdischarge_geoinfo(:,5))
disp('Min river discharge:')
min(HBCdischarge_geoinfo(:,5))
disp('HBC annual river discharge:')
sum(HBCdischarge_geoinfo(:,5))

DischargeCM=[0,0.025,0.05,0.075,0.1,0.5,1,2.5,5,7.5,10,15,20,40,70,100,130];
% Make colour scale
colourmatrix=[0,0,175;
    0,0,223;
    0,16,255; % royal blue
    0,48,255;
    0,159,255;
    0,207,255;
    0,255,255; % cyan/light blue
    96,255,159;
    143,255,111;
    191,255,64;
    255,223,0;
    255,175,0;
    255,128,0;
    255,32,0;
    239,0,0;
    191,0,0]; % dark red
colourmatrix=colourmatrix./255;
CBticks=[0:(1/16):1];
DischargeCB=num2str(DischargeCM(:));

figure
m_proj('lambert','long',[-96.25 -60],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
for c=2:length(DischargeCM)
    if c<=length(DischargeCM)
        ind=find(HBCdischarge_geoinfo(:,5)<=DischargeCM(c) & HBCdischarge_geoinfo(:,5)>DischargeCM(c-1));
    end
    m_line(HBCdischarge_geoinfo(ind,4),HBCdischarge_geoinfo(ind,3),'LineStyle','none','Marker','.','MarkerSize',10,'Color',colourmatrix(c-1,:));
    hold on
end
title('HYPE integrated discharge [km3/yr]')
colormap(jet(length(DischargeCM(1:end-1))))
colorbar('Ticks',CBticks,'TickLabels',DischargeCB(1:end,:))


