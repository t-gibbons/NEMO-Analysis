% Get HYPE integrated data ready to put into ANHA4 grid
% this script creates the mat file used in remapping the runoff.
close all
clear all
clc

% Load in Integrated HYPE data
% this file was modified to add decimal places to a column (379) because it was
% being problematic. No values were actually changed.
% The two values at the end say ignore the first 2 rows and first 2 columns.
HYPEdischarge_int = csvread('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBCdischarge_integratedcopy.csv',2,2); 

% Load in Nelson River data provided on April 29, 2016
% Load in the monthly mean file (also have daily and yearly) 
load('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/month_mm.nn')
nelsonriver = zeros(444,3);

i = 1;
for y = 1:length(month_mm)%month_mm(1,1):month_mm(end,1)
    for m = 1:12
         nelsonriver(i,1) = month_mm(y,1);
         nelsonriver(i,2) = m;
         nelsonriver(i,3) = month_mm(y,m+1);
         i                = i+1;
    end
end

% But I need to get the haroids in the HBCdischarge_integratedcopy.csv file

HYPEdischargeinfo=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBSdischarge_HYPE151030.csv'); 
%HYPEdischargeinfo.textdata(1,1:400)=rname;

% now with this data in this format, I want to replace the HYPEint Nelson
% river data with this one. ind 25 is Jan 1979 to ind 408 which is Dec
% 2010. In the hYPEint data, nelson river is column 4
HYPEdischarge_int(:,4)=nelsonriver(25:408,3);
HYPEdischarge_intNelson=HYPEdischarge_int;%*3600*24; % m3/day
% for n=1:398
% HYPEdischarge_intNelson(:,n)=HYPEdischarge_intNelson(:,n).*HYPEdischargeinfo.data(:,1); % multiply by numbers in month, so should have m3/month
% end
% HYPEdischarge_intNelson=HYPEdischarge_intNelson/1000000000; % now should have km3/month

% test=sum(HYPEdischarge_intNelson);
% test=test/32;
% sum(test)
szHYPE=size(HYPEdischarge_intNelson);

HYPEdischarge_intNelson(:,:,2)=zeros(szHYPE);
HYPEdischarge_intNelson(:,:,3)=zeros(szHYPE);

%save('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HYPEdischarge_intNelson.csv','HYPEdischarge_intNelson','-ascii') % Units of m3/s





haroidCOL=HYPEdischargeinfo.textdata(2,3:400); % HAROIDS

% get HYPE integrated into annual means 

HYPEsz=size(HYPEdischarge_int);
months=1:12;
for ny=1:(HYPEsz(1)/12)
temp=mean(HYPEdischarge_int(months,:));
temp=temp*3600*24*sum(HYPEdischargeinfo.data(months,1)); % km3/yr
HYPEdischarge_intYM(ny,:)=temp;
months=months+12;
end

HYPEdischarge_intYM=mean(HYPEdischarge_intYM);
%HYPEdischarge_intmean=mean(HYPEdischarge_int); % mean m3/s over the whole time series

load /mnt/storage0/xhu/PROGRAM/RUNOFF/runoff_HBC/HBSdischarge_geoinfo.csv

% Now I want to put the annual mean of each river into
% HBSdischarge_geoinfo.csv file by mathing the haroids

for i=1:length(HBSdischarge_geoinfo)
    
    
     temphid_char=haroidCOL{1,i}; % get cell haroid to str
     temphid_num=str2double(temphid_char(2:end)) % convert str haroid to number
     ind = find(HBSdischarge_geoinfo(:,1)==temphid_num);
     
     
     HBSdischarge_geoinfo(ind,5)= HYPEdischarge_intYM(1,i);
% also want to put lat and lon into intNelson matrix
    lat=HBSdischarge_geoinfo(ind,3); % latitude
    lon=HBSdischarge_geoinfo(ind,4); % longitude
    HYPEdischarge_intNelson(:,i,2)=lat; % put latitude in each row in the column
    HYPEdischarge_intNelson(:,i,3)=lon; % put longitude in each row in the column
end

save('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HYPEint_HBC/HYPEdischarge_intNelson.mat','HYPEdischarge_intNelson') % Units of m3/s


HBSdischarge_geoinfo(:,5)= HBSdischarge_geoinfo(:,5)/1000000000; % change form m3 to km3


 % saving HBSdischarge_geoinfo which has the haroid, km2 drainage area,
 % lat, lon, discharge km3/year
save('HBSdischarge_geoinfo.csv','HBSdischarge_geoinfo','-ascii')

disp('Check:')
disp('Max river discharge:')
max(HBSdischarge_geoinfo(:,5))
disp('Min river discharge:')
min(HBSdischarge_geoinfo(:,5))
disp('HBC annual river discharge:')
sum(HBSdischarge_geoinfo(:,5))

DischargeCM=[0,0.025,0.05,0.075,0.1,0.5,1,2.5,5,7.5,10,15,20,40,70,100,130];
% Make colour scale
colourmatrix=[0,0,175;
    0,0,223;
    0,16,255; % royal blue
    0,48,255;
  %  0,111,255;
    0,159,255;
    0,207,255;
    0,255,255; % cyan/light blue
  %  48,255,207;
    96,255,159;
    143,255,111;
    191,255,64;
%    239,255,16;
    255,223,0;
    255,175,0;
    255,128,0;
  %  255,80,0;
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
        ind=find(HBSdischarge_geoinfo(:,5)<=DischargeCM(c) & HBSdischarge_geoinfo(:,5)>DischargeCM(c-1));       
    end
    m_line(HBSdischarge_geoinfo(ind,4),HBSdischarge_geoinfo(ind,3),'LineStyle','none','Marker','.','MarkerSize',10,'Color',colourmatrix(c-1,:));
    hold on
end
title('HYPE integrated discharge')
colormap(jet(length(DischargeCM(1:end-1))))
colorbar('Ticks',CBticks,'TickLabels',DischargeCB(1:end,:))

