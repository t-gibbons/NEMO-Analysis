% COMPARING RUNOFF
% Currently (20/11/2015) NEMO is using Dai and Trenberth inter-annual
% fields for river runoff data. Greg McCullough has sent preliminary output
% from the HYPE continental scale model that BaySys Team 2 is using to
% prepare river discharge data for the Hudson Bay System Watershed. 

clear all
close all
clc

DTrunoff = GetNcVar('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/DaiTren/runoff.daitren.iaf.10FEB2011.nc','runoff'); % gives runoff at time, lat, lon
DTtime   = GetNcVar('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/DaiTren/runoff.daitren.iaf.10FEB2011.nc','Time'); % time, days since 1948-01-01 00:00:00
DTlon    = GetNcVar('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/DaiTren/runoff.daitren.iaf.10FEB2011.nc','xc'); % longitude, degrees east
DTlat    = GetNcVar('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/DaiTren/runoff.daitren.iaf.10FEB2011.nc','yc'); % latitude, degrees north
DTarea   = GetNcVar('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/DaiTren/runoff.daitren.iaf.10FEB2011.nc','area'); % area of grid cell,m2 

tavrun = squeeze(mean(DTrunoff(:,1:180,1:360))); % total time average of runoff

figDT = figure;
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
[h] = m_pcolor(DTlon,DTlat,tavrun);
shading flat
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
cmap = colormap(jet);
c = colorbar;
y = ylabel(c,'kg/s/m2','FontSize',14,'FontWeight','bold');
title('Dai and Trenberth runoff data presently using in NEMO')
print(figDT,'-dpng','DTrunoff_map')
 

% Separate the data in to different regions -- find the indices
% DTlat-- ind 141 to 163 is 50N to 72.5N
% DTlon--ind 84 is 96.5W, ind 117 is 63.5W

% Indices for HBC
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(141:152,85:115),DTlat(141:152,85:115),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
hold on
m_line(DTlon(153:155,85:111),DTlat(153:155,85:111),'LineStyle','none','Marker','.','MarkerSize',7,'Color','g');
hold on
m_line(DTlon(156:157,94:111),DTlat(156:157,94:111),'LineStyle','none','Marker','.','MarkerSize',7,'Color','y');
hold on
m_line(DTlon(158:160,97:111),DTlat(158:160,97:111),'LineStyle','none','Marker','.','MarkerSize',7,'Color','m');
title('DT data section for HBC')

DThbc1 = DTrunoff(361:720,141:152,85:115);
DThbc1 = squeeze(sum(DThbc1,2));
DThbc1 = sum(DThbc1,2);

DThbc2 = DTrunoff(361:720,153:155,85:111);
DThbc2 = squeeze(sum(DThbc2,2));
DThbc2 = sum(DThbc2,2);

DThbc3 = DTrunoff(361:720,156:157,94:111);
DThbc3 = squeeze(sum(DThbc2,2));
DThbc3 = sum(DThbc2,2);

DThbc4 = DTrunoff(361:720,158:160,97:111);
DThbc4 = squeeze(sum(DThbc2,2));
DThbc4 = sum(DThbc2,2);

DThbc=DThbc1+DThbc2 +DThbc3+DThbc4;


% Indices for James Bay
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(141:145,97:103),DTlat(141:145,97:103),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('DT data section for James Bay')

% Indices for Hudson Bay Proper
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(146:154,85:102),DTlat(146:154,85:102),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
hold on
m_line(DTlon(146:152,103:105),DTlat(146:152,103:105),'LineStyle','none','Marker','.','MarkerSize',7,'Color','b');
hold on
m_line(DTlon(155:155,85:97),DTlat(155:155,85:97),'LineStyle','none','Marker','.','MarkerSize',7,'Color','g');
hold on
m_line(DTlon(156:156,85:95),DTlat(156:156,85:95),'LineStyle','none','Marker','.','MarkerSize',7,'Color','y');
title('DT data section for Hudson Bay Proper')

% Indices for Hudson Bay Proper --EAST
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(146:153,100:102),DTlat(146:153,100:102),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
hold on
m_line(DTlon(146:152,103:105),DTlat(146:152,103:105),'LineStyle','none','Marker','.','MarkerSize',7,'Color','b');
title('DT data section for Hudson Bay Proper-EAST')

% Indices for Hudson Bay Proper -- WEST
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(146:154,85:99),DTlat(146:154,85:99),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
hold on
m_line(DTlon(155:155,85:97),DTlat(155:155,85:97),'LineStyle','none','Marker','.','MarkerSize',7,'Color','g');
hold on
m_line(DTlon(156:156,85:95),DTlat(156:156,85:95),'LineStyle','none','Marker','.','MarkerSize',7,'Color','y');
title('DT data section for Hudson Bay Proper-WEST')

% Indices for Hudson Strait and Ungava Bay
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(153:155,103:117),DTlat(153:155,103:117),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
hold on
m_line(DTlon(148:152,107:117),DTlat(148:152,107:117),'LineStyle','none','Marker','.','MarkerSize',7,'Color','b');
title('DT data section for Hudson Strait and Ungava Bay')


% Indices for Hudson Strait
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(153:155,103:117),DTlat(153:155,103:117),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
hold on
m_line(DTlon(152:152,107:117),DTlat(152:152,107:117),'LineStyle','none','Marker','.','MarkerSize',7,'Color','b');
title('DT data section for Hudson Strait')

% Indices for Ungava Bay
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(148:151,110:117),DTlat(148:151,110:117),'LineStyle','none','Marker','.','MarkerSize',7,'Color','b');
title('DT data section for Ungava Bay')


% Indices for Foxe Basin
figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.2 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(DTlon(158:161,96:109),DTlat(158:161,96:109),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
hold on
m_line(DTlon(157:157,93:109),DTlat(157:157,93:109),'LineStyle','none','Marker','.','MarkerSize',7,'Color','y');
hold on
m_line(DTlon(156:156,96:109),DTlat(156:156,96:109),'LineStyle','none','Marker','.','MarkerSize',7,'Color','g');
hold on
m_line(DTlon(155:155,98:102),DTlat(155:155,98:102),'LineStyle','none','Marker','.','MarkerSize',7,'Color','b');
title('DT data section for Foxe Basin')

% HAVE THE POINTS FOR EACH REGION, SO NOW AVERAGE OVER TIME AND GET A SUM
% OF ALL THE FRESHWATER INPUT

% DO TIME INTERVAL OF 1979(JAN) TO DEC 2010

% Get a date array
DTdate = zeros(720,1);
DTdate(:,1) = datenum(1949,01,01);
DTtime_r = reshape(DTtime,720,1);
DTdate(:,1) = DTdate(:,1)+DTtime_r(:,1);
DTdatev = datevec(DTdate);
% the dates I want to compare are indices: 361 to 720 (dec 2008)

% James Bay -- lat=141:145,lon=97:103
DTjb = DTrunoff(361:720,141:145,97:103);
DTjb = squeeze(sum(DTjb,2));
DTjb = sum(DTjb,2);

figure
plot(DTdate(361:720,1),DTjb(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for James Bay, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% plot time series for Hudson Strait and Ungava Bay
 
DThsub1 = DTrunoff(361:720,153:155,103:117);
DThsub1 = squeeze(sum(DThsub1,2));
DThsub1 = sum(DThsub1,2);

DThsub2 = DTrunoff(361:720,148:152,107:117);
DThsub2 = squeeze(sum(DThsub2,2));
DThsub2 = sum(DThsub2,2);

DThsub=DThsub1+DThsub2;

figure
plot(DTdate(361:720,1),DThsub(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Strait, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% plot time series for Hudson Strait
 
DThs1 = DTrunoff(361:720,153:155,103:117);
DThs1 = squeeze(sum(DThs1,2));
DThs1 = sum(DThs1,2);

DThs2 = DTrunoff(361:720,152:152,107:117);
DThs2 = squeeze(sum(DThs2,2));
DThs2 = sum(DThs2,2);

DThs=DThs1+DThs2;

figure
plot(DTdate(361:720,1),DThs(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Strait, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% plot time series for Ungava Bay
 
DTub1 = DTrunoff(361:720,148:151,110:117);
DTub1 = squeeze(sum(DTub1,2));
DTub = sum(DTub1,2);

figure
plot(DTdate(361:720,1),DTub(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Ungava Bay, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )



% Plot time series for Foxe Basin

DTfb1 = DTrunoff(361:720,158:161,96:109);
DTfb1 = squeeze(sum(DTfb1,2));
DTfb1 = sum(DTfb1,2);

DTfb2 = DTrunoff(361:720,157:157,93:109);
DTfb2 = squeeze(sum(DTfb2,2));
DTfb2 = sum(DTfb2,2);

DTfb3 = DTrunoff(361:720,156:156,96:109);
DTfb3 = squeeze(sum(DTfb3,2));
DTfb3 = sum(DTfb3,2);

DTfb4 = DTrunoff(361:720,155:155,98:102);
DTfb4 = squeeze(sum(DTfb4,2));
DTfb4 = sum(DTfb4,2);

DTfb=DTfb1+DTfb2+DTfb3+DTfb4;

figure
plot(DTdate(361:720,1),DTfb(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Foxe Basin, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Plot time series for Hudson Bay Proper

DThb1 = DTrunoff(361:720,146:154,85:102);
DThb1 = squeeze(sum(DThb1,2));
DThb1 = sum(DThb1,2);

DThb2 = DTrunoff(361:720,146:152,103:105);
DThb2 = squeeze(sum(DThb2,2));
DThb2 = sum(DThb2,2);

DThb3 = DTrunoff(361:720,155:155,85:97);
DThb3 = squeeze(sum(DThb3,2));
DThb3 = sum(DThb3,2);

DThb4 = DTrunoff(361:720,156:156,85:95);
DThb4 = squeeze(sum(DThb4,2));
DThb4 = sum(DThb4,2);

DThb = DThb1+DThb2+DThb3+DThb4;

figure
plot(DTdate(361:720,1),DThb(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Bay, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Plot time series for Hudson Bay Proper-EAST

DThb1E = DTrunoff(361:720,146:153,100:102);
DThb1E = squeeze(sum(DThb1E,2));
DThb1E = sum(DThb1E,2);

DThb2E = DTrunoff(361:720,146:152,103:105);
DThb2E = squeeze(sum(DThb2E,2));
DThb2E = sum(DThb2E,2);

DThbE = DThb1E+DThb2E;

figure
plot(DTdate(361:720,1),DThbE(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Bay EAST, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Plot time series for Hudson Bay Proper-WEST

DThb1W = DTrunoff(361:720,146:154,85:99);
DThb1W = squeeze(sum(DThb1W,2));
DThb1W = sum(DThb1W,2);

DThb2W = DTrunoff(361:720,155:155,85:97);
DThb2W = squeeze(sum(DThb2W,2));
DThb2W = sum(DThb2W,2);

DThb3W = DTrunoff(361:720,156:156,85:95);
DThb3W = squeeze(sum(DThb3W,2));
DThb3W = sum(DThb3W,2);

DThbW = DThb1W+DThb2W+DThb3W;

figure
plot(DTdate(361:720,1),DThbW(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Bay WEST, 1979-2008')
y=ylabel('Kg/(sm^2)','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% GET UNITS CONVERTED TO KM3/month  *** THIS WILL HAVE TO BE CHECKED *** 

% so to convert from kg/s/m2 to m3/s I need to divide by the density of
% water. Assume 1000kg/m3. And then I multiply the result (units of m/s) by the size of
% the grid cell (m2).
DTrunoff_vf =DTrunoff/1000; % divide by density
for t = 1:720
    tempDTrunoff = squeeze(DTrunoff_vf(t,:,:));
    DTrunoff_vf(t,1:180,1:360) = tempDTrunoff.*DTarea; 
end

load numdays_month % this from the HYPEdata from Jan 1979 to Dec 2010, and contains the number of days in each month
% The DT data does not do leap years, so need to change all 29 to 28
for n=1:length(numdays_month)
if numdays_month(n)==29
    numdays_month(n)=28;
end
end

%% 

DThbc1_vf = DTrunoff_vf(361:720,141:152,85:115);
DThbc1_vf = squeeze(sum(DThbc1_vf,2));
DThbc1_vf = sum(DThbc1_vf,2);

DThbc2_vf = DTrunoff_vf(361:720,153:155,85:111);
DThbc2_vf = squeeze(sum(DThbc2_vf,2));
DThbc2_vf = sum(DThbc2_vf,2);

DThbc3_vf = DTrunoff_vf(361:720,156:157,94:111);
DThbc3_vf = squeeze(sum(DThbc2_vf,2));
DThbc3_vf = sum(DThbc2_vf,2);

DThbc4_vf = DTrunoff_vf(361:720,158:160,97:111);
DThbc4_vf = squeeze(sum(DThbc2_vf,2));
DThbc4_vf = sum(DThbc2_vf,2);

DThbc_vf=DThbc1_vf+DThbc2_vf +DThbc3_vf+DThbc4_vf;


DThbc_vf=DThbc_vf*3600*24.*numdays_month(1:360);
DThbc_vf = DThbc_vf/1000000000; % to get km3


% DThbc_vf=DTrunoff_vf(361:720,:,:);
% DThbc_vf = squeeze(sum(DThbc_vf,2));
% DThbc_vf = sum(DThbc_vf,2);
% 
% DThbc_vf= DThbc_vf*3600*24.*numdays_month(1:360); % to get input per month
% DThbc_vf = DThbc_vf/1000000000; % to get km3
%%

% for n=361:720
% DTrunoff_vf = DTrunoff_vf(n,:,:)*3600*24.*numdays_month(1:360);
% DTrunoff_vf = DTrunoff_vf/1000000000; % to get km3
% end
% James Bay -- lat=141:145,lon=97:103
DTjb_vf = DTrunoff_vf(361:720,141:145,97:103);
DTjb_vf = squeeze(sum(DTjb_vf,2));
DTjb_vf = sum(DTjb_vf,2);

DTjb_vf= DTjb_vf*3600*24.*numdays_month(1:360); % to get input per month
DTjb_vf = DTjb_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DTjb_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for James Bay, 1979-2008')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% plot time series for Hudson Strait and Ungava Bay
 
DThsub1_vf = DTrunoff_vf(361:720,153:155,103:117);
DThsub1_vf = squeeze(sum(DThsub1_vf,2));
DThsub1_vf = sum(DThsub1_vf,2);

DThsub2_vf = DTrunoff_vf(361:720,148:152,107:117);
DThsub2_vf = squeeze(sum(DThsub2_vf,2));
DThsub2_vf = sum(DThsub2_vf,2);

DThsub_vf=DThsub1_vf+DThsub2_vf;

DThsub_vf= DThsub_vf*3600*24.*numdays_month(1:360); % to get input per month
DThsub_vf = DThsub_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DThsub_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Strait and Ungava Bay, 1979-2008')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot time series for Hudson Strait

DThs1_vf = DTrunoff_vf(361:720,153:155,103:117);
DThs1_vf = squeeze(sum(DThs1_vf,2));
DThs1_vf = sum(DThs1_vf,2);

DThs2_vf = DTrunoff_vf(361:720,152:152,107:117);
DThs2_vf = squeeze(sum(DThs2_vf,2));
DThs2_vf = sum(DThs2_vf,2);

DThs_vf=DThs1_vf+DThs2_vf;

DThs_vf= DThs_vf*3600*24.*numdays_month(1:360); % to get input per month
DThs_vf = DThs_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DThs_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Strait, 1979-2008')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% plot time series for Ungava Bay

DTub1_vf = DTrunoff_vf(361:720,148:151,110:117);
DTub1_vf = squeeze(sum(DTub1_vf,2));
DTub_vf = sum(DTub1_vf,2);

DTub_vf= DTub_vf*3600*24.*numdays_month(1:360); % to get input per month
DTub_vf = DTub_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DTub_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Ungava Bay, 1979-2008')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Plot time series for Foxe Basin

DTfb1_vf = DTrunoff_vf(361:720,158:161,96:109);
DTfb1_vf = squeeze(sum(DTfb1_vf,2));
DTfb1_vf = sum(DTfb1_vf,2);

DTfb2_vf = DTrunoff_vf(361:720,157:157,93:109);
DTfb2_vf = squeeze(sum(DTfb2_vf,2));
DTfb2_vf = sum(DTfb2_vf,2);

DTfb3_vf = DTrunoff_vf(361:720,156:156,96:109);
DTfb3_vf = squeeze(sum(DTfb3_vf,2));
DTfb3_vf = sum(DTfb3_vf,2);

DTfb4_vf = DTrunoff_vf(361:720,155:155,98:102);
DTfb4_vf = squeeze(sum(DTfb4_vf,2));
DTfb4_vf = sum(DTfb4_vf,2);

DTfb_vf=DTfb1_vf+DTfb2_vf+DTfb3_vf+DTfb4_vf;

DTfb_vf= DTfb_vf*3600*24.*numdays_month(1:360); % to get input per month
DTfb_vf = DTfb_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DTfb_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Foxe Basin, 1979-2008')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Plot time series for Hudson Bay Proper

DThb1_vf = DTrunoff_vf(361:720,146:154,85:102);
DThb1_vf = squeeze(sum(DThb1_vf,2));
DThb1_vf = sum(DThb1_vf,2);

DThb2_vf = DTrunoff_vf(361:720,146:152,103:105);
DThb2_vf = squeeze(sum(DThb2_vf,2));
DThb2_vf = sum(DThb2_vf,2);

DThb3_vf = DTrunoff_vf(361:720,155:155,85:97);
DThb3_vf = squeeze(sum(DThb3_vf,2));
DThb3_vf = sum(DThb3_vf,2);

DThb4_vf = DTrunoff_vf(361:720,156:156,85:95);
DThb4_vf = squeeze(sum(DThb4_vf,2));
DThb4_vf = sum(DThb4_vf,2);

DThb_vf = DThb1_vf+DThb2_vf+DThb3_vf+DThb4_vf;

DThb_vf= DThb_vf*3600*24.*numdays_month(1:360); % to get input per month
DThb_vf = DThb_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DThb_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Bay, 1979-2008')
y=ylabel('Km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Plot time series for Hudson Bay Proper-EAST

DThb1E_vf = DTrunoff_vf(361:720,146:153,100:102);
DThb1E_vf = squeeze(sum(DThb1E_vf,2));
DThb1E_vf = sum(DThb1E_vf,2);

DThb2E_vf = DTrunoff_vf(361:720,146:152,103:105);
DThb2E_vf = squeeze(sum(DThb2E_vf,2));
DThb2E_vf = sum(DThb2E_vf,2);

DThbE_vf = DThb1E_vf+DThb2E_vf;

DThbE_vf= DThbE_vf*3600*24.*numdays_month(1:360); % to get input per month
DThbE_vf = DThbE_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DThbE_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Bay EAST, 1979-2008')
y=ylabel('Km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Plot time series for Hudson Bay Proper-WEST

DThb1W_vf = DTrunoff_vf(361:720,146:154,85:99);
DThb1W_vf = squeeze(sum(DThb1W_vf,2));
DThb1W_vf = sum(DThb1W_vf,2);

DThb2W_vf = DTrunoff_vf(361:720,155:155,85:97);
DThb2W_vf = squeeze(sum(DThb2W_vf,2));
DThb2W_vf = sum(DThb2W_vf,2);

DThb3W_vf = DTrunoff_vf(361:720,156:156,85:95);
DThb3W_vf = squeeze(sum(DThb3W_vf,2));
DThb3W_vf = sum(DThb3W_vf,2);

DThbW_vf = DThb1W_vf+DThb2W_vf+DThb3W_vf;

DThbW_vf= DThbW_vf*3600*24.*numdays_month(1:360); % to get input per month
DThbW_vf = DThbW_vf/1000000000; % to get km3

figure
plot(DTdate(361:720,1),DThbW_vf(:,1),'LineWidth',2)
title('River input from Dai Trenberth data for Hudson Bay WEST, 1979-2008')
y=ylabel('Km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:6:360,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )




%**************************************************************************
% HYPE DATA

% Currently in this directory there are 2 csv files and one xlsk file.
% There are also 3 kmz (Google Earth files) and a pdf explaining the files.
HYPErunoff='/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE';

% GET FIRST ROW WITH RIVER NAMES
fid=fopen('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBSdischarge_integrated151116.csv');
data=textscan(fid, '%18q','Delimiter',',','EmptyValue',-Inf);
fclose(fid);

rnametemp=[data{1,1}];
rname(1,1:400)=rnametemp(1:400,1);

clear data rnametemp
rname_int=strrep(rname,'.','0');

%HYPEdischarge_int=importdata('/mnt/storage1/natasha/HBSdischarge_HYPE/HBSdischarge_integrated151116.csv'); % best one so far but how to get the river names in the right position
%HYPEdischarge_int=importdata('/mnt/storage1/natasha/HBSdischarge_HYPE/HBCdischarge_integratedcopy.csv'); % best one so far but how to get the river names in the right position

% this file was modified to add decimal places to a column (379) because it was
% being problematic. No values were actually changed.
% The two values at the end say ignore the first 2 rows and first 2 columns.
HYPEdischarge_int=csvread('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBCdischarge_integratedcopy.csv',2,2); 

HYPEdischarge_metadata=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBSdischarge_metadata.xlsx'); % best one so far but how to get the river names in the right position

Lat = HYPEdischarge_metadata.data(:,3);
Lon = HYPEdischarge_metadata.data(:,4);
haroid = HYPEdischarge_metadata.data(:,1);
basin = HYPEdischarge_metadata.textdata(:,5);

% Figure of locations of river input
fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon,Lat,'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids')
% how to get map to display haroid on point??
print(fig,'-dpng','HYPErunoff_map')


% GET FIRST ROW WITH RIVER NAMES
fid=fopen('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBSdischarge_HYPE151030.csv');
data=textscan(fid, '%18q','Delimiter',',','EmptyValue',-Inf);
fclose(fid);

rnametemp=[data{1,1}];
rname(1,1:400)=rnametemp(1:400,1);

clear data rnametemp
rname=strrep(rname,'.','0');

% load in HBSdischarge_HYPE151030.csv -- first column is the number of days
% in each month
HYPEdischargeinfo=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBSdischarge_HYPE151030.csv'); 

HYPEdischarge=csvread('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HBSdischarge_HYPE151030.csv',2,2);

HYPEdischargeinfo.textdata(1,1:400)=rname;
HYPEdischarge_date=datenum(HYPEdischargeinfo.textdata(3:end,1));

HYPEdischargehid = char(HYPEdischargeinfo.textdata(2,3:end)); % get all haroids 
HYPEdischargehid = str2num(HYPEdischargehid(:,2:end));
HYPEdischargehid = reshape(HYPEdischargehid,1,398);

clear HYPEdischargeinfo

% Load in Nelson River data provided on April 29, 2016
% Load in the monthly mean file (also have daily and yearly) 
load('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/month_mm.nn')
nelsonriver=zeros(444,3);
i=1;
for y=1:length(month_mm)%month_mm(1,1):month_mm(end,1)
    for m=1:12
         nelsonriver(i,1)=month_mm(y,1);
         nelsonriver(i,2)=m;
         nelsonriver(i,3)=month_mm(y,m+1);
         i=i+1;
    end
end
% now with this data in this format, I want to replace the HYPEint Nelson
% river data with this one. ind 25 is Jan 1979 to ind 408 which is Dec
% 2010. In the hYPEint data, nelson river is column 4
HYPEdischarge_int(:,4)=nelsonriver(25:408,3);

% Plot by region

basinregion = {'JB','HB','HS','FB'};
% not integrated data
JBInd =zeros(45,1);
HBInd=zeros(200,1);
HSInd = zeros(100,1);
FBInd = zeros(100,1);
% Integrated data
JBInd_int =zeros(45,1);
HBInd_int=zeros(200,1);
HSInd_int = zeros(100,1);
FBInd_int = zeros(100,1);

c = 1;
for n=1:4
    for m=1:length(HYPEdischarge)
        t = strcmp(basinregion(1,n),basin(m,1)); % compare strings
        if t && n==1
            JBInd(c,1)=m;
            JBInd_int(c,1)=m;
            c = c+1;
        end
        if t && n==2
            HBInd(c,1)=m;
            HBInd_int(c,1)=m;
            c = c+1;
        end
        if t && n==3
            HSUBInd(c,1)=m;
            HSUBInd_int(c,1)=m;
            c = c+1;
        end
        if t && n==4
            FBInd(c,1)=m;
            FBInd_int(c,1)=m;
            c = c+1;
        end
    end
    c =1;
end

% so now have indices of the haroids in each basin

% JAMES BAY
JBInd = JBInd(JBInd~=0); % remove zeros at end
JBInd=JBInd-1; % -1 to match the number for the Lon and lat
JBhid = haroid(JBInd(:,1),1); % get all haroids in James Bay

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(JBInd),Lat(JBInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in James Bay')

% now that I have the right haroids for each region, get timeseries data
JBrunoff = zeros(length(HYPEdischarge_date)+1,length(JBhid));
for n=1:length(JBhid)
    i = find(JBhid(n,1)==HYPEdischargehid(1,:));
    JBrunoff(1,n) = JBhid(n); % first row is haroid
    JBrunoff(2:end, n)=HYPEdischarge(1:end,i);
end

JBrunoff_int = zeros(length(HYPEdischarge_date)+1,length(JBhid));
for n=1:length(JBhid)
    i = find(JBhid(n,1)==HYPEdischargehid(1,:));
    JBrunoff_int(1,n) = JBhid(n); % first row is haroid
    JBrunoff_int(2:end,n)=HYPEdischarge_int(1:end,i);
end

% HUDSON BAY
HBInd = HBInd(HBInd~=0); % remove zeros at end
HBInd=HBInd-1; % -1 to match the number for the Lon and lat
HBhid = haroid(HBInd(:,1),1); % get all haroids in Hudson Bay

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(HBInd),Lat(HBInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Hudson Bay')

% now that I have the right haroids for each region, get timeseries data
HBrunoff = zeros(length(HYPEdischarge_date)+1,length(HBhid));
for n=1:length(HBhid)
    i = find(HBhid(n,1)==HYPEdischargehid(1,:));
    HBrunoff(1,n) = HBhid(n); % first row is haroid, REST IS TIME SERIES DATA
    HBrunoff(2:end, n)=HYPEdischarge(1:end,i);
end

HBrunoff_int = zeros(length(HYPEdischarge_date)+1,length(HBhid));
for n=1:length(HBhid)
    i = find(HBhid(n,1)==HYPEdischargehid(1,:));
    HBrunoff_int(1,n) = HBhid(n); % first row is haroid, REST IS TIME SERIES DATA
    HBrunoff_int(2:end, n)=HYPEdischarge_int(1:end,i);
end


% HUDSON BAY EAST

a=find(Lon>-80.7 & Lon < -75.7 );
b=find(Lat < 63 & Lat > 54.5);
HBEInd=intersect(a,b);

HBEhid=haroid(HBEInd(:,1),1);

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(HBEInd),Lat(HBEInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Hudson Bay EAST')

% now that I have the right haroids for each region, get timeseries data
HBErunoff = zeros(length(HYPEdischarge_date)+1,length(HBEhid));
for n=1:length(HBEhid)
    i = find(HBEhid(n,1)==HYPEdischargehid(1,:));
    HBErunoff(1,n) = HBEhid(n); % first row is haroid, REST IS TIME SERIES DATA
    HBErunoff(2:end, n)=HYPEdischarge(1:end,i);
end

HBErunoff_int = zeros(length(HYPEdischarge_date)+1,length(HBEhid));
for n=1:length(HBEhid)
    i = find(HBEhid(n,1)==HYPEdischargehid(1,:));
    HBErunoff_int(1,n) = HBEhid(n); % first row is haroid, REST IS TIME SERIES DATA
    HBErunoff_int(2:end,n)=HYPEdischarge_int(1:end,i);
end

% HUDSON BAY WEST

a=find(Lon<-80.3);
b=find(Lat < 66 & Lat > 54.5);
HBWInd=intersect(a,b);

t1=find(Lon<-80.3 & Lon > -85.5);
t2=find(Lat < 66 & Lat > 64.2);
notincl1=intersect(t1,t2);

t3=find(Lon<-80.3 & Lon > -82.8);
t4=find(Lat < 64.5 & Lat > 64);
notincl2=intersect(t3,t4);

notinc=vertcat(notincl1,notincl2);

[t5,hbi,nii]=intersect(HBWInd,notinc)

HBWInd(hbi)=0;
HBWInd = HBWInd(HBWInd~=0);

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(notinc),Lat(notinc),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Hudson Bay notinclude')

HBWhid=haroid(HBWInd(:,1),1);

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(HBWInd),Lat(HBWInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Hudson Bay WEST')

% now that I have the right haroids for each region, get timeseries data
HBWrunoff = zeros(length(HYPEdischarge_date)+1,length(HBWhid));
for n=1:length(HBWhid)
    i = find(HBWhid(n,1)==HYPEdischargehid(1,:));
    HBWrunoff(1,n) = HBWhid(n); % first row is haroid, REST IS TIME SERIES DATA
    HBWrunoff(2:end, n)=HYPEdischarge(1:end,i);
end

HBWrunoff_int = zeros(length(HYPEdischarge_date)+1,length(HBWhid));
for n=1:length(HBWhid)
    i = find(HBWhid(n,1)==HYPEdischargehid(1,:));
    HBWrunoff_int(1,n) = HBWhid(n); % first row is haroid, REST IS TIME SERIES DATA
    HBWrunoff_int(2:end,n)=HYPEdischarge_int(1:end,i);
end


% HUDSON STRAIT AND UNGAVA BAY
HSUBInd = HSUBInd(HSUBInd~=0); % remove zeros at end
HSUBInd=HSUBInd-1; % -1 to match the number for the Lon and lat
HSUBhid = haroid(HSUBInd(:,1),1); % get all haroids in Hudson Strait

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(HSUBInd),Lat(HSUBInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Hudson Strait and Ungava Bay')

HSUBrunoff = zeros(length(HYPEdischarge_date)+1,length(HSUBhid));
for n=1:length(HSUBhid)
    i = find(HSUBhid(n,1)==HYPEdischargehid(1,:));
    HSUBrunoff(1,n) = HSUBhid(n); % first row is haroid
    HSUBrunoff(2:end, n)=HYPEdischarge(1:end,i);
end

HSUBrunoff_int = zeros(length(HYPEdischarge_date)+1,length(HSUBhid));
for n=1:length(HSUBhid)
    i = find(HSUBhid(n,1)==HYPEdischargehid(1,:));
    HSUBrunoff_int(1,n) = HSUBhid(n); % first row is haroid
    HSUBrunoff_int(2:end, n)=HYPEdischarge_int(1:end,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Split Hudson Strait and Ungava Bay

% HUDSON STRAIT

a=find(Lon > -77.9 & Lon <= -64);
b=find(Lat < 64.6 & Lat > 60.7);
HSInd=intersect(a,b);

t3=find(Lon<-77.5 & Lon > -78);
t4=find(Lat < 62 & Lat > 60);
notinc=intersect(t3,t4);

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(notinc),Lat(notinc),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Hudson Strait notincluded')

[t5,hsi,nii]=intersect(HSInd,notinc)

HSInd(hsi)=0;
HSInd = HSInd(HSInd~=0);

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(HSInd),Lat(HSInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Hudson Strait')

HShid=haroid(HSInd(:,1),1);

HSrunoff = zeros(length(HYPEdischarge_date)+1,length(HShid));
for n=1:length(HShid)
    i = find(HShid(n,1)==HYPEdischargehid(1,:));
    HSrunoff(1,n) = HShid(n); % first row is haroid
    HSrunoff(2:end, n)=HYPEdischarge(1:end,i);
end

HSrunoff_int = zeros(length(HYPEdischarge_date)+1,length(HShid));
for n=1:length(HShid)
    i = find(HShid(n,1)==HYPEdischargehid(1,:));
    HSrunoff_int(1,n) = HShid(n); % first row is haroid
    HSrunoff_int(2:end, n)=HYPEdischarge_int(1:end,i);
end

% UNGAVA BAY

a=find(Lon > -70.5 & Lon <= -64);
b=find(Lat < 60.8 & Lat > 58);
UBInd=intersect(a,b);

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(UBInd),Lat(UBInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Ungava Bay')

UBhid=haroid(UBInd(:,1),1);

UBrunoff = zeros(length(HYPEdischarge_date)+1,length(UBhid));
for n=1:length(UBhid)
    i = find(UBhid(n,1)==HYPEdischargehid(1,:));
    UBrunoff(1,n) = UBhid(n); % first row is haroid
    UBrunoff(2:end, n)=HYPEdischarge(1:end,i); 
end

UBrunoff_int = zeros(length(HYPEdischarge_date)+1,length(UBhid));
for n=1:length(UBhid)
    i = find(UBhid(n,1)==HYPEdischargehid(1,:));
    UBrunoff_int(1,n) = UBhid(n); % first row is haroid
    UBrunoff_int(2:end,n)=HYPEdischarge_int(1:end,i);
end

% FOXE BASIN
FBInd = FBInd(FBInd~=0); % remove zeros at end
FBInd=FBInd-1; % -1 to match the number for the Lon and lat
FBhid = haroid(FBInd(:,1),1); % get all haroids in Foxe Basin

fig=figure
m_proj('lambert','long',[-96.25 -64],'lat',[50.3 71.5]);
hold on
m_grid('box','fancy','tickdir','in','fontweight','bold');
m_coast('patch',[.7 .7 .7]);
m_line(Lon(FBInd),Lat(FBInd),'LineStyle','none','Marker','.','MarkerSize',7,'Color','r');
title('Location of HYPE haroids in Foxe Basin')
close

FBrunoff = zeros(length(HYPEdischarge_date)+1,length(FBhid));
for n=1:length(FBhid)
    i = find(FBhid(n,1)==HYPEdischargehid(1,:));
    FBrunoff(1,n) = FBhid(n); % first row is haroid
    FBrunoff(2:end, n)=HYPEdischarge(1:end,i);
end

FBrunoff_int = zeros(length(HYPEdischarge_date)+1,length(FBhid));
for n=1:length(FBhid)
    i = find(FBhid(n,1)==HYPEdischargehid(1,:));
    FBrunoff_int(1,n) = FBhid(n); % first row is haroid
    FBrunoff_int(2:end, n)=HYPEdischarge_int(1:end,i);
end

% now have time series data of runoff for each region.
% Best way to compare is get sum of all river input and plot that over the
% years? will have monthly sums from January 1979 to December 2010

JBrunoff_sum=sum(JBrunoff,2); 

% Plot total runoff per month
JBrunoff_sum = JBrunoff_sum(2:end)*3600*24.*numdays_month;
JBrunoff_sum = JBrunoff_sum/1000000000;

figure
plot(HYPEdischarge_date, JBrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for James Bay')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% James Bay Integrated data
JBrunoff_intsum=sum(JBrunoff_int,2); 

% Plot total runoff per month
JBrunoff_intsum = JBrunoff_intsum(2:end)*3600*24.*numdays_month;
JBrunoff_intsum = JBrunoff_intsum/1000000000;

figure
plot(HYPEdischarge_date, JBrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for James Bay')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Hudson Bay
HBrunoff_sum=sum(HBrunoff,2);
% Plot total runoff per month
HBrunoff_sum = HBrunoff_sum(2:end)*3600*24.*numdays_month; % get to m3/month
HBrunoff_sum = HBrunoff_sum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HBrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for Hudson Bay Proper')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Hudson Bay integrated
HBrunoff_intsum=sum(HBrunoff_int,2);
% Plot total runoff per month
HBrunoff_intsum = HBrunoff_intsum(2:end)*3600*24.*numdays_month; % get to m3/month
HBrunoff_intsum = HBrunoff_intsum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HBrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for Hudson Bay Proper')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Hudson Bay East
HBErunoff_sum=sum(HBErunoff,2);
% Plot total runoff per month
HBErunoff_sum = HBErunoff_sum(2:end)*3600*24.*numdays_month; % get to m3/month
HBErunoff_sum = HBErunoff_sum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HBErunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for Hudson Bay EAST')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Hudson Bay East integrated
HBErunoff_intsum=sum(HBErunoff_int,2);
% Plot total runoff per month
HBErunoff_intsum = HBErunoff_intsum(2:end)*3600*24.*numdays_month; % get to m3/month
HBErunoff_intsum = HBErunoff_intsum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HBErunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for Hudson Bay EAST')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Hudson Bay West
HBWrunoff_sum=sum(HBWrunoff,2);
% Plot total runoff per month
HBWrunoff_sum = HBWrunoff_sum(2:end)*3600*24.*numdays_month; % get to m3/month
HBWrunoff_sum = HBWrunoff_sum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HBWrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for Hudson Bay WEST')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Hudson Bay West integrated
HBWrunoff_intsum=sum(HBWrunoff_int,2);
% Plot total runoff per month
HBWrunoff_intsum = HBWrunoff_intsum(2:end)*3600*24.*numdays_month; % get to m3/month
HBWrunoff_intsum = HBWrunoff_intsum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HBWrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for Hudson Bay WEST')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Hudson Strait and Ungava Bay
HSUBrunoff_sum=sum(HSUBrunoff,2);
HSUBrunoff_sum = HSUBrunoff_sum(2:end)*3600*24.*numdays_month; % get to m3/month
HSUBrunoff_sum = HSUBrunoff_sum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HSUBrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for Hudson Strait and Ungava Bay')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Hudson Strait and Ungava Bay integrated
HSUBrunoff_intsum=sum(HSUBrunoff_int,2);
HSUBrunoff_intsum = HSUBrunoff_intsum(2:end)*3600*24.*numdays_month; % get to m3/month
HSUBrunoff_intsum = HSUBrunoff_intsum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HSUBrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for Hudson Strait and Ungava Bay')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Hudson Strait
HSrunoff_sum=sum(HSrunoff,2);
% Plot total runoff per month
HSrunoff_sum = HSrunoff_sum(2:end)*3600*24.*numdays_month; % get to m3/month
HSrunoff_sum = HSrunoff_sum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HSrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for Hudson Strait')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Hudson Strait integrated
HSrunoff_intsum=sum(HSrunoff_int,2);
% Plot total runoff per month
HSrunoff_intsum = HSrunoff_intsum(2:end)*3600*24.*numdays_month; % get to m3/month
HSrunoff_intsum = HSrunoff_intsum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, HSrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for Hudson Strait')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Ungava Bay
UBrunoff_sum=sum(UBrunoff,2);
% Plot total runoff per month
UBrunoff_sum = UBrunoff_sum(2:end)*3600*24.*numdays_month; % get to m3/month
UBrunoff_sum = UBrunoff_sum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, UBrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for Ungava Bay')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Ungava Bay integrated
UBrunoff_intsum=sum(UBrunoff_int,2);
% Plot total runoff per month
UBrunoff_intsum = UBrunoff_intsum(2:end)*3600*24.*numdays_month; % get to m3/month
UBrunoff_intsum = UBrunoff_intsum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, UBrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for Ungava Bay')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


% Foxe Basin Runoff
FBrunoff_sum=sum(FBrunoff,2);
% Plot total runoff per month
FBrunoff_sum = FBrunoff_sum(2:end)*3600*24.*numdays_month; % get to m3/month
FBrunoff_sum = FBrunoff_sum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, FBrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for Foxe Basin')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )

% Foxe Basin Runoff integrated
FBrunoff_intsum=sum(FBrunoff_int,2);
% Plot total runoff per month
FBrunoff_intsum = FBrunoff_intsum(2:end)*3600*24.*numdays_month; % get to m3/month
FBrunoff_intsum = FBrunoff_intsum/1000000000; % convert to km3

figure
plot(HYPEdischarge_date, FBrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for Foxe Basin')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )


%Plot all HYPE regions on the same plot
figure
plot(HYPEdischarge_date, HBrunoff_sum,'LineWidth',2)
hold on
plot(HYPEdischarge_date, JBrunoff_sum,'LineWidth',2)
hold on
plot(HYPEdischarge_date, HSUBrunoff_sum,'LineWidth',2)
hold on
plot(HYPEdischarge_date, FBrunoff_sum,'LineWidth',2)
title('Monthly river runoff from HYPE for the Hudson Bay System')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
legend('Hudson Bay','James Bay','Hudson Strait','Foxe Basin','Location','northeast')


%Plot all integrated HYPE regions on the same plot
figure
plot(HYPEdischarge_date, HBrunoff_intsum,'LineWidth',2)
hold on
plot(HYPEdischarge_date, JBrunoff_intsum,'LineWidth',2)
hold on
plot(HYPEdischarge_date, HSUBrunoff_intsum,'LineWidth',2)
hold on
plot(HYPEdischarge_date, FBrunoff_intsum,'LineWidth',2)
title('Monthly river runoff from integrated HYPE for the Hudson Bay Complex')
y=ylabel('km^3/month','FontSize',14,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
xlabel('Date (month/year)','FontSize',14)
xlim([min(HYPEdischarge_date) max(HYPEdischarge_date)])
tick_locations = datenum(1979,1:6:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
legend('Hudson Bay','James Bay','Hudson Strait','Foxe Basin','Location','northeast')



%*************************************************************************
% plot to compare the two sets of data

% Compare Hudson Bay Runoff with not integrated HYPE
figHB=figure('name', 'Runoff comparison for Hudson Bay');
plot(DTdate(361:720,1), HBrunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThb_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HBrunoff_sum(1:360,1)-DThb_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHB, '-dpng','HB_runoff_comparison')

% Compare Hudson Bay Runoff with Integrated HYPE
figHB=figure('name', 'Runoff comparison for Hudson Bay');
plot(DTdate(361:720,1), HBrunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThb_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HBrunoff_intsum(1:360,1)-DThb_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHB, '-dpng','HB_intrunoff_comparison')

%Compare Hudson Bay East Runoff
figHBE=figure('name', 'Runoff comparison for Hudson Bay EAST');
plot(DTdate(361:720,1), HBErunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThbE_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HBErunoff_sum(1:360,1)-DThbE_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Bay East 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHBE, '-dpng','HBEast_runoff_comparison')

%Compare Hudson Bay East Runoff with integrated HYPE
figHBE=figure('name', 'Runoff comparison for Hudson Bay EAST');
plot(DTdate(361:720,1), HBErunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThbE_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HBErunoff_intsum(1:360,1)-DThbE_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Bay East 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHBE, '-dpng','HBEast_intrunoff_comparison')


% Compare Hudson Bay West runoff
figHBW=figure('name', 'Runoff comparison for Hudson Bay WEST');
plot(DTdate(361:720,1), HBWrunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThbW_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HBWrunoff_sum(1:360,1)-DThbW_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Bay West 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHBW, '-dpng','HBWest_runoff_comparison')

% Compare Hudson Bay West runoff with integrated HYPE
figHBW=figure('name', 'Runoff comparison for Hudson Bay WEST');
plot(DTdate(361:720,1), HBWrunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThbW_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HBWrunoff_intsum(1:360,1)-DThbW_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Bay West 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHBW, '-dpng','HBWest_intrunoff_comparison')


%Compare James Bay Runoff
figJB=figure('name', 'Runoff comparison for James Bay');
plot(DTdate(361:720,1), JBrunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DTjb_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), JBrunoff_sum(1:360,1)-DTjb_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for James Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figJB, '-dpng','JB_runoff_comparison')

%Compare James Bay Runoff with integrated HYPE
figJB=figure('name', 'Runoff comparison for James Bay');
plot(DTdate(361:720,1), JBrunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DTjb_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), JBrunoff_intsum(1:360,1)-DTjb_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for James Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figJB, '-dpng','JB_intrunoff_comparison')


% Compare Hudson Strait and Ungava Bay runoff
figHSUB=figure('name', 'Runoff comparison for Hudson Strait and Ungava Bay');
plot(DTdate(361:720,1), HSUBrunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThsub_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HSUBrunoff_sum(1:360,1)-DThsub_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Strait and Ungava Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHSUB,'-dpng','HSUB_runoff_comparison')

% Compare Hudson Strait and Ungava Bay runoff with integrated HYPE
figHSUB=figure('name', 'Runoff comparison for Hudson Strait and Ungava Bay');
plot(DTdate(361:720,1), HSUBrunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThsub_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HSUBrunoff_intsum(1:360,1)-DThsub_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Strait and Ungava Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHSUB,'-dpng','HSUB_intrunoff_comparison')


% Compare Hudson Strait runoff
figHS=figure('name', 'Runoff comparison for Hudson Strait');
plot(DTdate(361:720,1), HSrunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThs_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HSrunoff_sum(1:360,1)-DThs_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Strait 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHS,'-dpng','HS_runoff_comparison')

% Compare Hudson Strait runoff with integrated HYPE
figHS=figure('name', 'Runoff comparison for Hudson Strait');
plot(DTdate(361:720,1), HSrunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DThs_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), HSrunoff_intsum(1:360,1)-DThs_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Hudson Strait 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figHS,'-dpng','HS_intrunoff_comparison')


% Compare Ungava Bay runoff
figUB=figure('name', 'Runoff comparison for Ungava Bay');
plot(DTdate(361:720,1), UBrunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DTub_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), UBrunoff_sum(1:360,1)-DTub_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Ungava Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figUB,'-dpng','UB_runoff_comparison')

% Compare Ungava Bay runoff with integrated HYPE
figUB=figure('name', 'Runoff comparison for Ungava Bay');
plot(DTdate(361:720,1), UBrunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DTub_vf,'LineWidth',2)
hold on
plot(DTdate(361:720,1), UBrunoff_intsum(1:360,1)-DTub_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Ungava Bay 1979-2008')
y=ylabel('km^3/month','FontSize',14);%,'rotation',0);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figUB,'-dpng','UB_intrunoff_comparison')

% Compare Foxe Basin Runoff
figFB=figure('name', 'Runoff comparison for Foxe Basin');
%subplot(2,1,1)
plot(DTdate(361:720,1), FBrunoff_sum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DTfb_vf,'LineWidth',1.85)
hold on
plot(DTdate(361:720,1), FBrunoff_sum(1:360,1)-DTfb_vf,'LineWidth',1.7,'Color','k')
legend('HYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Foxe Basin 1979-2008')
y=ylabel('km^3/month','FontSize',14);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figFB,'-dpng','FB_runoff_comparison')

% Compare Foxe Basin Runoff with integrated HYPE
figFB=figure('name', 'Runoff comparison for Foxe Basin');
%subplot(2,1,1)
plot(DTdate(361:720,1), FBrunoff_intsum(1:360,1),'LineWidth',2)
hold on
plot(DTdate(361:720,1), DTfb_vf,'LineWidth',1.85)
hold on
plot(DTdate(361:720,1), FBrunoff_intsum(1:360,1)-DTfb_vf,'LineWidth',1.7,'Color','k')
legend('intHYPE','Dai Trenberth','HYPE-DT')
title('Monthly river runoff for Foxe Basin 1979-2008')
y=ylabel('km^3/month','FontSize',14);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); % moves label to the left
%xlabel('Date (month/year)','FontSize',14)
xlim([DTdate(361,1) DTdate(720,1)])
tick_locations = datenum(1979,1:12:384,1);
set(gca,'XTick',tick_locations)
datetick('x','mm/yy','keeplimits','keepticks')
grid on
rotateXLabels( gca(), 45 )
print(figFB,'-dpng','FB_intrunoff_comparison')

close all 

%% Calculate seasonal cycle of each basin and get the annual mean river runoff

% variables I work with are int, insum and vf and DTdate

% get mean HBC area
%DThbc_vf=squeeze(mean(DTrunoff_vf,2)); % time, lat, lon
%DThbc_vf=squeeze(mean(DThbc_vf,2)); % have spatial averages in time
%DThbc_vf=DThbc_vf(361:720,1);% now just haev values from 1979 to 2008

HBCrunoff=squeeze(mean(HYPEdischarge,2));
HBCrunoff=HBCrunoff(1:360,1);

HBCrunoff_int=squeeze(mean(HYPEdischarge_int,2));
HBCrunoff_int=HBCrunoff_int(1:360,1);

date79_08=DTdatev(361:720,1:3);

DTfb_vfseas=zeros(12,1);
for m=1:12
    moni=find(date79_08(:,2)==m);
    % Foxe Basin
    DTfb_vfseas(m,1)=mean(DTfb_vf(moni)); % units are km^3 per month
    FBrunoff_sumseas(m,1)=mean(FBrunoff_sum(moni));
    FBrunoff_intsumseas(m,1)=mean(FBrunoff_intsum(moni));
    % Hudson Bay East
    DThbE_vfseas(m,1)=mean(DThbE_vf(moni));
    HBErunoff_intsumseas(m,1)=mean(HBErunoff_intsum(moni));
    HBErunoff_sumseas(m,1)=mean(HBErunoff_sum(moni));
    % Hudson Bay west
    DThbW_vfseas(m,1)=mean(DThbW_vf(moni));
    HBWrunoff_intsumseas(m,1)=mean(HBWrunoff_intsum(moni));
    HBWrunoff_sumseas(m,1)=mean(HBWrunoff_sum(moni));
    % Hudson bay whole    
    DThb_vfseas(m,1)=mean(DThb_vf(moni));
    HBrunoff_sumseas(m,1)=mean(HBrunoff_sum(moni));
    HBrunoff_intsumseas(m,1)=mean(HBrunoff_intsum(moni));
    % Hudson Strait    
    DThs_vfseas(m,1)=mean(DThs_vf(moni));
    HSrunoff_intsumseas(m,1)=mean(HSrunoff_intsum(moni));
    HSrunoff_sumseas(m,1)=mean(HSrunoff_sum(moni));
    % Hudson Strait and Ungava bay    
    DThsub_vfseas(m,1)=mean(DThsub_vf(moni));
    HSUBrunoff_intsumseas(m,1)=mean(HSUBrunoff_intsum(moni));
    HSUBrunoff_sumseas(m,1)=mean(HSUBrunoff_sum(moni));
    % Ungava Bay
    DTrunoff_vfseas(m,1)=mean(DTub_vf(moni));
    UBrunoff_intsumseas(m,1)=mean(UBrunoff_intsum(moni));
    UBrunoff_sumseas(m,1)=mean(UBrunoff_sum(moni));
    % James Bay
    DTjb_vfseas(m,1)=mean(DTjb_vf(moni));
    JBrunoff_intsumseas(m,1)=mean(JBrunoff_intsum(moni));
    JBrunoff_sumseas(m,1)=mean(JBrunoff_sum(moni));
    % HBC   
    DThbc_vfseas(m,1)=mean(DThbc_vf(moni)); 
    HBCrunoff_seas(m,1)=mean(HBCrunoff(moni)); %time and river
    HBCrunoff_intseas(m,1)=mean(HBCrunoff_int(moni)); % isthis in km^3???
    
    
end

colourmatrix=[0.93 0.69 0.126; 0.49 0.18 0.56; 0.47 0.67 0.188]; % yellow, purple, green
fbmeans(1:12,1)=mean(DTfb_vfseas);
fbmeans(1:12,2)=mean(FBrunoff_sumseas);
fbmeans(1:12,3)=mean(FBrunoff_intsumseas);

% Foxe Basin
fig1=figure
plot(1:12, DTfb_vfseas,'LineWidth',1.75,'Color',colourmatrix(1,1:3))
hold on
plot(1:12, fbmeans(:,1),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(1,1:3))
hold on
plot(1:12, FBrunoff_sumseas,'LineWidth',2.5,'Color',colourmatrix(2,1:3))
hold on
plot(1:12, fbmeans(:,2),'LineWidth',2.5,'LineStyle','--','Color',colourmatrix(2,1:3))
hold on 
plot(1:12, FBrunoff_intsumseas,'LineWidth',1.75,'Color',colourmatrix(3,1:3))
hold on
plot(1:12, fbmeans(:,3),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(3,1:3))
legend('DaiTren','mean' ,'HYPE','mean','HYPE int','mean','Location','Northwest')
title('Seasonal river runoff for Foxe Basin 1979-2008')
ylabel('km^3','FontSize',14);
xlabel('Month','FontSize',14)
xlim([1 12])
set(gca,'XTick',1:1:12,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
grid on
rotateXLabels( gca(), 45 )
print(fig1,'-dpng','FB_3runoff_comparison')

AnnualRunoffFB_DT=sum(DTfb_vfseas)
AnnualRunoffFB_HYPE=sum(FBrunoff_sumseas)
AnnualRunoffFB_HYPEint=sum(FBrunoff_intsumseas)

hsubmeans(1:12,1)=mean(DThsub_vfseas);
hsubmeans(1:12,2)=mean(HSUBrunoff_sumseas);
hsubmeans(1:12,3)=mean(HSUBrunoff_intsumseas);
% Hudson Strait and Ungava Bay
fig1=figure
plot(1:12, DThsub_vfseas,'LineWidth',1.75,'Color',colourmatrix(1,1:3))
hold on
plot(1:12, hsubmeans(:,1),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(1,1:3))
hold on
plot(1:12, HSUBrunoff_sumseas,'LineWidth',1.75,'Color',colourmatrix(2,1:3))
hold on
plot(1:12, hsubmeans(:,2),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(2,1:3))
hold on 
plot(1:12, HSUBrunoff_intsumseas,'LineWidth',1.75,'Color',colourmatrix(3,1:3))
hold on
plot(1:12, hsubmeans(:,3),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(3,1:3))
legend('DaiTren','mean' ,'HYPE','mean','HYPE int','mean','Location','Northwest')
title('Seasonal river runoff for HSUB 1979-2008')
ylabel('km^3','FontSize',14);
xlabel('Month','FontSize',14)
xlim([1 12])
set(gca,'XTick',1:1:12,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
grid on
rotateXLabels( gca(), 45 )
print(fig1,'-dpng','HSUB_3runoff_comparison')

AnnualRunoffHSUB_DT=sum(DThsub_vfseas)
AnnualRunoffHSUB_HYPE=sum(HSUBrunoff_sumseas)
AnnualRunoffHSUB_HYPEint=sum(HSUBrunoff_intsumseas)


hbmeans(1:12,1)=mean(DThb_vfseas);
hbmeans(1:12,2)=mean(HBrunoff_sumseas);
hbmeans(1:12,3)=mean(HBrunoff_intsumseas);
% Hudson Bay 
fig1=figure
plot(1:12, DThb_vfseas,'LineWidth',1.75,'Color',colourmatrix(1,1:3))
hold on
plot(1:12, hbmeans(:,1),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(1,1:3))
hold on
plot(1:12, HBrunoff_sumseas,'LineWidth',1.75,'Color',colourmatrix(2,1:3))
hold on
plot(1:12, hbmeans(:,2),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(2,1:3))
hold on 
plot(1:12, HBrunoff_intsumseas,'LineWidth',1.75,'Color',colourmatrix(3,1:3))
hold on
plot(1:12, hbmeans(:,3),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(3,1:3))
legend('DaiTren','mean' ,'HYPE','mean','HYPE int','mean','Location','Northwest')
title('Seasonal river runoff for Hudson Bay 1979-2008')
ylabel('km^3','FontSize',14);
xlabel('Month','FontSize',14)
xlim([1 12])
set(gca,'XTick',1:1:12,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
grid on
rotateXLabels( gca(), 45 )
print(fig1,'-dpng','HB_3runoff_comparison')

AnnualRunoffHB_DT=sum(DThb_vfseas)
AnnualRunoffHB_HYPE=sum(HBrunoff_sumseas)
AnnualRunoffHB_HYPEint=sum(HBrunoff_intsumseas)


jbmeans(1:12,1)=mean(DTjb_vfseas);
jbmeans(1:12,2)=mean(JBrunoff_sumseas);
jbmeans(1:12,3)=mean(JBrunoff_intsumseas);
% James Bay 
fig1=figure
plot(1:12, DTjb_vfseas,'LineWidth',1.75,'Color',colourmatrix(1,1:3))
hold on
plot(1:12, jbmeans(:,1),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(1,1:3))
hold on
plot(1:12, JBrunoff_sumseas,'LineWidth',1.75,'Color',colourmatrix(2,1:3))
hold on
plot(1:12, jbmeans(:,2),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(2,1:3))
hold on 
plot(1:12, JBrunoff_intsumseas,'LineWidth',1.75,'Color',colourmatrix(3,1:3))
hold on
plot(1:12, jbmeans(:,3),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(3,1:3))
legend('DaiTren','mean' ,'HYPE','mean','HYPE int','mean','Location','Northwest')
title('Seasonal river runoff for James Bay 1979-2008')
ylabel('km^3','FontSize',14);
xlabel('Month','FontSize',14)
xlim([1 12])
set(gca,'XTick',1:1:12,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
grid on
rotateXLabels( gca(), 45 )
print(fig1,'-dpng','JB_3runoff_comparison')

AnnualRunoffJB_DT=sum(DTjb_vfseas)
AnnualRunoffJB_HYPE=sum(JBrunoff_sumseas)
AnnualRunoffJB_HYPEint=sum(JBrunoff_intsumseas)

hbcmeans(1:12,1)=mean(DThbc_vfseas);
hbcmeans(1:12,2)=mean(HBCrunoff_seas);
hbcmeans(1:12,3)=mean(HBCrunoff_intseas);
% HBC 
fig1=figure
plot(1:12, DThbc_vfseas,'LineWidth',1.75,'Color',colourmatrix(1,1:3))
hold on
plot(1:12, hbcmeans(:,1),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(1,1:3))
hold on
plot(1:12, HBCrunoff_seas,'LineWidth',1.75,'Color',colourmatrix(2,1:3))
hold on
plot(1:12, hbcmeans(:,2),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(2,1:3))
hold on 
plot(1:12, HBCrunoff_intseas,'LineWidth',1.75,'Color',colourmatrix(3,1:3))
hold on
plot(1:12, hbcmeans(:,3),'LineWidth',1.75,'LineStyle','--','Color',colourmatrix(3,1:3))
legend('DaiTren','mean' ,'HYPE','mean','HYPE int','mean','Location','Northwest')
title('Seasonal river runoff for HBC 1979-2008')
ylabel('km^3','FontSize',14);
xlabel('Month','FontSize',14)
xlim([1 12])
set(gca,'XTick',1:1:12,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
grid on
rotateXLabels( gca(), 45 )
print(fig1,'-dpng','HBC_3runoff_comparison')

AnnualRunoffHBC_DT=sum(DThbc_vfseas)
AnnualRunoffHBC_HYPE=sum(HBCrunoff_seas)
AnnualRunoffHBC_HYPEint=sum(HBCrunoff_intseas)















