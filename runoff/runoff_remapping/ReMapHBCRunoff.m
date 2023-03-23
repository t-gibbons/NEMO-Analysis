function ReMapHBCRunoff( baysys_code, baysys_conf, baysys_domn, YS, YE )
% based on a pre-created coast-buffer zone, 
% remap Bamber's 5x5 km resolution data to ANHA model grid
% history:
%          2014-09: first coded for Dai & Trenberth 1x1 data,  xianmin@ualberta.ca
%          2014-xx: modified for Greenland Runoff data, xianmin@ualberta.ca
%          2015-07: added detailed comments, xianmin@ualberta.ca
%          2016-04: modified for HBC case

%% STEP 1: 
%  load the pre-defined coastal-buffer-polygons 
%  coast_buffer is a cell array contains information of all coastal-buffer-polygons
%  each element of coast_buffer, i.e., coast_buffer{n}, is a structure variable that contains 
%  the longitude and latitude of every vertex of current polygon, coast_buffer{n}.lon and coast_buffer{n}.lat
load('Support_Files/HBCCoast_buffer.mat','coast_buffer');

%% STEP 2:
%  find which coast-buffer-zone does each of your individual river (runoff data point) belong to
%  this step is better to include all possible positions of your rivers (even no data some period)
%  because only the river locations are needed in this step
rnfFile0 = importdata(['GeoInfo_',baysys_conf,'_',baysys_code,'_',baysys_domn,'.mat']);% ** CHANGE
%Arcticind=find(rnfFile0(5,:)==0);
%rnfFile0(:,Arcticind)=[];

rnfInfo.runoff       = rnfFile0(4,:); % annual runoff for each river
rnfInfo.temp         = rnfFile0(5,:);   % annual water temp for each river, do I need this??
rnfInfo.lon          = rnfFile0(3,:);
rnfInfo.lat          = rnfFile0(2,:);
rnfInfo.coast_buffer = coast_buffer;
rnfInfo.savefile     = 'coast_buffer_HBC_runoff.mat';
rnfInfo.isGridded    = 0;

if exist(rnfInfo.savefile,'file')
   load(rnfInfo.savefile)
else
   myRunoffInfo = labelRunoff(rnfInfo);
   %eval(['save ',rnfInfo.savefile,' myRunoffInfo'])
end

%% STEP 3:
%  create the mask file for coast-buffer-zone under current model configuration
%  It is simply based on the geographic location of each coast-buffer-zone (longitudes, latitudes)
modelInfo.CF           = 'ANHA4';
maskfile               = 'Support_Files/ANHA4_mask.nc'; % model mask files
modelInfo.coast_buffer = coast_buffer;

NX = 544; %x-dimension
NY = 800; %y-dimension

modelInfo.lsmask = GetNcVar(maskfile,'tmask',[0 0 0 0],[NX NY 1 1]);
modelInfo.lon    = GetNcVar(maskfile,'nav_lon');
modelInfo.lat    = GetNcVar(maskfile,'nav_lat');

if exist(['Support_Files/HBC_buffer_mask_',modelInfo.CF,'_stp1.mat'],'file')
   load(['Support_Files/HBC_buffer_mask_',modelInfo.CF,'_stp1.mat'])
   modelInfo.bufferMask  = bufferMask;
   modelInfo.isFromCoast = [];
   % 1: distribute runoff from "model" coast (most cases, yes)
   % 0: distribute runoff at the closest (to data source)  water points in model 
else
   modelInfo.savefile   = ['Support_Files/HBC_buffer_mask_',modelInfo.CF,'.mat'];
   modelInfo.bufferMask = getCoastBufferModelMask(modelInfo);
   disp(['Edit the bufferMask if necessary and save as HBC_buffer_mask_',modelInfo.CF,'_stp1.mat']);
   % Need modification because runoff is accumulated in some long channel-like bays, which keep freshening the ocean
   % and eventually blow up the model due to negative salinity
   return
end

numP           = numel(coast_buffer);
modelInfo.nWPT = zeros(numP,1);

for np = 1:numP
    modelInfo.nWPT(np) = numel(find(modelInfo.bufferMask == np & modelInfo.lsmask == 1));
end
myRunoffInfo.extraPolygonExist = 1; % contains buffer-zones not in model domain (some rivers do not flow into model region)

%% STEP 4:
%  get model mesh information (area of each model cell, which does not vary in time) for later interpolation
meshhfile          = '/mnt/storage1/xhu/ANHA4-I/ANHA4_mesh_hgr.nc';     % model horizontal mesh file, which contains model mesh size
e1t                = GetNcVar(meshhfile,'e1t');
e2t                = GetNcVar(meshhfile,'e2t');
modelInfo.e1e2tInv = 1./(e1t.*e2t);
clear('e1t','e2t');

%% STEP 5: 
% create the output file and get information of river source data
% YS=1981;
% YE=2070; % ** CHANGE
mm=repmat((1:12)',1,YE-YS+1); 
yy=repmat(YS:YE,12,1);
timeCounter=datenum(yy(:),mm(:),15);
%timeCounter=1:12;

varList           = {'nav_lon','nav_lat','time_counter','socoefr','runoff'};
[ncfid,varIDList] = createNC([modelInfo.CF,'_',baysys_conf,'_',baysys_code,'_',baysys_domn,'_runoff_monthly_y',num2str(YS),'_y',num2str(YE),'.nc'],NY,NX,varList,['re-mapped from HYPE ',baysys_code,' data']);% ** CHANGE

netcdf.putVar(ncfid,varIDList(1),(modelInfo.lon)');
netcdf.putVar(ncfid,varIDList(2),(modelInfo.lat)');

netcdf.putVar(ncfid,varIDList(3),0,numel(timeCounter),timeCounter);
%runoffMonthlyFile=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HYPEcal_HBC/monthlyQ-20170410_WFDEI.csv'); %
%runoffMonthlyFile=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HYPEcal_HBC/monthlyQ-20170410_WFDEI.mat'); % ** CHANGE
runoffMonthlyFile = importdata(['MonthlyDischarge_',baysys_conf,'_',baysys_code,'_',baysys_domn,'.mat']);% ** CHANGE
runoffMonthlyFile = runoffMonthlyFile(:,3:end);
runoffMonthlyFile = runoffMonthlyFile(4:end,:);

tempMonthlyFile = importdata(['MonthlyTemp_',baysys_conf,'_',baysys_code,'_',baysys_domn,'.mat']);
tempMonthlyFile=tempMonthlyFile(:,3:end);
tempMonthlyFile=tempMonthlyFile(4:end,:);

NXGL    = 397; %x-dimesnion in runoffMonthlyFile -- number of discharge locations
NYGL    = size(runoffMonthlyFile,1); %y-dimension in runoffMonthlyFile -- months from 1979 to 2010 ** CHANGE
socoefr = modelInfo.lsmask*0;

%% STEP 6:
%  re-mapping source data to model grid
%  (assume: the larger runoff should be distributed onto more model cells (to avoid extreme values for big rivers)
%  classify simply based on the total volume (km^3) in each buffer-zone (see details in the loop below)
%  This is not necessary if the runoff is not very big (just set a constant value for nClosestPT, again see the loop below)

volLevPTs=[ 5000  20000 20; ...
            4000  5000 12; ...
            3000  4000  7; ...
           0.0   2000 2];
       
KM3MonToMS        = 1e9/(24*3600*(365/12));  % factor to convert km^3/month to m^3/s --- not using
myRunoffInfo.avgN = 25;               % max number of cells used to distribute the runoff for one "river"

for nmon = 1:numel(timeCounter)
    nClosestPT = zeros(size(myRunoffInfo.isInsideDomain))+1; % myrunoffInfo is the annual mean geo info spreadsheet
    % read original runoff -- this would be my large HYPE csv file?
    % and I think I want to read 1 time step at a time (so read each row)
    cRunoff    = runoffMonthlyFile(nmon,:); % use nmon???? not sure about this. Data is m/s. cRunoff has runoff, lat and lon
    cTemp      = tempMonthlyFile(nmon,:);
    
    % what information do I need for this loop
   %*******
    % If 
    totalRunoffInPolygon = getHBCRunoffInEachPolygon(cRunoff,myRunoffInfo);  % m3/s ** CHANGE IN FUNCTION
    for nlev = 1:size(volLevPTs,1)
        indCP = find(totalRunoffInPolygon >= volLevPTs(nlev,1) & totalRunoffInPolygon < volLevPTs(nlev,2));
        if ~isempty(indCP)
            for np=1:numel(indCP)
                nClosestPT(myRunoffInfo.isInsideDomain==indCP(np))=volLevPTs(nlev,3);
            end
        end
    end

%%
    %[modelRunoff,myRunoffInfo]=computeModelRunoff(cRunoff,myRunoffInfo,modelInfo,nClosestPT);
    [modelRunoff,modelTemp,myRunoffInfo]=computeModelRunoff(cRunoff,cTemp,myRunoffInfo,modelInfo,nClosestPT);
    modelRunoff(modelRunoff<0)=0;  % No negative runoff here
    modelTemp(modelTemp<0)=0;  % No negative runoff here
    netcdf.putVar(ncfid,varIDList(5),[0 0 nmon-1],[NX NY 1],modelTemp');
    socoefr(modelRunoff>1e-7)=0.5; % where NEMO will do enhanced mixing because the existence of runoff
end
netcdf.putVar(ncfid,varIDList(4),socoefr');
netcdf.close(ncfid);

%-------------------------------------------------------------------------
function [ncfid,varIDList]=createNC(ncfileName,ny,nx,varList,srcTag)
ncfid=netcdf.create(ncfileName,'NC_CLOBBER');
% dimension
mydimID(1,1)=netcdf.defDim(ncfid,'x',nx);
mydimID(1,2)=netcdf.defDim(ncfid,'y',ny);
mydimID(1,3)=netcdf.defDim(ncfid,'time_counter',netcdf.getConstant('NC_UNLIMITED'));

%variable
varIDList=zeros(numel(varList),1);
varIDList(1)=netcdf.defVar(ncfid,'nav_lon','float',mydimID(1,1:2));
varIDList(2)=netcdf.defVar(ncfid,'nav_lat','float',mydimID(1,1:2));
varIDList(3)=netcdf.defVar(ncfid,'time_counter','double',mydimID(3));
varIDList(4)=netcdf.defVar(ncfid,varList{4},'float',mydimID(1,1:2));
varIDList(5)=netcdf.defVar(ncfid,varList{5},'float',mydimID);
netcdf.putAtt(ncfid,netcdf.getConstant('NC_GLOBAL'),'source',srcTag);
netcdf.endDef(ncfid);
