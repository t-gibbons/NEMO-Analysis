function CreateGeoInfoMatrix( baysys_code, baysys_conf, baysys_vers, baysys_domn, YS, YE )
% Create geoinfo file
% baysys code indicates the codification fo the .txt files (i.e., MR3-R41 for MRI-CGCM3 RCP4.5)

runoffData = importdata(['NEMO_HYPE_Results_v',baysys_vers,'_',baysys_domn,'_',baysys_conf,'_',num2str(YS),'_',num2str(YE),'_',baysys_code,'.txt']); % In this format, data is input as a structure, data units are m3/s
tempData = importdata(['NEMO_HYPE_Results_v',baysys_vers,'_',baysys_domn,'_',baysys_conf,'_temp_',num2str(YS),'_',num2str(YE),'_',baysys_code,'.txt']);
saveFile   = ['GeoInfo_',baysys_conf,'_',baysys_code,'_',baysys_domn];

riverID        = runoffData.data(1,:);
riverLat       = runoffData.data(2,:);
riverLon       = runoffData.data(3,:);
riverDischarge = runoffData.data(4:end,:);
riverTemp      = tempData.data(4:end,:);

% need mean annual discharge in m3/s
riverDischarge = nanmean(riverDischarge,1);
riverTemp      = nanmean(riverTemp,1);

% combine ID, lat, lon, annual river discharge for each region into one
% file
runoffFile(1,:) = riverID;
runoffFile(2,:) = riverLat;
runoffFile(3,:) = riverLon;
runoffFile(4,:) = riverDischarge;
runoffFile(5,:) = riverTemp;

if ~exist(saveFile,'file')
    save(saveFile,'runoffFile')
else
    disp('file already exists')
end


