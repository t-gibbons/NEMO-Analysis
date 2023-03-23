function CreateMonthlyMeanMatrix( baysys_code, baysys_conf, baysys_vers, baysys_domn, YS, YE )

runoffData = importdata(['NEMO_HYPE_Results_v',baysys_vers,'_',baysys_domn,'_',baysys_conf,'_',num2str(YS),'_',num2str(YE),'_',baysys_code,'.txt']); % In this format, data is input as a structure, data units are degrees celcius
tempData = importdata(['NEMO_HYPE_Results_v',baysys_vers,'_',baysys_domn,'_',baysys_conf,'_temp_',num2str(YS),'_',num2str(YE),'_',baysys_code,'.txt']);
saveFile   = ['MonthlyDischarge_',baysys_conf,'_',baysys_code,'_',baysys_domn];
saveFileTemp = ['MonthlyTemp_',baysys_conf,'_',baysys_code,'_',baysys_domn];

dates       = tempData.textdata(4:end);
dates       = datevec(dates);
dailyRunoff = runoffData.data(4:end,:);
dailyTemp   = tempData.data(4:end,:);

monthlyRunoff = zeros(((dates(end,1)-dates(1,1))+1)*12, size(dailyRunoff,2)+2 );
monthlyTemp   = zeros(((dates(end,1)-dates(1,1))+1)*12, size(dailyRunoff,2)+2 );
i = 1;
for yr = dates(1,1):dates(end,1) %dates(end,1)
    for mon = 1:12
        ind = find(dates(:,1) == yr & dates(:,2) == mon);
        
        monthlyRunoff(i,3:end) = nanmean(dailyRunoff(ind,:));
        monthlyRunoff(i,1)     = yr;
        monthlyRunoff(i,2)     = mon;
        
        monthlyTemp(i,3:end) = nanmean(dailyTemp(ind,:));
        monthlyTemp(i,1)     = yr;
        monthlyTemp(i,2)     = mon;
        
        i = i+1;
    end
end

riverID  = [0,0,runoffData.data(1,:)];
riverLat = [0,0,runoffData.data(2,:)];
riverLon = [0,0,runoffData.data(3,:)];

monthlyRunoff = vertcat(riverID,riverLat,riverLon,monthlyRunoff);
monthlyTemp = vertcat(riverID,riverLat,riverLon,monthlyTemp);

if ~exist(saveFile,'file')
    save(saveFile,'monthlyRunoff')
else
    disp('file already exists')
end

if ~exist(saveFileTemp,'file')
    save(saveFileTemp,'monthlyTemp')
else
    disp('file already exists')
end