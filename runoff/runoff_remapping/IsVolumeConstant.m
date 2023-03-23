function IsVolumeConstant( geoinf_code )
% look at ratio between pre and post remapping --  specifically for
% remapping that was done when I was away
% clear all
close all

%opengl software % for matlab2014

solid = 0;
YS    = 1981;
YE    = 2070;

if solid == 0
RemappedHBC = GetNcVar(['Monthly_All_Data_nc_Files/ANHA4_',geoinf_code,'_regulated_runoff_monthly_y' num2str(YS) '_y' num2str(YE) '.nc'],'runoff'); % kg/m2/s ** CHANGE
else
Remapped = GetNcVar(['ANHA4_RCP85_RACMO_solidFWrunoff_monthly_y' num2str(YS) '_y' num2str(YE) '.nc'],'runoff'); % 1999 - 1788-1800
end

directory = '/mnt/storage2/xhu/NEMO/ANHA4-EXH001/'; % ** All ANHA4 runs use this directory for the meshfiles and stuff

meshfile=[directory 'mask.nc'];
meshfileh=[directory 'mesh_hgr.nc'];
meshfilez=[directory 'mesh_zgr.nc'];

e1t = GetNcVar(meshfileh,'e1t'); % Dy
e2t = GetNcVar(meshfileh,'e2t'); % Dx

gridArea=e1t.*e2t;

% for i=1:size(RemappedA,1)
%     RemappedA(i,:,:)=squeeze(RemappedA(i,:,:)).*gridArea; % kg/s
% end
% 
% test=RemappedA(5,:,:);
% SumRemappedA = nansum(squeeze(nansum(RemappedA,2)),2);

for i=1:size(RemappedHBC,1)
    RemappedHBC(i,:,:)=squeeze(RemappedHBC(i,:,:)).*gridArea; % kg/s
end

test=RemappedHBC(5,:,:);
SumRemappedHBC = nansum(squeeze(nansum(RemappedHBC,2)),2);

SumRemapped=(SumRemappedHBC)/1000; % m3/s


% load in pre remapped runoff
load(['Monthly_and_GeoInfo_mat_Files/MonthlyDischarge_REG_',geoinf_code,'_HBC.mat'],'monthlyRunoff'); % m3/s

%Arcticind=find(resultmatrix(3,:)==0);
%resultmatrix(:,Arcticind)=[];
monthlyRunoff=monthlyRunoff(4:end,3:end);

SumOrig=sum(monthlyRunoff,2);

figure;
plot(SumOrig,'LineWidth',2);hold on; plot(SumRemapped);legend('orig','remapped')
print(['VolumeCheck1_',geoinf_code],'-djpeg');

ratio=SumRemapped./SumOrig;
figure;
plot(ratio)
print(['VolumeCheck2_',geoinf_code],'-djpeg');

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


runoffFile='/mnt/storage4/laura/RESEARCH/REMAP_HYPE/ResultDateLatLonMatrix_noleap_monthly_nodates_inDomain_lc.mat';
rnfFile0=GetNcVar(runoffFile,'FW_liquid'); rnfFile0(rnfFile0==0)=nan; rnfFile0(rnfFile0>1000000)=nan;
rnfInfo.runoff=squeeze(nanmean(rnfFile0,1)); % annual runoff for each river
lon=GetNcVar(runoffFile,'lon');
if max(max(lon))>180
    temp=find(lon>180);
    lon(temp)=lon(temp)-360;
    rnfInfo.lon=GetNcVar(runoffFile,'lon');
else
    rnfInfo.lon=GetNcVar(runoffFile,'lon');
end
lat=GetNcVar(runoffFile,'lat');


for l=1:size(lon,2)
    RACMO_mesh_y(:,l)=m_lldist(lon(:,l),lat(:,l));
    % RACMO_mesh_yT(:,l)=myll2dist(lon(:,l),lat(:,l));
end
RACMO_mesh_y(384,:)=RACMO_mesh_y(1,:);
for l=1:size(lon,1)
    RACMO_mesh_x(l,:)=m_lldist(lon(l,:),lat(l,:));
    % RACMO_mesh_xT(:,l)=myll2dist(lon(:,l),lat(:,l));
end
RACMO_mesh_x(:,320)=RACMO_mesh_x(:,1);
RACMOgridArea=RACMO_mesh_x.*RACMO_mesh_y*1000*1000;% m2


i=1;
for y=YS:YE
   % if YS==2006
        temp=GetNcVar(['/mnt/storage3/natasha/RUNOFF/RACMO/RCP2.6/FW_GrIS_AIS_' num2str(y) '_RCP2.6.nc'],'FW_liquid');
   % else
   %     temp=GetNcVar(['/mnt/storage3/natasha/RUNOFF/RACMO/hist/FW_GrIS_AIS_' num2str(y) '.nc'],'FW_liquid');
   % end
    
    for n=1:12
        OriginalRACMO(n,:,:)=squeeze(temp(n,:,:)).*RACMOgridArea;
        if solid==1
           % Remove Antarctic 
           OriginalRACMO(n,1:187,:)=nan;
        end
    end
    SumOrigRACMO(:,i)=nansum(squeeze(nansum(OriginalRACMO,2)),2);
    i=i+1;
end

SumOrigRACMO=reshape(SumOrigRACMO,size(SumOrigRACMO,1)*size(SumOrigRACMO,2),1);

figure
plot(SumOrigRACMO,'LineWidth',2)
hold on
plot(SumRemappedRACMO)
legend('orig','remapped')

factor2=SumRemappedRACMO./SumOrigRACMO;
max(factor2)
min(factor2)

figure
plot(factor2)


%% Check combined files in IA directory
GLmask=GetNcVar('ANHA4_RCP26_RACMO_liquidFWrunoff_monthly_y2006_y2017.nc','socoefr');
GLmask(GLmask>0)=1;

i=1;
for y=YS:YE
    if solid==0
        temp=GetNcVar(['/mnt/storage3/natasha/RUNOFF/RACMO/Remap2ANHAgrid/IA/ANHA4_runoff_monthly_combined_Dai_Trenberth_RCP26_RACMO_y' num2str(y) '.nc'],'runoff');
    else
        temp=GetNcVar(['/mnt/storage3/natasha/RUNOFF/RACMO/Remap2ANHAgrid/IA/ANHA4_solidrunoff_monthly_RCP26_RACMO_y' num2str(y) '.nc'],'runoff');
    end
    for n=1:12
        CombinedRACMO(n,:,:)=squeeze(temp(n,:,:)).*GLmask.*gridArea;
    end
    SumCombRACMO(:,i)=nansum(squeeze(nansum(CombinedRACMO,2)),2);
    i=i+1;
end

%
% GLmask(GLmask==0)=nan;
% test=squeeze(nansum(CombinedRACMO,1));
% test(test==0)=nan;
% 
% temp=CombinedRACMO(8,:,:);
% temp(temp>0)=2;
% temp(temp==0)=nan;
% figure; 
% mypcolor(temp)%2
% hold on
% mypcolor(GLmask) % 1




SumCombRACMO=reshape(SumCombRACMO,size(SumCombRACMO,1)*size(SumCombRACMO,2),1);

figure
plot(SumCombRACMO,'LineWidth',2)
hold on
plot(SumRemappedRACMO)
legend('merged','not merged')

rat=SumCombRACMO./SumRemappedRACMO;
figure
plot(rat)

%% 
if ~isempty(find(rat<0.9))
    
    GLmask=GetNcVar('GLmask.nc','tmask');   
    GLmask(GLmask==0)=nan;
    
    RemappedRACMOmask=GetNcVar('ANHA4_RACMO_liquidFWrunoff_monthly_y1958_y2005.nc','socoefr');
    RemappedRACMOmask(RemappedRACMOmask==0)=nan;
    
    figure
    mypcolor(RemappedRACMOmask)
    hold on
    mypcolor(GLmask)
    
   
    
    
end
