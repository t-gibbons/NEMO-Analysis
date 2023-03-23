function createFinalCombinedRunoff_IA( CF, YS, YE, geoinf_code, baysys_code, baysys_conf )
% combine inter-annual Dai & Trenberth's data and Bamber's data
% usage:
%      createFinalCombinedRunoff_IA(CF,YS,YE)
%          CF: 'ANHA4' or 'ANHA12'
%          YS: start year (Bamber's data starts from 1958)
%          YE:   end year (Bamber's data ends by 2010 but Dai & Trenberth's ends in 2007) 
%
%  note: 
%      the ratio (totalFW2RunoffRatio) used to estimate the runoff into the ocean from total freshwater (Bamber's data)
%      is set to 0.46.

% history:
%          2015-01: xianmin@ualberta.ca
%          2015-12: ste CF, YS, YE as function inputs

% if nargin==2
%    YE=YS;
% elseif nargin==1
%    if ischar(CF)
%       YS=1958; YE=2010;
%    elseif isnumeric(CF)
%       YS=CF; YE=YS; CF='ANHA4';
%    else
%       help createFinalCombinedRunoff_IA
%       return
%    end
% elseif nargin~=3
%    help createFinalCombinedRunoff_IA
%    return
% end

rnfDTpath='/mnt/storage5/natasha/FINAL_HYPE/Naturalized/'; % original ANHA4 runoff path
rnfGL=['Monthly_All_Data_nc_Files/ANHA4_',geoinf_code,'_',baysys_conf,'_runoff_monthly_y',num2str(YS),'_y',num2str(YE),'.nc'];
GLmask=load('Support_Files/HBCRegion_ANHA4.mat');                               % a mask file of HBC region ** CHANGE

[NY,NX]=size(GLmask.tmask);

totalFW2RunoffRatio=1; % ratio used to estimate the runoff into the ocean from total freshwater (Bamber's data)
%ratioStr=strrep(num2str(totalFW2RunoffRatio),'.','p');         % used to create the output file name

varList={'nav_lon','nav_lat','time_counter','socoefr','runoff'};
GLclim=zeros(NY,NX);
for ny=YS:YE
    tTag=['y',num2str(ny)]; % used to create the output file name
    %if ny<=2010
        rnfCLIM=[rnfDTpath,CF,'_',geoinf_code,'_monthly_final_',tTag,'.nc'];
   % else
   %     rnfCLIM=[rnfDTpath,CF,'_runoff_monthly_combined_Dai_Trenberth_Bamber_y2010m00_r1.nc'];
   % end

    %% create the output file
    modLon=GetNcVar(rnfCLIM,'nav_lon');
    modLat=GetNcVar(rnfCLIM,'nav_lat');
    %saveFileName=['IA/',CF,'_runoff_monthly_combined_Dai_Trenberth_Bamber_HYPEcal_MI5_',tTag,'m00_r1.nc']; % ** CHANGE
    saveFileName=['Yearly_nc_Files_FINAL_INPUT/',baysys_code,'/',CF,'_',baysys_conf,'_',geoinf_code,'_monthly_final_',tTag,'.nc'];
    [ncfid,varIDList]=createNC(saveFileName,NY,NX,varList,['HYPE ',baysys_conf,' runoff for ',geoinf_code,' BaySys run for ',num2str(YS),' to ',num2str(YE),' (',date,')']);
    disp(['Saved File - ',saveFileName]);
    netcdf.putVar(ncfid,varIDList(1),(modLon)');
    netcdf.putVar(ncfid,varIDList(2),(modLat)');
    netcdf.putVar(ncfid,varIDList(3),0,12,GetNcVar(rnfCLIM,'time_counter'));

    socoefrCLIM=GetNcVar(rnfCLIM,'socoefr');
    socoefrGL=GLmask.tmask*0;

    for nmon=1:12
        if ny==YS % change 1979 to YS
           % take it from 1958
           indC=nmon-1;
        else
           indC=(ny-YS)*12+nmon-1;
        end
        GLclim(:,:)=GetNcVar(rnfGL,'runoff',[0 0 indC],[NX NY 1])*totalFW2RunoffRatio;
        RunoffCLIM=GetNcVar(rnfCLIM,'runoff',[0 0 nmon-1],[NX NY 1]);
        RunoffCLIM(GLmask.tmask==1)=GLclim(GLmask.tmask==1);
        RunoffCLIM(RunoffCLIM<0)=0;
        netcdf.putVar(ncfid,varIDList(5),[0 0 nmon-1],[NX NY 1],RunoffCLIM');
        socoefrGL(GLclim>1e-7)=0.5;
    end % loop of month
    socoefrCLIM(GLmask.tmask==1)=socoefrGL(GLmask.tmask==1);
    netcdf.putVar(ncfid,varIDList(4),socoefrCLIM');
    netcdf.close(ncfid);
end % loop of year

function [ncfid,varIDList]=createNC(ncfileName,ny,nx,varList,srcTag)
if exist(ncfileName,'file')
    eval(['!rm ',ncfileName])
end
ncfid=netcdf.create(ncfileName,'NETCDF4');
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
netcdf.defVarDeflate(ncfid,varIDList(1),true,true,9);
netcdf.defVarDeflate(ncfid,varIDList(2),true,true,9);
netcdf.defVarDeflate(ncfid,varIDList(4),true,true,9);
netcdf.defVarDeflate(ncfid,varIDList(5),true,true,9);

netcdf.putAtt(ncfid,varIDList(3),'units','days since Jan-1-0000 00:00:00');
netcdf.putAtt(ncfid,varIDList(3),'calendar','noleap');
netcdf.putAtt(ncfid,varIDList(5),'units','kg/s/m^2');
netcdf.putAtt(ncfid,netcdf.getConstant('NC_GLOBAL'),'source',srcTag);
netcdf.endDef(ncfid);
