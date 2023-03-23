function updateRunoffCOEF(isRemove)
% add more area with enhanced mixing for runoff
% update the socoefr in runoff file (0.5: enhanced mixing, 0: rest region) 
if nargin==0
   isRemove=0;
end

clc;
global NX NY tmask socoefr xxLog yyLog isRemoveGlo
xxLog=[]; yyLog=[]; isRemoveGlo=isRemove;

maskfile='ANHA4_mask.nc';
runofffile='Copy_of_HBCRiverMouthMask.nc';
NX=GetNcDimLen(maskfile,'x');
NY=GetNcDimLen(maskfile,'y');

% load surface mask  and socoefr
tmask=GetNcVar(maskfile,'tmask',[0 0 0 0],[NX NY 1 1]);
socoefr=GetNcVar(runofffile,'socoefr');
rcoef0=GetNcVar('Copy_of_ANHA4_HBC_HYPEcal_WFDEI_runoff_monthly_y1979_y2013TEST.nc','socoefr');
socoefr(rcoef0==0.5)=0.5;  % change with updated runoff

figure;
set(gcf,'render','zbuffer','color','w','CloseRequestFcn',@clearMyGlobal);
imagesc(1:NX,1:NY,socoefr); set(gca,'ydir','normal','tickdir','out');
hold on;
contour(1:NX,1:NY,tmask,[0.5 0.5],'k-')
caxis([0 1])
axis equal; axis tight;

% datacursor
mycobj=datacursormode(gcf);
datacursormode on;
set(mycobj,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on','DisplayStyle','window');

end

%---------------------------------------------------------------
function clearMyGlobal(varargin)
 % save some variables
 global xxLog yyLog NX NY tmask socoefr isRemoveGlo
 % making the mask
 bwmask=poly2mask(xxLog,yyLog,NY,NX);
 if isRemoveGlo==0
    socoefr(bwmask==1&tmask==1)=0.5;
    newRunoff='river_mouth_updated.nc';
 else
    socoefr(bwmask==1&tmask==1)=0;
    newRunoff='river_mouth_degrade.nc';
 end
 ncfid=netcdf.create(newRunoff,'NC_CLOBBER');
 % dimension
 mydimID(1,1)=netcdf.defDim(ncfid,'x',NX);
 mydimID(1,2)=netcdf.defDim(ncfid,'y',NY);

 varIDList(1)=netcdf.defVar(ncfid,'socoefr','float',mydimID);
 netcdf.putAtt(ncfid,netcdf.getConstant('NC_GLOBAL'),'history','updated to solve negative SSS in NENMO');
 netcdf.endDef(ncfid);
 netcdf.putVar(ncfid,varIDList(1),socoefr');
 netcdf.close(ncfid);
 
 % clear global variables when close the window
 clearvars -global *
 delete(gcf);
end

function output_txt = myupdatefcn(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

global preXX preYY xxLog yyLog
pos = get(event_obj,'Position');
if ~ishold, hold on; end
plot(pos(1),pos(2),'ko');  % label the previous selections

if ~isempty(preXX)
    line([preXX pos(1)],[preYY pos(2)])
end
% update previous location
preXX=pos(1); preYY=pos(2);
xxLog=[xxLog; pos(1)]; yyLog=[yyLog;pos(2)];
output_txt = {['x: ',num2str(preXX,4)],[' y: ',num2str(preYY,4)]};
end
