function drawMask()
% make a tmask file by drawing a polygon

clc;
global NX NY tmask
xxLog=[]; yyLog=[];

CF='ANHA4';

if strcmpi(CF,'anha4')
   NX=544; NY=800;
   maskfile='/mnt/storage0/xhu/CREG025-I/CREG025-CICEREF_tmask.nc';
elseif strcmpi(CF,'anha12')
   NX=1632; NY=2400;
   maskfile='/mnt/storage0/xhu/CREG012-I/2014/mask.nc';
elseif strcmpi(CF,'anha2')
   NX=272; NY=400;
   maskfile='/mnt/storage1/xhu/ANHA2-I/ANHA2_mesh_mask.nc';
elseif strcmpi(CF,'creg12')
   NX=1580; NY=1817;
   maskfile='/home/xhu/RUNOFF/creg12/CREG12.L75_tmask.nc';
end
% load surface mask file
tmask=squeeze(GetNcVar(maskfile,'tmask',[0 0 0 0],[NX NY 1 1]));

HBCRiver=GetNcVar('ANHA4_HBC_HYPEint_runoff_monthly_y1979_y2010.nc','runoff');  % better to loop all years
[xx,yy]=meshgrid(1:NX,1:NY);

% sum average distribution for all years -- find non0 points
HBCRiver=squeeze(nansum(HBCRiver,1));


figure;
set(gcf,'render','zbuffer','color','w','CloseRequestFcn',@clearMyGlobal);
imagesc(1:NX,1:NY,tmask); 
hold on;
%% plot the location of all non-zero HBC runoff points
plot(xx(HBCRiver>0),yy(HBCRiver>0),'g+');
set(gca,'ydir','normal','tickdir','out');
caxis([0 2])
axis equal; axis tight;

% datacursor
mycobj=datacursormode(gcf);
datacursormode on;
set(mycobj,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on','DisplayStyle','window');

end

%---------------------------------------------------------------
function clearMyGlobal(varargin)
 % save some variables
 global xxLog yyLog NX NY tmask
 % making the mask
 bwmask=poly2mask(xxLog,yyLog,NY,NX);
 tmask(bwmask==0)=0;
 assignin('base','xxLog',xxLog);
 assignin('base','yyLog',yyLog);
 assignin('base','tmask',tmask);
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
