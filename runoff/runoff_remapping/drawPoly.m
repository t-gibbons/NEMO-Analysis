function varargout = drawPoly(varargin)
% DRAWPOLY M-file for drawPoly.fig
%      DRAWPOLY, by itself, creates a new DRAWPOLY or raises the existing
%      singleton*.
%
%      H = DRAWPOLY returns the handle to a new DRAWPOLY or the handle to
%      the existing singleton*.
%
%      DRAWPOLY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRAWPOLY.M with the given input arguments.
%
%      DRAWPOLY('Property','Value',...) creates a new DRAWPOLY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drawPoly_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drawPoly_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help drawPoly

% History:
%          2014-2015 : xianmin@ualberta.ca
%        Apr 21 2016 : modified for HBC case
%                      add model grid for a better interface

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drawPoly_OpeningFcn, ...
                   'gui_OutputFcn',  @drawPoly_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before drawPoly is made visible.
function drawPoly_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drawPoly (see VARARGIN)

% Choose default command line output for drawPoly
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes drawPoly wait for user response (see UIRESUME)
% uiwait(handles.drawPolygon);
set(gcf, 'units','normalized');
set(gca,'units','normalized','position',[0.13 0.11 0.8 0.8]);

global lonLog latLog NX NY myPolyRec
lonLog=[]; latLog=[];

% load original runoff data (modify it to load runoff information)
isPointData=1;
if isPointData==0
   % the case of Dai & Trenberth data and Greenland melt runoff data
   runoffFile='/mnt/storage0/myers/DATA/RACMO/FW_GrIS_1850_RACMO2_RU_basins+D.nc';
   runoffLon=GetNcVar(runoffFile,'lon');
   runoffLat=GetNcVar(runoffFile,'lat');
   if isvector(runoffLon)
      [runoffLon,runoffLat]=meshgrid(runoffLon,runoffLat);
   end
   [NY,NX]=size(runoffLon);
   runoff=GetNcVar(runoffFile,'FW_calving'); runoff(runoff==0)=nan;
else
   %runoffFile='HBC_runoff_climatology_data.dat'; % assume to have three columns: lon, lat, volume-flux
   %runoffData=load(runoffFile);
   %runoffLon=runoffData(:,1); 
   %runoffLat=runoffData(:,2); 
   %runoff=runoffData(:,3); clear runoffData

   runoffData=load('/mnt/storage4/natasha/REMAP_HYPE/HBC/GeoInfo_GFDL-CM3_RCP45_HBC.mat');
   runoffLat=runoffData.runoffFile(2,:);
   runoffLon=runoffData.runoffFile(3,:);
   runoffVol=runoffData.runoffFile(4,:); % m3/s
   % convert to km3/yr
   runoffVol=runoffVol*((3600*24*365)/1000000000);
   
end

% modified by NR
% set up colormatrix and  volume for colormatrix
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

% load model tmask
modelMaskFile='/mnt/storage1/xhu/ANHA4-I/ANHA4_mask.nc';
tmask=GetNcVar(modelMaskFile,'tmask',[0 0 0 0],[544 800 1 1]);
modLon=GetNcVar(modelMaskFile,'nav_lon');
modLat=GetNcVar(modelMaskFile,'nav_lat');
modLon(modLon>0)=nan; % not covered by the HBC

if exist('HBCRiverMouthMask.nc')
   runoffmask=GetNcVar('HBCRiverMouthMask.nc','socoefr',[0 0],[544 800]);
   runoffmask(runoffmask==0)=nan;
end

set(gcf,'color','w','Pointer','cross','CloseRequestFcn',@clearMyGlobal);
%m_proj('Azimuthal Equal-area','lat',72.5,'long',-45,'radius',22,'rect','on') % Greenland
m_proj('lambert','long',[-100 -50],'lat',[50 76]);
if isPointData==0
   indRunoff=find(runoff>0);
   m_plot(runoffLon(indRunoff),runoffLat(indRunoff),'r.','markersize',8)
else
    
   % better to plot the river location with different color or markers to classify rivers with different discharge
   for c=2:length(DischargeCM)
    if c<=length(DischargeCM)
        ind=find(runoffVol(:)<=DischargeCM(c) & runoffVol(:)>DischargeCM(c-1));       
    end
    m_line(runoffLon(ind),runoffLat(ind),'LineStyle','none','Marker','.','markersize',10,'Color',colourmatrix(c-1,:))
  %  m_line(HBSdischarge_geoinfo(ind,4),HBSdischarge_geoinfo(ind,3),'LineStyle','none','Marker','.','MarkerSize',10,'Color',colourmatrix(c-1,:));
    hold on
   end
end
hold on;
m_ShowGridLine(modLon,modLat,1,[0.6 0.6 0.6])
tmask(tmask==0)=nan;
m_pcolor(modLon,modLat,tmask); %set(findobj(gca,'tag','m_pcolor'),'linestyle','none'); % may need this line
hold on
m_pcolor(modLon,modLat,runoffmask); set(findobj(gca,'tag','m_pcolor'),'linestyle','none');
m_gshhs_i('color',[0.7 0.7 0.7])
m_grid; hold on;
cleanmap;

% global variables
if exist('HBC_buffer.mat','file')
    load HBC_buffer.mat coast_buffer
    pID=numel(coast_buffer);
    for np=1:pID
        hline=m_plot(coast_buffer{np}.lon,coast_buffer{np}.lat,'-o'); set(hline,'color',[0.3 0.3 0.3],'markeredgecolor','b','markersize',8)
    end
    set(handles.bt_next,'UserData',pID+1);  % polygon ID
else
   set(handles.bt_next,'UserData',1);  % polygon ID
   coast_buffer={};
end
% datacursor
mycobj=datacursormode(gcf);
datacursormode on;
set(mycobj,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','off','Enable','on','DisplayStyle','window');
title('Click Polygon Points','fontweight','bold')

% --- Outputs from this function are returned to the command line.
function varargout = drawPoly_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in bt_load.
function bt_load_Callback(hObject, eventdata, handles)
% hObject    handle to bt_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oriNcFile=uigetfile('*.nc', 'Select a mask file');
if ~exist(oriNcFile,'file')
  disp('Please load the mask file')
  return
end
ncfid=netcdf.open(oriNcFile,'NC_NOWRITE');
idX=netcdf.inqDimID(ncfid,'x'); [~,NX]=netcdf.inqDim(ncfid,idX);
idY=netcdf.inqDimID(ncfid,'y'); [~,NY]=netcdf.inqDim(ncfid,idY);
netcdf.close(ncfid);
fileInfo.name=oriNcFile;
fileInfo.NX=NX;
fileInfo.NY=NY;
set(handles.bt_load,'UserData',fileInfo);

% --- Executes on button press in bt_next.
function bt_next_Callback(hObject, eventdata, handles)
% hObject    handle to bt_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global lonLog latLog coast_buffer preXX preYY
pID=get(handles.bt_next,'UserData');
if ~isempty(lonLog)
   [lonLog,latLog]=myUnique(lonLog,latLog); 
   coast_buffer{pID}.lon=lonLog;
   coast_buffer{pID}.lat=latLog;
   lonLog=[]; latLog=[]; preXX=[]; preYY=[];
   set(handles.bt_next,'UserData',pID+1);
end

% --- Executes on button press in bt_save.
function bt_save_Callback(hObject, eventdata, handles)
% hObject    handle to bt_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global coast_buffer
save HBCCoast_buffer.mat coast_buffer


function clearMyGlobal(varargin)
clearvars -global *
delete(gcf);
 
function output_txt = myupdatefcn(~,~)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

 global preXX preYY lonLog latLog
 pID=get(findobj(gcf,'tag','bt_next'),'UserData');
 tmpH=datacursormode(gcf);
 if strcmpi(tmpH.Enable,'on') 
    pos = get(gca,'CurrentPoint');     
    [clon,clat]=m_xy2ll(pos(1,1),pos(1,2));
    if (clon<-180), clon=clon+360; end
    if ~ishold, hold on; end
    plot(pos(1,1),pos(1,2),'ko');  % label the previous selections
    lonLog=[lonLog clon];
    latLog=[latLog clat];

    if ~isempty(preXX)
       plot([preXX pos(1,1)],[preYY pos(1,2)]);
    end
    % update previous location
    preXX=pos(1,1); preYY=pos(1,2);
    output_txt = {['Lon: ',num2str(clon,4)],[' Lat: ',num2str(clat,4)],[' id: ',num2str(pID)]};
 end

% --- Executes when drawPolygon is resized.
function drawPolygon_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to drawPolygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gca,'position',[0.13 0.11 0.8 0.8]);


% --- Executes on button press in bt_cancel.
function bt_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to bt_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global lonLog latLog preXX preYY
lonLog=[]; latLog=[]; preXX=[]; preYY=[];
