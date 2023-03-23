function editModelBufferMaskById(polyID,varargin)
% edit model bufferzone mask  by polygon id
% usage:
%       editModelBufferMaskById(polyID,numNB,varargin)
%       varargin:
%          matfile: the mat file name for model buffer-mask data
%         maskfile: model mask file contains tmask for land-ocean information
%            nloop: number of loop to edit the polygon
%             numb: number of neighbour polygons to show
%               pp: the IDs of (extra) polygons shown as a line feature  
%               ii: x-limit
%               jj: y-limit
%  history:
%       Jan 2015: xianmin@ualberta.ca
opengl software
if nargin==0
   help editModelBufferMaskById
   return
end

if numel(polyID)>1
   for np=1:numel(polyID)
       editModelBufferMaskById(polyID(np),varargin);
   end
   return
end

clc;
numNB=6; 
isCloseFig=1;
isIIRange=0;
isJJRange=0;
isExtraPoly=0;
nSub=1;
matfile='HBC_buffer_mask_ANHA4_stp1.mat';
maskfile='/mnt/storage1/xhu/ANHA4-I/ANHA4_mask.nc';
meshfile='/mnt/storage1/xhu/ANHA4-I/ANHA4_mesh_hgr.nc';

rnffile='coast_buffer_HBC_runoff.mat';
isProj=1;

while(size(varargin,2)>0)
   switch lower(varargin{1})
     case {'fig'}
          isCloseFig=0;
          varargin(1)=[];
     case {'noproj','nomap'}
          isProj=0;
          varargin(1)=[];
     case {'proj','map'}
          isProj=1;
          varargin(1)=[];      
     case {'matfile','bufferfile'}
          matfile=varargin{2};
          varargin(1:2)=[];
     case {'maskfile','mask'}
          maskfile=varargin{2};
          varargin(1:2)=[];     
     case {'mesh','meshfile'}
          meshfile=varargin{2};
          varargin(1:2)=[];
     case {'rnf','rnffile','runoff'}
          rnffile=varargin{2};
          varargin(1:2)=[];      
     case {'nb','numb','neighbour'}
          numNB=varargin{2};
          varargin(1:2)=[];
     case {'loop','nloop'}
          nSub=varargin{2};
          varargin(1:2)=[];
     case {'ii','irange','xx','xlim'}
          isIIRange=1;
          iiRange=varargin{2};
          varargin(1:2)=[];
     case {'jj','jrange','yy','ylim'}
          isJJRange=1;
          jjRange=varargin{2};
          varargin(1:2)=[];
     case {'more','extra','plotpoly','pp'}
          isExtraPoly=1;
          extraPolyID=varargin{2};
          varargin(1:2)=[];
     otherwise
          error('unkown input')
   end
end
if isCloseFig==1
    close all;
end

if ~exist(matfile,'file')
   error([matfile, ' is not found!'])
end
eval(['load ',matfile]);

%% myRunoffInfo
if ~exist(rnffile,'file')
   error([rnffile, ' is not found!'])
end
eval(['load ',rnffile]);

numP=max(bufferMask(:)); %#ok<NODEF>
if polyID>numP || polyID<=0
   error(['polyID must be positive and less than ',num2str(numP)]);
end
modBufferMask=bufferMask;

[NY,NX]=size(bufferMask); [ii,jj]=meshgrid(1:NX,1:NY);
if ischar(maskfile)
   if ~exist(maskfile,'file'), error([maskfile, ' is not found!']); end
   lsmask=GetNcVar(maskfile,'tmask',[0 0  0 0],[NX NY 1 1]);
elseif isnumeric(maskfile)
   lsmask=maskfile; clear maskfile;
end
modBufferMask(bufferMask==polyID&lsmask==0)=0;
modBufferMask(bufferMask==polyID&lsmask==1)=-1;

xx=GetNcVar(meshfile,'nav_lon');
yy=GetNcVar(meshfile,'nav_lat');

indSub=find(bufferMask>=max(1,polyID-numNB) & bufferMask<=polyID+numNB);
bufferMask(bufferMask<polyID-numNB | bufferMask>polyID+numNB)=nan;

figure;
if isProj==0
   [jjSub,iiSub]=ind2sub([NY NX],indSub);
   if isIIRange==0
      myXLIM=[max(1,min(iiSub)-5) min(NX,max(iiSub)+5)];
   else
      myXLIM=iiRange;
   end
   if isJJRange==0
      myYLIM=[max(1,min(jjSub)-5) min(NY,max(jjSub)+5)];
   else
      myYLIM=jjRange;
   end
   mypcolor_imagesc(bufferMask);
   if ~ishold, hold on; end
   contour(lsmask,[0.5 0.5],'k-');
   plot(ii(bufferMask==polyID),jj(bufferMask==polyID),'k+');
else
    xx(xx>0)=nan; % not covered by the HBC
   %m_proj('stereographic','lat',73,'long',-40,'radius',14,'rotate',0);
   m_proj('lambert','long',[-100 -50],'lat',[50 76]);
   m_pcolor(xx,yy,bufferMask); set(findobj('tag','m_pcolor'),'linestyle','none')
   if ~ishold, hold on; end
   m_contour(xx,yy,lsmask,[0.5 0.5],'k-');
   m_plot(xx(bufferMask==polyID),yy(bufferMask==polyID),'k+');
   
   indP=find(myRunoffInfo.isInsideDomain==polyID);
   
   
   runoffData=load('/mnt/storage4/natasha/REMAP_HYPE/HBC/GeoInfo_GFDL-CM3_RCP45_HBC.mat'); 
   runoffLat=runoffData.runoffFile(2,:);
   runoffLon=runoffData.runoffFile(3,:);
   m_plot(runoffLon,runoffLat,'r.')
   % plot rivers where isInsideDomain==0
   ind=find(myRunoffInfo.isInsideDomain==0);
   for t=1:length(ind)
   m_plot(runoffLon(ind(t)),runoffLat(ind(t)),'cp')
   end
   if ~isempty(indP)
      m_plot(myRunoffInfo.runoffLon(indP),myRunoffInfo.runoffLat(indP),'mp')
   end
end
caxis([polyID-numNB-0.5 polyID+numNB+0.5]); 

%% add extra features if requested
if exist('ijLog','var')
   if numel(ijLog)>=polyID %#ok<NODEF>
      if isfield(ijLog{polyID},'ii')
         if isProj==0
            plot(ijLog{polyID}.ii,ijLog{polyID}.jj,'k-o');
         else
            tmpInd=sub2ind([NY NX],fix(ijLog{polyID}.jj),fix(ijLog{polyID}.ii));
            m_plot(xx(tmpInd),yy(tmpInd),'k-o');
         end
      end
   end
   % tricks
   if isExtraPoly==1
      for np=1:numel(extraPolyID)
         if numel(ijLog)>=extraPolyID(np)
            if isfield(ijLog{extraPolyID(np)},'ii')
               if isProj==0
                  plot(ijLog{extraPolyID(np)}.ii,ijLog{extraPolyID(np)}.jj,'w-o');
               else
                  tmpInd=sub2ind([NY NX],fix(ijLog{extraPolyID(np)}.jj),fix(ijLog{extraPolyID(np)}.ii));
                  m_plot(xx(tmpInd),yy(tmpInd),'m-o');
               end
            end
         end
      end
   end
end
colorbar;
cmap=ColorN(numNB*2+1,jet(64)); colormap(cmap);
if isProj==0
   figurenicer;
   axis equal;
   xlim(myXLIM); ylim(myYLIM);
else
   m_grid; cleanmap;
   if ~ishold, hold on; end
   set(gcf,'color','w');
end

title('Click the new polygon (green) points [hit key ENTER when finish]', ...
      'fontweight','bold','color','r');

for nloop=1:nSub
  myCh=input('room into to ROI','s');
  if ~(isempty(myCh) || strcmpi(myCh,'yes'))
     return
  end
  if nloop==1
     [newII,newJJ]=ginput();
     plot(newII,newJJ,'rs-');
  else
     [cII,cJJ]=ginput();
     newII=[newII;cII]; %#ok<*AGROW>
     newJJ=[newJJ;cJJ];
     plot(cII,cJJ,'rs-');
  end
end

if isProj==0
   tmpMask=poly2mask(newII,newJJ,NY,NX);
   ijLog{polyID}.ii=newII; ijLog{polyID}.jj=newJJ;
else
    [pLon,pLat]=m_xy2ll(newII,newJJ);
    indP=zeros(size(pLon));
    for np=1:numel(indP)
        [~,indP(np)]=min(myll2dist(pLon(np),pLat(np),xx(:),yy(:)));
    end
    tmpMask=poly2mask(ii(indP),jj(indP),NY,NX);
    ijLog{polyID}.ii=ii(indP); ijLog{polyID}.jj=jj(indP);
end
modBufferMask(tmpMask==1)=polyID;

%%  manually relocate runoff within current polygon #polyID
myCh=input(['Relocate all runoff in Polygon #',num2str(polyID),' [yes|no]'],'s');
if strcmpi(myCh,'yes')
   title('Click on the new river mouth location','fontweight','bold','fontsize',18)
   [rnfXX,rnfYY]=ginput(1);
   if isProj==1
       [rnfLon,rnfLat]=m_xy2ll(rnfXX,rnfYY);
       ijLog{polyID}.newRunoffLon=rnfLon;
       ijLog{polyID}.newRunoffLat=rnfLat;
   else
       ijLog{polyID}.newRunoffLon=xx(fix(rnfYY),fix(rnfXX));
       ijLog{polyID}.newRunoffLat=yy(fix(rnfYY),fix(rnfXX));
   end
else
   disp('No relocation')
end

%% is from the coast?
myCh=input('Distribute runoff from the coast [yes|no]','s');
if isempty(myCh) || strcmpi(myCh,'yes')
   disp('Will distribute from the coast')
   isFromCoast{polyID}=1; %#ok<*NASGU>
else
   disp('Will distribute from the closest water points')
   isFromCoast{polyID}=0;
end


%% save the file
myCh=input(['save to ',matfile,' [yes|no]'],'s');
if isempty(myCh) || strcmpi(myCh,'yes')
   disp('save the changes')
   bufferMask=modBufferMask;
   eval(['save ',matfile,' bufferMask ijLog isFromCoast']);  
else
   disp('changes are discarded.')
end
