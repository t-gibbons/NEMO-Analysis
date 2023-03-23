function [tarRunoff,myRunoffInfo]=computeModelRunoff(oriRunoff,myRunoffInfo,modelInfo,nClosestPT)
%% interpolate the source runoff onto model grid
%  usage:
%       tarRunoff=computeModelRunoff(oriRunoff,myRunoffInfo,modelInfo)
%          oriRunoff: runoff flux at each data point (make sure same order as the runoff field used in labelRunoff.m
%                     expected unit: m^3/s here
%          myRunoffInfo: (input)
%                      .indRunoff  (index of non-zero runoff points)
%                      .isInsideDomain (0: not belong to any pre-defined coast-buffer polygon)
%                      .runoffLon  (longitude of non-zero runoff points)
%                      .runoffLat  (latitude of non-zero runoff points)
%                      .runoffArea (area of non-zero runoff data point, m^2)             
%                      .avgN (maximum number of points used to place each individual runoff, default:1)
%                      .rnfWGT.rnf{nr}.WPT{n} (weight structure array for the case placing runoff on n points, n could be 1 to avgN)
%                      .rnfWGT.isComputed(nr,n) (WPT is computed already or not)
%                      .rnfWGT.rnf{nr}.WPT{n}.wgt (weight on each  model cell to place the #nr runoff)
%                      .rnfWGT.rnf{nr}.WPT{n}.index (global index of the model cells to place the #nr runoff)
%                      .rnfWGT.rnf{nr}.WPT{n}.dist  (the distances used for weight calculation)
%            modelInfo :
%                      .e1e2tInv (1/e1t/e2t)
%                      .CF    (configuration string)
%                      .isFromCoast (is the runoff is distributed from coast in this polygon (index is the polygon id)
%  see also 
%          labelRunoff, getRunoffArea, getCoastBufferModel, getRunoffIDW
% history:
%          2014-09: xianmin@ualberta.ca

if nargin<3
   help computeModelRunoff
   return
elseif nargin==3
   nClosestPT=3;
   disp('nClosestPT is set to 3: will place runoff over maximum 3 water points!')
end

if ~isstruct(myRunoffInfo)
   error('myRunoffInfo should be a structure');
end
if ~isstruct(modelInfo)
   error('modelInfo should be a structure')
end

%% read in runoffInfo
if ~isfield(myRunoffInfo,'rnfWGT')
   disp('rnfWGT field is not found in myRunoffInfo, computing');
   rnfWGT.isComputed=zeros(myRunoffInfo.avgN,numel(myRunoffInfo.runoffLat));
else
   rnfWGT=myRunoffInfo.rnfWGT;
end
if ~isfield(myRunoffInfo,'indRunoff')
   error('indRunoff field is not found in myRunoffInfo')
else
   indRunoff=myRunoffInfo.indRunoff; 
end
if ~isfield(myRunoffInfo,'isInsideDomain')
   error('isInsideDomain field is not found in myRunoffInfo')
else
   isInsideDomain=myRunoffInfo.isInsideDomain; 
end

%% read in modelInfo
if isfield(modelInfo,'e1e2tInv')
    e1e2tInv=modelInfo.e1e2tInv;
elseif isfield(modelInfo,'e1e2t')
    e1e2tInv=1./modelInfo.e1e2t;
elseif isfield(modelInfo,'e1t') && isfield(modelInfo,'e2t')
    e1e2tInv=1./(e1t.*e2t);
end

tarRunoff=zeros(size(e1e2tInv));
modelInfo.NX=size(tarRunoff,2);
modelInfo.NY=size(tarRunoff,1);

%% fix nClosestPT based on water points in current model mask
if numel(nClosestPT)==1
   nClosestPT=zeros(size(isInsideDomain))+nClosestPT;
end
for nr=1:numel(indRunoff)
    if (myRunoffInfo.isInsideDomain(nr)>0)
       if modelInfo.nWPT(myRunoffInfo.isInsideDomain(nr))==0
          if myRunoffInfo.extraPolygonExist==0
             error(['polyon #',num2str(myRunoffInfo.isInsideDomain(nr)),' is out of model domain'])
          end
          nClosestPT(nr)=0;
       else
          nClosestPT(nr)=min(nClosestPT(nr),modelInfo.nWPT(myRunoffInfo.isInsideDomain(nr)));
       end
    else
       nClosestPT(nr)=0;
    end
end

%% interpolate the runoff onto model grid based on the weights
for nr=1:numel(indRunoff)
    if (isInsideDomain(nr)>0)
       if nClosestPT(nr)==0
          disp(['river #',num2str(nr),' is not in model domain within polygon #',num2str(isInsideDomain(nr))])
       else
          % convert runoff data to volume equivalent if necessary 
          rnfVol=oriRunoff(indRunoff(nr))*1000;           % m^3/s * kg/m^3  ==> kg/s
          
          % spread onto the closest water point in model
          if (rnfWGT.isComputed(nClosestPT(nr),nr)==0)
             [wgtC,distC,indexC]=computeWGT(myRunoffInfo,modelInfo,nr,nClosestPT(nr));
             rnfWGT.isComputed(nClosestPT(nr),nr)=1;
             rnfWGT.rnf{nr}.WPT{nClosestPT(nr)}.wgt=wgtC;
             rnfWGT.rnf{nr}.WPT{nClosestPT(nr)}.dist=distC;
             rnfWGT.rnf{nr}.WPT{nClosestPT(nr)}.index=indexC;
          else
             wgtC=rnfWGT.rnf{nr}.WPT{nClosestPT(nr)}.wgt;
             indexC=rnfWGT.rnf{nr}.WPT{nClosestPT(nr)}.index;
          end
          if nClosestPT(nr)==1
             tarRunoff(indexC)=tarRunoff(indexC)+rnfVol*e1e2tInv(indexC);
          else
             tarRunoff(indexC)=tarRunoff(indexC) + rnfVol.*wgtC.*e1e2tInv(indexC);
          end
       end
    %else
    %   disp(['river #',num2str(nr),' does not belong to any polygon'])
    end
end
% save for later re-mapping
myRunoffInfo.rnfWGT=rnfWGT;

function [wgt,dist,index]=computeWGT(myRunoffInfo,modelInfo,nr,avgN)
  indModelPWTC=find(modelInfo.bufferMask==myRunoffInfo.isInsideDomain(nr) & modelInfo.lsmask==1);
  if ~isfield(modelInfo,'isFromCoast')
     isFromCoast=1;
  else
     if numel(modelInfo.isFromCoast)<myRunoffInfo.isInsideDomain(nr)
        isFromCoast=1;
     else
        isFromCoast=modelInfo.isFromCoast{myRunoffInfo.isInsideDomain(nr)};
        if isempty(isFromCoast), isFromCoast=1; end
     end
  end

  nValidWPT=numel(indModelPWTC);
  NX=modelInfo.NX;
  NY=modelInfo.NY;
  if nValidWPT<avgN
     error('computeWGT: should not see this msg')
  elseif avgN==1
      wgt=1; dist=0;
      ppdist=myll2dist(myRunoffInfo.runoffLon(nr),myRunoffInfo.runoffLat(nr),modelInfo.lon(indModelPWTC),modelInfo.lat(indModelPWTC));
      [~,indexSort]=sort(ppdist);
      n=1; 
      if isFromCoast==1
         % make sure find the coastal point as river mouth
         [jj,ii]=ind2sub([NY NX],indModelPWTC(indexSort));
         tmpSumMask=modelInfo.lsmask(max(1,jj(n)-1):min(NY,jj(n)+1),max(1,ii(n)-1):min(NX,ii(n)+1));
         while min(tmpSumMask(:))==1 && n<nValidWPT % all water points
           n=n+1;
           tmpSumMask=modelInfo.lsmask(max(1,jj(n)-1):min(NY,jj(n)+1),max(1,ii(n)-1):min(NX,ii(n)));
         end
      end
      index=indModelPWTC(indexSort(n));
      return
  else
     % first search the closest water point (might be okay to use a simplified distance calculation?)
     ppdist=myll2dist(myRunoffInfo.runoffLon(nr),myRunoffInfo.runoffLat(nr),modelInfo.lon(indModelPWTC),modelInfo.lat(indModelPWTC));
     
     if isFromCoast==1
        [~,indexSort]=sort(ppdist);
         
        % make sure find the coastal point as river mouth
        [jj,ii]=ind2sub([NY NX],indModelPWTC(indexSort));
        n=1;
        tmpSumMask=modelInfo.lsmask(max(1,jj(n)-1):min(NY,jj(n)+1),max(1,ii(n)-1):min(NX,ii(n)+1));
        while min(tmpSumMask(:))==1 && n<nValidWPT % all water points
              n=n+1;
              tmpSumMask=modelInfo.lsmask(max(1,jj(n)-1):min(NY,jj(n)+1),max(1,ii(n)-1):min(NX,ii(n)));
        end
         
        % reset the center location
        tmpLonC=modelInfo.lon(indModelPWTC(indexSort(n)));  tmpLatC=modelInfo.lat(indModelPWTC(indexSort(n)));
        ppdist=myll2dist(tmpLonC,tmpLatC,modelInfo.lon(indModelPWTC),modelInfo.lat(indModelPWTC));
        [distSort,indexSort]=sort(ppdist);
     else
        [distSort,indexSort]=sort(ppdist);
     end
     index(1,1:avgN)=indModelPWTC(indexSort(1:avgN));
     dist(1,1:avgN)=[distSort(2) reshape(distSort(2:avgN),1,[])];
     
     if avgN<4
        wgt(1,1:avgN)=1/avgN;
     else
        tmpwgt=cos(distSort(1:avgN)*pi*0.25/distSort(avgN));
        wgt=reshape(tmpwgt/sum(tmpwgt),1,[]);
    end
  end

function dist=myll2dist(lon0,lat0,lon,lat)
% origins from m_lldist (m_map package)
% compute the distance from (lon0,lat0)
% lon0, lat0: start point(s)
% lon, lat:   end point(s)
pi180=pi/180;
earth_radius=6378.137;
lat0=lat0*pi180;
lat=lat*pi180;
if numel(lon0)~=1
   np0=numel(lon0);
   np=numel(lon);
   lon0=repmat(lon0(:),1,np);
   lat0=repmat(lat0(:),1,np);
   lon=repmat(lon(:)',np0,1);
   lat=repmat(lat(:)',np0,1);
end
dlon =(lon-lon0)*pi180;
dlat = lat-lat0;
a = (sin(dlat/2)).^2 + cos(lat0)* cos(lat) .* (sin(dlon/2)).^2;
angles = 2 * atan2( sqrt(a), sqrt(1-a) );
dist = earth_radius * angles;
