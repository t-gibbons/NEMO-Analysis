function bufferMask = getCoastBufferModelMask(modelInfo)
%% create the mask for coast buffer polygons
%  usage:
%       bufferMask=getCoastBufferModelMask(modelInfo)
%       modelInfo:
%                .lsmask (model surface, 1st level, land-ocean mask, 1 for water and 0 for land)
%                .CF    (configuration string)
%                .lon   (longitude of model t-point, -180 --180, 2D)
%                .lat   (latitude  of model t-point, -90  -- 90, 2D)
%                .coast_buffer (pre-defined coast-buffer polygons)
%                       if not given, will try read it from a file named coast_buffer.mat
%                .savefile (file name to save bufferMask)
%      bufferMask: 
%                 0: land 
%                -1: ocean but not belong to any pre-defined coast-buffer polygon
%                >0: ocean/land within pre-defined polygons
%            nWPT:
%                  number of water points in each polygon
% history:
%          2014-09: xianmin@ualberta.ca

clc;
if nargin==0
   help getCoastBufferModelMask
   return
end

%% read in model
CF = 'model';
if isstruct(modelInfo)
   if ~isfield(modelInfo,'lsmask')
      error('lsmask field is not found in input modelInfo!')
   end
   
   if ~isfield(modelInfo,'lon')
      error('lon field is not found in input modelInfo!')
   end
   
   if ~isfield(modelInfo,'lat')
      error('lat field is not found in input modelInfo!')
   end
   
   if isvector(modelInfo.lon)
       error('lon field is supposed to be a 2d array');
   end
   
   if isvector(modelInfo.lat)
       error('lon field is supposed to be a 2d array');
   end
   
   bufferMask = modelInfo.lsmask;
   navLon     = modelInfo.lon;
   navLat     = modelInfo.lat;
   
   if isfield(modelInfo,'coast_buffer')
      coast_buffer=modelInfo.coast_buffer;
   end
   
   if isfield(modelInfo,'savefile')
      saveMatFile=modelInfo.savefile;
   end

else
   help getCoastBufferModelMask
   return
end

bufferMask(bufferMask == 1) = -1; % set water-point to -1
[NY,NX]                     = size(bufferMask);

%% load coast buffer zone
if ~exist('coast_buffer','var')
   if ~exist('coast_buffer.mat','file')
      error('need coast_buffer.mat')
   end
   load coast_buffer
end

numP = numel(coast_buffer);

for np = 1:numP
    tmpLon = coast_buffer{np}.lon;
    tmpLat = coast_buffer{np}.lat;
    tmpInd = zeros(size(tmpLon));% find closest i,j
    
    for nn = 1:numel(tmpLat)
        [~,tmpInd(nn)] = min((navLon(:)-tmpLon(nn)).^2 + (navLat(:)-tmpLat(nn)).^2);
    end
    
    [tmpJJ,tmpII]           = ind2sub([NY NX],tmpInd);
    bwmask                  = poly2mask(tmpII,tmpJJ,NY,NX); % polygon 2 mask
    bufferMask(bwmask == 1) = np;
end

%% check the output
if exist('saveMatFile','var')
   if length(saveMatFile)>4
      eval(['save ',saveMatFile,' bufferMask'])
   end
else
   if nargout == 0
      eval(['save coast_buffer_mask_',CF,'.mat bufferMask'])
   end
end
