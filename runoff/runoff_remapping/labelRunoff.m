function [ myRunoffInfo ] = labelRunoff( rnfInfo )
%% find the polygon that each runoff data point belong to
% usage:
%       myRunoffInfo=labelRunoff(rnfInfo)
%
%       rnfInfo:
%              .runoff  (including zeros)
%              .lon     (longitude of input runoff data, -180 -- 180)
%              .lat     (latitude of input runoff data, -90 -- 90)
%              .coast_buffer (pre-defined coast-buffer polygons)
%                       if not given, will try read it from a file named coast_buffer.mat
%              .savefile (file name to save the result)
%       myRunoffInfo:
%              .indRunoff  (index of non-zero runoff points)
%              .isInsideDomain (0: not belong to any pre-defined coast-buffer polygon)
%              .runoffLon  (longitude of non-zero runoff points)
%              .runoffLat  (latitude of non-zero runoff points)
% history:
%       2014-09: xianmin@ualberta.ca

clc;
if nargin == 0
   help labelRunoff
   return
end

%% read in the runoff data
if ischar(rnfInfo) 
   if strcmpi(rnfInfo,'default')
      runoffFile      = 'runoff_ave.nc';
      navLon          = GetNcVar(runoffFile,'LON');
      navLat          = GetNcVar(runoffFile,'LAT');
      [navLon,navLat] = meshgrid(navLon,navLat);
      runoff          = GetNcVar(runoffFile,'runoff');
   else
       help labelRunoff
       return
   end
elseif isstruct(rnfInfo)
   if ~isfield(rnfInfo,'runoff')
      error('runoff field is not found in input rnfInfo!')
   end
   
   if ~isfield(rnfInfo,'lon')
      error('lon field is not found in input rnfInfo!')
   end
   
   if ~isfield(rnfInfo,'lat')
      error('lat field is not found in input rnfInfo!')
   end
   
   if ~isfield(rnfInfo,'isGridded')
       isGridded = 0;
   else
       isGridded = rnfInfo.isGridded;
   end
   
   if isvector(rnfInfo.lon) && isGridded == 1
      [navLon,navLat] = meshgrid(rnfInfo.lon,rnfInfo.lat);
   else
      navLon = rnfInfo.lon;
      navLat = rnfInfo.lat;
   end
   
   runoff = rnfInfo.runoff;
   
   if isfield(rnfInfo,'coast_buffer')
      coast_buffer = rnfInfo.coast_buffer;
   end
   
   if isfield(rnfInfo,'savefile')
      saveMatFile = rnfInfo.savefile;
   end
else
   help labelRunoff
   return
end

%% find the non-zero value points
indRunoff      = find(runoff>0);
runoffLon      = navLon(indRunoff);
runoffLat      = navLat(indRunoff);
isInsideDomain = zeros(size(indRunoff));

%% load coast buffer zone
if ~exist('coast_buffer','var')
   error('need coast_buffer')
end

latMax = 90;
latMin = -90;
lonMax = 180;
%lonMin = -80;
lonMin = -180;


for np=1:numel(coast_buffer)
    
    latMax = max([latMax,max(coast_buffer{np}.lat(:))]);
    latMin = min([latMin,min(coast_buffer{np}.lat(:))]);
    lonMax = max([lonMax,max(coast_buffer{np}.lon)]);
    lonMin = min([lonMin,min(coast_buffer{np}.lon)]);
    
end

%% find the polygon that each river is located in
for nr = 1:length(indRunoff)
    
    tmpLon = runoffLon(nr);
    tmpLat = runoffLat(nr);
    
    if tmpLat < latMin || tmpLat > latMax || tmpLon < lonMin || tmpLon > lonMax
       isInsideDomain(nr) = 0;
    else
        
       for np = 1:numel(coast_buffer)
           isInThisPoly = inpolygon(tmpLon,tmpLat,coast_buffer{np}.lon,coast_buffer{np}.lat);
           
           if isInThisPoly == 1
              isInsideDomain(nr) = np;
              break
           end 
       end
       %if mouth isn't in polygon but in domain, need to find nearest polygon to assign runoff to
       % only exception will be for those outside the Bering Strait part of the domain
       if abs(tmpLon) > 150 & tmpLat < 65
	  isInsideDomain(nr) = 0;
       else
	  dist_min=100000000.0;
	  iwhich=0;
	  for np = 1:numel(coast_buffer)
  		dist=myll2dist(tmpLon,tmpLat,coast_buffer{np}.lon(1),coast_buffer{np}.lat(1));
	 	if dist < dist_min
		   dist_min=dist;
		   iwhich=np;
 		end
          end
	  isInsideDomain(nr) = iwhich;
       end
    end
end

%% check the output
myRunoffInfo.isInsideDomain = isInsideDomain;
myRunoffInfo.indRunoff      = indRunoff;
myRunoffInfo.runoffLon      = runoffLon;
myRunoffInfo.runoffLat      = runoffLat;
myRunoffInfo.numRiver       = numel(runoffLat);

if exist('saveMatFile','var')
   if length(saveMatFile)>4
      eval(['save ',saveMatFile,' myRunoffInfo'])
   end
else
   if nargout == 0
      save runoff_labeled.mat myRunoffInfo
   end
end
