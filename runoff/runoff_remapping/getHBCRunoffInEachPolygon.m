function totalRunoff=getHBCRunoffInEachPolygon(oriRunoff,myRunoffInfo)
%% sum up the runoff flowing into each polygon
%  usage:
%       totalRunoff=getHBCRunoffInEachPolygon(oriRunoff,myRunoffInfo)
%             oriRunoff: km^3/month
%          myRunoffInfo: (input)
%                      .indRunoff  (index of non-zero runoff points)
%                      .isInsideDomain (0: not belong to any pre-defined coast-buffer polygon)
%          totalRunoff : m3/s  %km^3 / month
% history:
%          2014-09: xianmin@ualberta.ca

if nargin==1
   if ischar(oriRunoff) 
      if strcmpi(oriRunoff,'test')
         %runoffMonthlyFile='HYPEdischarge_intNelson.csv';
      %   runoffMonthlyFile=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HBSdischarge_HYPE/HYPEdischarge_intNelson.mat');
      %   nRec=1;
      elseif strcmpi(oriRunoff,'all')
         runoffMonthlyFile=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HYPEcal_HBC/monthlyQ-20170410_GF3.mat'); % this has m3/s
         nRec=430; % WHAT IS THIS NUMBER? number of records?
      else
         error(['unknown input',oriRunoff])
      end
      load coast_buffer_HBC_runoff.mat
      
      for nmon=1:nRec
           cRunoff=squeeze(runoffMonthlyFile(nmon,:,1)); 
           totalRunoff(nmon,:)=getHBCRunoffInEachPolygon(cRunoff,myRunoffInfo);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%% OLD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%       NX=GetNcDimLen(runoffMonthlyFile,'x'); %%% CHANGE
%       NY=GetNcDimLen(runoffMonthlyFile,'y');%%% CHANGE
%       for nmon=1:nRec
%            cRunoff=GetNcVar(runoffMonthlyFile,'runoff',[0 0 nmon-1],[NX NY 1]); % CHANGE
%            totalRunoff(:,nmon)=getHBCRunoffInEachPolygon(cRunoff,myRunoffInfo);
%       end
      %%%%%%%%%%%%%%%%%%%%%%% END OF OLD CODE %%%%%%%%%%%%%%%%%%%%%
       return
   else
      help getHBCRunoffInEachPolygon
      return
   end
elseif nargin~=2
   help getHBCRunoffInEachPolygon
   return
end

if ~isstruct(myRunoffInfo)
   error('myRunoffInfo should be a structure');
end

%% read in runoffInfo
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

totalRunoff=zeros(max(isInsideDomain),1);

for nr=1:numel(indRunoff)
    polyID=isInsideDomain(nr);

    if (polyID>0)
        totalRunoff(polyID)=totalRunoff(polyID)+oriRunoff(indRunoff(nr));                      % m^3/s
    end
end
