% Compare HYPE runoff on original and ANHA grid


runoffMonthlyFile=importdata('/mnt/storage1/natasha/DATA/RIVER_RUNOFF/HYPEcal_HBC/monthlyQ-20170410_WFDEI.mat'); % ** CHANGE
runoffMonthlyFile=squeeze(runoffMonthlyFile(:,:,1));


NX=544; %x-dimension
NY=800; %y-dimension


FWANHA=GetNcVar('ANHA4_HBC_HYPEcal_WFDEI_runoff_monthly_y1979_y2013.nc','runoff',[0 0 0],[NX NY 420]);%*****
e1t=GetNcVar('/mnt/storage1/xhu/ANHA4-I/ANHA4_mesh_hgr.nc','e1t');
e2t=GetNcVar('/mnt/storage1/xhu/ANHA4-I/ANHA4_mesh_hgr.nc','e2t');
nav_lon=GetNcVar('/mnt/storage1/xhu/ANHA4-I/ANHA4_mesh_hgr.nc','nav_lon');
nav_lat=GetNcVar('/mnt/storage1/xhu/ANHA4-I/ANHA4_mesh_hgr.nc','nav_lat');


for y=1:420
    monthly_ori(y,1)=sum(runoffMonthlyFile(y,:));
    monthly_ANHA(y,1)=sum(sum((squeeze(FWANHA(y,:,:)).*e1t.*e2t)/1000));    % kg/m2/s ==> m3/s
end

test=monthly_ANHA-monthly_ori;
figure;
plot(monthly_ori)
hold on
plot(monthly_ANHA)
hold on
plot(test,'Color','k')



