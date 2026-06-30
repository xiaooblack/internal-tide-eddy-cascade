clear
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
grdname='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc.2';
h=ncread(grdname,'h');f=ncread(grdname,'f');
lon_rho= ncread(grdname,'lon_rho',[300 1],[1700 1300]);   lat_rho= ncread(grdname,'lat_rho',[300 1],[1700 1300]);   % X/Y rho  :km
x_rho= ncread(grdname,'x_rho',[300 1],[1700 1300]);   y_rho= ncread(grdname,'y_rho',[300 1],[1700 1300]);   % X/Y rho  :km


% create nc file
nccreate('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','u','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','v','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','zeta','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');

step1=0;
for ii=1:695
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    u=ncread('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide.nc','u',[300 1 ii],[1700 1300 1]);
    v=ncread('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide.nc','v',[300 1 ii],[1700 1300 1]);
    [dudx,dudy]=model_gradient(x_rho,y_rho,u);
    [dvdx,dvdy]=model_gradient(x_rho,y_rho,v);
    zeta=(dvdx-dudy);


% write nc file
    ncwrite('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','u',single(u),[1 1 step1]);
    ncwrite('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','v',single(v),[1 1 step1]);
    ncwrite('I:\NCS_wu_0.5km_notide_15m\NCS_wu_0.5km_15m_notide_origin.nc','zeta',single(zeta),[1 1 step1]);
end




clear
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
grdname='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc.2';
h=ncread(grdname,'h');f=ncread(grdname,'f');
lon_rho= ncread(grdname,'lon_rho',[300 1],[1700 1300]);   lat_rho= ncread(grdname,'lat_rho',[300 1],[1700 1300]);   % X/Y rho  :km
x_rho= ncread(grdname,'x_rho',[300 1],[1700 1300]);   y_rho= ncread(grdname,'y_rho',[300 1],[1700 1300]);   % X/Y rho  :km


% create nc file
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','x_rho',single(x_rho),[1 1]);
ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','y_rho',single(y_rho),[1 1]);
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','u','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','v','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','zeta','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');

step1=0;
for ii=1:695
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    u=ncread('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','u',[300 1 ii],[1700 1300 1]);
    v=ncread('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','v',[300 1 ii],[1700 1300 1]);
    [dudx,dudy]=model_gradient(x_rho,y_rho,u);
    [dvdx,dvdy]=model_gradient(x_rho,y_rho,v);
    zeta=(dvdx-dudy);


% write nc file
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','u',single(u),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','v',single(v),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc','zeta',single(zeta),[1 1 step1]);
end
