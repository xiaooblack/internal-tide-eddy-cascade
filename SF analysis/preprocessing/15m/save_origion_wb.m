%save_origion w&b
%%%% tide %%%%
clear
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
grdname='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc.2';
h=ncread(grdname,'h');f=ncread(grdname,'f');
lon_rho= ncread(grdname,'lon_rho',[449 100],[1200 1200]);   lat_rho= ncread(grdname,'lat_rho',[449 100],[1200 1200]);   % X/Y rho  :km
x_rho= ncread(grdname,'x_rho',[449 100],[1200 1200]);   y_rho= ncread(grdname,'y_rho',[449 100],[1200 1200]);   % X/Y rho  :km


N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='F:\NCS_wu_0.5km_tide_hourly_output';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

% create nc file
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','y_rho',single(y_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','lat_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','w','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','b','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');


step1=0;
for ii=1:m1
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[filelist(ii).folder,'\',filelist(ii).name];  

    temp=ncread(filename,'temp',[449 100 1 1],[1200 1200 inf inf]);
    salt=ncread(filename,'salt',[449 100 1 1],[1200 1200 inf inf]);
    w0=ncread(filename,'w',[449 100 1 1],[1200 1200 inf inf]);
    rho=sw_dens0(salt,temp);
    b0=-g.*rho./rho_r;
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(41:70,449:1648,100:1299);
    z_rho=permute(z_rho,[2 3 1]);


    w=interp_model_z(w0,x_rho,y_rho,z_rho,z_need);
    b=interp_model_z(b0,x_rho,y_rho,z_rho,z_need);
    % 
    % write nc file
    ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','w',single(w),[1 1 step1]);
    ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc','b',single(b),[1 1 step1]);

end

%%%% notide %%%%

clear
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
grdname='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc.2';
h=ncread(grdname,'h');f=ncread(grdname,'f');
lon_rho= ncread(grdname,'lon_rho',[449 100],[1200 1200]);   lat_rho= ncread(grdname,'lat_rho',[449 100],[1200 1200]);   % X/Y rho  :km
x_rho= ncread(grdname,'x_rho',[449 100],[1200 1200]);   y_rho= ncread(grdname,'y_rho',[449 100],[1200 1200]);   % X/Y rho  :km


N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='F:\NCS_wu_0.5km_notide_hourly_output';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

% create nc file
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','y_rho',single(y_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','x_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','y_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','w','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','b','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');


step1=0;
for ii=1:m1
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[filelist(ii).folder,'\',filelist(ii).name];  

    temp=ncread(filename,'temp',[449 100 1 1],[1200 1200 inf inf]);
    salt=ncread(filename,'salt',[449 100 1 1],[1200 1200 inf inf]);
    w0=ncread(filename,'w',[449 100 1 1],[1200 1200 inf inf]);
    rho=sw_dens0(salt,temp);
    b0=-g.*rho./rho_r;
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(41:70,449:1648,100:1299);
    z_rho=permute(z_rho,[2 3 1]);


    w=interp_model_z(w0,x_rho,y_rho,z_rho,z_need);
    b=interp_model_z(b0,x_rho,y_rho,z_rho,z_need);
    % 
    % write nc file
    ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','w',single(w),[1 1 step1]);
    ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc','b',single(b),[1 1 step1]);

end

