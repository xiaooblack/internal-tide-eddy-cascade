% 取模式15m数据存到NC文件里 u v temp salt w dudz(单层还真没法算这个)
clear
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
grdname='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc.2';
h=ncread(grdname,'h');f=ncread(grdname,'f');
lon_rho= ncread(grdname,'lon_rho');   lat_rho= ncread(grdname,'lat_rho');   % X/Y rho  :km
lon_u  = ncread(grdname,'lon_u');     lat_u  = ncread(grdname,'lat_u');     % X/Y u    :km
lon_v  = ncread(grdname,'lon_v');     lat_v  = ncread(grdname,'lat_v');     % X/Y v    :km
x_rho= ncread(grdname,'x_rho');   y_rho= ncread(grdname,'y_rho');   % X/Y rho  :km
x_u  = ncread(grdname,'x_u');     y_u  = ncread(grdname,'y_u');     % X/Y u    :km
x_v  = ncread(grdname,'x_v');     y_v  = ncread(grdname,'y_v');     % X/Y v    :km

N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='F:\NCS_wu_0.5km_tide_hourly_output';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

% create nc file
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','u','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','v','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','w','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','temp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','salt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');

step1=0;
for ii=1:m1
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    u_u=ncread(filename,'u');
    v_v=ncread(filename,'v');
    w_m=ncread(filename,'w');
    temp_m=ncread(filename,'temp');
    salt_m=ncread(filename,'salt');
    for iii=1:size(u_u,3)
        u_m(:,:,iii)=v2rho_2d(u_u(:,:,iii));
        v_m(:,:,iii)=u2rho_2d(v_v(:,:,iii));
    end
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(41:70,:,:);
    z_rho=permute(z_rho,[2 3 1]);

    u=interp_model_z(u_m,x_rho,y_rho,z_rho,z_need);
    v=interp_model_z(v_m,x_rho,y_rho,z_rho,z_need);
    w=interp_model_z(w_m,x_rho,y_rho,z_rho,z_need);
    temp=interp_model_z(temp_m,x_rho,y_rho,z_rho,z_need);
    salt=interp_model_z(salt_m,x_rho,y_rho,z_rho,z_need);

    % write nc file
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','u',single(u),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','v',single(v),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','w',single(w),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','temp',single(temp),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide.nc','salt',single(salt),[1 1 step1]);
end

%%
clear
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
grdname='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc.2';
h=ncread(grdname,'h');f=ncread(grdname,'f');
lon_rho= ncread(grdname,'lon_rho');   lat_rho= ncread(grdname,'lat_rho');   % X/Y rho  :km
lon_u  = ncread(grdname,'lon_u');     lat_u  = ncread(grdname,'lat_u');     % X/Y u    :km
lon_v  = ncread(grdname,'lon_v');     lat_v  = ncread(grdname,'lat_v');     % X/Y v    :km
x_rho= ncread(grdname,'x_rho');   y_rho= ncread(grdname,'y_rho');   % X/Y rho  :km
x_u  = ncread(grdname,'x_u');     y_u  = ncread(grdname,'y_u');     % X/Y u    :km
x_v  = ncread(grdname,'x_v');     y_v  = ncread(grdname,'y_v');     % X/Y v    :km

N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='F:\NCS_wu_0.5km_notide_hourly_output';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

% create nc file
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','u','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','v','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','w','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','temp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');
nccreate('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','salt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',m1},'Datatype','single');

step1=0;
for ii=1:m1
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    u_u=ncread(filename,'u');
    v_v=ncread(filename,'v');
    w_m=ncread(filename,'w');
    temp_m=ncread(filename,'temp');
    salt_m=ncread(filename,'salt');
    for iii=1:size(u_u,3)
        u_m(:,:,iii)=v2rho_2d(u_u(:,:,iii));
        v_m(:,:,iii)=u2rho_2d(v_v(:,:,iii));
    end
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(41:70,:,:);
    z_rho=permute(z_rho,[2 3 1]);

    u=interp_model_z(u_m,x_rho,y_rho,z_rho,z_need);
    v=interp_model_z(v_m,x_rho,y_rho,z_rho,z_need);
    w=interp_model_z(w_m,x_rho,y_rho,z_rho,z_need);
    temp=interp_model_z(temp_m,x_rho,y_rho,z_rho,z_need);
    salt=interp_model_z(salt_m,x_rho,y_rho,z_rho,z_need);

    % write nc file
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','u',single(u),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','v',single(v),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','w',single(w),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','temp',single(temp),[1 1 step1]);
    ncwrite('F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide.nc','salt',single(salt),[1 1 step1]);
end
