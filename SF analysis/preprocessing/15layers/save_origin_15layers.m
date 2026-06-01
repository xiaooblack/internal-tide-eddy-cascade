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

rpath1='F:\NCS_wu_0.5km_tide_uv';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-210 -195 -180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 -5];

% create nc file
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'x_rho',size(x_rho,2)},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','depth','Dimensions',{'depth',length(z_need)},'Datatype','single');
% 
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','x_rho',single(x_rho),[1 1]);
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','y_rho',single(y_rho),[1 1]);
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','lon_rho',single(lon_rho),[1 1]);
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','lat_rho',single(lat_rho),[1 1]);
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','depth',single(z_need),[1]);
% 
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','u','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(z_need),'t',m1},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','v','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(z_need),'t',m1},'Datatype','single');

step1=460;
for ii=461:m1
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    u_u=ncread(filename,'u');
    v_v=ncread(filename,'v');
    for iii=1:size(u_u,3)
        u_m1(:,:,iii)=v2rho_2d(u_u(:,:,iii));
        v_m1(:,:,iii)=u2rho_2d(v_v(:,:,iii));
    end
    u_m=u_m1(449:1648,100:1299,:);v_m=v_m1(449:1648,100:1299,:);
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(:,449:1648,100:1299);
    z_rho=permute(z_rho,[2 3 1]);

    u=interp_model_z(u_m,x_rho,y_rho,z_rho,z_need);
    v=interp_model_z(v_m,x_rho,y_rho,z_rho,z_need);

    % write nc file
    ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','u',single(u),[1 1 1 step1]);
    ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc','v',single(v),[1 1 1 step1]);

end



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

rpath1='F:\NCS_wu_0.5km_notide_uv';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-105 -90 -75 -60 -45 -30 -15 -5];

% create nc file
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'x_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','depth','Dimensions',{'depth',length(z_need)},'Datatype','single');

ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','y_rho',single(y_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','lat_rho',single(lat_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','depth',single(z_need),[1]);

nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','u','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(z_need),'t',m1},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','v','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(z_need),'t',m1},'Datatype','single');

step1=0;
for ii=1:m1
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    u_u=ncread(filename,'u');
    v_v=ncread(filename,'v');
    for iii=1:size(u_u,3)
        u_m1(:,:,iii)=v2rho_2d(u_u(:,:,iii));
        v_m1(:,:,iii)=u2rho_2d(v_v(:,:,iii));
    end
    u_m=u_m1(449:1648,100:1299,:);v_m=v_m1(449:1648,100:1299,:);
    clear u_m1 v_m1
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(:,449:1648,100:1299);
    z_rho=permute(z_rho,[2 3 1]);

    u=interp_model_z(u_m,x_rho,y_rho,z_rho,z_need);
    v=interp_model_z(v_m,x_rho,y_rho,z_rho,z_need);

    % write nc file
    ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','u',single(u),[1 1 1 step1]);
    ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc','v',single(v),[1 1 1 step1]);

end

%%
clear
filename='I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');
lon_rho=ncread(filename,'lon_rho');lat_rho=ncread(filename,'lat_rho');
depth=ncread(filename,'depth');

nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lat_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','u_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','v_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','u_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','v_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');

for kk = 1:length(depth)
    for ii = 1:5
        disp(ii);
        disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        for jj = 1:5
            ist = 240*(ii-1) + 1;
            jst = 240*(jj-1) + 1;
            
            u_use = squeeze(ncread(filename, 'u', [ist jst kk 1], [240 240 1 inf]));
            v_use = squeeze(ncread(filename, 'v', [ist jst kk 1], [240 240 1 inf]));
            nt = size(u_use, 3);
            
            % 预分配（用 NaN，这样陆地点自动就是 NaN）
            u_lp = NaN(240, 240, nt, 'single');
            v_lp = NaN(240, 240, nt, 'single');
            
            for mm = 1:240
                parfor nn = 1:240
                    u_ts = squeeze(u_use(mm, nn, :));
                    v_ts = squeeze(v_use(mm, nn, :));
                    
                    if all(isnan(u_ts))
                        % 陆地点，跳过，保持 NaN
                        u_lp(mm, nn, :) = NaN;
                        v_lp(mm, nn, :) = NaN;
                    else
                        u_lp(mm, nn, :) = filt1('lp', u_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        v_lp(mm, nn, :) = filt1('lp', v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                    end
                end
            end
            
            u_hp = single(u_use) - u_lp;
            v_hp = single(v_use) - v_lp;
            
            % reshape 回 4D 再写入
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc', 'u_lp', reshape(u_lp, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc', 'v_lp', reshape(v_lp, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc', 'u_hp', reshape(u_hp, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc', 'v_hp', reshape(v_hp, [240,240,1,nt]), [ist jst kk 1]);
        end
    end
end

%%

clear
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
filename='I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');
lon_rho=ncread(filename,'lon_rho');lat_rho=ncread(filename,'lat_rho');
depth=ncread(filename,'depth');

nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','lat_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','u_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','v_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','u_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc','v_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');

for kk = 1:length(depth)
    for ii = 1:5
        disp(ii);
        disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        for jj = 1:5
            ist = 240*(ii-1) + 1;
            jst = 240*(jj-1) + 1;
            
            u_use = squeeze(ncread(filename, 'u', [ist jst kk 1], [240 240 1 inf]));
            v_use = squeeze(ncread(filename, 'v', [ist jst kk 1], [240 240 1 inf]));
            nt = size(u_use, 3);
            
            % 预分配（用 NaN，这样陆地点自动就是 NaN）
            u_lp = NaN(240, 240, nt, 'single');
            v_lp = NaN(240, 240, nt, 'single');
            
            for mm = 1:240
                parfor nn = 1:240
                    u_ts = squeeze(u_use(mm, nn, :));
                    v_ts = squeeze(v_use(mm, nn, :));
                    
                    if all(isnan(u_ts))
                        % 陆地点，跳过，保持 NaN
                        u_lp(mm, nn, :) = NaN;
                        v_lp(mm, nn, :) = NaN;
                    else
                        u_lp(mm, nn, :) = filt1('lp', u_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        v_lp(mm, nn, :) = filt1('lp', v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                    end
                end
            end
            
            u_hp = single(u_use) - u_lp;
            v_hp = single(v_use) - v_lp;
            
            % reshape 回 4D 再写入
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc', 'u_lp', reshape(u_lp, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc', 'v_lp', reshape(v_lp, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc', 'u_hp', reshape(u_hp, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc', 'v_hp', reshape(v_hp, [240,240,1,nt]), [ist jst kk 1]);
        end
    end
end