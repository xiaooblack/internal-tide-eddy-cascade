%%
clear
filename='I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');
lon_rho=ncread(filename,'lon_rho');lat_rho=ncread(filename,'lat_rho');
depth=ncread(filename,'depth');
% 
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','x_rho',single(x_rho),[1 1]);
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','y_rho',single(y_rho),[1 1]);
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lon_rho',single(lon_rho),[1 1]);
% ncwrite('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','lat_rho',single(lat_rho),[1 1]);
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','u_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','v_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','u_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
% nccreate('I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc','v_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');

for kk = 11:length(depth)
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
%%
clear
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
filename_tide='I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_tide_filt_order4.nc';
filename_notide='I:\NCS_wu_0.5km_15layer_all\NCS_wu_0.5km_15layers_notide_filt_order4.nc';
lon_rho= ncread(filename_tide,'lon_rho');   lat_rho= ncread(filename_tide,'lat_rho');   % X/Y rho  :km
x_rho= ncread(filename_tide,'x_rho');   y_rho= ncread(filename_tide,'y_rho');   % X/Y rho  :km
f=sw_f(18)

u_lp_tide=ncread(filename_tide,'u_lp',[1 1 5 300],[inf inf 1 1]);
v_lp_tide=ncread(filename_tide,'v_lp',[1 1 5 300],[inf inf 1 1]);

u_hp_tide=ncread(filename_tide,'u_hp',[1 1 5 300],[inf inf 1 1]);
v_hp_tide=ncread(filename_tide,'v_hp',[1 1 5 300],[inf inf 1 1]);

u_lp_notide=ncread(filename_notide,'u_lp',[1 1 5 300],[inf inf 1 1]);
v_lp_notide=ncread(filename_notide,'v_lp',[1 1 5 300],[inf inf 1 1]);

u_hp_notide=ncread(filename_notide,'u_hp',[1 1 5 300],[inf inf 1 1]);
v_hp_notide=ncread(filename_notide,'v_hp',[1 1 5 300],[inf inf 1 1]);


[dudx,dudy]=model_gradient(x_rho,y_rho,u_lp_tide);
[dvdx,dvdy]=model_gradient(x_rho,y_rho,v_lp_tide);
zeta_lp_tide=(dvdx-dudy);
div_lp_tide=(dudx+dvdy);

[dudx,dudy]=model_gradient(x_rho,y_rho,u_hp_tide);
[dvdx,dvdy]=model_gradient(x_rho,y_rho,v_hp_tide);
zeta_hp_tide=(dvdx-dudy);
div_hp_tide=(dudx+dvdy);

[dudx,dudy]=model_gradient(x_rho,y_rho,u_lp_notide);
[dvdx,dvdy]=model_gradient(x_rho,y_rho,v_lp_notide);
zeta_lp_notide=(dvdx-dudy);
div_lp_notide=(dudx+dvdy);

[dudx,dudy]=model_gradient(x_rho,y_rho,u_hp_notide);
[dvdx,dvdy]=model_gradient(x_rho,y_rho,v_hp_notide);
zeta_hp_notide=(dvdx-dudy);
div_hp_notide=(dudx+dvdy);

figure
pcolor(lon_rho,lat_rho,zeta_lp_tide./f);shading interp
caxis([-1 1])
colortable = textread('NCV_blu_red.txt');
colormap(colortable);
hhh = colorbar;
set(hhh,'fontsize',11);
title('tide lp zeta')

figure
pcolor(lon_rho,lat_rho,zeta_hp_tide./f);shading interp
caxis([-1 1])
colortable = textread('NCV_blu_red.txt');
colormap(colortable);
hhh = colorbar;
set(hhh,'fontsize',11);
title('tide hp zeta')

figure
pcolor(lon_rho,lat_rho,zeta_lp_notide./f);shading interp
caxis([-1 1])
colortable = textread('NCV_blu_red.txt');
colormap(colortable);
hhh = colorbar;
set(hhh,'fontsize',11);
title('notide lp zeta')

figure
pcolor(lon_rho,lat_rho,zeta_hp_notide./f);shading interp
caxis([-1 1])
colortable = textread('NCV_blu_red.txt');
colormap(colortable);
hhh = colorbar;
set(hhh,'fontsize',11);
title('notide hp zeta')