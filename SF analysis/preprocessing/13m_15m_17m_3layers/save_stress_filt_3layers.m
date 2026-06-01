clear
filename='I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_origin.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');
lon_rho=ncread(filename,'lon_rho');lat_rho=ncread(filename,'lat_rho');
depth=ncread(filename,'depth');

% 
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','lat_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','depth','Dimensions',{'depth',length(depth)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','depth',single(depth),[1]);
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','uu_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','uv_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','vv_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','uw_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc','vw_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');

for kk = 1:length(depth)
    disp(kk);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    for ii = 1:5
        for jj = 1:5
            ist = 240*(ii-1) + 1;
            jst = 240*(jj-1) + 1;
            
            u_use = squeeze(ncread(filename, 'u', [ist jst kk 1], [240 240 1 inf]));
            v_use = squeeze(ncread(filename, 'v', [ist jst kk 1], [240 240 1 inf]));
            w_use = squeeze(ncread(filename, 'w', [ist jst kk 1], [240 240 1 inf]));
            nt = size(u_use, 3);
            
            % 预分配（用 NaN，这样陆地点自动就是 NaN）
            uu_filt = NaN(240, 240, nt, 'single');
            uv_filt = NaN(240, 240, nt, 'single');
            vv_filt = NaN(240, 240, nt, 'single');
            uw_filt = NaN(240, 240, nt, 'single');
            vw_filt = NaN(240, 240, nt, 'single');

            for mm = 1:240
                parfor nn = 1:240
                    u_ts = squeeze(u_use(mm, nn, :));
                    v_ts = squeeze(v_use(mm, nn, :));
                    w_ts = squeeze(w_use(mm, nn, :));
                    
                    if all(isnan(u_ts))
                        % 陆地点，跳过，保持 NaN
                        uu_filt(mm, nn, :) = NaN;
                        uv_filt(mm, nn, :) = NaN;
                        vv_filt(mm, nn, :) = NaN;
                        uw_filt(mm, nn, :) = NaN;
                        vw_filt(mm, nn, :) = NaN;
                    else
                        uu_filt(mm, nn, :) = filt1('lp', u_ts.*u_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        uv_filt(mm, nn, :) = filt1('lp', u_ts.*v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        vv_filt(mm, nn, :) = filt1('lp', v_ts.*v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        uw_filt(mm, nn, :) = filt1('lp', w_ts.*u_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        vw_filt(mm, nn, :) = filt1('lp', w_ts.*v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                    end
                end
            end
            

            % reshape 回 4D 再写入
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc', 'uu_filt', reshape(uu_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc', 'uv_filt', reshape(uv_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc', 'vv_filt', reshape(vv_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc', 'uw_filt', reshape(uw_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc', 'vw_filt', reshape(vw_filt, [240,240,1,nt]), [ist jst kk 1]);
        end
    end
end



%%%

clear
filename='I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_origin.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');
lon_rho=ncread(filename,'lon_rho');lat_rho=ncread(filename,'lat_rho');
depth=ncread(filename,'depth');

% 
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','lat_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','depth','Dimensions',{'depth',length(depth)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','depth',single(depth),[1]);
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','uu_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','uv_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','vv_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','uw_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc','vw_filt','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'depth',length(depth),'t',695},'Datatype','single');

for kk = 1:length(depth)
    disp(kk);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    for ii = 1:5
        for jj = 1:5
            ist = 240*(ii-1) + 1;
            jst = 240*(jj-1) + 1;
            
            u_use = squeeze(ncread(filename, 'u', [ist jst kk 1], [240 240 1 inf]));
            v_use = squeeze(ncread(filename, 'v', [ist jst kk 1], [240 240 1 inf]));
            w_use = squeeze(ncread(filename, 'w', [ist jst kk 1], [240 240 1 inf]));
            nt = size(u_use, 3);
            
            % 预分配（用 NaN，这样陆地点自动就是 NaN）
            uu_filt = NaN(240, 240, nt, 'single');
            uv_filt = NaN(240, 240, nt, 'single');
            vv_filt = NaN(240, 240, nt, 'single');
            uw_filt = NaN(240, 240, nt, 'single');
            vw_filt = NaN(240, 240, nt, 'single');

            for mm = 1:240
                parfor nn = 1:240
                    u_ts = squeeze(u_use(mm, nn, :));
                    v_ts = squeeze(v_use(mm, nn, :));
                    w_ts = squeeze(w_use(mm, nn, :));
                    
                    if all(isnan(u_ts))
                        % 陆地点，跳过，保持 NaN
                        uu_filt(mm, nn, :) = NaN;
                        uv_filt(mm, nn, :) = NaN;
                        vv_filt(mm, nn, :) = NaN;
                        uw_filt(mm, nn, :) = NaN;
                        vw_filt(mm, nn, :) = NaN;
                    else
                        uu_filt(mm, nn, :) = filt1('lp', u_ts.*u_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        uv_filt(mm, nn, :) = filt1('lp', u_ts.*v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        vv_filt(mm, nn, :) = filt1('lp', v_ts.*v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        uw_filt(mm, nn, :) = filt1('lp', w_ts.*u_ts, 'fs', 1, 'fc', 1/48,'order',4);
                        vw_filt(mm, nn, :) = filt1('lp', w_ts.*v_ts, 'fs', 1, 'fc', 1/48,'order',4);
                    end
                end
            end
            

            % reshape 回 4D 再写入
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc', 'uu_filt', reshape(uu_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc', 'uv_filt', reshape(uv_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc', 'vv_filt', reshape(vv_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc', 'uw_filt', reshape(uw_filt, [240,240,1,nt]), [ist jst kk 1]);
            ncwrite('I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc', 'vw_filt', reshape(vw_filt, [240,240,1,nt]), [ist jst kk 1]);
        end
    end
end

