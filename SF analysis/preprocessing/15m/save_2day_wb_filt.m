%% 2 day filt w b
clear
filename='I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_origin.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');
lon_rho=ncread(filename,'lon_rho');lat_rho=ncread(filename,'lat_rho');

% 
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','lat_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','w_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','b_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','w_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc','b_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');

for ii = 1:5
    disp(ii);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    for jj = 1:5
        ist = 240*(ii-1) + 1;
        jst = 240*(jj-1) + 1;
        
        b_use = squeeze(ncread(filename, 'b', [ist jst 1], [240 240 inf]));
        w_use = squeeze(ncread(filename, 'w', [ist jst 1], [240 240 inf]));
        nt = size(b_use, 3);
        
        % 预分配（用 NaN，这样陆地点自动就是 NaN）
        w_lp = NaN(240, 240, nt, 'single');
        b_lp = NaN(240, 240, nt, 'single');
        
        for mm = 1:240
            parfor nn = 1:240
                w_ts = squeeze(w_use(mm, nn, :));
                b_ts = squeeze(b_use(mm, nn, :));
                
                if all(isnan(w_ts))
                    % 陆地点，跳过，保持 NaN
                    w_lp(mm, nn, :) = NaN;
                    b_lp(mm, nn, :) = NaN;
                else
                    w_lp(mm, nn, :) = filt1('lp', w_ts, 'fs', 1, 'fc', 1/48,'order',4);
                    b_lp(mm, nn, :) = filt1('lp', b_ts, 'fs', 1, 'fc', 1/48,'order',4);
                end
            end
        end
        
        w_hp = single(w_use) - w_lp;
        b_hp = single(b_use) - b_lp;
        
        % reshape 回 3D 再写入
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc', 'w_lp', reshape(w_lp, [240,240,nt]), [ist jst 1]);
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc', 'b_lp', reshape(b_lp, [240,240,nt]), [ist jst 1]);
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc', 'w_hp', reshape(w_hp, [240,240,nt]), [ist jst 1]);
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc', 'b_hp', reshape(b_hp, [240,240,nt]), [ist jst 1]);
    end
end

%%%%
clear
filename='I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_origin.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');
lon_rho=ncread(filename,'lon_rho');lat_rho=ncread(filename,'lat_rho');

% 
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','lon_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','lat_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','lon_rho',single(lon_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','lat_rho',single(lat_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','w_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','b_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','w_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc','b_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');

for ii = 1:5
    disp(ii);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    for jj = 1:5
        ist = 240*(ii-1) + 1;
        jst = 240*(jj-1) + 1;
        
        b_use = squeeze(ncread(filename, 'b', [ist jst 1], [240 240 inf]));
        w_use = squeeze(ncread(filename, 'w', [ist jst 1], [240 240 inf]));
        nt = size(b_use, 3);
        
        % 预分配（用 NaN，这样陆地点自动就是 NaN）
        w_lp = NaN(240, 240, nt, 'single');
        b_lp = NaN(240, 240, nt, 'single');
        
        for mm = 1:240
            parfor nn = 1:240
                w_ts = squeeze(w_use(mm, nn, :));
                b_ts = squeeze(b_use(mm, nn, :));
                
                if all(isnan(w_ts))
                    % 陆地点，跳过，保持 NaN
                    w_lp(mm, nn, :) = NaN;
                    b_lp(mm, nn, :) = NaN;
                else
                    w_lp(mm, nn, :) = filt1('lp', w_ts, 'fs', 1, 'fc', 1/48,'order',4);
                    b_lp(mm, nn, :) = filt1('lp', b_ts, 'fs', 1, 'fc', 1/48,'order',4);
                end
            end
        end
        
        w_hp = single(w_use) - w_lp;
        b_hp = single(b_use) - b_lp;
        
        % reshape 回 3D 再写入
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc', 'w_lp', reshape(w_lp, [240,240,nt]), [ist jst 1]);
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc', 'b_lp', reshape(b_lp, [240,240,nt]), [ist jst 1]);
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc', 'w_hp', reshape(w_hp, [240,240,nt]), [ist jst 1]);
        ncwrite('I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc', 'b_hp', reshape(b_hp, [240,240,nt]), [ist jst 1]);
    end
end
