clear
cd('C:\Users\Lenovo\Desktop\结构函数')

tide_lp_file     = 'I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_filt_order4.nc';
tide_stress_file = 'I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc';

x_rho = ncread(tide_lp_file,'x_rho');
y_rho = ncread(tide_lp_file,'y_rho');
depth = ncread(tide_lp_file,'depth');

time_idx = 49:647;
nt = numel(time_idx);

% 这里一定要按你的坐标单位来
dx = 500;   % 千万注意单位是m
dy = 500;   % 千万注意单位是m
numbin=50;

step1 = 0;

for it = 1:nt
    ii = time_idx(it);
    step1 = step1 + 1;

    disp(['step = ', num2str(step1), ' / ', num2str(nt)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    u_filt = ncread(tide_lp_file,'u_lp',[1 1 1 ii],[inf inf inf 1]);
    v_filt = ncread(tide_lp_file,'v_lp',[1 1 1 ii],[inf inf inf 1]);
    w_filt = ncread(tide_lp_file,'w_lp',[1 1 1 ii],[inf inf inf 1]);

    uu_filt = squeeze(ncread(tide_stress_file,'uu_filt',[1 1 2 ii],[inf inf 1 1]));
    vv_filt = squeeze(ncread(tide_stress_file,'vv_filt',[1 1 2 ii],[inf inf 1 1]));
    uv_filt = squeeze(ncread(tide_stress_file,'uv_filt',[1 1 2 ii],[inf inf 1 1]));
    uw_filt = squeeze(ncread(tide_stress_file,'uw_filt',[1 1 2 ii],[inf inf 1 1]));
    vw_filt = squeeze(ncread(tide_stress_file,'vw_filt',[1 1 2 ii],[inf inf 1 1]));

    u2 = squeeze(u_filt(:,:,2));
    v2 = squeeze(v_filt(:,:,2));
    w2 = squeeze(w_filt(:,:,2));

    tao_xx = uu_filt - u2 .* u2;
    tao_yy = vv_filt - v2 .* v2;
    tao_xy = uv_filt - u2 .* v2;
    tao_xz = uw_filt - u2 .* w2;
    tao_yz = vw_filt - v2 .* w2;

    [du_filtdx, du_filtdy] = model_gradient(x_rho, y_rho, u2);
    [dv_filtdx, dv_filtdy] = model_gradient(x_rho, y_rho, v2);

    du_filtdz = (squeeze(u_filt(:,:,3)) - squeeze(u_filt(:,:,1))) ./ (depth(3) - depth(1));
    dv_filtdz = (squeeze(v_filt(:,:,3)) - squeeze(v_filt(:,:,1))) ./ (depth(3) - depth(1));

    shear_xx = du_filtdx;
    shear_xy = du_filtdy + dv_filtdx;
    shear_yy = dv_filtdy;
    shear_xz = du_filtdz;
    shear_yz = dv_filtdz;

    [k_tide, Ck_xx, dk_tide1] = cospec2d_1d(tao_xx, shear_xx, dx, dy, numbin);
    [~,      Ck_xy, ~]       = cospec2d_1d(tao_xy, shear_xy, dx, dy, numbin);
    [~,      Ck_yy, ~]       = cospec2d_1d(tao_yy, shear_yy, dx, dy, numbin);
    [~,      Ck_xz, ~]       = cospec2d_1d(tao_xz, shear_xz, dx, dy, numbin);
    [~,      Ck_yz, ~]       = cospec2d_1d(tao_yz, shear_yz, dx, dy, numbin);

    Ck_all = Ck_xx + Ck_xy + Ck_yy + Ck_xz + Ck_yz;

    if step1 == 1
        k_iso = k_tide.*1000;%cpm-->cpkm
        dk_tide = dk_tide1;
        nbin = numel(k_iso);
        Ck_all_tide_sum = nan(nbin, nt);
    end

    Ck_all_tide_sum(:, step1) = Ck_all;
end

Ck_all_tide_mean = mean(Ck_all_tide_sum, 2, 'omitnan');


save('C:\Users\Lenovo\Desktop\结构函数\pai_time_cospec.mat','k_tide','Ck_all_tide_mean','dk_tide','Ck_all_tide_sum')
%
clear
cd('C:\Users\Lenovo\Desktop\结构函数')

notide_lp_file     = 'I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_filt_order4.nc';
notide_stress_file = 'I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc';

x_rho = ncread(notide_lp_file,'x_rho');
y_rho = ncread(notide_lp_file,'y_rho');
depth = ncread(notide_lp_file,'depth');

time_idx = 49:647;
nt = numel(time_idx);

% 这里一定要按你的坐标单位来
dx = 500;   % 千万注意单位是m
dy = 500;   % 千万注意单位是m
numbin=50;

step1 = 0;

for it = 1:nt
    ii = time_idx(it);
    step1 = step1 + 1;

    disp(['step = ', num2str(step1), ' / ', num2str(nt)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    u_filt = ncread(notide_lp_file,'u_lp',[1 1 1 ii],[inf inf inf 1]);
    v_filt = ncread(notide_lp_file,'v_lp',[1 1 1 ii],[inf inf inf 1]);
    w_filt = ncread(notide_lp_file,'w_lp',[1 1 1 ii],[inf inf inf 1]);

    uu_filt = squeeze(ncread(notide_stress_file,'uu_filt',[1 1 2 ii],[inf inf 1 1]));
    vv_filt = squeeze(ncread(notide_stress_file,'vv_filt',[1 1 2 ii],[inf inf 1 1]));
    uv_filt = squeeze(ncread(notide_stress_file,'uv_filt',[1 1 2 ii],[inf inf 1 1]));
    uw_filt = squeeze(ncread(notide_stress_file,'uw_filt',[1 1 2 ii],[inf inf 1 1]));
    vw_filt = squeeze(ncread(notide_stress_file,'vw_filt',[1 1 2 ii],[inf inf 1 1]));

    u2 = squeeze(u_filt(:,:,2));
    v2 = squeeze(v_filt(:,:,2));
    w2 = squeeze(w_filt(:,:,2));

    tao_xx = uu_filt - u2 .* u2;
    tao_yy = vv_filt - v2 .* v2;
    tao_xy = uv_filt - u2 .* v2;
    tao_xz = uw_filt - u2 .* w2;
    tao_yz = vw_filt - v2 .* w2;

    [du_filtdx, du_filtdy] = model_gradient(x_rho, y_rho, u2);
    [dv_filtdx, dv_filtdy] = model_gradient(x_rho, y_rho, v2);

    du_filtdz = (squeeze(u_filt(:,:,3)) - squeeze(u_filt(:,:,1))) ./ (depth(3) - depth(1));
    dv_filtdz = (squeeze(v_filt(:,:,3)) - squeeze(v_filt(:,:,1))) ./ (depth(3) - depth(1));

    shear_xx = du_filtdx;
    shear_xy = du_filtdy + dv_filtdx;
    shear_yy = dv_filtdy;
    shear_xz = du_filtdz;
    shear_yz = dv_filtdz;

    [k_notide, Ck_xx, dk_notide1] = cospec2d_1d(tao_xx, shear_xx, dx, dy, numbin);
    [~,      Ck_xy, ~]       = cospec2d_1d(tao_xy, shear_xy, dx, dy, numbin);
    [~,      Ck_yy, ~]       = cospec2d_1d(tao_yy, shear_yy, dx, dy, numbin);
    [~,      Ck_xz, ~]       = cospec2d_1d(tao_xz, shear_xz, dx, dy, numbin);
    [~,      Ck_yz, ~]       = cospec2d_1d(tao_yz, shear_yz, dx, dy, numbin);

    Ck_all = Ck_xx + Ck_xy + Ck_yy + Ck_xz + Ck_yz;

    if step1 == 1
        k_iso = k_notide.*1000;%cpm-->cpkm
        dk_notide = dk_notide1;
        nbin = numel(k_iso);
        Ck_all_notide_sum = nan(nbin, nt);
    end

    Ck_all_notide_sum(:, step1) = Ck_all;
end

Ck_all_notide_mean = mean(Ck_all_notide_sum, 2, 'omitnan');


save('C:\Users\Lenovo\Desktop\结构函数\pai_time_cospec.mat','-append','k_notide','Ck_all_notide_mean','dk_notide','Ck_all_notide_sum')



clear
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask'); % 使用原ncread函数
;
tide_file   = 'I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_wb_filt_order4.nc';
notide_file = 'I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_wb_filt_order4.nc';

dx = 500;   % 千万注意单位是m
dy = 500;   % 千万注意单位是m

numbin = 50;
time_idx = 49:647;
nt = numel(time_idx);

k_iso = [];
dk = [];

wb_tide_sum   = [];
wb_notide_sum = [];

check_tide_phys   = nan(nt,1);
check_tide_spec   = nan(nt,1);
check_notide_phys = nan(nt,1);
check_notide_spec = nan(nt,1);

for it = 1:nt
    ii = time_idx(it);

    disp(['step = ', num2str(it), ' / ', num2str(nt), ', file index = ', num2str(ii)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    %---------------- read data ----------------
    b_tide   = ncread(tide_file,   'b_lp', [1 1 ii], [inf inf 1]);
    w_tide   = ncread(tide_file,   'w_lp', [1 1 ii], [inf inf 1]);
    b_notide = ncread(notide_file, 'b_lp', [1 1 ii], [inf inf 1]);
    w_notide = ncread(notide_file, 'w_lp', [1 1 ii], [inf inf 1]);

    %---------------- directly compute 1D isotropic cospectrum ----------------
    [k_now, Ck_tide, dk_now, check_tide_phys(it), check_tide_spec(it)] = ...
        cospec2d_1d(w_tide, b_tide, dx, dy, numbin);

    [~, Ck_notide, ~, check_notide_phys(it), check_notide_spec(it)] = ...
        cospec2d_1d(w_notide, b_notide, dx, dy, numbin);

    %---------------- allocate on first step ----------------
    if it == 1
        k_iso = k_now.*1000;%cpm-->cpkm
        dk_wb    = dk_now;

        nbin = numel(k_iso);
        wb_tide_sum   = nan(nbin, nt);
        wb_notide_sum = nan(nbin, nt);
    end

    %---------------- save ----------------
    wb_tide_sum(:, it)   = Ck_tide;
    wb_notide_sum(:, it) = Ck_notide;
end

%---------------- time average ----------------
C1d_tide_sum_timeavg   = mean(wb_tide_sum,   2, 'omitnan');
C1d_notide_sum_timeavg = mean(wb_notide_sum, 2, 'omitnan');

%---------------- consistency check ----------------
check_tide_phys_mean   = mean(check_tide_phys,   'omitnan');
check_tide_spec_mean   = mean(check_tide_spec,   'omitnan');
check_notide_phys_mean = mean(check_notide_phys, 'omitnan');
check_notide_spec_mean = mean(check_notide_spec, 'omitnan');

fprintf('\n');
fprintf('TIDE   : mean(check_phys) = %.6e, mean(check_spec) = %.6e\n', ...
    check_tide_phys_mean, check_tide_spec_mean);
fprintf('NOTIDE : mean(check_phys) = %.6e, mean(check_spec) = %.6e\n', ...
    check_notide_phys_mean, check_notide_spec_mean);


wb_spec_tide=C1d_tide_sum_timeavg;
wb_spec_notide=C1d_notide_sum_timeavg;

k_center=k_now;
save('wb_cospec.mat','k_center','wb_spec_tide','wb_spec_notide','dk_wb','wb_tide_sum','wb_notide_sum');