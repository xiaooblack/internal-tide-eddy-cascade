clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数

r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc','r')';
sum_D3_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_D3_lp', [1 1 240], [inf inf 408]));
sum_DLL_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_DLL_lp', [1 1 240], [inf inf 408]));
sum_DTT_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_DTT_lp', [1 1 240], [inf inf 408]));
counts_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'counts_lp', [1 1 240], [inf inf 408]));

depth=[-210 -195 -180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 -5];
S3_tide_lp = sum(sum_D3_tide_lp, 3) ./ sum(counts_tide_lp, 3);
S2_tide_lp=sum(sum_DLL_tide_lp, 3) ./ sum(counts_tide_lp, 3)+sum(sum_DTT_tide_lp, 3) ./ sum(counts_tide_lp, 3);

% --- 数据参与拟合的范围（尽量用全部数据）---
r_min_data = 1;     % km  所有 >= 1 km 的数据都参与拟合
r_max_data = 500;   % km  所有 <= 100 km 的数据都参与拟合

% --- 你关心的物理波数范围（对应目标尺度 5~90 km）---
r_min_target = 1;   % km  → k_max = 1/5 = 0.2  1/km
r_max_target = 99;  % km  → k_min = 1/90 ≈ 0.011 1/km

smooth_win  = 8;
lambda_reg  = 1e-2;
Nk          = 100;

for ii=1:15
    [k_tide_lp,R_tide_lp(:,ii),F_k_tide_lp(:,ii),Pidk_tide_lp(:,ii),S3_recon_tide_lp(:,ii)]=calc_epsilon_Fk(r,S3_tide_lp(:,ii),r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);

end

r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc','r')';
sum_D3_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_D3_lp', [1 1 240], [inf inf 408]));
sum_DLL_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_DLL_lp', [1 1 240], [inf inf 408]));
sum_DTT_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_DTT_lp', [1 1 240], [inf inf 408]));
counts_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'counts_lp', [1 1 240], [inf inf 408]));

depth=[-210 -195 -180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 -5];
S3_notide_lp = sum(sum_D3_notide_lp, 3) ./ sum(counts_notide_lp, 3);
S2_notide_lp=sum(sum_DLL_notide_lp, 3) ./ sum(counts_notide_lp, 3)+sum(sum_DTT_notide_lp, 3) ./ sum(counts_notide_lp, 3);

% --- 数据参与拟合的范围（尽量用全部数据）---
r_min_data = 1;     % km  所有 >= 1 km 的数据都参与拟合
r_max_data = 500;   % km  所有 <= 100 km 的数据都参与拟合

% --- 你关心的物理波数范围（对应目标尺度 5~90 km）---
r_min_target = 1;   % km  → k_max = 1/5 = 0.2  1/km
r_max_target = 99;  % km  → k_min = 1/90 ≈ 0.011 1/km

smooth_win  = 8;
lambda_reg  = 1e-2;
Nk          = 100;

for ii=1:15
    [k_notide_lp,R_notide_lp(:,ii),F_k_notide_lp(:,ii),Pidk_notide_lp(:,ii),S3_recon_notide_lp(:,ii)]=calc_epsilon_Fk(r,S3_notide_lp(:,ii),r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);

end
%% figure S11

[depthh,rr]=meshgrid(depth,r);
figure
set(gcf, 'unit', 'centimeters', 'position', [1.5 2 27 10]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 9 6]);
pcolor(rr,depthh,S3_tide_lp./repmat(r',[1 15]));shading interp
caxis([-2e-5 2e-5])
colortable = textread('NCV_blue_red.txt');
colormap(colortable);
set(gca,'xscale','log');
xlim([1 100])
ylim([-200 0])
set(gca,'layer','top');
set(gca,'fontsize',12,'fontname','Arial','linewidth',1);
ylabel('depth (m)','fontsize',13,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',13,'fontname','Arial')
tit=title('\bf(a) \rmVertical structure of velocity SF3 (Tide lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.45 6.1 0])

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [14 2 9 6]);
pcolor(rr,depthh,S3_notide_lp./repmat(r',[1 15]));shading interp
caxis([-2e-5 2e-5])
colortable = textread('NCV_blue_red.txt');
colormap(colortable);
set(gca,'xscale','log');
xlim([1 100])
ylim([-200 0])
set(gca,'layer','top');
set(gca,'fontsize',12,'fontname','Arial','linewidth',1);
ylabel('depth (m)','fontsize',13,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',13,'fontname','Arial')
tit=title('\bf(b) \rmVertical structure of velocity SF3 (Notide lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.7 6.1 0])
hhh = colorbar;
set(hhh, 'unit', 'centimeters','Position', [24.5 2.5 0.5 5]);
tith=title(hhh,'D3(r)/r','fontsize',11,'fontname','Arial')
set(tith,'Units','points','Position',[-5 145.5 0])
print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S11.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi


%% figure s11
[depthh,kk_tide_lp]=meshgrid(depth,k_tide_lp);
[depthh,kk_notide_lp]=meshgrid(depth,k_notide_lp);

figure
set(gcf, 'unit', 'centimeters', 'position', [1.5 2 26 19]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 11 9 6]);
pcolor(kk_tide_lp,depthh,F_k_tide_lp);shading interp
caxis([-1.5e-5 1.5e-5])
colortable = textread('NCV_blue_red.txt');
colormap(colortable);
set(gca,'xscale','log');
xlim([1e-2 1])
ylim([-200 0])
set(gca,'layer','top');
set(gca,'fontsize',12,'fontname','Arial','linewidth',1);
ylabel('depth (m)','fontsize',13,'fontname','Arial');
xlabel('\itk \rm(1/km)','fontsize',13,'fontname','Arial')
tit=title('\bf(a) \rmVertical structure of \itF(k)\rm (Tide lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.65 6.1 0])

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 11 9 6]);
pcolor(kk_tide_lp,depthh,F_k_notide_lp);shading interp
caxis([-1.5e-5 1.5e-5])
colortable = textread('NCV_blue_red.txt');
colormap(colortable);
set(gca,'xscale','log');
xlim([1e-2 1])
ylim([-200 0])
set(gca,'layer','top');
set(gca,'fontsize',12,'fontname','Arial','linewidth',1);
ylabel('depth (m)','fontsize',13,'fontname','Arial');
xlabel('\itk \rm(1/km)','fontsize',13,'fontname','Arial')
tit=title('\bf(b) \rmVertical structure of \itF(k)\rm (Notide lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.85 6.1 0])
hhh = colorbar;
set(hhh, 'unit', 'centimeters','Position', [23.5 11.5 0.5 5]);
tith=title(hhh,'\itF(k)','fontsize',11,'fontname','Arial')
set(tith,'Units','points','Position',[4 144 0])

AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 9 6]);
pcolor(kk_tide_lp,depthh,Pidk_tide_lp.*kk_tide_lp);shading interp
caxis([-1.5e-5 1.5e-5])
colortable = textread('NCV_blue_red.txt');
colormap(colortable);
set(gca,'xscale','log');
xlim([1e-2 1])
ylim([-200 0])
set(gca,'layer','top');
set(gca,'fontsize',12,'fontname','Arial','linewidth',1);
ylabel('depth (m)','fontsize',13,'fontname','Arial');
xlabel('\itk \rm(1/km)','fontsize',13,'fontname','Arial')
tit=title('\bf(c) \rmVertical structure of \it\epsilon_i·k\rm (Tide lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.6 6.1 0])

AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 2 9 6]);
pcolor(kk_tide_lp,depthh,Pidk_notide_lp.*kk_notide_lp);shading interp
caxis([-1.5e-5 1.5e-5])
colortable = textread('NCV_blue_red.txt');
colormap(colortable);
set(gca,'xscale','log');
xlim([1e-2 1])
ylim([-200 0])
set(gca,'layer','top');
set(gca,'fontsize',12,'fontname','Arial','linewidth',1);
ylabel('depth (m)','fontsize',13,'fontname','Arial');
xlabel('\itk \rm(1/km)','fontsize',13,'fontname','Arial')
tit=title('\bf(d) \rmVertical structure of \it\epsilon_i·k\rm (Notide lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.8 6.1 0])
hhh = colorbar;
set(hhh, 'unit', 'centimeters','Position', [23.5 2.5 0.5 5]);
tith=title(hhh,'\it\epsilon_i·k','fontsize',11,'fontname','Arial')
set(tith,'Units','points','Position',[4 144 0])
print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S12.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi
