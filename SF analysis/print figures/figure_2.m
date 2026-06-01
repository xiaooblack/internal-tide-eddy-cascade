    clear
    cd 'C:\Users\Lenovo\Desktop\结构函数'
    warning off all;  
    %read grid
    rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
    x_rho=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt.nc','x_rho');
    y_rho=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt.nc','y_rho');
    load('I:\NCS_wu_0.5km_2day_filt\lonlat_rho_origin.mat')
    u_lp_tide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','u_lp',[1 1 590],[inf inf 1]);
    v_lp_tide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','v_lp',[1 1 590],[inf inf 1]);
    
    u_hp_tide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','u_hp',[1 1 590],[inf inf 1]);
    v_hp_tide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','v_hp',[1 1 590],[inf inf 1]);
    
    u_lp_notide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','u_lp',[1 1 590],[inf inf 1]);
    v_lp_notide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','v_lp',[1 1 590],[inf inf 1]);
    
    u_hp_notide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','u_hp',[1 1 590],[inf inf 1]);
    v_hp_notide=ncread('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','v_hp',[1 1 590],[inf inf 1]);
    
    [dudx,dudy]=model_gradient(x_rho,y_rho,u_lp_tide);
    [dvdx,dvdy]=model_gradient(x_rho,y_rho,v_lp_tide);
    zeta_lp_tide=(dvdx-dudy)./sw_f(18);
    div_lp_tide=(dudx+dvdy)./sw_f(18);
    
    [dudx,dudy]=model_gradient(x_rho,y_rho,u_lp_notide);
    [dvdx,dvdy]=model_gradient(x_rho,y_rho,v_lp_notide);
    zeta_lp_notide=(dvdx-dudy)./sw_f(18);
    div_lp_notide=(dudx+dvdy)./sw_f(18);
    
    zeta_lp_tide(zeta_lp_tide>10)=nan;
    zeta_lp_tide(zeta_lp_tide<-10)=nan;
    zeta_lp_tide=inpaint_nans(zeta_lp_tide,1);
    %%% SF3 & F(k)
    % ===================== 1. 核心可调参数 =====================
    % --- 数据参与拟合的范围（尽量用全部数据）---
    r_min_data = 1;     % km  所有 >= 1 km 的数据都参与拟合
    r_max_data = 500;   % km  所有 <= 100 km 的数据都参与拟合
    
    % --- 你关心的物理波数范围（对应目标尺度 5~90 km）---
    r_min_target = 1;   % km  → k_max = 1/5 = 0.2  1/km
    r_max_target = 99;  % km  → k_min = 1/90 ≈ 0.011 1/km
    
    smooth_win  = 8;
    lambda_reg  = 1e-2;
    Nk          = 100;
    
    % ===================== 2. 读取 & 预处理数据 =====================
    r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc','r');
    sum_D3_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
        'sum_D3_lp', [1 14 49], [inf 1 599]));
    sum_DLL_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
        'sum_DLL_lp', [1 14 49], [inf 1 599]));
    sum_DTT_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
        'sum_DTT_lp', [1 14 49], [inf 1 599]));
    sum_D3_zeta_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
        'sum_D3_zeta_lp', [1 14 49], [inf 1 599]));
    sum_D2_zeta_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
        'sum_D2_zeta_lp', [1 14 49], [inf 1 599]));
    counts_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
        'counts_lp', [1 14 49], [inf 1 599]));
    
    S3_tide_lp = sum(sum_D3_tide_lp, 2) ./ sum(counts_tide_lp, 2);
    S3_zeta_tide_lp = sum(sum_D3_zeta_tide_lp, 2) ./ sum(counts_tide_lp, 2);
    S2_zeta_tide_lp = sum(sum_D2_zeta_tide_lp, 2) ./ sum(counts_tide_lp, 2);
    S2_tide_lp=sum(sum_DLL_tide_lp, 2) ./ sum(counts_tide_lp, 2)+sum(sum_DTT_tide_lp, 2) ./ sum(counts_tide_lp, 2);
    [k_tide_lp,R_tide_lp,F_k_tide_lp,Pidk_tide_lp,S3_recon_tide_lp]=calc_epsilon_Fk(r,S3_tide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
    [k_zeta_tide_lp,R_zeta_tide_lp,F_k_zeta_tide_lp,Pidk_zeta_tide_lp,S3_zeta_recon_tide_lp]=calc_epsilon_Fk_zeta(r,S3_zeta_tide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
    
    r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc','r');
    sum_D3_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
        'sum_D3_lp', [1 14 49], [inf 1 599]));
    sum_DLL_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
        'sum_DLL_lp', [1 14 49], [inf 1 599]));
    sum_DTT_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
        'sum_DTT_lp', [1 14 49], [inf 1 599]));
    sum_D3_zeta_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
        'sum_D3_zeta_lp', [1 14 49], [inf 1 599]));
    sum_D2_zeta_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
        'sum_D2_zeta_lp', [1 14 49], [inf 1 599]));
    counts_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
        'counts_lp', [1 14 49], [inf 1 599]));
    
    S3_notide_lp = sum(sum_D3_notide_lp, 2) ./ sum(counts_notide_lp, 2);
    S3_zeta_notide_lp = sum(sum_D3_zeta_notide_lp, 2) ./ sum(counts_notide_lp, 2);
    S2_zeta_notide_lp = sum(sum_D2_zeta_notide_lp, 2) ./ sum(counts_notide_lp, 2);
    S2_notide_lp=sum(sum_DLL_notide_lp, 2) ./ sum(counts_notide_lp, 2)+sum(sum_DTT_notide_lp, 2) ./ sum(counts_notide_lp, 2);
    [k_notide_lp,R_notide_lp,F_k_notide_lp,Pidk_notide_lp,S3_recon_notide_lp]=calc_epsilon_Fk(r,S3_notide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
    [k_zeta_notide_lp,R_zeta_notide_lp,F_k_zeta_notide_lp,Pidk_zeta_notide_lp,S3_zeta_recon_notide_lp]=calc_epsilon_Fk_zeta(r,S3_zeta_notide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);

%% Vorticity & Velocity SF2\SF3
color_red=[229,76,94]/255;
color_blue=[72,116,203]/255;
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 26 16]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [4.5 8.5 7.5 6]);
pcolor(lon_rho,lat_rho,zeta_lp_tide);shading interp
rectangle('Position', [116, 18, 1, 0.8], 'EdgeColor', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
caxis([-1 1])
colortable = textread('MPL_bwr.txt');
colormap(colortable);
xlim([113 119]);ylim([16.7 21.5]);box on
yticks([17:2:21]);xticks([113:2:119]);
yticklabels({'17°N','19°N','21°N'});
xticklabels({'113°E','115°E','117°E','119°E'});
set(AH1,'fontsize',11,'fontname','Arial');
set(AH1,'layer','top');
set(AH1,'linewidth', 1)
tit=title('\bf(a) \rm2-day Low-pass Ro (Tide)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.8 6 0])

AH1b= axes('Parent', gcf, 'unit', 'centimeters', 'Position', [8.74 11.88 3.25 2.6]);
pcolor(lon_rho,lat_rho,zeta_lp_tide);shading flat;hold on
contour(lon_rho,lat_rho,zeta_lp_tide,[0.2,0.2],'LineStyle','-','linecolor',[0.3 0.3 0.3],'linewidth',1);
caxis([-1 1])
colortable = textread('MPL_bwr.txt');
colormap(colortable);
xlim([116 117]);ylim([18 18.8]);box on
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);box on
set(gca, 'LineWidth', 1.2);
set(gca,'layer','top');set(gca,'ticklength',[0 0]);


AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [15 8.5 7.5 6]);
pcolor(lon_rho,lat_rho,zeta_lp_notide);shading flat
rectangle('Position', [115.3, 18, 1, 0.8], 'EdgeColor', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
caxis([-1 1])
colortable = textread('MPL_bwr.txt');
colormap(colortable);
xlim([112 118]);ylim([16.7 21.5]);box on
yticks([17:2:21]);xticks([112:2:118]);
yticklabels({'17°N','19°N','21°N'});
xticklabels({'112°E','114°E','116°E','118°E'});
set(AH2,'fontsize',11,'fontname','Arial');
set(AH2,'linewidth', 1)
set(AH2,'layer','top');
hhh = colorbar;
set(hhh, 'unit', 'centimeters','Position', [23.5 9.5 0.4 4]);
title(hhh,'\zeta/\itf');
tit=title('\bf(b) \rm2-day Low-pass Ro (Notide)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.2 6 0])

AH2b= axes('Parent', gcf, 'unit', 'centimeters', 'Position', [19.25 11.88 3.25 2.6]);
pcolor(lon_rho,lat_rho,zeta_lp_notide);shading interp;hold on
contour(lon_rho,lat_rho,zeta_lp_notide,[0.2,0.2],'LineStyle','-','linecolor',[0.3 0.3 0.3],'linewidth',1);
caxis([-1 1])
colortable = textread('MPL_bwr.txt');
colormap(colortable);
xlim([115.3 116.3]);ylim([18 18.8]);box on
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);box on
set(gca, 'LineWidth', 1.2);
set(gca,'layer','top');set(gca,'ticklength',[0 0]);

AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 6 4.8]);
plot(r, S2_zeta_tide_lp./r.*1e10,'-','color',color_red, 'LineWidth', 1.2, 'DisplayName', 'Tide lp'); hold on
plot(r, S2_zeta_notide_lp./r.*1e10, '-','color',color_blue, 'LineWidth', 1.2, 'DisplayName', 'Notide lp'); 
set(gca, 'XScale', 'log');
% ylim([0 7])
xlim([1 1e2])
set(gca,'FontName', 'Arial','FontSize',11)
set(gca,'linewidth', 1)
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
ylabel('D2_\omega(r)/r (s^{-2}km^{-1}×10^{-10})','fontsize',12,'fontname','Arial');
tit=title('\bf(c) \rmVorticity SF2','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [1.55 4.9 0])
set(AH3,'layer','top');
legend('show', 'Location', 'southwest','FontSize', 12);box on
set(gca,'layer','top');


AH6 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [10.5 2 6 4.8]);
loglog(r, S3_zeta_tide_lp,'-','color',color_red, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); hold on
loglog(r, -S3_zeta_tide_lp, '--','color',color_red, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); 
loglog(r, S3_zeta_notide_lp,'-','color',color_blue, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); hold on
loglog(r, -S3_zeta_notide_lp, '--','color',color_blue, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); 
text(3,0.5e-13,'Solid: Positive Values','fontsize',11,'fontname','Arial');
text(1.85,0.3e-13,'Dashed: Negative Values','fontsize',11,'fontname','Arial');
ylim([2e-14 4e-12])
xlim([1 1e2])
set(gca,'FontName', 'Arial','FontSize',11)
set(gca,'linewidth', 1)
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
ylabel('D3_\omega(r) (ms^{-3})','fontsize',12,'fontname','Arial');
tit=title('\bf(d) \rmVorticity SF3','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [1.55 4.9 0])
set(gca,'layer','top');

AH7 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [19 2 6 4.8]);
loglog(r, S3_tide_lp,'-','color',color_red, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); hold on
loglog(r, -S3_tide_lp, '--','color',color_red, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); 
loglog(r, S3_notide_lp,'-','color',color_blue, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); hold on
loglog(r, -S3_notide_lp, '--','color',color_blue, 'LineWidth', 1.2, 'DisplayName', 'tide lp'); 
text(5,0.0015e-5,'Solid: Positive Values','fontsize',11,'fontname','Arial');
text(3.05,0.0003e-5,'Dashed: Negative Values','fontsize',11,'fontname','Arial');
ylim([1e-9 5e-3])
xlim([1 1e2])
set(gca,'FontName', 'Arial','FontSize',11)
set(gca,'linewidth', 1)
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
ylabel('D3(r) (m^3s^{-3})','fontsize',12,'fontname','Arial');
tit=title('\bf(e) \rmVelocity SF3','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [1.55 4.9 0])
set(gca,'layer','top');

print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig 2.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi