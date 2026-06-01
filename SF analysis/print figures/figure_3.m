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
load('SF_lp_model_witherrobar.mat')

[k_tide_lp,R_tide_lp,F_k_tide_lp,Pidk_tide_lp,S3_recon_tide_lp]=calc_epsilon_Fk(r,S3_tide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_tide_lp,R_zeta_tide_lp,F_k_zeta_tide_lp,Pidk_zeta_tide_lp,S3_zeta_recon_tide_lp]=calc_epsilon_Fk_zeta(r,S3_zeta_tide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_tide_lp_low,R_tide_lp_low,F_k_tide_lp_low,Pidk_tide_lp_low,S3_recon_tide_lp_low]=calc_epsilon_Fk(r,S3_tide_lp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_tide_lp_high,R_tide_lp_high,F_k_tide_lp_high,Pidk_tide_lp_high,S3_recon_tide_lp_high]=calc_epsilon_Fk(r,S3_tide_lp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_tide_lp_low,R_zeta_tide_lp_low,F_k_zeta_tide_lp_low,Pidk_zeta_tide_lp_low,S3_zeta_recon_tide_lp_low]=calc_epsilon_Fk_zeta(r,S3_zeta_tide_lp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_tide_lp_high,R_zeta_tide_lp_high,F_k_zeta_tide_lp_high,Pidk_zeta_tide_lp_high,S3_zeta_recon_tide_lp_high]=calc_epsilon_Fk_zeta(r,S3_zeta_tide_lp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);

[k_notide_lp,R_notide_lp,F_k_notide_lp,Pidk_notide_lp,S3_recon_notide_lp]=calc_epsilon_Fk(r,S3_notide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_notide_lp,R_zeta_notide_lp,F_k_zeta_notide_lp,Pidk_zeta_notide_lp,S3_zeta_recon_notide_lp]=calc_epsilon_Fk_zeta(r,S3_zeta_notide_lp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_notide_lp_low,R_notide_lp_low,F_k_notide_lp_low,Pidk_notide_lp_low,S3_recon_notide_lp_low]=calc_epsilon_Fk(r,S3_notide_lp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_notide_lp_high,R_notide_lp_high,F_k_notide_lp_high,Pidk_notide_lp_high,S3_recon_notide_lp_high]=calc_epsilon_Fk(r,S3_notide_lp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_notide_lp_low,R_zeta_notide_lp_low,F_k_zeta_notide_lp_low,Pidk_zeta_notide_lp_low,S3_zeta_recon_notide_lp_low]=calc_epsilon_Fk_zeta(r,S3_zeta_notide_lp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_notide_lp_high,R_zeta_notide_lp_high,F_k_zeta_notide_lp_high,Pidk_zeta_notide_lp_high,S3_zeta_recon_notide_lp_high]=calc_epsilon_Fk_zeta(r,S3_zeta_notide_lp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);

%%
color_red=[229, 76, 94]/255;
color_blue=[72, 116, 203]/255;
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 26 12]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [14.5 2 10 8]);
plot(k_tide_lp, F_k_tide_lp.*1e6,'-','color',color_red, 'LineWidth', 1.2,'DisplayName', 'Tide lp'); hold on;
fill([k_tide_lp  fliplr(k_tide_lp)], [F_k_tide_lp_high.*1e6 ; flipud(F_k_tide_lp_low.*1e6)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(k_notide_lp, F_k_notide_lp.*1e6,'-','color',color_blue, 'LineWidth', 1.2,'DisplayName', 'Notide lp'); hold on;
fill([k_notide_lp  fliplr(k_notide_lp)], [F_k_notide_lp_high.*1e6 ; flipud(F_k_notide_lp_low.*1e6)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
xline(0.12, 'LineStyle','--','Color',[1 0 1], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
text(0.055,-15.5,'8.3km','fontsize',15,'fontname','Arial','Color',[1 0 1]);
text(0.0105,0.9,'Downscale transfer (+)','fontsize',11,'fontname','Arial','Color',[0.3 0.3 0.3]);
text(0.0105,-0.8,'Upscale transfer (-)','fontsize',11,'fontname','Arial','Color',[0.3 0.3 0.3]);
xlabel('\itk \rm(1/km)'); ylabel('F(\itk\rm) (m^3s^{-3}km^{-1} ×10^{-6})');
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12) 
tit=title('\bf(b) \rmSpectral Energy flux','FontName', 'Arial','FontSize',13); 
set(tit, 'unit', 'centimeters','Position', [2.5 8 0])
ylim([-17 2])
xlim([1e-2 1]);
legend('show', 'Location', 'southeast','FontSize', 13);box on

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 10 8]);
plot(k_zeta_tide_lp, F_k_zeta_tide_lp.*1e14,'-','color',color_red, 'LineWidth', 1.2,'DisplayName', 'Tide lp'); hold on;
fill([k_zeta_tide_lp  fliplr(k_zeta_tide_lp)], [F_k_zeta_tide_lp_high.*1e14 ; flipud(F_k_zeta_tide_lp_low.*1e14)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(k_zeta_notide_lp, F_k_zeta_notide_lp.*1e14,'-','color',color_blue, 'LineWidth', 1.2,'DisplayName', 'Notide lp'); hold on;
fill([k_zeta_notide_lp  fliplr(k_zeta_notide_lp)], [F_k_zeta_notide_lp_high.*1e14 ; flipud(F_k_zeta_notide_lp_low.*1e14)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
yline(8.2, 'LineStyle','--','Color',color_red, 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
text(1.1e-2,8.7,'8.2×10^{-14}','fontsize',12,'fontname','Arial','Color',color_red);
yline(7.2, 'LineStyle','--','Color',color_blue, 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
text(1.1e-2,6.7,'7.2×10^{-14}','fontsize',12,'fontname','Arial','Color',color_blue);
xlabel('\itk \rm(1/km)','FontSize',13); ylabel('F_\omega(\itk\rm) (ms^{-3}km^{-1} ×10^{-14})','FontSize',13);
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12)
tit=title('\bf(a) \rmSpectral Enstrophy flux','FontName', 'Arial','FontSize',13); 
set(tit, 'unit', 'centimeters','Position', [2.85 8 0])
ylim([-0.5 10])
xlim([1e-2 1]);
legend('show', 'Location', 'southeast','FontSize', 13);box on
 
print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig 3.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi