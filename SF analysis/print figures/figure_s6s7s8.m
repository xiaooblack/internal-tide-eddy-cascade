%%% figure S6 model hp errorbar
clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数

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
load('SF_hp_model_witherrobar.mat')

[k_tide_hp,R_tide_hp,F_k_tide_hp,Pidk_tide_hp,S3_recon_tide_hp]=calc_epsilon_Fk(r,S3_tide_hp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_tide_hp,R_zeta_tide_hp,F_k_zeta_tide_hp,Pidk_zeta_tide_hp,S3_zeta_recon_tide_hp]=calc_epsilon_Fk_zeta(r,S3_zeta_tide_hp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_tide_hp_low,R_tide_hp_low,F_k_tide_hp_low,Pidk_tide_hp_low,S3_recon_tide_hp_low]=calc_epsilon_Fk(r,S3_tide_hp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_tide_hp_high,R_tide_hp_high,F_k_tide_hp_high,Pidk_tide_hp_high,S3_recon_tide_hp_high]=calc_epsilon_Fk(r,S3_tide_hp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_tide_hp_low,R_zeta_tide_hp_low,F_k_zeta_tide_hp_low,Pidk_zeta_tide_hp_low,S3_zeta_recon_tide_hp_low]=calc_epsilon_Fk_zeta(r,S3_zeta_tide_hp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_tide_hp_high,R_zeta_tide_hp_high,F_k_zeta_tide_hp_high,Pidk_zeta_tide_hp_high,S3_zeta_recon_tide_hp_high]=calc_epsilon_Fk_zeta(r,S3_zeta_tide_hp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);

[k_notide_hp,R_notide_hp,F_k_notide_hp,Pidk_notide_hp,S3_recon_notide_hp]=calc_epsilon_Fk(r,S3_notide_hp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_notide_hp,R_zeta_notide_hp,F_k_zeta_notide_hp,Pidk_zeta_notide_hp,S3_zeta_recon_notide_hp]=calc_epsilon_Fk_zeta(r,S3_zeta_notide_hp,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_notide_hp_low,R_notide_hp_low,F_k_notide_hp_low,Pidk_notide_hp_low,S3_recon_notide_hp_low]=calc_epsilon_Fk(r,S3_notide_hp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_notide_hp_high,R_notide_hp_high,F_k_notide_hp_high,Pidk_notide_hp_high,S3_recon_notide_hp_high]=calc_epsilon_Fk(r,S3_notide_hp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_notide_hp_low,R_zeta_notide_hp_low,F_k_zeta_notide_hp_low,Pidk_zeta_notide_hp_low,S3_zeta_recon_notide_hp_low]=calc_epsilon_Fk_zeta(r,S3_zeta_notide_hp_low,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);
[k_zeta_notide_hp_high,R_zeta_notide_hp_high,F_k_zeta_notide_hp_high,Pidk_zeta_notide_hp_high,S3_zeta_recon_notide_hp_high]=calc_epsilon_Fk_zeta(r,S3_zeta_notide_hp_high,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk);

%%
color_red=[229, 76, 94]/255;
color_blue=[72, 116, 203]/255;
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 25 20]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 12 9 6]);
loglog(r,S2_tide_hp, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide hp');hold on
fill([r  fliplr(r)], [S2_tide_hp_high   fliplr(S2_tide_hp_low)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S2_notide_hp, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide hp');hold on
fill([r  fliplr(r)], [S2_notide_hp_high  fliplr(S2_notide_hp_low)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2(r) (m^2s^{-2})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH1,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(a) \rmVelocity SF2 with errorbar (Model hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 12 9 6]);
plot(r,S2_tide_hp./r.*1e3, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide hp');hold on
fill([r  fliplr(r)], [S2_tide_hp_high./r.*1e3  fliplr(S2_tide_hp_low./r.*1e3)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,S2_notide_hp./r.*1e3, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide hp');hold on
fill([r  fliplr(r)], [S2_notide_hp_high./r.*1e3  fliplr(S2_notide_hp_low./r.*1e3)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2(r)/r (m^2s^{-2}km^{-1} ×10^{-3})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH2,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(b) \rmVelocity SF2/r with errorbar (Model hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.15 6.1 0])
lll=legend('show', 'Location', 'northwest','FontSize', 11);

AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 3.5 9 6]);
loglog(r,S2_zeta_tide_hp, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide hp');hold on
fill([r  fliplr(r)], [S2_zeta_tide_hp_high  fliplr(S2_zeta_tide_hp_low)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S2_zeta_notide_hp, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide hp');hold on
fill([r  fliplr(r)], [S2_zeta_notide_hp_high  fliplr(S2_zeta_notide_hp_low)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2_\omega(r) (s^{-2})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH3,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(c) \rmVorticity SF2 with errorbar (Model hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 3.5 9 6]);
plot(r,S2_zeta_tide_hp./r.*1e11, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide hp');hold on
fill([r  fliplr(r)], [S2_zeta_tide_hp_high./r.*1e11  fliplr(S2_zeta_tide_hp_low./r.*1e11)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,S2_zeta_notide_hp./r.*1e11, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide hp');hold on
fill([r  fliplr(r)], [S2_zeta_notide_hp_high./r.*1e11  fliplr(S2_zeta_notide_hp_low./r.*1e11)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2_\omega(r)/r (s^{-2}km^{-1} ×10^{-11})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH4,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(d) \rmVorticity SF2/r with errorbar (Model hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.15 6.1 0])
lll=legend('show', 'Location', 'northeast','FontSize', 11);
color_red=[229, 76, 94]/255;
color_blue=[72, 116, 203]/255;
print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S6.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi

%% SF3
color_red=[229, 76, 94]/255;
color_blue=[72, 116, 203]/255;

figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 23 19]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 10.5 8 6]);
plot(r,S3_zeta_tide_hp./r.*1e13, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide hp');hold on
fill([r  fliplr(r)], [S3_zeta_tide_hp_high./r.*1e13  fliplr(S3_zeta_tide_hp_low./r.*1e13)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
ylim([-15 2]);
xlim([1 100])
ylabel('D3_\omega(r)/r (ms^{-3}km^{-1} ×10^{-13})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH1,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(a) \rmVorticity SF3/r with errorbar (Tide hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13 10.5 8 6]);
plot(r,S3_tide_hp./r.*1e5, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide hp');hold on
fill([r  fliplr(r)], [S3_tide_hp_high./r.*1e5  fliplr(S3_tide_hp_low./r.*1e5)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
ylim([-16 1]);
xlim([1 100])
set(gca,'xscale','log');
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-5})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH2,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(b) \rmVelocity SF3/r with errorbar (Tide hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.02 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);


AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 8 6]);
plot(r,S3_zeta_notide_hp./r.*1e13, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide hp');hold on
fill([r  fliplr(r)], [S3_zeta_notide_hp_high./r.*1e13  fliplr(S3_zeta_notide_hp_low./r.*1e13)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
ylim([-1.2 0.2]);
xlim([1 100])
ylabel('D3_\omega(r)/r (ms^{-3}km^{-1} ×10^{-13})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH3,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(c) \rmVorticity SF3/r with errorbar (Notide hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.25 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13 2 8 6]);
plot(r,S3_notide_hp./r.*1e5, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide hp');hold on
fill([r  fliplr(r)], [S3_notide_hp_high./r.*1e5  fliplr(S3_notide_hp_low./r.*1e5)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
ylim([-0.04 0.02]);
xlim([1 100])
set(gca,'xscale','log');
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-5})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH4,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(d) \rmVelocity SF3/r with errorbar (Notide hp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.25 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);
print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S7.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi

%% F(k) hp
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 23 19]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 10.5 8 6]);
plot(k_zeta_tide_hp, F_k_zeta_tide_hp.*1e14,'-','color',color_red, 'LineWidth', 1.2,'DisplayName', 'Tide hp'); hold on;
fill([k_zeta_tide_hp  fliplr(k_zeta_tide_hp)], [F_k_zeta_tide_hp_high.*1e14 ; flipud(F_k_zeta_tide_hp_low.*1e14)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
xlabel('\itk \rm(1/km)','FontSize',13); ylabel('F_\omega(\itk\rm) (ms^{-3}km^{-1} ×10^{-14})','FontSize',13);
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12)
tit=title('\bf(a) \rmSpectral Enstrophy flux (Tide hp)','FontName', 'Arial','FontSize',13); 
set(tit, 'unit', 'centimeters','Position', [3.9 6 0])
ylim([-5 100])
xlim([1e-2 1]);
legend('show', 'Location', 'northwest','FontSize', 13);box on

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13 10.5 8 6]);
plot(k_tide_hp, F_k_tide_hp.*1e6,'-','color',color_red, 'LineWidth', 1.2,'DisplayName', 'Tide hp'); hold on;
fill([k_tide_hp  fliplr(k_tide_hp)], [F_k_tide_hp_high.*1e6 ; flipud(F_k_tide_hp_low.*1e6)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
xlabel('\itk \rm(1/km)'); ylabel('F(\itk\rm) (m^3s^{-3}km^{-1} ×10^{-6})');
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12) 
tit=title('\bf(b) \rmSpectral Energy flux (Tide hp)','FontName', 'Arial','FontSize',13); 
set(tit, 'unit', 'centimeters','Position', [3.6 6 0])
ylim([-5 80])
xlim([1e-2 1]);
legend('show', 'Location', 'northwest','FontSize', 13);box on


AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 8 6]);
plot(k_zeta_notide_hp, F_k_zeta_notide_hp.*1e14,'-','color',color_blue, 'LineWidth', 1.2,'DisplayName', 'Notide hp'); hold on;
fill([k_zeta_notide_hp  fliplr(k_zeta_notide_hp)], [F_k_zeta_notide_hp_high.*1e14 ; flipud(F_k_zeta_notide_hp_low.*1e14)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
xlabel('\itk \rm(1/km)','FontSize',13); ylabel('F_\omega(\itk\rm) (ms^{-3}km^{-1} ×10^{-14})','FontSize',13);
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12)
tit=title('\bf(c) \rmSpectral Enstrophy flux (Notide hp)','FontName', 'Arial','FontSize',13); 
set(tit, 'unit', 'centimeters','Position', [4.1 6 0])
ylim([-1 4])
xlim([1e-2 1]);
legend('show', 'Location', 'northwest','FontSize', 13);box on

AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13 2 8 6]);
plot(k_notide_hp, F_k_notide_hp.*1e6,'-','color',color_blue, 'LineWidth', 1.2,'DisplayName', 'Notide hp'); hold on;
fill([k_notide_hp  fliplr(k_notide_hp)], [F_k_notide_hp_high.*1e6 ; flipud(F_k_notide_hp_low.*1e6)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
xlabel('\itk \rm(1/km)'); ylabel('F(\itk\rm) (m^3s^{-3}km^{-1} ×10^{-6})');
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12) 
tit=title('\bf(d) \rmSpectral Energy flux (Notide hp)','FontName', 'Arial','FontSize',13); 
set(tit, 'unit', 'centimeters','Position', [3.85 6 0])
ylim([-0.2 0.2])
xlim([1e-2 1]);
legend('show', 'Location', 'northwest','FontSize', 13);box on
print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S8.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi
