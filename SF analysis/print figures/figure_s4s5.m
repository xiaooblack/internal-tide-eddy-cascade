%%% figure S3 model lp errorbar
clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
load('SF_lp_model_witherrobar.mat')

%% figure
color_red=[229, 76, 94]/255;
color_blue=[72, 116, 203]/255;
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 25 20]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 12 9 6]);
loglog(r,S2_tide_lp, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide lp');hold on
fill([r  fliplr(r)], [S2_tide_lp_high   fliplr(S2_tide_lp_low)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S2_notide_lp, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide lp');hold on
fill([r  fliplr(r)], [S2_notide_lp_high  fliplr(S2_notide_lp_low)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2(r) (m^2s^{-2})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH1,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(a) \rmVelocity SF2 with errorbar (Model lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 12 9 6]);
plot(r,S2_tide_lp./r.*1e3, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide lp');hold on
fill([r  fliplr(r)], [S2_tide_lp_high./r.*1e3  fliplr(S2_tide_lp_low./r.*1e3)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,S2_notide_lp./r.*1e3, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide lp');hold on
fill([r  fliplr(r)], [S2_notide_lp_high./r.*1e3  fliplr(S2_notide_lp_low./r.*1e3)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2(r)/r (m^2s^{-2}km^{-1} ×10^{-3})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH2,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(b) \rmVelocity SF2/r with errorbar (Model lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 3.5 9 6]);
loglog(r,S2_zeta_tide_lp, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide lp');hold on
fill([r  fliplr(r)], [S2_zeta_tide_lp_high  fliplr(S2_zeta_tide_lp_low)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S2_zeta_notide_lp, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide lp');hold on
fill([r  fliplr(r)], [S2_zeta_notide_lp_high  fliplr(S2_zeta_notide_lp_low)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2_\omega(r) (s^{-2})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH3,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(c) \rmVorticity SF2 with errorbar (Model lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 3.5 9 6]);
plot(r,S2_zeta_tide_lp./r.*1e11, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide lp');hold on
fill([r  fliplr(r)], [S2_zeta_tide_lp_high./r.*1e11  fliplr(S2_zeta_tide_lp_low./r.*1e11)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,S2_zeta_notide_lp./r.*1e11, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide lp');hold on
fill([r  fliplr(r)], [S2_zeta_notide_lp_high./r.*1e11  fliplr(S2_zeta_notide_lp_low./r.*1e11)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% ylim([-4 1]);
xlim([1 100])
ylabel('D2_\omega(r)/r (s^{-2}km^{-1} ×10^{-11})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH4,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(d) \rmVorticity SF2/r with errorbar (Model lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.15 6.1 0])
lll=legend('show', 'Location', 'southwest','FontSize', 11);

print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S4.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi

%% SF3
color_red=[229, 76, 94]/255;
color_blue=[72, 116, 203]/255;

figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 23 10]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 8 6]);
plot(r,S3_zeta_tide_lp./r.*1e13, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide lp');hold on
fill([r  fliplr(r)], [S3_zeta_tide_lp_high./r.*1e13  fliplr(S3_zeta_tide_lp_low./r.*1e13)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,S3_zeta_notide_lp./r.*1e13, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide lp');hold on
fill([r  fliplr(r)], [S3_zeta_notide_lp_high./r.*1e13  fliplr(S3_zeta_notide_lp_low./r.*1e13)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
% ylim([-4 1]);
xlim([1 100])
ylabel('D3_\omega(r)/r (ms^{-3}km^{-1} ×10^{-13})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH1,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(a) \rmVorticity SF3/r with errorbar (Model lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.1 6.1 0])
lll=legend('show', 'Location', 'northwest','FontSize', 11);

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13 2 8 6]);
plot(r,S3_tide_lp./r.*1e5, 'Color', color_red,'LineWidth', 1.2, 'DisplayName', 'Tide lp');hold on
fill([r  fliplr(r)], [S3_tide_lp_high./r.*1e5  fliplr(S3_tide_lp_low./r.*1e5)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S3_notide_lp./r.*1e5, 'Color', color_blue,'LineWidth', 1.2, 'DisplayName', 'Notide lp');hold on
fill([r  fliplr(r)], [S3_notide_lp_high./r.*1e5  fliplr(S3_notide_lp_low./r.*1e5)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
ylim([-0.5 3]);
xlim([1 100])
set(gca,'xscale','log');
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-5})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH2,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(b) \rmVelocity SF3/r with errorbar (Model lp)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [4.1 6.1 0])
lll=legend('show', 'Location', 'northwest','FontSize', 11);
print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S5.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi
