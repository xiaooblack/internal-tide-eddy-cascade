%%% figure S1 N pair & errorbar
clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
load('SF_origion_model_witherrobar.mat');
load('SF_origion_drifters_witherrobar.mat');
load('SF_origion_ADCP_witherrobar.mat');


%% figure
colorlist_SF(1,:)=[0.85 0.1 0.1];
colorlist_SF(2,:)=[0.7 0.7 0.7];
colorlist_SF(3,:)=[0, 0, 0];
colorlist_SF(4,:)=[0 0.45 0.74];
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 35 20]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 12 9 6]);
loglog(r,N_drifters, 'Color', colorlist_SF(3,:),'LineWidth', 1.5, 'DisplayName', 'Drifters');hold on
loglog(r,N_ADCP, 'Color', colorlist_SF(4,:),'LineWidth', 1.5, 'DisplayName', 'ADCP');hold on
ylim([1e3 5e4]);
xlim([1 100])
ylabel('Number of pairs','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH1,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(a) \rmnumber of pairs in observation','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.35 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 12 9 6]);
loglog(r,N_tide, 'Color', colorlist_SF(1,:),'LineWidth', 1.5, 'DisplayName', 'Model');hold on
% ylim([1e3 5e4]);
xlim([1 100])
ylabel('Number of pairs','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH2,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(b) \rmnumber of pairs in model','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.85 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 3.5 9 6]);
loglog(r,S_total_drifters, 'Color', colorlist_SF(3,:),'LineWidth', 1.2, 'DisplayName', 'drifters');hold on
fill([r; flipud(r)], [S_total_drifters_high  fliplr(S_total_drifters_low)], colorlist_SF(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S_total_ADCP, 'Color', colorlist_SF(4,:),'LineWidth', 1.2, 'DisplayName', 'ADCP');hold on
fill([r; flipud(r)], [S_total_ADCP_high  fliplr(S_total_ADCP_low)], colorlist_SF(4,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S_total_tide, 'Color', colorlist_SF(1,:),'LineWidth', 1.2, 'DisplayName', 'Tide Model');hold on
fill([r; flipud(r)], [S_total_tide_high  fliplr(S_total_tide_low)], colorlist_SF(1,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
loglog(r,S_total_notide, 'Color', colorlist_SF(2,:),'LineWidth', 1.2, 'DisplayName', 'Notide Model');hold on
fill([r; flipud(r)], [S_total_notide_high  fliplr(S_total_notide_low)], colorlist_SF(2,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
ylim([1e-4 3])
xlim([1 100])
ylabel('D2(r) (m^2s^{-2})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH3,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(c) \rmSF2 with errorbar','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.05 6.1 0])
lll=legend('show', 'Location', 'southeast','FontSize', 11);

AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 3.5 9 6]);
plot(r,D3_drifters./r'.*1e3, 'Color', colorlist_SF(3,:),'LineWidth', 1.2, 'DisplayName', 'drifters');hold on
fill([r; flipud(r)], [D3_drifters_high./r'.*1e3  fliplr(D3_drifters_low./r'.*1e3)], colorlist_SF(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,D3_ADCP./r'.*1e3, 'Color', colorlist_SF(4,:),'LineWidth', 1.2, 'DisplayName', 'ADCP');hold on
fill([r; flipud(r)], [D3_ADCP_high./r'.*1e3  fliplr(D3_ADCP_low./r'.*1e3)], colorlist_SF(4,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,D3_tide./r'.*1e3, 'Color', colorlist_SF(1,:),'LineWidth', 1.2, 'DisplayName', 'Tide Model');hold on
fill([r; flipud(r)], [D3_tide_high./r'.*1e3  fliplr(D3_tide_low./r'.*1e3)], colorlist_SF(1,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,D3_notide./r'.*1e3, 'Color', colorlist_SF(2,:),'LineWidth', 1.2, 'DisplayName', 'Notide Model');hold on
fill([r; flipud(r)], [D3_notide_high./r'.*1e3  fliplr(D3_notide_low./r'.*1e3)], colorlist_SF(2,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
x_box = [1, 10, 10, 1, 1]; % 首尾重合形成闭合
y_box = [-0.8, -0.8, 0.2, 0.2, -0.8];
% 用m_line绘制，确保线条遵循地图投影
line(x_box, y_box, ...
    'Color', 'k', ...
    'LineWidth', 1.5, ...
    'LineStyle', '-'); 
ylim([-4 1])
xlim([1 100])
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-6})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(gca,'xscale','log');
set(AH4,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(d) \rmSF3/r with errorbar','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.2 6.1 0])
lll=legend('show', 'Location', 'southwest','FontSize', 11);


AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [25.5 7.5 9 7]);
plot(r,D3_drifters./r'.*1e3, 'Color', colorlist_SF(3,:),'LineWidth', 1.2, 'DisplayName', 'drifters');hold on
fill([r; flipud(r)], [D3_drifters_high./r'.*1e3  fliplr(D3_drifters_low./r'.*1e3)], colorlist_SF(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,D3_ADCP./r'.*1e3, 'Color', colorlist_SF(4,:),'LineWidth', 1.2, 'DisplayName', 'ADCP');hold on
fill([r; flipud(r)], [D3_ADCP_high./r'.*1e3  fliplr(D3_ADCP_low./r'.*1e3)], colorlist_SF(4,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,D3_tide./r'.*1e3, 'Color', colorlist_SF(1,:),'LineWidth', 1.2, 'DisplayName', 'Tide Model');hold on
fill([r; flipud(r)], [D3_tide_high./r'.*1e3  fliplr(D3_tide_low./r'.*1e3)], colorlist_SF(1,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(r,D3_notide./r'.*1e3, 'Color', colorlist_SF(2,:),'LineWidth', 1.2, 'DisplayName', 'Notide Model');hold on
fill([r; flipud(r)], [D3_notide_high./r'.*1e3  fliplr(D3_notide_low./r'.*1e3)], colorlist_SF(2,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
annotation(gcf,'arrow',[0.514047619047619 0.678089285714286],...
    [0.426017405063291 0.546000582944703],'LineWidth',1.5,'HeadWidth',12,...
    'HeadLength',12);
ylim([-0.8 0.2])
xlim([1 10])
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-6})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(gca,'xscale','log');
set(AH4,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(e) \rm Zoom of panel (d)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.3 7.1 0])
% lll=legend('show', 'Location', 'southwest','FontSize', 11);

print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S2.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi
%%
colorlist_SF(1,:)=[0.85 0.1 0.1];
colorlist_SF(2,:)=[0.7 0.7 0.7];
colorlist_SF(3,:)=[0, 0, 0];
colorlist_SF(4,:)=[0 0.45 0.74];
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 25 20]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 12 9 6]);
plot(r,D3_drifters./r'.*1e3, 'Color', colorlist_SF(3,:),'LineWidth', 1.2, 'DisplayName', 'drifters');hold on
fill([r; flipud(r)], [D3_drifters_high./r'.*1e3  fliplr(D3_drifters_low./r'.*1e3)], colorlist_SF(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
ylim([-4 1]);
xlim([1 100])
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-3})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(gca,'xscale','log');
set(AH1,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(a) \rmSF3/r with errorbar (drifters)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.1 6.1 0])
lll=legend('show', 'Location', 'southwest','FontSize', 11);

AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 12 9 6]);
plot(r,D3_ADCP./r'.*1e3, 'Color', colorlist_SF(4,:),'LineWidth', 1.2, 'DisplayName', 'ADCP');hold on
fill([r; flipud(r)], [D3_ADCP_high./r'.*1e3  fliplr(D3_ADCP_low./r'.*1e3)], colorlist_SF(4,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
ylim([-4 1]);
xlim([1 100])
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-3})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH2,'fontsize',12,'fontname','Arial','linewidth',1);
set(gca,'xscale','log');
tit=title('\bf(b) \rmSF3/r with errorbar (ADCP)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.05 6.1 0])
lll=legend('show', 'Location', 'southwest','FontSize', 11);

AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 3.5 9 6]);
plot(r,D3_tide./r'.*1e5, 'Color', colorlist_SF(1,:),'LineWidth', 1.2, 'DisplayName', 'Tide Model');hold on
fill([r; flipud(r)], [D3_tide_high./r'.*1e5  fliplr(D3_tide_low./r'.*1e5)], colorlist_SF(1,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
set(gca,'xscale','log');
ylim([-20 5])
xlim([1 100])
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-5})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(AH3,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(c) \rmSF3/r with errorbar (Tide Model)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.55 6.1 0])
lll=legend('show', 'Location', 'southwest','FontSize', 11);

AH4 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [13.5 3.5 9 6]);
plot(r,D3_notide./r'.*1e5, 'Color', colorlist_SF(2,:),'LineWidth', 1.2, 'DisplayName', 'Notide Model');hold on
fill([r; flipud(r)], [D3_notide_high./r'.*1e5  fliplr(D3_notide_low./r'.*1e5)], colorlist_SF(2,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
ylim([-20 5])
xlim([1 100])
ylabel('D3(r)/r (m^3s^{-3}km^{-1} ×10^{-5})','fontsize',11,'fontname','Arial');
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
set(gca,'xscale','log');
set(AH4,'fontsize',12,'fontname','Arial','linewidth',1);
tit=title('\bf(d) \rmSF3/r with errorbar (Notide Model)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.8 6.1 0])
lll=legend('show', 'Location', 'southwest','FontSize', 11);

print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S3.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi