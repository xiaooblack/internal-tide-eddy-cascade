%% figure 4 energy inject & R(w*b) & R(tao*shear)
clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
warning off all;  
%read grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
% SF3 & F(k)
% ===================== 1. 核心可调参数 =====================
% --- 数据参与拟合的范围（尽量用全部数据）---
r_min_data = 1;     % km  所有 >= 1 km 的数据都参与拟合
r_max_data = 100;   % km  所有 <= 100 km 的数据都参与拟合

% --- 你关心的物理波数范围（对应目标尺度 5~90 km）---
r_min_target = 1;   % km  → k_max = 1/5 = 0.2  1/km
r_max_target = 500;  % km  → k_min = 1/90 ≈ 0.011 1/km

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


load('C:\Users\Lenovo\Desktop\结构函数\pai_time_cospec.mat');
load('C:\Users\Lenovo\Desktop\结构函数\wb_cospec.mat');

Ck_all_tide_sum1 = movmean(Ck_all_tide_sum,6,1);
Ck_all_tide_sum1 = movmean(Ck_all_tide_sum1,3,1)';

N = 48;% 每组元素个数
n = ceil(size(Ck_all_tide_sum1,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
Ck_all_tide_sum1_block = zeros(n, length(k_center));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(Ck_all_tide_sum1,1));
    Ck_all_tide_sum1_block(i,:) = sum(Ck_all_tide_sum1(startIdx:endIdx,:))./(endIdx-startIdx+1);
end

[Ck_all_tide_low,Ck_all_tide_high]=block_bootstrap_cospec(Ck_all_tide_sum1_block,1000);


Ck_all_notide_sum1 = movmean(Ck_all_notide_sum,6,1);
Ck_all_notide_sum1 = movmean(Ck_all_notide_sum1,3,1)';

N = 48;% 每组元素个数
n = ceil(size(Ck_all_notide_sum1,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
Ck_all_notide_sum1_block = zeros(n, length(k_center));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(Ck_all_notide_sum1,1));
    Ck_all_notide_sum1_block(i,:) = sum(Ck_all_notide_sum1(startIdx:endIdx,:))./(endIdx-startIdx+1);
end

[Ck_all_notide_low,Ck_all_notide_high]=block_bootstrap_cospec(Ck_all_notide_sum1_block,1000);

wb_tide_sum1 = movmean(wb_tide_sum,6,1);
wb_tide_sum1 = movmean(wb_tide_sum1,3,1)';
N = 6;% 每组元素个数
n = ceil(size(wb_tide_sum1,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
wb_tide_sum1_block = zeros(n, length(k_center));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(wb_tide_sum1,1));
    wb_tide_sum1_block(i,:) = sum(wb_tide_sum1(startIdx:endIdx,:))./(endIdx-startIdx+1);
end
[wb_tide_low,wb_tide_high]=block_bootstr    ap_cospec(wb_tide_sum1_block,1000);


wb_notide_sum1 = movmean(wb_notide_sum,6,1);
wb_notide_sum1 = movmean(wb_notide_sum1,3,1)';
N = 48;% 每组元素个数
n = ceil(size(wb_notide_sum1,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
wb_notide_sum1_block = zeros(n, length(k_center));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(wb_notide_sum1,1));
    wb_notide_sum1_block(i,:) = sum(wb_notide_sum1(startIdx:endIdx,:))./(endIdx-startIdx+1);
end
[wb_notide_low,wb_notide_high]=block_bootstrap_cospec(wb_notide_sum1_block,1000);
%%
color_red=[229, 76, 94]/255;
color_blue=[72, 116, 203]/255;
figure
set(gcf, 'unit', 'centimeters', 'position', [1 2 22 18]);% [left, bottom, width, height]
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 10.5 18 6]);
plot(k_tide_lp, Pidk_tide_lp.*k_tide_lp'.*1e6','-','color',color_red, 'LineWidth', 1.2,'DisplayName', 'Tide lp'); hold on;
fill([k_tide_lp  fliplr(k_tide_lp)], [Pidk_tide_lp_high.*k_tide_lp'.*1e6 ; flipud(Pidk_tide_lp_low.*k_tide_lp'.*1e6)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(k_notide_lp, Pidk_notide_lp.*k_notide_lp'.*1e6','-','color',color_blue, 'LineWidth', 1.2,'DisplayName', 'Notide lp'); hold on;
fill([k_notide_lp  fliplr(k_notide_lp)], [Pidk_notide_lp_high.*k_notide_lp'.*1e6 ; flipud(Pidk_notide_lp_low.*k_notide_lp'.*1e6)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
xline(1/14,'LineStyle','--','Color',[1 0 1], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
text(0.085,16,'ML Radius ~14km','fontsize',13,'fontname','Arial','Color',[1 0 1]);
annotation(gcf,'doublearrow',[0.400482954545455,0.400300000000001],...
    [0.858425925925926,0.799991203703704],...
    'Color',[0.184313725490196 0.745098039215686 0.937254901960784],...
    'LineWidth',1.5);
text(0.02,14,'~37.4% larger','fontsize',13,'fontname','Arial','Color',[0.184313725490196 0.745098039215686 0.937254901960784]);
xlabel('\itk \rm(1/km)'); ylabel('\it\epsilon_i·k\rm (m^3s^{-3}km^{-1} ×10^{-6})');
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12)
tit=title('\bf(a) \rmEnergy injection rate \epsilon_i','FontName', 'Arial');
set(tit, 'unit', 'centimeters','Position', [2.85 6 0])
xlim([1e-2 1]);
ylim([-2 20]);
legend('show', 'Location', 'northeast','FontSize', 13,'FontName', 'Arial');box on

wb_spec_tide1 = movmean(wb_spec_tide,6);
wb_spec_tide1 = movmean(wb_spec_tide1,3);
wb_spec_notide1 = movmean(wb_spec_notide,6);
wb_spec_notide1 = movmean(wb_spec_notide1,3);
AH2 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 2 8 6]);
plot(k_center*2400, wb_spec_tide1.*k_center*1000.*1e6,'-','color',color_red, 'LineWidth', 1.2,'DisplayName', 'Tide lp'); hold on;
fill([k_center*2400 ; flipud(k_center*2400)], [wb_tide_high.*k_center'*1000.*1e6  fliplr(wb_tide_low.*k_center'*1000.*1e6)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(k_center*2400, wb_spec_notide1.*k_center*1000.*1e6,'-','color',color_blue, 'LineWidth', 1.2,'DisplayName', 'Notide lp'); hold on;
fill([k_center*2400 ; flipud(k_center*2400)], [wb_notide_high.*k_center'*1000.*1e6  fliplr(wb_notide_low.*k_center'*1000.*1e6)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
annotation(gcf,'arrow',[0.274612187367359 0.235963064700842],...
    [0.393935185185185 0.371886574074074],...
    'Color',[0.184313725490196 0.745098039215686 0.937254901960784],...
    'LineWidth',2);
text(0.28, 11.3, {'Baroclinic energy', 'release'},'FontSize', 13,'FontName', 'Arial','Color', [0.184313725490196 0.745098039215686 0.937254901960784], ...
    'HorizontalAlignment', 'center')  
xlabel('\itk \rm(1/km)'); ylabel('\itE_{VBF}(k)·k\rm (m^3s^{-3}km^{-1} ×10^{-6})')
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12)
tit=title('\bf(b) \rmCospectrum of buoyancy flux','FontName', 'Arial');
set(tit, 'unit', 'centimeters','Position', [3.45 6 0])
xlim([1e-2 1])
ylim([0 13])
% legend('show', 'Location', 'northeast','FontSize', 13,'FontName', 'Arial');box on


AH3 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [12 2 8 6]);
Ck_all_tide_mean1 = movmean(Ck_all_tide_mean,6);
Ck_all_tide_mean1 = movmean(Ck_all_tide_mean1,3);
Ck_all_notide_mean1 = movmean(Ck_all_notide_mean,6);
Ck_all_notide_mean1 = movmean(Ck_all_notide_mean1,3);
plot(k_tide*2400, Ck_all_tide_mean1.*k_tide*1000.*1e6,'-','color',color_red, 'LineWidth', 1.2,'DisplayName', 'Tide lp'); hold on;
fill([k_center*2400 ; flipud(k_center*2400)], [Ck_all_tide_high.*k_center'*1000.*1e6  fliplr(Ck_all_tide_low.*k_center'*1000.*1e6)], color_red, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot(k_notide*2400, Ck_all_notide_mean1.*k_notide*1000.*1e6,'-','color',color_blue, 'LineWidth', 1.2,'DisplayName', 'Notide lp'); hold on;
fill([k_center*2400 ; flipud(k_center*2400)], [Ck_all_notide_high.*k_center'*1000.*1e6  fliplr(Ck_all_notide_low.*k_center'*1000.*1e6)], color_blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
yline(0, 'LineStyle','--','Color',[0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility','off','Layer','bottom');hold on;
annotation(gcf,'arrow',[0.756612804660523 0.720388257575758],...
    [0.395405092592592 0.374826388888889],...
    'Color',[0.184313725490196 0.745098039215686 0.937254901960784],...
    'LineWidth',2);
text(0.33, 6, {'kinetic energy', 'transfer'},'FontSize', 13,'FontName', 'Arial','Color', [0.184313725490196 0.745098039215686 0.937254901960784], ...
    'HorizontalAlignment', 'center')  
xlabel('\itk \rm(1/km)'); 
ylabel('\itE_{\Pi}(k)·k\rm (m^3s^{-3}km^{-1} ×10^{-6})')
set(gca, 'XScale', 'log');set(gca,'linewidth', 1)
set(gca,'FontName', 'Arial','FontSize',12)
tit=title('\bf(c) \rmCospectrum of KE flux','FontName', 'Arial');
set(tit, 'unit', 'centimeters','Position', [2.8 6 0])
xlim([1e-2 1])
% legend('show', 'Location', 'northeast','FontSize', 13,'FontName', 'Arial');box on
ylim([-1.5 7.5])

print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig 4.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi