clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
% read outputfile
% drifters 轨迹
% 筛选 797594400至808664400
% 第一批
rpath1='C:\Users\Lenovo\Desktop\结构函数\表漂观测一批15min插值\15min插值后';
list1=dir(fullfile(rpath1,'*.nc'));
m1=length(list1);
step1=0;
% start 797594400
% end 808664400
for ii=1:m1
    step1=step1+1;
    disp(['stepnum: ',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[list1(ii).folder,'\',list1(ii).name];
    drifter1(ii).time=ncread(filename,'time');
    drifter1(ii).lon=ncread(filename,'lon');
    drifter1(ii).lat=ncread(filename,'lat');
    drifter1(ii).u=ncread(filename,'u');
    drifter1(ii).v=ncread(filename,'v');

    id_time= drifter1(ii).time>=795906000 & drifter1(ii).time<=800380800;
    drifter1(ii).lon_use=drifter1(ii).lon(id_time);
    drifter1(ii).lat_use=drifter1(ii).lat(id_time);
end


% 第二批
rpath1='C:\Users\Lenovo\Desktop\结构函数\表漂观测二批15min插值\15min插值后';
list1=dir(fullfile(rpath1,'*.nc'));
m1=length(list1);

step1=0;
% start 797594400
% end 808664400
for ii=1:m1
    step1=step1+1;
    disp(['stepnum: ',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[list1(ii).folder,'\',list1(ii).name];
    drifter2(ii).time=ncread(filename,'time');
    drifter2(ii).lon=ncread(filename,'lon');
    drifter2(ii).lat=ncread(filename,'lat');
    drifter2(ii).u=ncread(filename,'u');
    drifter2(ii).v=ncread(filename,'v');

    id_time= drifter2(ii).time>=797202000 & drifter2(ii).time<=800380800;
    drifter2(ii).lon_use=drifter2(ii).lon(id_time);
    drifter2(ii).lat_use=drifter2(ii).lat(id_time);
end

% 第三批
rpath1='C:\Users\Lenovo\Desktop\结构函数\表漂观测三批15min插值\15min插值后';
list1=dir(fullfile(rpath1,'*.nc'));
m1=length(list1);

step1=0;
% start 797594400
% end 808664400
for ii=1:m1
    step1=step1+1;
    disp(['stepnum: ',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    filename=[list1(ii).folder,'\',list1(ii).name];
    drifter3(ii).time=ncread(filename,'time');
    drifter3(ii).lon=ncread(filename,'lon');
    drifter3(ii).lat=ncread(filename,'lat');
    drifter3(ii).u=ncread(filename,'u');
    drifter3(ii).v=ncread(filename,'v');
    
    id_time= drifter3(ii).time>=797202000 & drifter3(ii).time<=800380800;
    drifter3(ii).lon_use=drifter3(ii).lon(id_time);
    drifter3(ii).lat_use=drifter3(ii).lat(id_time);
end

%% SLA uvgo
% 2025-04-05 至 2025-05-16
lon_sla=ncread('C:\Users\Lenovo\Desktop\结构函数\cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D_1761221629064.nc','longitude',[1],[inf]);
lat_sla=ncread('C:\Users\Lenovo\Desktop\结构函数\cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D_1761221629064.nc','latitude',[1],[inf]);
ugos=ncread('C:\Users\Lenovo\Desktop\结构函数\cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D_1761221629064.nc','ugos',[1 1 1],[inf inf inf]);
vgos=ncread('C:\Users\Lenovo\Desktop\结构函数\cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D_1761221629064.nc','vgos',[1 1 1],[inf inf inf]);
sla=ncread('C:\Users\Lenovo\Desktop\结构函数\cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D_1761221629064.nc','sla',[1 1 1],[inf inf inf]);
time=ncread('C:\Users\Lenovo\Desktop\结构函数\cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D_1761221629064.nc','time');
KEos=(ugos.^2+vgos.^2)./2;
%% FSLE
fslefile='I:\FSLE\nrt_global_allsat_madt_fsle_20250407_20250413.nc';
[lonf,latf]=meshgrid(double(ncread(fslefile,'lon',[1],[inf])),double(ncread(fslefile,'lat',[1],[inf])));
bfsle=ncread(fslefile,'fsle_max',[1 1 1],[inf inf 1])';
%% K-F spectrum
load('K-F spec.mat')
%% 地形 
topo=ncread('E:\CROCO_external_datasets\DATASETS_CROCOTOOLS\Topo\etopo2.nc','topo',[8301 2701],[1101 1051]);
lon_topo=ncread('E:\CROCO_external_datasets\DATASETS_CROCOTOOLS\Topo\etopo2.nc','lon',[8301],[1101]);
lat_topo=ncread('E:\CROCO_external_datasets\DATASETS_CROCOTOOLS\Topo\etopo2.nc','lat',[2701],[1051]);
[lat_topo lon_topo]=meshgrid(lat_topo,lon_topo);
%% SF2 & SF3
load('SF_origion_drifters_witherrobar.mat');
load('SF_origion_ADCP_witherrobar.mat');
sum_DLL_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','sum_DLL');
sum_DTT_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','sum_DTT');
sum_D3_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','sum_D3');
counts_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','counts');

S_long_tide=sum(sum_DLL_tide,1)./sum(counts_tide,1);
S_trans_tide=sum(sum_DTT_tide,1)./sum(counts_tide,1);
S_total_tide=S_long_tide+S_trans_tide;
D3_tide=sum(sum_D3_tide,1)./sum(counts_tide,1);

sum_DLL_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','sum_DLL');
sum_DTT_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','sum_DTT');
sum_D3_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','sum_D3');
counts_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','counts');

S_long_notide=sum(sum_DLL_notide,1)./sum(counts_notide,1);
S_trans_notide=sum(sum_DTT_notide,1)./sum(counts_notide,1);
S_total_notide=S_long_notide+S_trans_notide;
D3_notide=sum(sum_D3_notide,1)./sum(counts_notide,1);

%% figure

[lat1 lon1]=meshgrid(lat_sla,lon_sla);
sla_use=sla(:,:,19);
figure
set(gcf, 'unit', 'centimeters','position', [1 1 23 23]);
%%%%% 第一个大子图
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [2 16 18 6]);
lon_mint=108;lon_maxt=122;
lat_mint=16.5;lat_maxt=21.575;
m_proj('miller','lon',[lon_mint lon_maxt],'lat',[lat_mint lat_maxt]);
m_pcolor(lon1,lat1,sla_use);shading interp;hold on
% FSLE
m_scatter(lonf(bfsle<-0.15),latf(bfsle<-0.15),1.5,[.7 .7 .7],'filled');
% 拖体线路
rpath1='I:\ADCP_ship\my process all depth';
list1=dir(fullfile(rpath1,'*.nc'));
m1=length(list1);
for ii=[2 4 6 7 9 10 11 12 13 15 16]
    filename_ADCP=[list1(ii).folder,'\',list1(ii).name];
    lon_ADCP=ncread(filename_ADCP,'lon_ADCP',[1],[inf]);
    lat_ADCP=ncread(filename_ADCP,'lat_ADCP',[1],[inf]);
    depth_ADCP=ncread(filename_ADCP,'depth',[1],[1]);
    depth_ADCP=repmat(depth_ADCP,length(lat_ADCP),1);
    p=m_plot(lon_ADCP,lat_ADCP,'-k','LineWidth', 1.2);hold on
    uistack(p, 'top');
end
% 绘制起始点
for i = 1:length(drifter1)
    m_plot(drifter1(i).lon_use(1), drifter1(i).lat_use(1), "pentagram", ...
    'MarkerEdgeColor', 'y', ...
    'MarkerFaceColor', 'y', ...
    'MarkerSize', 3);
end
for i = 1:length(drifter2)
    m_plot(drifter2(i).lon_use(1), drifter2(i).lat_use(1), "pentagram", ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 3);
end
for i = 1:length(drifter3)
    m_plot(drifter3(i).lon_use(1), drifter3(i).lat_use(1), "pentagram", ...
    'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', 'b', ...
    'MarkerSize', 3);
end
% 陆地
m_gshhs_h('patch', [222 204 176]/255, 'edgecolor', 'none'); 
hold on;
m_gshhs_h('color', 'k', 'linewidth', 0.5); % 用黑色('k')或深灰([0.2 0.2 0.2])，线宽0.5或0.8
m_grid('box','on','linestyle','none','xtick',lon_mint+1:2:lon_maxt,'ytick',lat_mint+0.5:1:lat_maxt,'fontsize',11,'fontname','Arial','linewidth',1);  % 'XaxisLocation','top',
tit=title('\bf(a)\rm Observational Site','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.1 6 0])
caxis([-0.3 0.3])
colortable = flipud(textread('MPL_RdBu.txt'));
colormap(colortable);
hhh = colorbar;
set(hhh,'fontsize',11,'fontname','Arial');
set(hhh, 'unit', 'centimeters','Position', [20 17 0.4 4]);
set(hhh,'ticks',[-0.3 0 0.3])
title(hhh,'SLA','fontsize',11,'fontname','Arial')

%%%%% 圆的那个 用etop2填色地形

AH2=axes('Parent',gcf,'unit', 'centimeters','Position',[15.7 18.7 3.2 3.2]);
m_proj('ortho','radius',15,'lat',18,'long',115,'rot',0);
m_contourf(lon_topo,lat_topo,topo,[-14000:100:0 100:100:7000],'edgecolor','none');
colormap(AH2,[m_colmap('blues',70);m_colmap('gland',70)]);  
caxis(AH2,[-7000 7000]);
m_grid('axes',AH2,'linest','none','xticklabels',[],'yticklabels',[],'linewidth',1.3);
m_line([108 122 122 108 108],[21.575 21.575 16.5 16.5 21.575],'linewidth',1.5,'color','r');

%%%%% 第二张子图 drifts
color_drifts=textread('WhiteBlueGreenYellowRed.txt');
color_drifts=color_drifts(10:200,:);
AH3=axes('Parent',gcf,'unit', 'centimeters','Position',[3 9 7 5.5]);
% 背景填KE
KEos_use=nanmean(KEos(:,:,27:42),3);hold on
% contourf(lon1,lat1,KEos_use,[0:0.02:0.08]);
pcolor(lon1,lat1,KEos_use);shading interp
caxis([0 0.1])
colortable = textread('MPL_PuRd.txt');
colortable=colortable(1:60,:);
colormap(AH3,colortable);
hhh = colorbar;
set(hhh,'fontsize',11,'fontname','Arial','location','eastoutside');
set(hhh, 'unit', 'centimeters','Position', [19.5 9.5 0.4 4]);
set(hhh,'ticks',[0 0.05 0.1])
title(hhh,'KE_g','fontsize',11,'fontname','Arial');

% drifts 轨迹
stepcolor=0;
for i = 1:length(drifter1)
    stepcolor=stepcolor+1;
    % 使用 m_line 绘制地理坐标系下的线，并赋予随机颜色
    p=plot(drifter1(i).lon_use, drifter1(i).lat_use, 'LineWidth', 1,'linestyle','-','Color',squeeze(color_drifts(2*stepcolor-1,:)));hold on
    
end

for i = 1:length(drifter2)
    stepcolor=stepcolor+1;
    % 使用 m_line 绘制地理坐标系下的线，并赋予随机颜色
    p=plot(drifter2(i).lon_use, drifter2(i).lat_use, 'LineWidth', 1,'linestyle','-','Color',squeeze(color_drifts(2*stepcolor-1,:)));hold on
    
end

for i = 1:length(drifter3)
    stepcolor=stepcolor+1;
    % 使用 m_line 绘制地理坐标系下的线，并赋予随机颜色
    p=plot(drifter3(i).lon_use, drifter3(i).lat_use, 'LineWidth', 1,'linestyle','-','Color',squeeze(color_drifts(2*stepcolor-1,:)));hold on
    
end
% 绘制起始点
for i = 1:length(drifter1)
    plot(drifter1(i).lon_use(1), drifter1(i).lat_use(1), "pentagram", ...
    'MarkerEdgeColor', 'y', ...
    'MarkerFaceColor', 'y', ...
    'MarkerSize', 3);
end
for i = 1:length(drifter2)
    plot(drifter2(i).lon_use(1), drifter2(i).lat_use(1), "pentagram", ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 3);
end
for i = 1:length(drifter3)
    plot(drifter3(i).lon_use(1), drifter3(i).lat_use(1), "pentagram", ...
    'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', 'b', ...
    'MarkerSize', 3);
end

tit=title('\bf(b) \rmTrajectory of Drifters','fontsize',12,'fontname','Arial');
    set(tit, 'unit', 'centimeters','Position', [2.35 5.5 0])
box on
yticks([16:2:21]);xticks([111:2:117]);
yticklabels({'16°N','18°N','20°N'});
xticklabels({'111°E','113°E','115°E','117°E'});
ylim([15.5 21]);
xlim([110.5 117.5])
set(AH3,'fontsize',11,'fontname','Arial','linewidth',1);
set(AH3,'layer','top');

%%%%%% K-F谱
AH4=axes('Parent',gcf,'unit', 'centimeters','Position',[12 9 7 5.5]);
pcolor(k_plot_tide, f_plot_tide, log10(S_plot_tide.*k_plot_tide.*f_plot_tide));
shading flat
set(gca, 'XScale', 'log', 'YScale', 'log','linewidth',1);
xlabel('\itk \rm(cpkm)')
ylabel('\itf \rm(cph)')
tit=title('\bf(c) \rm\itk-f \rmSpectrum of KE (Tide)','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.95 5.5 0])

colortable = textread('MPL_rainbow.txt');
colormap(AH4,colortable);
hhh = colorbar;
hhh.Label.String = 'log_{10}[KE×\omega×\itk\rm]';
set(hhh,'fontsize',11,'fontname','Arial','location','eastoutside');
set(hhh, 'unit', 'centimeters','Position', [21 9.5 0.4 4]);
set(hhh,'ticks',[-6:2:-2])


set(AH4,'fontsize',11,'fontname','Arial');
axis tight
xlim([min(k_plot_tide(1,:)), max(k_plot_tide(1,:))])
ylim([min(f_plot_tide(:,1)), max(f_plot_tide(:,1))])
caxis([-6.5 -2])
hold on


% -----------------------------
% 叠加内波模态理论曲线
% ω^2 = f^2 + c_n^2 k^2
% 注意单位换算：
% k [cpkm] -> rad/m: k_rad = 2*pi*k/1000
% omega [rad/s]
% f_plot_theory [cph] = omega/(2*pi)*3600
% -----------------------------
mode_num = [1 3 10];
cn = [2.8, 1.1, 0.32];   % [m/s]

f_cor = sw_f(18);
f0_cph = abs(f_cor) / (2*pi) * 3600;

k_theory = logspace(log10(min(k_plot_notide(1,:))), log10(max(k_plot_notide(1,:))), 300);

for i = 1:length(cn)
    omega = sqrt(f_cor^2 + (cn(i) * (2*pi*k_theory/1000)).^2);
    f_theory = omega / (2*pi) * 3600;

    plot(k_theory, f_theory, 'w-', 'LineWidth', 1.3);
end
text(0.025,0.0023,'Vertical Mode: 1 3 10','FontSize',12,'Color','w','fontname','Arial');

yline(f0_cph, 'w--', '\itf', 'LineWidth', 1.3,'fontname','Arial','FontSize',14);
box on
set(AH4,'layer','top');

%%%% SF2
AH5=axes('Parent',gcf,'unit', 'centimeters','Position',[3 1.5 7 5.5]);

colorlist_SF(1,:)=[0.85 0.1 0.1];
colorlist_SF(2,:)=[0.7 0.7 0.7];
colorlist_SF(3,:)=[0, 0, 0];
colorlist_SF(4,:)=[0 0.45 0.74];


p1 = loglog(r, S_total_drifters, 'Color', colorlist_SF(3,:), 'LineWidth', 1.5, 'DisplayName', 'Drifters'); hold on
p2 = loglog(r, S_total_ADCP, 'Color', colorlist_SF(4,:), 'LineWidth', 1.2, 'DisplayName', 'ADCP'); 
p3 = loglog(r, S_total_tide, 'Color',colorlist_SF(1,:) , 'LineWidth', 1.5, 'DisplayName', 'Tide Model'); 
p4 = loglog(r, S_total_notide, 'Color', colorlist_SF(2,:), 'LineWidth', 1.2, 'DisplayName', 'Notide Model'); 
x_l = logspace(log10(min(r)), log10(max(r)), 100);
C1 = 0.2e-4; % 调整这个常数
C3 = 1.2e-2;
l_r2 = C1 * x_l.^2;
l_r1 = C3 * x_l.^(1);
loglog(x_l, l_r2, '--k', 'handlevisibility','off', 'LineWidth', 1);
loglog(x_l, l_r1, '--k', 'handlevisibility','off', 'LineWidth', 1);
ylim([1e-4 1])
xlim([1 1e2])
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
ylabel('D2(r) (m^2s^{-2})','fontsize',11,'fontname','Arial');
tit=title('\bf(d) \rmSecond-Order Structure Function','fontsize',12,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.55 5.5 0])
set(AH5,'fontsize',11,'fontname','Arial','linewidth',1);
text(22,0.5,'r^1','fontsize',13,'fontname','Arial')
text(56,0.03,'r^2','fontsize',13,'fontname','Arial')
lll=legend('show', 'Location', 'northwest','FontSize', 11);
set(lll, 'unit', 'centimeters','Position', [19.2 3 3.5 2])

AH6=axes('Parent',gcf,'unit', 'centimeters','Position',[12 1.5 7 5.5]);
loglog(r, D3_drifters, 'Color', colorlist_SF(3,:),'LineStyle','-' , 'LineWidth', 1.2); hold on
loglog(r, -D3_drifters, 'Color', colorlist_SF(3,:),'LineStyle','--' , 'LineWidth', 1.2); 
loglog(r, D3_ADCP, 'Color', colorlist_SF(4,:),'LineStyle','-' , 'LineWidth', 1.2); 
loglog(r, -D3_ADCP, 'Color', colorlist_SF(4,:),'LineStyle','--' , 'LineWidth', 1.2); 
loglog(r, D3_tide, 'Color',colorlist_SF(1,:),'LineStyle','-' , 'LineWidth', 1.5); 
loglog(r, -D3_tide, 'Color',colorlist_SF(1,:),'LineStyle','--' , 'LineWidth', 1.5); 
loglog(r, D3_notide, 'Color', colorlist_SF(2,:),'LineStyle','-' , 'LineWidth', 1.2); 
loglog(r, -D3_notide, 'Color', colorlist_SF(2,:),'LineStyle','--' , 'LineWidth', 1.2); 

% x_l = logspace(log10(min(r)), log10(max(r)), 100);
% C1 = 0.2e-4; % 调整这个常数
% C3 = 1.2e-2;
% l_r2 = C1 * x_l.^2;
% l_r1 = C3 * x_l.^(1);
% loglog(x_l, l_r2, '--k', 'handlevisibility','off', 'LineWidth', 1);
% loglog(x_l, l_r1, '--k', 'handlevisibility','off', 'LineWidth', 1);
text(8,0.1e-5,'Solid: Positive Values','fontsize',11,'fontname','Arial');
text(5.2,0.03e-5,'Dashed: Negative Values','fontsize',11,'fontname','Arial');
ylim([1e-7 1.5e-1])
xlim([1 1e2])
xlabel('\itr \rm(km)','fontsize',11,'fontname','Arial')
ylabel('D3(r) (m^3s^{-3})','fontsize',11,'fontname','Arial');
tit=title('\bf(e) \rmThird-Order Structure Function','fontsize',11,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [3.3 5.5 0])
set(AH6,'fontsize',11,'fontname','Arial','linewidth',1);
% text(22,0.5,'r^1','fontsize',13,'fontname','Arial')
% text(56,0.03,'r^2','fontsize',13,'fontname','Arial')

print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig 1.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi