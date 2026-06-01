%%% Fig S1 model config
clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
%read parent grid
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
grdname_parent='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc';
h_parent=ncread(grdname_parent,'h');f_parent=ncread(grdname_parent,'f');
lon_rho_parent= ncread(grdname_parent,'lon_rho');   lat_rho_parent= ncread(grdname_parent,'lat_rho');   % X/Y rho  :km
x_rho_parent= ncread(grdname_parent,'x_rho');   y_rho_parent= ncread(grdname_parent,'y_rho');   % X/Y rho  :km


N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='F:\NCS_wu_7.5km\7.5km output_70layers_from2.10';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

for ii=57
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    u_u=ncread(filename,'u');
    v_v=ncread(filename,'v');
    for iii=1:size(u_u,3)
        u_m(:,:,iii)=v2rho_2d(u_u(:,:,iii));
        v_m(:,:,iii)=u2rho_2d(v_v(:,:,iii));
    end
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h_parent,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(:,:,:);
    z_rho=permute(z_rho,[2 3 1]);

    u_parent=interp_model_z(u_m,x_rho_parent,y_rho_parent,z_rho,z_need);
    v_parent=interp_model_z(v_m,x_rho_parent,y_rho_parent,z_rho,z_need);
end

[dudx,dudy]=model_gradient(x_rho_parent,y_rho_parent,u_parent);
[dvdx,dvdy]=model_gradient(x_rho_parent,y_rho_parent,v_parent);
Ro_parent=(dvdx-dudy)./f_parent;


clear u_u v_v u_m v_m

grdname_child1='H:\NCS_wu_1.5km\croco_grd.nc.1';
h_child1=ncread(grdname_child1,'h');f_child1=ncread(grdname_child1,'f');
lon_rho_child1= ncread(grdname_child1,'lon_rho');   lat_rho_child1= ncread(grdname_child1,'lat_rho');   % X/Y rho  :km
x_rho_child1= ncread(grdname_child1,'x_rho');   y_rho_child1= ncread(grdname_child1,'y_rho');   % X/Y rho  :km


N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='H:\NCS_wu_1.5km\tide';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

for ii=4
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    u_u=ncread(filename,'u');
    v_v=ncread(filename,'v');
    for iii=1:size(u_u,3)
        u_m(:,:,iii)=v2rho_2d(u_u(:,:,iii));
        v_m(:,:,iii)=u2rho_2d(v_v(:,:,iii));
    end
    w=ncread(filename,'w');
    zeta_child1_tide=ncread(filename,'zeta');
    z_rho1=zlevs(h_child1,zeta_child1_tide,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(:,:,:);
    z_rho=permute(z_rho,[2 3 1]);

    u_child1=interp_model_z(u_m,x_rho_child1,y_rho_child1,z_rho,z_need);
    v_child1=interp_model_z(v_m,x_rho_child1,y_rho_child1,z_rho,z_need);
    w_child1_tide=interp_model_z(w,x_rho_child1,y_rho_child1,z_rho,z_need);
end
KE_child1_tide=(u_child1.^2+v_child1.^2)/2;

[dudx,dudy]=model_gradient(x_rho_child1,y_rho_child1,u_child1);
[dvdx,dvdy]=model_gradient(x_rho_child1,y_rho_child1,v_child1);
Ro_child1_tide=(dvdx-dudy)./f_child1;

clear u_u v_v u_m v_m

rpath1='H:\NCS_wu_1.5km\notide';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

for ii=4
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    u_u=ncread(filename,'u');
    v_v=ncread(filename,'v');
    for iii=1:size(u_u,3)
        u_m(:,:,iii)=v2rho_2d(u_u(:,:,iii));
        v_m(:,:,iii)=u2rho_2d(v_v(:,:,iii));
    end
    w=ncread(filename,'w');
    zeta_child1_notide=ncread(filename,'zeta');
    z_rho1=zlevs(h_child1,zeta_child1_notide,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(:,:,:);
    z_rho=permute(z_rho,[2 3 1]);

    u_child1=interp_model_z(u_m,x_rho_child1,y_rho_child1,z_rho,z_need);
    v_child1=interp_model_z(v_m,x_rho_child1,y_rho_child1,z_rho,z_need);
    w_child1_notide=interp_model_z(w,x_rho_child1,y_rho_child1,z_rho,z_need);

end

KE_child1_notide=(u_child1.^2+v_child1.^2)/2;
[dudx,dudy]=model_gradient(x_rho_child1,y_rho_child1,u_child1);
[dvdx,dvdy]=model_gradient(x_rho_child1,y_rho_child1,v_child1);
Ro_child1_notide=(dvdx-dudy)./f_child1;

clear u_u v_v u_m v_m


grdname_child2='F:\NCS_wu_7.5km\croco file_70layers\CROCO_FILES\croco_grd.nc.2';
h_child2=ncread(grdname_child2,'h');f_child2=ncread(grdname_child2,'f');
lon_rho_child2= ncread(grdname_child2,'lon_rho');   lat_rho_child2= ncread(grdname_child2,'lat_rho');   % X/Y rho  :km
x_rho_child2= ncread(grdname_child2,'x_rho');   y_rho_child2= ncread(grdname_child2,'y_rho');   % X/Y rho  :km

N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='F:\NCS_wu_0.5km_tide_hourly_output';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

for ii=180
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    w_m=ncread(filename,'w');
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h_child2,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(41:70,:,:);
    z_rho=permute(z_rho,[2 3 1]);
    w_child2_tide=interp_model_z(w_m,x_rho_child2,y_rho_child2,z_rho,z_need);
end


N= 70; theta_s= 5; theta_b= 1; hc= 50;vtransform= 2.;
rho_r=1025;g=9.8;

rpath1='F:\NCS_wu_0.5km_notide_hourly_output';
filelist=dir(fullfile(rpath1,'*avg*.nc'));
m1=length(filelist);
z_need=[-15];

for ii=180
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    w_m=ncread(filename,'w');
    zeta=ncread(filename,'zeta');
    z_rho1=zlevs(h_child2,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    z_rho=z_rho1(41:70,:,:);
    z_rho=permute(z_rho,[2 3 1]);
    w_child2_notide=interp_model_z(w_m,x_rho_child2,y_rho_child2,z_rho,z_need);
end

%% figure
figure
set(gcf, 'unit', 'centimeters','position', [1 1 42 23]);
%%%%% 第一个大子图
AH1 = axes('Parent', gcf, 'unit', 'centimeters', 'Position', [1 5 13 13]);
lon_mint=105;lon_maxt=130;
lat_mint=5;lat_maxt=30;
m_proj('miller','lon',[lon_mint lon_maxt],'lat',[lat_mint lat_maxt]);
m_pcolor(lon_rho_parent,lat_rho_parent,Ro_parent);shading interp;hold on
m_gshhs_h('patch', [222 204 176]/255, 'edgecolor', 'none'); 
hold on;
m_gshhs_h('color', 'k', 'linewidth', 0.5); 
m_grid('box','on','linestyle','none','xtick',lon_mint+5:5:lon_maxt,'ytick',lat_mint:5:lat_maxt,'fontsize',11,'fontname','Arial','linewidth',1);  % 'XaxisLocation','top',
lon_box = [108, 125.6, 125.6, 108, 108]; % 首尾重合形成闭合
lat_box = [8.3, 8.3, 25.3, 25.3, 8.3];
% 用m_line绘制，确保线条遵循地图投影
m_line(lon_box, lat_box, ...
    'Color', 'k', ...
    'LineWidth', 2, ...
    'LineStyle', '-'); 
m_text(105.25,29.25,'Parent grid','fontsize',15,'fontname','Arial','BackgroundColor', 'w','EdgeColor', 'none')
caxis([-1 1])
colortable = textread('NCV_blue_red.txt');
colormap(AH1,colortable);
hhh = colorbar;
set(hhh,'fontsize',11,'fontname','Arial','location','southoutside');
set(hhh, 'unit', 'centimeters','Position', [4.5 2 5 0.4]);
set(hhh,'ticks',[-1 0 1])
title(hhh,'Ro','fontsize',11,'fontname','Arial')
annotation(gcf,'arrow',[0.275922619047619 0.365376984126984],...
    [0.671961956521739 0.89268115942029],'LineWidth',2,'HeadStyle','deltoid');
annotation(gcf,'arrow',[0.275922619047619 0.360337301587302],...
    [0.288891304347826 0.103532608695652],'LineWidth',2,'HeadStyle','deltoid');
tit=title('\bf(a)\rm Rossby number','fontsize',15,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.3 13.1 0])

AH2= axes('Parent', gcf, 'unit', 'centimeters', 'Position', [16 13 10 9]);
lon_mint=108;lon_maxt=125.6;
lat_mint=8.3;lat_maxt=25.3;
m_proj('miller','lon',[lon_mint lon_maxt],'lat',[lat_mint lat_maxt]);
m_pcolor(lon_rho_child1,lat_rho_child1,Ro_child1_tide);shading interp;hold on
m_gshhs_h('patch', [222 204 176]/255, 'edgecolor', 'none'); 
hold on;
m_gshhs_h('color', 'k', 'linewidth', 0.5); 
m_grid('box','on','linestyle','none','xtick',lon_mint:5:lon_maxt,'ytick',lat_mint+0.7:5:lat_maxt,'fontsize',11,'fontname','Arial','linewidth',1);  % 'XaxisLocation','top',
lon_box = [110.6, 124, 124, 110.6, 110.6]; % 首尾重合形成闭合
lat_box = [16.1, 16.1, 23, 23,16.1];
% 用m_line绘制，确保线条遵循地图投影
m_line(lon_box, lat_box, ...
    'Color', 'k', ...
    'LineWidth', 2, ...
    'LineStyle', '-'); 
m_text(108.25,24.55,'Child1 grid Tide model','fontsize',15,'fontname','Arial','BackgroundColor', 'w','EdgeColor', 'none')
caxis([-1 1])
colortable = textread('NCV_blue_red.txt');
colormap(AH2,colortable);
tit=title('\bf(b)\rm Rossby number','fontsize',15,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.3 9.1 0])

% hhh = colorbar;
% set(hhh,'fontsize',11,'fontname','Arial','location','southoutside');
% set(hhh, 'unit', 'centimeters','Position', [4.5 2 5 0.4]);
% set(hhh,'ticks',[-1 0 1])
% title(hhh,'W','fontsize',11,'fontname','Arial')

AH3= axes('Parent', gcf, 'unit', 'centimeters', 'Position', [16 2 10 9]);
lon_mint=108;lon_maxt=125.6;
lat_mint=8.3;lat_maxt=25.3;
m_proj('miller','lon',[lon_mint lon_maxt],'lat',[lat_mint lat_maxt]);
m_pcolor(lon_rho_child1,lat_rho_child1,Ro_child1_notide);shading interp;hold on
m_gshhs_h('patch', [222 204 176]/255, 'edgecolor', 'none'); 
hold on;
m_gshhs_h('color', 'k', 'linewidth', 0.5); 
m_grid('box','on','linestyle','none','xtick',lon_mint:5:lon_maxt,'ytick',lat_mint+0.7:5:lat_maxt,'fontsize',11,'fontname','Arial','linewidth',1);  % 'XaxisLocation','top',
lon_box = [110.6, 124, 124, 110.6, 110.6]; % 首尾重合形成闭合
lat_box = [16.1, 16.1, 23, 23,16.1];
% 用m_line绘制，确保线条遵循地图投影
m_line(lon_box, lat_box, ...
    'Color', 'k', ...
    'LineWidth', 2, ...
    'LineStyle', '-'); 
m_text(108.25,24.55,'Child1 grid Notide model','fontsize',15,'fontname','Arial','BackgroundColor', 'w','EdgeColor', 'none')
caxis([-1 1])
colortable = textread('NCV_blue_red.txt');
colormap(AH3,colortable);
% hhh = colorbar;
% set(hhh,'fontsize',11,'fontname','Arial','location','southoutside');
% set(hhh, 'unit', 'centimeters','Position', [4.5 2 5 0.4]);
% set(hhh,'ticks',[-1 0 1])
% title(hhh,'W','fontsize',11,'fontname','Arial')
tit=title('\bf(c)\rm Rossby number','fontsize',15,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.3 9.1 0])

AH4= axes('Parent', gcf, 'unit', 'centimeters', 'Position', [29 13 12 9]);
lon_mint=110.6;lon_maxt=124;
lat_mint=16.1;lat_maxt=23;
m_proj('miller','lon',[lon_mint lon_maxt],'lat',[lat_mint lat_maxt]);
m_pcolor(lon_rho_child2,lat_rho_child2,w_child2_tide);shading interp;hold on
m_gshhs_h('patch', [222 204 176]/255, 'edgecolor', 'none'); 
hold on;
m_gshhs_h('color', 'k', 'linewidth', 0.5); 
m_grid('box','on','linestyle','none','xtick',lon_mint+0.4:3:lon_maxt,'ytick',lat_mint+0.9:3:lat_maxt,'fontsize',11,'fontname','Arial','linewidth',1);  % 'XaxisLocation','top',
annotation(gcf,'arrow',[0.622400793650794 0.657678571428571],...
    [0.764990942028985 0.764990942028985],'LineWidth',10,'HeadWidth',25,...
    'HeadStyle','deltoid',...
    'HeadLength',25,'color',[0.18,0.75,0.94]);
m_text(110.75,22.57,'Child2 grid Tide model','fontsize',15,'fontname','Arial','BackgroundColor', 'w','EdgeColor', 'none')
caxis([-8e-4 8e-4])
colortable = textread('cmp_b2r.txt');
colormap(AH4,colortable);
hhh = colorbar;
set(hhh,'fontsize',13,'fontname','Arial','location','southoutside');
set(hhh, 'unit', 'centimeters','Position', [32.5 12 5 0.4]);
% set(hhh,'ticks',[-1 0 1])
title(hhh,'W (m/s)','fontsize',13,'fontname','Arial')
tit=title('\bf(d)\rm Vertical velocity','fontsize',15,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.3 6.6 0])



AH5= axes('Parent', gcf, 'unit', 'centimeters', 'Position', [29 2 12 9]);
lon_mint=110.6;lon_maxt=124;
lat_mint=16.1;lat_maxt=23;
m_proj('miller','lon',[lon_mint lon_maxt],'lat',[lat_mint lat_maxt]);
m_pcolor(lon_rho_child2,lat_rho_child2,w_child2_notide);shading interp;hold on
m_gshhs_h('patch', [222 204 176]/255, 'edgecolor', 'none'); 
hold on;
m_gshhs_h('color', 'k', 'linewidth', 0.5); 
m_grid('box','on','linestyle','none','xtick',lon_mint+0.4:3:lon_maxt,'ytick',lat_mint+0.9:3:lat_maxt,'fontsize',11,'fontname','Arial','linewidth',1);  % 'XaxisLocation','top',
annotation(gcf,'arrow',[0.624290674603175 0.659568452380952],...
    [0.286440217391304 0.286440217391304],'LineWidth',10,'HeadWidth',25,...
    'HeadStyle','deltoid',...
    'HeadLength',25,'color',[0.18,0.75,0.94]);
m_text(110.75,22.57,'Child2 grid Notide model','fontsize',15,'fontname','Arial','BackgroundColor', 'w','EdgeColor', 'none')
caxis([-8e-4 8e-4])
colortable = textread('cmp_b2r.txt');
colormap(AH5,colortable);
tit=title('\bf(e)\rm Vertical velocity','fontsize',15,'fontname','Arial');
set(tit, 'unit', 'centimeters','Position', [2.3 6.6 0])


print(gcf, 'D:\学花招的素材\初稿\wave-flow interation influence flow energy cascade\figure\Fig S1.tiff', '-dtiff', '-r600')   % TIFF, 600 dpi