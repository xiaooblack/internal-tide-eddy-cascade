%% model with errobar
clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
r=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','r');
sum_DLL_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','sum_DLL',[240 1], [408 inf]);
sum_DTT_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','sum_DTT',[240 1], [408 inf]);
sum_D3_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','sum_D3',[240 1], [408 inf]);
counts_tide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_tide_origin.nc','counts',[240 1], [408 inf]);

S_long_tide=sum(sum_DLL_tide,1)./sum(counts_tide,1);
S_trans_tide=sum(sum_DTT_tide,1)./sum(counts_tide,1);
S_total_tide=S_long_tide+S_trans_tide;
D3_tide=sum(sum_D3_tide,1)./sum(counts_tide,1);
N_tide=sum(counts_tide,1);

N = 48;% 每组元素个数
n = ceil(size(sum_DLL_tide,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
sum_DLL_tide_block = zeros(n, length(r));
sum_DTT_tide_block = zeros(n, length(r));
sum_D3_tide_block = zeros(n, length(r));
counts_tide_block = zeros(n, length(r));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(sum_DLL_tide,1));
    sum_DLL_tide_block(i,:) = sum(sum_DLL_tide(startIdx:endIdx,:));
    sum_DTT_tide_block(i,:) = sum(sum_DTT_tide(startIdx:endIdx,:));
    sum_D3_tide_block(i,:) = sum(sum_D3_tide(startIdx:endIdx,:));
    counts_tide_block(i,:) = sum(counts_tide(startIdx:endIdx,:));
end


sum_DLL_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','sum_DLL',[240 1], [408 inf]);
sum_DTT_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','sum_DTT',[240 1], [408 inf]);
sum_D3_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','sum_D3',[240 1], [408 inf]);
counts_notide=ncread('I:\NCS_wu_0.5km_15m_all\SF_NCS_wu_0.5km_15m_notide_origin.nc','counts',[240 1], [408 inf]);

S_long_notide=sum(sum_DLL_notide,1)./sum(counts_notide,1);
S_trans_notide=sum(sum_DTT_notide,1)./sum(counts_notide,1);
S_total_notide=S_long_notide+S_trans_notide;
D3_notide=sum(sum_D3_notide,1)./sum(counts_notide,1);
N_notide=sum(counts_notide,1);

% 24h block 每24个时刻合成一个

N = 48;% 每组元素个数
n = ceil(size(sum_DLL_notide,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
sum_DLL_notide_block = zeros(n, length(r));
sum_DTT_notide_block = zeros(n, length(r));
sum_D3_notide_block = zeros(n, length(r));
counts_notide_block = zeros(n, length(r));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(sum_DLL_notide,1));
    sum_DLL_notide_block(i,:) = sum(sum_DLL_notide(startIdx:endIdx,:));
    sum_DTT_notide_block(i,:) = sum(sum_DTT_notide(startIdx:endIdx,:));
    sum_D3_notide_block(i,:) = sum(sum_D3_notide(startIdx:endIdx,:));
    counts_notide_block(i,:) = sum(counts_notide(startIdx:endIdx,:));
end


[S_long_tide_low,S_long_tide_high]=block_bootstrap(sum_DLL_tide_block,counts_tide_block,1000);
[S_trans_tide_low,S_trans_tide_high]=block_bootstrap(sum_DTT_tide_block,counts_tide_block,1000);
S_total_tide_low=S_long_tide_low+S_trans_tide_low;
S_total_tide_high=S_long_tide_high+S_trans_tide_high;
[D3_tide_low,D3_tide_high]=block_bootstrap(sum_D3_tide_block,counts_tide_block,1000);

[S_long_notide_low,S_long_notide_high]=block_bootstrap(sum_DLL_notide_block,counts_notide_block,1000);
[S_trans_notide_low,S_trans_notide_high]=block_bootstrap(sum_DTT_notide_block,counts_notide_block,1000);
S_total_notide_low=S_long_notide_low+S_trans_notide_low;
S_total_notide_high=S_long_notide_high+S_trans_notide_high;
[D3_notide_low,D3_notide_high]=block_bootstrap(sum_D3_notide_block,counts_notide_block,1000);

figure
plot(r,S_total_tide,'-k');hold on
plot(r,S_total_tide_low,'-r');hold on
plot(r,S_total_tide_high,'-b');hold on

figure
plot(r,S_total_notide,'-k');hold on
plot(r,S_total_notide_low,'-r');hold on
plot(r,S_total_notide_high,'-b');hold on

figure
plot(r,D3_tide,'-k');hold on
plot(r,D3_tide_low,'-r');hold on
plot(r,D3_tide_high,'-b');hold on

figure
plot(r,D3_notide,'-k');hold on
plot(r,D3_notide_low,'-r');hold on
plot(r,D3_notide_high,'-b');hold on


save('SF_origion_model_witherrobar.mat','r','S_total_tide','D3_tide','S_total_notide','D3_notide','S_total_tide_low','S_total_tide_high','D3_tide_low','D3_tide_high',...
    'S_total_notide_low','S_total_notide_high','D3_notide_low','D3_notide_high','N_tide','N_notide');
%% drifters with errobar
% clear
% cd 'C:\Users\Lenovo\Desktop\结构函数'
% rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
% read outputfile
% 先试试0415 15m
rpath1='C:\Users\Lenovo\Desktop\结构函数\表漂观测一二三批15min插值';
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
    drifter(ii).time=ncread(filename,'time');
    drifter(ii).lon=ncread(filename,'lon');
    drifter(ii).lat=ncread(filename,'lat');
    drifter(ii).u=ncread(filename,'u');
    drifter(ii).v=ncread(filename,'v');
    drifter(ii).u=filt1('lp',drifter(ii).u,'Ts',15,'Tc',120);
    drifter(ii).v=filt1('lp',drifter(ii).v,'Ts',15,'Tc',120);
    time_start_list(ii)=drifter(ii).time(1);
    time_end_list(ii)=drifter(ii).time(end);
end

time_all=[798163200:900:799448400];
R=6371393;
for ii=1:length(time_all)
    aaa=0;
    for jj=1:length(drifter)
        id_t=drifter(jj).time==time_all(ii);
        if sum(id_t)~=0
            aaa=aaa+1;
            tdrifter(ii).lon(aaa)=drifter(jj).lon(id_t);
            tdrifter(ii).lat(aaa)=drifter(jj).lat(id_t);
            tdrifter(ii).u(aaa)=drifter(jj).u(id_t);
            tdrifter(ii).v(aaa)=drifter(jj).v(id_t);
            tdrifter(ii).x(aaa)=deg2rad(drifter(jj).lon(id_t)).*cosd(drifter(jj).lat(id_t)).*R;
            tdrifter(ii).y(aaa)=deg2rad(drifter(jj).lat(id_t)).*R;
        else
        end
    end
end

for ii=1:length(tdrifter)
    [r, sum_DLL_drifters(ii,:), sum_DTT_drifters(ii,:), sum_D3_drifters(ii,:),counts_drifters(ii,:)] = calc_binned_accum_SF_parallel_single(tdrifter(ii).u, tdrifter(ii).v, tdrifter(ii).x./1000, tdrifter(ii).y./1000,1,100);
end

S_long_drifters=sum(sum_DLL_drifters,1)./sum(counts_drifters,1);
S_trans_drifters=sum(sum_DTT_drifters,1)./sum(counts_drifters,1);
S_total_drifters=S_long_drifters+S_trans_drifters;
D3_drifters=sum(sum_D3_drifters,1)./sum(counts_drifters,1);
N_drifters=sum(counts_drifters,1);

% 24h block 每96个时刻合成一个

N = 48;% 每组元素个数
n = ceil(size(sum_DLL_drifters,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
sum_DLL_block = zeros(n, length(r));
sum_DTT_block = zeros(n, length(r));
sum_D3_block = zeros(n, length(r));
counts_block = zeros(n, length(r));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(sum_DLL_drifters,1));
    sum_DLL_block(i,:) = sum(sum_DLL_drifters(startIdx:endIdx,:));
    sum_DTT_block(i,:) = sum(sum_DTT_drifters(startIdx:endIdx,:));
    sum_D3_block(i,:) = sum(sum_D3_drifters(startIdx:endIdx,:));
    counts_block(i,:) = sum(counts_drifters(startIdx:endIdx,:));
end

[S_long_drifters_low,S_long_drifters_high]=block_bootstrap(sum_DLL_block,counts_block,1000);
[S_trans_drifters_low,S_trans_drifters_high]=block_bootstrap(sum_DTT_block,counts_block,1000);
S_total_drifters_low=S_long_drifters_low+S_trans_drifters_low;
S_total_drifters_high=S_long_drifters_high+S_trans_drifters_high;
[D3_drifters_low,D3_drifters_high]=block_bootstrap(sum_D3_block,counts_block,1000);

figure
plot(r,S_long_drifters,'-k');hold on
plot(r,S_long_drifters_low,'-r');hold on
plot(r,S_long_drifters_high,'-b');hold on

figure
plot(r,D3_drifters,'-k');hold on
plot(r,D3_drifters_low,'-r');hold on
plot(r,D3_drifters_high,'-b');hold on

save('SF_origion_drifters_witherrobar.mat','r','S_total_drifters','D3_drifters',...
    'S_total_drifters_low','S_total_drifters_high','D3_drifters_low','D3_drifters_high','N_drifters');

%% ADCP with errorbar
% clear
% cd 'C:\Users\Lenovo\Desktop\结构函数'
% rmpath('E:\Roms_CROCO_source\croco_tools-v1.3\croco_tools-v1.3\UTILITIES\mask');% 使用原ncread函数
% 300k
rpath1='I:\ADCP_ship\my process\15m\等间距插值';
filelist=dir(fullfile(rpath1,'*_300k_equidistant_interpolation.nc'));
m1=length(filelist);
step1=0;

% [1 2 3 4 6 9 13 14 16] 成功
for ii=[2 4 6 7 9 10 11 12 13 15 16]
    step1=step1+1;
    filename=[filelist(ii).folder,'\',filelist(ii).name];  
    lon=ncread(filename,'lon_ADCP');
    lat=ncread(filename,'lat_ADCP');
    u300k=ncread(filename,'u300k');
    v300k=ncread(filename,'v300k');
    
    R=6371393;
    x=deg2rad(lon).*cosd(lat).*R;
    y=deg2rad(lat).*R;  

    % % 去除停船点 
    % dx_ADCP=zeros(length(x),1);dy_ADCP=zeros(length(x),1);
    % dx_ADCP(1:end-1)=x(2:end)-x(1:end-1);
    % dy_ADCP(1:end-1)=y(2:end)-y(1:end-1);
    % dx_ADCP(end)=dx_ADCP(end-1);dy_ADCP(end)=dy_ADCP(end-1);
    % 
    % v_ship=sqrt(dx_ADCP.^2+dy_ADCP.^2)/180;
    % x(v_ship<2)=[];
    % y(v_ship<2)=[];
    % u300k(v_ship<2)=[];
    % v300k(v_ship<2)=[];

    [r, sum_DLL_ADCP(step1,:), sum_DTT_ADCP(step1,:), sum_D3_ADCP(step1,:),counts_ADCP(step1,:)] = calc_binned_accum_SF_parallel(u300k, v300k, x./1000, y./1000, 1,100);
end

S_long_ADCP=sum(sum_DLL_ADCP,1)./sum(counts_ADCP,1);
S_trans_ADCP=sum(sum_DTT_ADCP,1)./sum(counts_ADCP,1);
S_total_ADCP=S_long_ADCP+S_trans_ADCP;
D3_ADCP=sum(sum_D3_ADCP,1)./sum(counts_ADCP,1);
N_ADCP=sum(counts_ADCP,1);

% 由于本来就是按照作业断面分的，不需要对块做合成
[S_long_ADCP_low,S_long_ADCP_high]=block_bootstrap(sum_DLL_ADCP,counts_ADCP,1000);
[S_trans_ADCP_low,S_trans_ADCP_high]=block_bootstrap(sum_DTT_ADCP,counts_ADCP,1000);
S_total_ADCP_low=S_long_ADCP_low+S_trans_ADCP_low;
S_total_ADCP_high=S_long_ADCP_high+S_trans_ADCP_high;
[D3_ADCP_low,D3_ADCP_high]=block_bootstrap(sum_D3_ADCP,counts_ADCP,1000);

figure
plot(r,S_long_ADCP,'-k');hold on
plot(r,S_long_ADCP_low,'-r');hold on
plot(r,S_long_ADCP_high,'-b');hold on

figure
plot(r,D3_ADCP,'-k');hold on
plot(r,D3_ADCP_low,'-r');hold on
plot(r,D3_ADCP_high-4e-4,'-b');hold on


figure
loglog(r,N_drifters);hold on
loglog(r,N_ADCP);
loglog(r,N_tide);
loglog(r,N_notide);

save('SF_origion_ADCP_witherrobar.mat','r','S_total_ADCP','D3_ADCP',...
    'S_total_ADCP_low','S_total_ADCP_high','D3_ADCP_low','D3_ADCP_high','N_ADCP');

%% lp model errorbar
% 读取的时候转置，为了匹配block的维度
r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc','r')';
sum_D3_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_D3_lp', [1 14 240], [inf 1 408]))';
sum_DLL_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_DLL_lp', [1 14 240], [inf 1 408]))';
sum_DTT_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_DTT_lp', [1 14 240], [inf 1 408]))';
sum_D3_zeta_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_D3_zeta_lp', [1 14 240], [inf 1 408]))';
sum_D2_zeta_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_D2_zeta_lp', [1 14 240], [inf 1 408]))';
counts_tide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'counts_lp', [1 14 240], [inf 1 408]))';


S3_tide_lp = sum(sum_D3_tide_lp, 1) ./ sum(counts_tide_lp, 1);
S3_zeta_tide_lp = sum(sum_D3_zeta_tide_lp, 1) ./ sum(counts_tide_lp, 1);
S2_zeta_tide_lp = sum(sum_D2_zeta_tide_lp, 1) ./ sum(counts_tide_lp, 1);
S2_tide_lp=sum(sum_DLL_tide_lp, 1) ./ sum(counts_tide_lp, 1)+sum(sum_DTT_tide_lp, 1) ./ sum(counts_tide_lp, 1);

N = 48;% 每组元素个数
n = ceil(size(sum_DLL_tide_lp,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
sum_DLL_tide_lp_block = zeros(n, length(r));
sum_DTT_tide_lp_block = zeros(n, length(r));
sum_D3_tide_lp_block = zeros(n, length(r));
sum_D2_zeta_tide_lp_block = zeros(n, length(r));
sum_D3_zeta_tide_lp_block = zeros(n, length(r));
counts_tide_lp_block = zeros(n, length(r));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(sum_DLL_tide_lp,1));
    sum_DLL_tide_lp_block(i,:) = sum(sum_DLL_tide_lp(startIdx:endIdx,:));
    sum_DTT_tide_lp_block(i,:) = sum(sum_DTT_tide_lp(startIdx:endIdx,:));
    sum_D3_tide_lp_block(i,:) = sum(sum_D3_tide_lp(startIdx:endIdx,:));
    sum_D2_zeta_tide_lp_block(i,:) = sum(sum_D2_zeta_tide_lp(startIdx:endIdx,:));
    sum_D3_zeta_tide_lp_block(i,:) = sum(sum_D3_zeta_tide_lp(startIdx:endIdx,:));
    counts_tide_lp_block(i,:) = sum(counts_tide_lp(startIdx:endIdx,:));
end

[S_long_tide_low,S_long_tide_high]=block_bootstrap(sum_DLL_tide_lp_block,counts_tide_lp_block,1000);
[S_trans_tide_low,S_trans_tide_high]=block_bootstrap(sum_DTT_tide_lp_block,counts_tide_lp_block,1000);
S2_tide_lp_low=S_long_tide_low+S_trans_tide_low;
S2_tide_lp_high=S_long_tide_high+S_trans_tide_high;
[S3_tide_lp_low,S3_tide_lp_high]=block_bootstrap(sum_D3_tide_lp_block,counts_tide_lp_block,1000);
[S2_zeta_tide_lp_low,S2_zeta_tide_lp_high]=block_bootstrap(sum_D2_zeta_tide_lp_block,counts_tide_lp_block,1000);
[S3_zeta_tide_lp_low,S3_zeta_tide_lp_high]=block_bootstrap(sum_D3_zeta_tide_lp_block,counts_tide_lp_block,1000);

r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc','r')';
sum_D3_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_D3_lp', [1 14 240], [inf 1 408]))';
sum_DLL_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_DLL_lp', [1 14 240], [inf 1 408]))';
sum_DTT_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_DTT_lp', [1 14 240], [inf 1 408]))';
sum_D3_zeta_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_D3_zeta_lp', [1 14 240], [inf 1 408]))';
sum_D2_zeta_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_D2_zeta_lp', [1 14 240], [inf 1 408]))';
counts_notide_lp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'counts_lp', [1 14 240], [inf 1 408]))';

S3_notide_lp = sum(sum_D3_notide_lp, 1) ./ sum(counts_notide_lp, 1);
S3_zeta_notide_lp = sum(sum_D3_zeta_notide_lp, 1) ./ sum(counts_notide_lp, 1);
S2_zeta_notide_lp = sum(sum_D2_zeta_notide_lp, 1) ./ sum(counts_notide_lp, 1);
S2_notide_lp=sum(sum_DLL_notide_lp, 1) ./ sum(counts_notide_lp, 1)+sum(sum_DTT_notide_lp, 1) ./ sum(counts_notide_lp, 1);

N = 48;% 每组元素个数
n = ceil(size(sum_DLL_notide_lp,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
sum_DLL_notide_lp_block = zeros(n, length(r));
sum_DTT_notide_lp_block = zeros(n, length(r));
sum_D3_notide_lp_block = zeros(n, length(r));
sum_D2_zeta_notide_lp_block = zeros(n, length(r));
sum_D3_zeta_notide_lp_block = zeros(n, length(r));
counts_notide_lp_block = zeros(n, length(r));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(sum_DLL_notide_lp,1));
    sum_DLL_notide_lp_block(i,:) = sum(sum_DLL_notide_lp(startIdx:endIdx,:));
    sum_DTT_notide_lp_block(i,:) = sum(sum_DTT_notide_lp(startIdx:endIdx,:));
    sum_D3_notide_lp_block(i,:) = sum(sum_D3_notide_lp(startIdx:endIdx,:));
    sum_D2_zeta_notide_lp_block(i,:) = sum(sum_D2_zeta_notide_lp(startIdx:endIdx,:));
    sum_D3_zeta_notide_lp_block(i,:) = sum(sum_D3_zeta_notide_lp(startIdx:endIdx,:));
    counts_notide_lp_block(i,:) = sum(counts_notide_lp(startIdx:endIdx,:));
end

[S_long_notide_low,S_long_notide_high]=block_bootstrap(sum_DLL_notide_lp_block,counts_notide_lp_block,100);
[S_trans_notide_low,S_trans_notide_high]=block_bootstrap(sum_DTT_notide_lp_block,counts_notide_lp_block,100);
S2_notide_lp_low=S_long_notide_low+S_trans_notide_low;
S2_notide_lp_high=S_long_notide_high+S_trans_notide_high;
[S3_notide_lp_low,S3_notide_lp_high]=block_bootstrap(sum_D3_notide_lp_block,counts_notide_lp_block,100);
[S2_zeta_notide_lp_low,S2_zeta_notide_lp_high]=block_bootstrap(sum_D2_zeta_notide_lp_block,counts_notide_lp_block,100);
[S3_zeta_notide_lp_low,S3_zeta_notide_lp_high]=block_bootstrap(sum_D3_zeta_notide_lp_block,counts_notide_lp_block,100);



save('SF_lp_model_witherrobar.mat','r','S2_tide_lp','S2_tide_lp_low','S2_tide_lp_high',...
    'S3_tide_lp','S3_tide_lp_low','S3_tide_lp_high',...
    'S2_zeta_tide_lp','S2_zeta_tide_lp_low','S2_zeta_tide_lp_high',...
    'S3_zeta_tide_lp','S3_zeta_tide_lp_low','S3_zeta_tide_lp_high',...
    'S2_notide_lp','S2_notide_lp_low','S2_notide_lp_high',...
    'S3_notide_lp','S3_notide_lp_low','S3_notide_lp_high',...
    'S2_zeta_notide_lp','S2_zeta_notide_lp_low','S2_zeta_notide_lp_high',...
    'S3_zeta_notide_lp','S3_zeta_notide_lp_low','S3_zeta_notide_lp_high');
%% hp model errorbar
% 读取的时候转置，为了匹配block的维度
r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc','r')';
sum_D3_tide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_D3_hp', [1 14 240], [inf 1 408]))';
sum_DLL_tide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_DLL_hp', [1 14 240], [inf 1 408]))';
sum_DTT_tide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_DTT_hp', [1 14 240], [inf 1 408]))';
sum_D3_zeta_tide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_D3_zeta_hp', [1 14 240], [inf 1 408]))';
sum_D2_zeta_tide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'sum_D2_zeta_hp', [1 14 240], [inf 1 408]))';
counts_tide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_tide_filt_order4_all.nc', ...
    'counts_hp', [1 14 240], [inf 1 408]))';


S3_tide_hp = sum(sum_D3_tide_hp, 1) ./ sum(counts_tide_hp, 1);
S3_zeta_tide_hp = sum(sum_D3_zeta_tide_hp, 1) ./ sum(counts_tide_hp, 1);
S2_zeta_tide_hp = sum(sum_D2_zeta_tide_hp, 1) ./ sum(counts_tide_hp, 1);
S2_tide_hp=sum(sum_DLL_tide_hp, 1) ./ sum(counts_tide_hp, 1)+sum(sum_DTT_tide_hp, 1) ./ sum(counts_tide_hp, 1);

N = 48;% 每组元素个数
n = ceil(size(sum_DLL_tide_hp,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
sum_DLL_tide_hp_block = zeros(n, length(r));
sum_DTT_tide_hp_block = zeros(n, length(r));
sum_D3_tide_hp_block = zeros(n, length(r));
sum_D2_zeta_tide_hp_block = zeros(n, length(r));
sum_D3_zeta_tide_hp_block = zeros(n, length(r));
counts_tide_hp_block = zeros(n, length(r));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(sum_DLL_tide_hp,1));
    sum_DLL_tide_hp_block(i,:) = sum(sum_DLL_tide_hp(startIdx:endIdx,:));
    sum_DTT_tide_hp_block(i,:) = sum(sum_DTT_tide_hp(startIdx:endIdx,:));
    sum_D3_tide_hp_block(i,:) = sum(sum_D3_tide_hp(startIdx:endIdx,:));
    sum_D2_zeta_tide_hp_block(i,:) = sum(sum_D2_zeta_tide_hp(startIdx:endIdx,:));
    sum_D3_zeta_tide_hp_block(i,:) = sum(sum_D3_zeta_tide_hp(startIdx:endIdx,:));
    counts_tide_hp_block(i,:) = sum(counts_tide_hp(startIdx:endIdx,:));
end

[S_long_tide_low,S_long_tide_high]=block_bootstrap(sum_DLL_tide_hp_block,counts_tide_hp_block,1000);
[S_trans_tide_low,S_trans_tide_high]=block_bootstrap(sum_DTT_tide_hp_block,counts_tide_hp_block,1000);
S2_tide_hp_low=S_long_tide_low+S_trans_tide_low;
S2_tide_hp_high=S_long_tide_high+S_trans_tide_high;
[S3_tide_hp_low,S3_tide_hp_high]=block_bootstrap(sum_D3_tide_hp_block,counts_tide_hp_block,1000);
[S2_zeta_tide_hp_low,S2_zeta_tide_hp_high]=block_bootstrap(sum_D2_zeta_tide_hp_block,counts_tide_hp_block,1000);
[S3_zeta_tide_hp_low,S3_zeta_tide_hp_high]=block_bootstrap(sum_D3_zeta_tide_hp_block,counts_tide_hp_block,1000);

r = ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc','r')';
sum_D3_notide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_D3_hp', [1 14 240], [inf 1 408]))';
sum_DLL_notide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_DLL_hp', [1 14 240], [inf 1 408]))';
sum_DTT_notide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_DTT_hp', [1 14 240], [inf 1 408]))';
sum_D3_zeta_notide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_D3_zeta_hp', [1 14 240], [inf 1 408]))';
sum_D2_zeta_notide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'sum_D2_zeta_hp', [1 14 240], [inf 1 408]))';
counts_notide_hp = squeeze(ncread('I:/NCS_wu_0.5km_15layer_all/SF_NCS_wu_0.5km_15layers_notide_filt_order4_all.nc', ...
    'counts_hp', [1 14 240], [inf 1 408]))';

S3_notide_hp = sum(sum_D3_notide_hp, 1) ./ sum(counts_notide_hp, 1);
S3_zeta_notide_hp = sum(sum_D3_zeta_notide_hp, 1) ./ sum(counts_notide_hp, 1);
S2_zeta_notide_hp = sum(sum_D2_zeta_notide_hp, 1) ./ sum(counts_notide_hp, 1);
S2_notide_hp=sum(sum_DLL_notide_hp, 1) ./ sum(counts_notide_hp, 1)+sum(sum_DTT_notide_hp, 1) ./ sum(counts_notide_hp, 1);

N = 48;% 每组元素个数
n = ceil(size(sum_DLL_notide_hp,1) / N);% 总组数（包括可能的不完整组）
% 预分配存储各组和
sum_DLL_notide_hp_block = zeros(n, length(r));
sum_DTT_notide_hp_block = zeros(n, length(r));
sum_D3_notide_hp_block = zeros(n, length(r));
sum_D2_zeta_notide_hp_block = zeros(n, length(r));
sum_D3_zeta_notide_hp_block = zeros(n, length(r));
counts_notide_hp_block = zeros(n, length(r));

for i = 1:n
    startIdx = (i-1)*N + 1;
    endIdx   = min(i*N, size(sum_DLL_notide_hp,1));
    sum_DLL_notide_hp_block(i,:) = sum(sum_DLL_notide_hp(startIdx:endIdx,:));
    sum_DTT_notide_hp_block(i,:) = sum(sum_DTT_notide_hp(startIdx:endIdx,:));
    sum_D3_notide_hp_block(i,:) = sum(sum_D3_notide_hp(startIdx:endIdx,:));
    sum_D2_zeta_notide_hp_block(i,:) = sum(sum_D2_zeta_notide_hp(startIdx:endIdx,:));
    sum_D3_zeta_notide_hp_block(i,:) = sum(sum_D3_zeta_notide_hp(startIdx:endIdx,:));
    counts_notide_hp_block(i,:) = sum(counts_notide_hp(startIdx:endIdx,:));
end

[S_long_notide_low,S_long_notide_high]=block_bootstrap(sum_DLL_notide_hp_block,counts_notide_hp_block,100);
[S_trans_notide_low,S_trans_notide_high]=block_bootstrap(sum_DTT_notide_hp_block,counts_notide_hp_block,100);
S2_notide_hp_low=S_long_notide_low+S_trans_notide_low;
S2_notide_hp_high=S_long_notide_high+S_trans_notide_high;
[S3_notide_hp_low,S3_notide_hp_high]=block_bootstrap(sum_D3_notide_hp_block,counts_notide_hp_block,100);
[S2_zeta_notide_hp_low,S2_zeta_notide_hp_high]=block_bootstrap(sum_D2_zeta_notide_hp_block,counts_notide_hp_block,100);
[S3_zeta_notide_hp_low,S3_zeta_notide_hp_high]=block_bootstrap(sum_D3_zeta_notide_hp_block,counts_notide_hp_block,100);



save('SF_hp_model_witherrobar.mat','r','S2_tide_hp','S2_tide_hp_low','S2_tide_hp_high',...
    'S3_tide_hp','S3_tide_hp_low','S3_tide_hp_high',...
    'S2_zeta_tide_hp','S2_zeta_tide_hp_low','S2_zeta_tide_hp_high',...
    'S3_zeta_tide_hp','S3_zeta_tide_hp_low','S3_zeta_tide_hp_high',...
    'S2_notide_hp','S2_notide_hp_low','S2_notide_hp_high',...
    'S3_notide_hp','S3_notide_hp_low','S3_notide_hp_high',...
    'S2_zeta_notide_hp','S2_zeta_notide_hp_low','S2_zeta_notide_hp_high',...
    'S3_zeta_notide_hp','S3_zeta_notide_hp_low','S3_zeta_notide_hp_high');