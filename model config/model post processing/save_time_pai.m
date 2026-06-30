
clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
tide_lp_file='I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_filt_order4.nc';
tide_stress_file='I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_tide_3layers_stress_filt_order4.nc';
x_rho=ncread(tide_lp_file,'x_rho');y_rho=ncread(tide_lp_file,'y_rho');
depth=ncread(tide_lp_file,'depth');

filt_threshold=[1 2 5:5:200];

step1=0;
for ii=240:647
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    u_filt=ncread(tide_lp_file,'u_lp',[1 1 1 ii],[inf inf inf 1]);
    v_filt=ncread(tide_lp_file,'v_lp',[1 1 1 ii],[inf inf inf 1]);
    w_filt=ncread(tide_lp_file,'w_lp',[1 1 1 ii],[inf inf inf 1]);

    uu_filt=squeeze(ncread(tide_stress_file,'uu_filt',[1 1 2 ii],[inf inf 1 1]));
    vv_filt=squeeze(ncread(tide_stress_file,'vv_filt',[1 1 2 ii],[inf inf 1 1]));
    uv_filt=squeeze(ncread(tide_stress_file,'uv_filt',[1 1 2 ii],[inf inf 1 1]));
    uw_filt=squeeze(ncread(tide_stress_file,'uw_filt',[1 1 2 ii],[inf inf 1 1]));
    vw_filt=squeeze(ncread(tide_stress_file,'vw_filt',[1 1 2 ii],[inf inf 1 1]));

    tao_xx=uu_filt-squeeze(u_filt(:,:,2)).*squeeze(u_filt(:,:,2));
    tao_yy=vv_filt-squeeze(v_filt(:,:,2)).*squeeze(v_filt(:,:,2));  
    tao_xy=uv_filt-squeeze(u_filt(:,:,2)).*squeeze(v_filt(:,:,2));
    tao_xz=uw_filt-squeeze(u_filt(:,:,2)).*squeeze(w_filt(:,:,2));
    tao_yz=vw_filt-squeeze(v_filt(:,:,2)).*squeeze(w_filt(:,:,2));

    [du_filtdx,du_filtdy]=model_gradient(x_rho,y_rho,u_filt(:,:,2));
    [dv_filtdx,dv_filtdy]=model_gradient(x_rho,y_rho,v_filt(:,:,2));

    du_filtdz=(u_filt(:,:,3)-u_filt(:,:,1))./(depth(3)-depth(1));
    dv_filtdz=(v_filt(:,:,3)-v_filt(:,:,1))./(depth(3)-depth(1));

    % 注意我们关注的是低频流动获得能量，负号取消！
    pai_tide_15m_V= (tao_xz.*du_filtdz+tao_yz.*dv_filtdz);
    pai_tide_15m_H= (tao_xx.*du_filtdx+tao_xy.*(du_filtdy+dv_filtdx)+tao_yy.*dv_filtdy);
    pai_tide_15m_all= pai_tide_15m_V+pai_tide_15m_H;

    parfor jj=1:length(filt_threshold)-1
        pai_tide_15m_V_bp=filt2(pai_tide_15m_V,0.5,[filt_threshold(jj) filt_threshold(jj+1)],'bp');
        pai_tide_15m_H_bp=filt2(pai_tide_15m_H,0.5,[filt_threshold(jj) filt_threshold(jj+1)],'bp');
        pai_tide_15m_all_bp=filt2(pai_tide_15m_all,0.5,[filt_threshold(jj) filt_threshold(jj+1)],'bp');
        pai_tide_15m_V_line(step1,jj)=nanmean(pai_tide_15m_V_bp,'all');
        pai_tide_15m_H_line(step1,jj)=nanmean(pai_tide_15m_H_bp,'all');
        pai_tide_15m_all_line(step1,jj)=nanmean(pai_tide_15m_all_bp,'all');
    end
end
band_center = 0.5*(filt_threshold(1:end-1)+filt_threshold(2:end));
save('C:\Users\Lenovo\Desktop\结构函数\pai_time.mat','band_center','pai_tide_15m_V_line','pai_tide_15m_H_line','pai_tide_15m_all_line')



clear
cd 'C:\Users\Lenovo\Desktop\结构函数'
notide_lp_file='I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_filt_order4.nc';
notide_stress_file='I:\NCS_wu_0.5km_13m_15m_17m_all\NCS_wu_0.5km_15m_notide_3layers_stress_filt_order4.nc';
x_rho=ncread(notide_lp_file,'x_rho');y_rho=ncread(notide_lp_file,'y_rho');
depth=ncread(notide_lp_file,'depth');

filt_threshold=[1 2 5:5:200];

step1=0;
for ii=240:647
    step1=step1+1;
    disp(['step =',num2str(step1)]);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    u_filt=ncread(notide_lp_file,'u_lp',[1 1 1 ii],[inf inf inf 1]);
    v_filt=ncread(notide_lp_file,'v_lp',[1 1 1 ii],[inf inf inf 1]);
    w_filt=ncread(notide_lp_file,'w_lp',[1 1 1 ii],[inf inf inf 1]);

    uu_filt=squeeze(ncread(notide_stress_file,'uu_filt',[1 1 2 ii],[inf inf 1 1]));
    vv_filt=squeeze(ncread(notide_stress_file,'vv_filt',[1 1 2 ii],[inf inf 1 1]));
    uv_filt=squeeze(ncread(notide_stress_file,'uv_filt',[1 1 2 ii],[inf inf 1 1]));
    uw_filt=squeeze(ncread(notide_stress_file,'uw_filt',[1 1 2 ii],[inf inf 1 1]));
    vw_filt=squeeze(ncread(notide_stress_file,'vw_filt',[1 1 2 ii],[inf inf 1 1]));

    tao_xx=uu_filt-squeeze(u_filt(:,:,2)).*squeeze(u_filt(:,:,2));
    tao_yy=vv_filt-squeeze(v_filt(:,:,2)).*squeeze(v_filt(:,:,2));  
    tao_xy=uv_filt-squeeze(u_filt(:,:,2)).*squeeze(v_filt(:,:,2));
    tao_xz=uw_filt-squeeze(u_filt(:,:,2)).*squeeze(w_filt(:,:,2));
    tao_yz=vw_filt-squeeze(v_filt(:,:,2)).*squeeze(w_filt(:,:,2));

    [du_filtdx,du_filtdy]=model_gradient(x_rho,y_rho,u_filt(:,:,2));
    [dv_filtdx,dv_filtdy]=model_gradient(x_rho,y_rho,v_filt(:,:,2));

    du_filtdz=(u_filt(:,:,3)-u_filt(:,:,1))./(depth(3)-depth(1));
    dv_filtdz=(v_filt(:,:,3)-v_filt(:,:,1))./(depth(3)-depth(1));

    % 注意我们关注的是低频流动获得能量，负号取消！
    pai_notide_15m_V= (tao_xz.*du_filtdz+tao_yz.*dv_filtdz);
    pai_notide_15m_H= (tao_xx.*du_filtdx+tao_xy.*(du_filtdy+dv_filtdx)+tao_yy.*dv_filtdy);
    pai_notide_15m_all= pai_notide_15m_V+pai_notide_15m_H;

    parfor jj=1:length(filt_threshold)-1
        pai_notide_15m_V_bp=filt2(pai_notide_15m_V,0.5,[filt_threshold(jj) filt_threshold(jj+1)],'bp');
        pai_notide_15m_H_bp=filt2(pai_notide_15m_H,0.5,[filt_threshold(jj) filt_threshold(jj+1)],'bp');
        pai_notide_15m_all_bp=filt2(pai_notide_15m_all,0.5,[filt_threshold(jj) filt_threshold(jj+1)],'bp');
        pai_notide_15m_V_line(step1,jj)=nanmean(pai_notide_15m_V_bp,'all');
        pai_notide_15m_H_line(step1,jj)=nanmean(pai_notide_15m_H_bp,'all');
        pai_notide_15m_all_line(step1,jj)=nanmean(pai_notide_15m_all_bp,'all');
    end
end
band_center = 0.5*(filt_threshold(1:end-1)+filt_threshold(2:end));

save('C:\Users\Lenovo\Desktop\结构函数\pai_time.mat','-append','pai_notide_15m_V_line','pai_notide_15m_H_line','pai_notide_15m_all_line')

%%
load('C:\Users\Lenovo\Desktop\结构函数\pai_time.mat')
pai_tide_15m_all_line_mean   = nanmean(pai_tide_15m_all_line,1);
pai_notide_15m_all_line_mean = nanmean(pai_notide_15m_all_line,1);
figure
plot(1./band_center,pai_tide_15m_all_line_mean, 'r-', 'LineWidth', 1.5);hold on
plot(1./band_center,pai_notide_15m_all_line_mean, 'b-', 'LineWidth', 1.5);hold on
set(gca,'xscale','log')
grid on;
xlim([1e-2 1]);
ylim([-3e-11 3e-11])

load('C:\Users\Lenovo\Desktop\结构函数\wb_spec.mat')
figure;
plot(k_center, wb_spec_tide, 'r-', 'LineWidth', 1.5); hold on;
plot(k_center, wb_spec_notide, 'b-', 'LineWidth', 1.5);
xlabel('Wavenumber (cycles/km)');
ylabel('1D cospectrum');
legend('tide','notide');
grid on;
set(gca,'XScale','log');
xlim([1e-2 1])