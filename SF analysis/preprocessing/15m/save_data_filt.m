clear
filename='F:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_tide_origin.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');


nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','u_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','v_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','u_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','v_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
for ii=1:5
    disp(ii);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    for jj=1:5
        u_use=squeeze(ncread(filename,'u',[340*ii-339 260*jj-259 1],[340 260 inf]));
        v_use=squeeze(ncread(filename,'v',[340*ii-339 260*jj-259 1],[340 260 inf]));
        for mm=1:340
            parfor nn=1:260
            u_lp(mm,nn,:)=filt1('lp',squeeze(u_use(mm,nn,:)),'fs',1,'fc',1/48,'order',4);
            v_lp(mm,nn,:)=filt1('lp',squeeze(v_use(mm,nn,:)),'fs',1,'fc',1/48,'order',4);
            end
        end

        u_hp=u_use-u_lp;v_hp=v_use-v_lp;

        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','u_lp',single(u_lp),[340*ii-339 260*jj-259 1]);
        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','v_lp',single(v_lp),[340*ii-339 260*jj-259 1]);
        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','u_hp',single(u_hp),[340*ii-339 260*jj-259 1]);
        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_tide_filt_order4.nc','v_hp',single(v_hp),[340*ii-339 260*jj-259 1]);
       
    end
end

%%%%%%%%%
clear
filename='I:\NCS_wu_0.5km_15m_all\NCS_wu_0.5km_15m_notide_origin.nc';
x_rho=ncread(filename,'x_rho');y_rho=ncread(filename,'y_rho');


nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','x_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','y_rho','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2)},'Datatype','single');
ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','x_rho',single(x_rho),[1 1]);
ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','y_rho',single(y_rho),[1 1]);
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','u_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','v_lp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','u_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
nccreate('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','v_hp','Dimensions',{'x_rho',size(x_rho,1),'y_rho',size(x_rho,2),'t',695},'Datatype','single');
for ii=1:5
    disp(ii);
    disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    for jj=1:5
        u_use=squeeze(ncread(filename,'u',[340*ii-339 260*jj-259 1],[340 260 inf]));
        v_use=squeeze(ncread(filename,'v',[340*ii-339 260*jj-259 1],[340 260 inf]));
        for mm=1:340
            parfor nn=1:260
            u_lp(mm,nn,:)=filt1('lp',squeeze(u_use(mm,nn,:)),'fs',1,'fc',1/48,'order',4);
            v_lp(mm,nn,:)=filt1('lp',squeeze(v_use(mm,nn,:)),'fs',1,'fc',1/48,'order',4);
            end
        end

        u_hp=u_use-u_lp;v_hp=v_use-v_lp;

        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','u_lp',single(u_lp),[340*ii-339 260*jj-259 1]);
        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','v_lp',single(v_lp),[340*ii-339 260*jj-259 1]);
        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','u_hp',single(u_hp),[340*ii-339 260*jj-259 1]);
        ncwrite('I:\NCS_wu_0.5km_2day_filt\NCS_wu_0.5km_15m_notide_filt_order4.nc','v_hp',single(v_hp),[340*ii-339 260*jj-259 1]);
       
    end
end