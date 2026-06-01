function [s_output] = interp_model_z(s,x_rho,y_rho,z_rho,z_need)
    %%%%%%%%%%%%%
    % s 为需要插值的变量 [lon, lat, layers]
    % x_rho,y_rho 为网格，二维 [lon lat]
    % z_rho为三维深度 [lon, lat, layers]，用z_levs生成 注意z_levs的维度与x_rho、y_rho一致
    % z_need 为需要获取的深度
    %%%%%%%%%%%%%
    % 将垂向S坐标插值到z坐标
    [n_lon, n_lat, ~] = size(s);
    n_z_new = length(z_need);
    
    % 预分配输出变量内存 (Parfor 必须步骤)
    % 输出维度: [lon, lat, new_layers]
    s_output = zeros(n_lon, n_lat, n_z_new); 
    
    % 并行循环
    parfor k = 1:n_lat
        % 预分配切片内存 (临时变量)
        % 维度: [lon, 1, new_layers] 以匹配 s_output(:, k, :)
        sSlice = nan(n_lon, 1, n_z_new); 
        for m = 1:n_lon
            % 提取当前点的深度和变量剖面
            z_profile = squeeze(z_rho(m,k,:));
            s_profile = squeeze(s(m,k,:));
            % 执行插值
            val = interp1(z_profile, s_profile, z_need, 'linear','extrap');
            sSlice(m, 1, :) = val; 
        end    
        % 将切片赋值回主变量
        s_output(:, k, :) = sSlice;
    end   
end