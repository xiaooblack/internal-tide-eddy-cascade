function [r_out, sum_DLL, sum_DTT, sum_D3, counts] = calc_binned_accum_SF_parallel_single(u, v, X, Y, dr, max_r)
% CALC_BINNED_ACCUM_SF_PARALLEL_SINGLE 
% 强制单精度计算版本，优化内存带宽和计算速度。

    %% 1. 数据预处理与排序 (强制转为 single)
    % 显式转换为 single，确保后续算术运算在单精度下进行
    u_vec = single(u(:));
    v_vec = single(v(:));
    X_vec = single(X(:));
    Y_vec = single(Y(:));
    
    % 排序 (排序索引操作不受精度影响)
    [X_vec, sort_idx] = sort(X_vec);
    Y_vec = Y_vec(sort_idx);
    u_vec = u_vec(sort_idx);
    v_vec = v_vec(sort_idx);
    
    N_points = length(u_vec); 
    
    if isempty(gcp('nocreate'))
        parpool; 
    end

    %% 2. 初始化 Bins (关键修改：声明 'single')
    % 确保 bins 的边缘也是 single，防止比较时发生类型提升
    r_edges = single(0 : dr : max_r); 
    num_bins = length(r_edges) - 1;
    max_r = single(max_r); % 确保判断阈值也是 single
    dr = single(dr);
    
    % 初始化全局累加变量为 single
    % 如果不加 'single'，MATLAB 默认创建 double，会导致计算回退
    sum_DLL = zeros(num_bins, 1, 'single'); 
    sum_DTT = zeros(num_bins, 1, 'single'); 
    sum_D3  = zeros(num_bins, 1, 'single'); 
    counts  = zeros(num_bins, 1, 'single'); 
    
    %% 3. 并行循环遍历
    parfor i = 1 : N_points
        
        % 关键修改：局部累加器必须初始化为 single
        local_sum_DLL = zeros(num_bins, 1, 'single');
        local_sum_DTT = zeros(num_bins, 1, 'single');
        local_sum_D3  = zeros(num_bins, 1, 'single');
        local_counts  = zeros(num_bins, 1, 'single');
        
        ui = u_vec(i); vi = v_vec(i);
        xi = X_vec(i); yi = Y_vec(i);
        
        if ~isnan(ui) && ~isnan(vi)
            
            for j = (i + 1) : N_points
                
                xj = X_vec(j);
                
                % 距离判断
                dx_val = xj - xi;
                if dx_val > max_r
                    break; 
                end
                
                yj = Y_vec(j);
                dy_val = yj - yi; 
                
                % 这里的 abs() 对 single 输入会返回 single
                if abs(dy_val) > max_r
                    continue; 
                end

                uj = u_vec(j); vj = v_vec(j);
                if isnan(uj) || isnan(vj), continue; end
                
                % --- 计算实际距离 ---
                % sqrt 对 single 输入会返回 single，速度更快
                dist = sqrt(dx_val^2 + dy_val^2);
                
                if dist > max_r || dist < 1e-10, continue; end
                
                % --- B. 确定 Bin 索引 ---
                bin_idx = floor(dist / dr) + 1;
                
                if bin_idx > num_bins, continue; end
                
                % --- C. 计算投影单位向量 ---
                % 所有变量此刻都是 single，除法也是单精度除法
                cos_theta = dx_val / dist;
                sin_theta = dy_val / dist;
                
                % --- D. 计算速度差 ---
                du = uj - ui;
                dv = vj - vi;
                
                % --- E. 投影 ---
                du_L = du * cos_theta + dv * sin_theta;
                du_T = -du * sin_theta + dv * cos_theta;
                
                % --- F. 累加 ---
                du_L2 = du_L^2;
                du_T2 = du_T^2;
                
                % 因为 local_sum_DLL 是 single 类型，这里不会发生类型转换
                local_sum_DLL(bin_idx) = local_sum_DLL(bin_idx) + du_L2;
                local_sum_DTT(bin_idx) = local_sum_DTT(bin_idx) + du_T2;
                local_sum_D3(bin_idx)  = local_sum_D3(bin_idx)  + (du_L2 + du_T2) * du_L;
                local_counts(bin_idx)  = local_counts(bin_idx)  + 1;
            end
        end
        
        % --- G. 归约 ---
        sum_DLL = sum_DLL + local_sum_DLL;
        sum_DTT = sum_DTT + local_sum_DTT;
        sum_D3  = sum_D3  + local_sum_D3;
        counts  = counts  + local_counts;
    end
    
    r_out = (r_edges(1:end-1) + r_edges(2:end))' / 2;
end
