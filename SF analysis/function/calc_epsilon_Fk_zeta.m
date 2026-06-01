function [k,R,F_k,g_j,S3_recon] = calc_epsilon_Fk_zeta(r,S3,r_min_data,r_max_data,r_min_target,r_max_target,smooth_win,lambda_reg,Nk)

    %%%%%%%%%%%%%%%
    % 输入
    % r 原始数据分离尺度
    % S3 这里传入拟能对应的 mixed structure function:
    % Vw(r) = <delta u_L * (delta omega)^2>
    % r_min_data,r_max_data 数据参与拟合的范围（尽量用全部数据）
    % r_min_target,r_max_target 拟合的波数范围（以防边界效应截取中间段）
    % smooth_win 平滑窗口，几个点平滑
    % lambda_reg Tikhonov 正则化求解参数
    % Nk r_min_target,r_max_target之间分为多少个点
    %
    % 输出
    % k 拟合的波数
    % R 拟合使用的分离尺度
    % F_k 拟能谱通量 Fw(k)
    % g_j 每个波数处的局地拟能注入率密度
    % S3_recon 重构的 Vw(r)

    S3 = S3(:);
    r  = r(:);

    % 用 r_min_data ~ r_max_data 筛选参与拟合的数据点
    idx = r >= r_min_data & r <= r_max_data;
    R      = r(idx);
    S3_raw = S3(idx);
    Nr     = length(R);

    S3_smooth = smoothdata(S3_raw, 'movmean', smooth_win);

    fprintf('拟合数据范围: [%.1f, %.1f] km | %d 个r点\n', r_min_data, r_max_data, Nr);
    fprintf('波数反演范围: k ∈ [%.4f, %.4f] 1/km  ↔  λ ∈ [%.1f, %.1f] km | %d 个k点\n', ...
        1/r_max_target, 1/r_min_target, r_min_target, r_max_target, Nk);

    % ===================== 3. 构建波数轴（由目标尺度决定）=====================
    k_min = 1 / r_max_target;
    k_max = 1 / r_min_target;

    k_edges = linspace(k_min, k_max, Nk + 1);
    dk_val  = k_edges(2) - k_edges(1);
    k       = 0.5 * (k_edges(1:end-1) + k_edges(2:end));

    % ===================== 4. 构建反演矩阵 A =====================
    % 对应推导:
    % Vw(r) = - sum_j [ 4*g_j/k_j * J1(k_j r) * dk ]
    A = zeros(Nr, Nk);
    for j = 1:Nk
        A(:, j) = -4 / k(j) * besselj(1, k(j) * R) * dk_val;
    end

    % 1/R 归一化（保留原脚本写法，尽量少改）
    b = S3_smooth ./ abs(R);
    for i = 1:Nr
        A(i, :) = A(i, :) / abs(R(i));
    end

    % ===================== 5. Tikhonov 正则化求解 =====================
    L = eye(Nk);
    g_j = (A' * A + lambda_reg * L) \ (A' * b);

    % ===================== 6. 累积积分 → Fw(k) =====================
    % Fw(k) = sum_{ki <= k} g_i * dk
    F_k = cumsum(g_j) * dk_val;

    % ===================== 7. 重构 Vw(r) =====================
    S3_recon = zeros(size(R));
    for j = 1:Nk
        S3_recon = S3_recon - 4 * g_j(j) / k(j) * besselj(1, k(j) * R) * dk_val;
    end
end
