function [S_low,S_high]=block_bootstrap(A_blocks,N_blocks,B);
%%% A_blocks 分块结构函数总和 M*K,M是分块数量，K是分离尺度维度
%%% N_blocks 对应的每个bins的个数 M*K,M是分块数量，K是分离尺度维度
%%% B 总共多少次bootstrap
    [M, K] = size(A_blocks);
    
    % 原始结构函数
    S = nansum(A_blocks, 1) ./ nansum(N_blocks, 1);
    
    S_boot = nan(B, K);
    
    for b = 1:B
    
        % 随机抽 M 个 block，允许重复
        idx = randi(M, [M, 1]);
    
        % 合并抽中的 block
        A_star = nansum(A_blocks(idx, :), 1);
        N_star = nansum(N_blocks(idx, :), 1);
    
        % 重新算结构函数
        S_boot(b, :) = A_star ./ N_star;
    
    end
    
    % 95% bootstrap interval
    S_low  = prctile(S_boot, 2.5, 1);
    S_high = prctile(S_boot, 97.5, 1);
    
end