function [S_low,S_high]=block_bootstrap_cospec(A_blocks,B);
%%% A_blocks 分块结构函数总和 M*K,M是分块数量，K是分离尺度维度
%%% N_blocks 对应的每个bins的个数 M*K,M是分块数量，K是分离尺度维度
%%% B 总共多少次bootstrap
    [M, K] = size(A_blocks);
       
    S_boot = nan(B, K);
    
    for b = 1:B
        % 随机抽 M 个 block，允许重复
        idx = randi(M, [M, 1]);
        % 对时间样本有放回抽样
        A_star(b,:) = mean(A_blocks(idx,:), 1, 'omitnan');
    end
    
    % 95% bootstrap interval
    S_low  = prctile(A_star, 2.5, 1);
    S_high = prctile(A_star, 97.5, 1);
    
end