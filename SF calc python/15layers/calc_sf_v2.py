import numpy as np
from numba import njit, prange

@njit(parallel=True, fastmath=False)
def calc_binned_accum_SF_zeta_uv_parallel(u, v, zeta, X, Y, dr, max_r):
    # 1. 展平并强制转 float32
    u_vec = u.ravel().astype(np.float32)
    v_vec = v.ravel().astype(np.float32)
    zeta_vec = zeta.ravel().astype(np.float32)
    X_vec = X.ravel().astype(np.float32)
    Y_vec = Y.ravel().astype(np.float32)

    # 2. 保留 Python 原版本的有限值过滤优化：同时过滤 nan 和 inf
    valid_mask = (
        np.isfinite(u_vec) &
        np.isfinite(v_vec) &
        np.isfinite(zeta_vec) &
        np.isfinite(X_vec) &
        np.isfinite(Y_vec)
    )

    u_vec = u_vec[valid_mask]
    v_vec = v_vec[valid_mask]
    zeta_vec = zeta_vec[valid_mask]
    X_vec = X_vec[valid_mask]
    Y_vec = Y_vec[valid_mask]

    # 3. 按 X 排序，保留原有加速策略
    sort_idx = np.argsort(X_vec)
    X_vec = X_vec[sort_idx]
    Y_vec = Y_vec[sort_idx]
    u_vec = u_vec[sort_idx]
    v_vec = v_vec[sort_idx]
    zeta_vec = zeta_vec[sort_idx]

    N_points = len(u_vec)

    # 4. bin 参数
    dr = np.float32(dr)
    max_r = np.float32(max_r)
    num_bins = int(np.floor(max_r / dr))

    # 5. 全局累加器
    # 累加建议仍用 float32
    sum_DLL_uv  = np.zeros(num_bins, dtype=np.float32)
    sum_DTT_uv  = np.zeros(num_bins, dtype=np.float32)
    sum_D3_uv   = np.zeros(num_bins, dtype=np.float32)
    sum_D2_zeta = np.zeros(num_bins, dtype=np.float32)
    sum_D3_zeta = np.zeros(num_bins, dtype=np.float32)
    counts      = np.zeros(num_bins, dtype=np.int32)

    # 6. 并行主循环
    for i in prange(N_points):
        local_sum_DLL_uv  = np.zeros(num_bins, dtype=np.float32)
        local_sum_DTT_uv  = np.zeros(num_bins, dtype=np.float32)
        local_sum_D3_uv   = np.zeros(num_bins, dtype=np.float32)
        local_sum_D2_zeta = np.zeros(num_bins, dtype=np.float32)
        local_sum_D3_zeta = np.zeros(num_bins, dtype=np.float32)
        local_counts      = np.zeros(num_bins, dtype=np.int32)

        ui = u_vec[i]
        vi = v_vec[i]
        zetai = zeta_vec[i]
        xi = X_vec[i]
        yi = Y_vec[i]

        for j in range(i + 1, N_points):
            xj = X_vec[j]
            dx_val = xj - xi

            if dx_val > max_r:
                break

            yj = Y_vec[j]
            dy_val = yj - yi

            if abs(dy_val) > max_r:
                continue

            uj = u_vec[j]
            vj = v_vec[j]
            zetaj = zeta_vec[j]

            if not (np.isfinite(uj) and np.isfinite(vj) and np.isfinite(zetaj)):
                continue

            dist = np.sqrt(dx_val * dx_val + dy_val * dy_val)

            if (not np.isfinite(dist)) or dist > max_r or dist < 1e-10:
                continue

            bin_idx = int(np.floor(dist / dr))
            if bin_idx >= num_bins:
                continue

            cos_theta = dx_val / dist
            sin_theta = dy_val / dist

            du = uj - ui
            dv = vj - vi
            dzeta = zetaj - zetai

            du_L = du * cos_theta + dv * sin_theta
            du_T = -du * sin_theta + dv * cos_theta

            if not (np.isfinite(du_L) and np.isfinite(du_T) and np.isfinite(dzeta)):
                continue

            du_L2 = du_L * du_L
            du_T2 = du_T * du_T
            dzeta_2 = dzeta * dzeta

            local_sum_DLL_uv[bin_idx]  += du_L2
            local_sum_DTT_uv[bin_idx]  += du_T2
            local_sum_D3_uv[bin_idx]   += (du_L2 + du_T2) * du_L
            local_sum_D2_zeta[bin_idx] += dzeta_2
            local_sum_D3_zeta[bin_idx] += dzeta_2 * du_L
            local_counts[bin_idx]      += 1

        sum_DLL_uv  += local_sum_DLL_uv
        sum_DTT_uv  += local_sum_DTT_uv
        sum_D3_uv   += local_sum_D3_uv
        sum_D2_zeta += local_sum_D2_zeta
        sum_D3_zeta += local_sum_D3_zeta
        counts      += local_counts

    # 7. 输出 bin 中心
    # r_out = (np.arange(num_bins, dtype=np.float32) + 0.5) * dr

    return sum_DLL_uv, sum_DTT_uv, sum_D3_uv, sum_D2_zeta, sum_D3_zeta, counts