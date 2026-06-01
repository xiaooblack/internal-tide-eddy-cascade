#!/usr/bin/env python3
"""
结构函数计算 - 超算版本（多层全域版，时间分片）
用法: python main.py <t_start> <t_end>
例如: python main.py 0 174      → 计算时间步 0~173
      python main.py 174 348    → 计算时间步 174~347
      python main.py 348 522    → 计算时间步 348~521
      python main.py 522 695    → 计算时间步 522~694
注意: t_end 是开区间，不包含在内（Python惯例）
"""

import sys
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import os

from calc_sf_v2 import calc_binned_accum_SF_zeta_uv_parallel
from model_gradient import model_gradient


def main():
    # ============ 命令行参数：手动输入起止 ============
    t_start = int(sys.argv[1])
    t_end   = int(sys.argv[2])
    num_times_local = t_end - t_start
    print(f"时间范围: [{t_start}, {t_end})，共 {num_times_local} 个时间步")

    # ============ 参数设置 ============
    input_file = '/XYFS01/HDD_POOL/ac_scsio_qiyiquan/ac_scsio_qiyiquanxy_3/wuby/python_project/calc_SF_15layers_2dayfilt_notide_part1/NCS_wu_0.5km_15layers_notide_filt_order4.nc'
    output_file = f'/XYFS01/HDD_POOL/ac_scsio_qiyiquan/ac_scsio_qiyiquanxy_3/wuby/python_project/calc_SF_15layers_2dayfilt_notide_part1/SF_NCS_wu_0.5km_15layers_notide_filt_order4_t{t_start}_{t_end}.nc'

    dr = 1.0      # km
    max_r = 100.0 # km
    r_edges = np.arange(0, max_r + dr, dr).astype(np.float32)
    num_bins = len(r_edges) - 1

    # ============ 打开输入文件 ============
    print(f"打开输入文件: {input_file}")
    nc_in = Dataset(input_file, 'r')

    u_lp_shape = nc_in.variables['u_lp'].shape
    num_layers = u_lp_shape[1]
    print(f"u_lp 形状: {u_lp_shape}")
    print(f"层数: {num_layers}")

    # 全域坐标
    # 梯度计算使用米
    x_m = np.asarray(nc_in.variables['x_rho'][:, :], dtype=np.float64)
    y_m = np.asarray(nc_in.variables['y_rho'][:, :], dtype=np.float64)

    # 结构函数距离计算使用km
    x_km = (x_m / 1000.0).astype(np.float32)
    y_km = (y_m / 1000.0).astype(np.float32)

    print(f"坐标形状: {x_km.shape}")
    print(f"x范围: [{np.nanmin(x_km):.1f}, {np.nanmax(x_km):.1f}] km, NaN数: {np.isnan(x_km).sum()}")
    print(f"y范围: [{np.nanmin(y_km):.1f}, {np.nanmax(y_km):.1f}] km, NaN数: {np.isnan(y_km).sum()}")
    print(f"x跨度: {np.nanmax(x_km)-np.nanmin(x_km):.1f} km")
    print(f"y跨度: {np.nanmax(y_km)-np.nanmin(y_km):.1f} km")

    # ============ 创建输出文件 ============
    print(f"创建输出文件: {output_file}")
    if os.path.exists(output_file):
        os.remove(output_file)

    nc_out = Dataset(output_file, 'w', format='NETCDF4')
    nc_out.createDimension('r', num_bins)
    nc_out.createDimension('depth', num_layers)
    nc_out.createDimension('t', num_times_local)

    var_r = nc_out.createVariable('r', 'f4', ('r',))
    var_r[:] = (r_edges[:-1] + r_edges[1:]) / 2
    var_r.units = 'km'

    var_t = nc_out.createVariable('time_index', 'i4', ('t',))
    var_t[:] = np.arange(t_start, t_end)
    var_t.long_name = 'global time index in original file'

    # ---------- LP ----------
    var_DLL_lp = nc_out.createVariable('sum_DLL_lp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_DTT_lp = nc_out.createVariable('sum_DTT_lp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_D3_lp  = nc_out.createVariable('sum_D3_lp',  'f4', ('r', 'depth', 't'), zlib=True)
    var_D2z_lp = nc_out.createVariable('sum_D2_zeta_lp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_D3z_lp = nc_out.createVariable('sum_D3_zeta_lp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_cnt_lp = nc_out.createVariable('counts_lp',  'i4', ('r', 'depth', 't'), zlib=True)

    # ---------- HP ----------
    var_DLL_hp = nc_out.createVariable('sum_DLL_hp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_DTT_hp = nc_out.createVariable('sum_DTT_hp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_D3_hp  = nc_out.createVariable('sum_D3_hp',  'f4', ('r', 'depth', 't'), zlib=True)
    var_D2z_hp = nc_out.createVariable('sum_D2_zeta_hp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_D3z_hp = nc_out.createVariable('sum_D3_zeta_hp', 'f4', ('r', 'depth', 't'), zlib=True)
    var_cnt_hp = nc_out.createVariable('counts_hp',  'i4', ('r', 'depth', 't'), zlib=True)

    # ============ 预热 Numba ============
    print("\n预热 Numba...")
    dummy = np.random.rand(30, 30).astype(np.float32)
    dummy_zeta = np.random.rand(30, 30).astype(np.float32)
    _ = calc_binned_accum_SF_zeta_uv_parallel(dummy, dummy, dummy_zeta, dummy, dummy, dr, max_r)
    print("JIT 完成！")

    # ============ 主循环 ============
    print(f"\n开始处理: {num_times_local} 时间步 × {num_layers} 层")
    total_start = datetime.now()

    for ii_local in range(num_times_local):
        ii_global = t_start + ii_local

        print(f"\n{'='*60}")
        print(f"本地 {ii_local+1}/{num_times_local} | 全局时间步 {ii_global}")
        print(f"开始: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        for kk in range(num_layers):

            # ===================== LP =====================
            u_lp_raw = nc_in.variables['u_lp'][ii_global, kk, :, :]
            v_lp_raw = nc_in.variables['v_lp'][ii_global, kk, :, :]

            u_lp = np.ma.filled(u_lp_raw, np.nan).astype(np.float64)
            v_lp = np.ma.filled(v_lp_raw, np.nan).astype(np.float64)

            # 计算梯度和涡度（坐标用米）
            dudx_lp, dudy_lp = model_gradient(x_m, y_m, u_lp)
            dvdx_lp, dvdy_lp = model_gradient(x_m, y_m, v_lp)
            zeta_lp = dvdx_lp - dudy_lp
            zeta_lp[np.abs(zeta_lp) > 100] = np.nan

            # 结构函数计算（传入 float32 更匹配 numba 内部）
            s_DLL_lp, s_DTT_lp, s_D3_lp, s_D2z_lp, s_D3z_lp, c_lp = \
                calc_binned_accum_SF_zeta_uv_parallel(
                    u_lp.astype(np.float32),
                    v_lp.astype(np.float32),
                    zeta_lp.astype(np.float32),
                    x_km,
                    y_km,
                    dr,
                    max_r
                )

            var_DLL_lp[:, kk, ii_local] = s_DLL_lp
            var_DTT_lp[:, kk, ii_local] = s_DTT_lp
            var_D3_lp[:, kk, ii_local]  = s_D3_lp
            var_D2z_lp[:, kk, ii_local] = s_D2z_lp
            var_D3z_lp[:, kk, ii_local] = s_D3z_lp
            var_cnt_lp[:, kk, ii_local] = c_lp

            # ===================== HP =====================
            u_hp_raw = nc_in.variables['u_hp'][ii_global, kk, :, :]
            v_hp_raw = nc_in.variables['v_hp'][ii_global, kk, :, :]

            u_hp = np.ma.filled(u_hp_raw, np.nan).astype(np.float64)
            v_hp = np.ma.filled(v_hp_raw, np.nan).astype(np.float64)

            dudx_hp, dudy_hp = model_gradient(x_m, y_m, u_hp)
            dvdx_hp, dvdy_hp = model_gradient(x_m, y_m, v_hp)
            zeta_hp = dvdx_hp - dudy_hp
            zeta_hp[np.abs(zeta_hp) > 100] = np.nan

            s_DLL_hp, s_DTT_hp, s_D3_hp, s_D2z_hp, s_D3z_hp, c_hp = \
                calc_binned_accum_SF_zeta_uv_parallel(
                    u_hp.astype(np.float32),
                    v_hp.astype(np.float32),
                    zeta_hp.astype(np.float32),
                    x_km,
                    y_km,
                    dr,
                    max_r
                )

            var_DLL_hp[:, kk, ii_local] = s_DLL_hp
            var_DTT_hp[:, kk, ii_local] = s_DTT_hp
            var_D3_hp[:, kk, ii_local]  = s_D3_hp
            var_D2z_hp[:, kk, ii_local] = s_D2z_hp
            var_D3z_hp[:, kk, ii_local] = s_D3z_hp
            var_cnt_hp[:, kk, ii_local] = c_hp

            # ===================== 诊断输出 =====================
            print(f"  层 {kk+1}/{num_layers} LP诊断:")
            print(f"    u_lp NaN数: {np.isnan(u_lp).sum()}/{u_lp.size}")
            print(f"    zeta_lp NaN数: {np.isnan(zeta_lp).sum()}/{zeta_lp.size}")
            print(f"    有效点对总数: {c_lp.sum():.2e}")
            print(f"    DLL 前5个bin: {s_DLL_lp[:5]}")
            print(f"    D2_zeta 前5个bin: {s_D2z_lp[:5]}")
            print(f"    D3_zeta 前5个bin: {s_D3z_lp[:5]}")
            print(f"    counts 前5个bin: {c_lp[:5]}")

            print(f"  层 {kk+1}/{num_layers} 完成 | lp点对: {c_lp.sum():.2e}, hp点对: {c_hp.sum():.2e}")

            nc_out.sync()

    nc_in.close()
    nc_out.close()

    total_elapsed = (datetime.now() - total_start).total_seconds()
    print(f"\n{'='*60}")
    print(f"完成！总耗时: {total_elapsed/3600:.2f} 小时")
    print(f"结果保存在: {output_file}")


if __name__ == '__main__':
    main()