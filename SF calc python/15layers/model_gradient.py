import numpy as np


def model_gradient(x, y, s):
    """
    计算不均匀二维水平网格上的梯度
    适用于 s 的最后两维是 (y, x) 的情况

    参数
    ----------
    x : ndarray, shape (ny, nx)
        x 坐标，单位 m
    y : ndarray, shape (ny, nx)
        y 坐标，单位 m
    s : ndarray
        待求梯度变量，最后两维必须是 (ny, nx)
        例如:
        - (ny, nx)
        - (nz, ny, nx)
        - (nt, nz, ny, nx)

    返回
    -------
    dsdx : ndarray
        x方向梯度，shape 与 s 相同
    dsdy : ndarray
        y方向梯度，shape 与 s 相同
    """

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    s = np.asarray(s, dtype=np.float64)

    dsdx = np.zeros_like(s, dtype=np.float64)
    dsdy = np.zeros_like(s, dtype=np.float64)

    # ---------- y方向梯度: 沿倒数第2维 ----------
    # 中心差分
    dsdy[..., 1:-1, :] = (
        s[..., 2:, :] - s[..., :-2, :]
    ) / (y[2:, :] - y[:-2, :])

    # 边界单侧差分
    dsdy[..., 0, :] = (
        s[..., 1, :] - s[..., 0, :]
    ) / (y[1, :] - y[0, :])

    dsdy[..., -1, :] = (
        s[..., -1, :] - s[..., -2, :]
    ) / (y[-1, :] - y[-2, :])

    # ---------- x方向梯度: 沿倒数第1维 ----------
    # 中心差分
    dsdx[..., :, 1:-1] = (
        s[..., :, 2:] - s[..., :, :-2]
    ) / (x[:, 2:] - x[:, :-2])

    # 边界单侧差分
    dsdx[..., :, 0] = (
        s[..., :, 1] - s[..., :, 0]
    ) / (x[:, 1] - x[:, 0])

    dsdx[..., :, -1] = (
        s[..., :, -1] - s[..., :, -2]
    ) / (x[:, -1] - x[:, -2])

    return dsdx, dsdy
