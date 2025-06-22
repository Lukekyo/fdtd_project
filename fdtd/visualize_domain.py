import matplotlib.pyplot as plt
import numpy as np

# relatvie
from .grid import Grid
from .backend import backend as bd
from .waveforms import *
from .backend import backend as bd
from .detectors import CurrentDetector
from .objects import get_n_map, get_k_map, get_eps_map



def plot_simulation_domain(grid, 
                           plane='xz', 
                           mode='n', 
                           vmin=None, 
                           vmax=None, 
                           show_grid=True,
                           x_range = None,
                           z_range = None
                        ):
    """在幾何建立後繪製材料分布圖"""
    if mode == 'n':
        data = get_n_map(grid)
        label = 'n'
        cmap = 'viridis'
    elif mode == 'k':
        data = get_k_map(grid)
        label = 'k'
        cmap = 'magma'
    elif mode == 'eps':
        data = bd.real(get_eps_map(grid))
        label = 'Re(ε)'
        cmap = 'plasma'
    else:
        raise ValueError("mode must be 'n', 'k' or 'eps'")

    dx = grid.grid_spacing
    Nx, Ny, Nz = grid.shape
    
    # 2. 裁切 index 範圍（μm → grid index）
    def um_to_idx(um_value): return int(um_value * 1e-6 / dx)

    if x_range:
        x_start_idx, x_end_idx = um_to_idx(x_range[0]), um_to_idx(x_range[1])
    else:
        x_start_idx, x_end_idx = 0, Nx

    if z_range:
        z_start_idx, z_end_idx = um_to_idx(z_range[0]), um_to_idx(z_range[1])
    else:
        z_start_idx, z_end_idx = 0, Nz
    
    # 3. 擷取資料範圍
    slice_data = data[x_start_idx:x_end_idx, Ny // 2, z_start_idx:z_end_idx].T

    # 4. 計算實體位置
    extent = [
        x_start_idx * dx * 1e6,
        x_end_idx * dx * 1e6,
        z_start_idx * dx * 1e6,
        z_end_idx * dx * 1e6,
    ]
    xlabel, ylabel = "x (μm)", "z (μm)"
    
    # 5. 繪圖
    plt.figure(figsize=(5, 4))
    im = plt.imshow(
        slice_data,
        extent=extent,
        origin='lower',
        cmap=cmap,
        aspect='auto',
        vmin=vmin,
        vmax=vmax
    )
    plt.colorbar(im, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{label} Distribution ({plane.upper()} plane)")

    # 6. 疊加 grid（如果需要）
    if show_grid:
        x_ticks = np.arange(x_start_idx, x_end_idx + 1) * dx * 1e6
        z_ticks = np.arange(z_start_idx, z_end_idx + 1) * dx * 1e6
        for xi in x_ticks:
            plt.axvline(x=xi, color='white', linewidth=0.3, linestyle='--', alpha=0.3)
        for zi in z_ticks:
            plt.axhline(y=zi, color='white', linewidth=0.3, linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.show()