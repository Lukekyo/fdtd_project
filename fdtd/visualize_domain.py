import matplotlib.pyplot as plt
import numpy as np

# relatvie
from .grid import Grid
from .backend import backend as bd
from .waveforms import *
from .backend import backend as bd
from .detectors import CurrentDetector
from .objects import get_n_map, get_k_map, get_eps_map



def plot_simulation_domain(grid, plane='xz', mode='n'):
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

    # 擷取要畫的切面 (預設 XZ)
    if plane == 'xz':
        slice_data = data[:, grid.shape[1] // 2, :].T
        xlabel, ylabel = 'x', 'z'
    elif plane == 'xy':
        slice_data = data[:, :, grid.shape[2] // 2].T
        xlabel, ylabel = 'x', 'y'
    elif plane == 'yz':
        slice_data = data[grid.shape[0] // 2, :, :].T
        xlabel, ylabel = 'y', 'z'
    else:
        raise ValueError("Invalid plane")

    plt.figure(figsize=(5, 4))
    plt.imshow(slice_data, origin='lower', cmap=cmap, aspect='auto')
    plt.colorbar(label=label)
    plt.title(f"{label} Distribution ({plane.upper()} plane)")
    plt.xlabel(f"{xlabel} (grid index)")
    plt.ylabel(f"{ylabel} (grid index)")
    plt.tight_layout()
    plt.show()