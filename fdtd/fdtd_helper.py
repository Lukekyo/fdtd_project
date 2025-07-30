# fdtd_helper.py - 完整整理版本
import numpy as np
from .backend import backend as bd
import matplotlib.pyplot as plt

def um(x: float) -> float:
    """Convert micrometer to meters"""
    return x * 1e-6

def nm(x: float) -> float:
    """Convert nanometer to meters"""
    return x * 1e-9

def to_grid(length: float, grid_spacing: float) -> int:
    """Convert physical length to grid index"""
    return int(length / grid_spacing + 0.5)

def from_grid(index: int, grid_spacing: float) -> float:
    """Convert grid index to physical length"""
    return index * grid_spacing
# FDTD週期邊界條件診斷工具

# 在 fdtd_helper.py 末尾添加這些函數
def get_monitor_power_at_wavelength(grid, monitor_name, wavelength=1550):
    """
    簡單獲取Monitor在特定波長的功率
    
    Args:
        grid: FDTD網格
        monitor_name: Monitor名稱
        wavelength_nm: 波長(nm)
    
    Returns:
        float: 該波長的功率
    """
    
    # 獲取數據
    monitor = getattr(grid, monitor_name)
    if hasattr(monitor, 'S'):
        data = bd.array(monitor.S)
    elif hasattr(monitor, 'monitor_data'):
        data = bd.array(monitor.monitor_data['E_field'])
    else:
        return 0
    print(data)
    # FFT
    fft_result = bd.fft(data)
    freqs = bd.fftfreq(len(data), grid.time_step)
    
    # 找目標頻率
    target_freq = bd.c0 / (wavelength * 1e-9)
    freq_idx = bd.argmin(bd.abs(freqs - target_freq))
    
    # 返回功率
    return bd.abs(fft_result[freq_idx])**2


def calculate_T(grid, wavelength=1550):
    """
    超簡單T/R計算
    
    Returns:
        tuple: (T, R)
    """
    
    # 獲取各Monitor的功率
    P_source = get_monitor_power_at_wavelength(grid, 'source', wavelength)
    P_T = get_monitor_power_at_wavelength(grid, 'T', wavelength) 
    # P_R = get_monitor_power_at_wavelength(grid, 'R', wavelength)
    
    # 計算T/R
    T = P_T / P_source if P_source > 0 else 0
    # R = P_R / P_source if P_source > 0 else 0
    
    print(f"λ={wavelength}nm: T={T:.3f}")
    
    return T