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

def get_data(grid, monitor_name, wavelength=1550e-9, prop_direction='z'):
    """
    簡單獲取Monitor在特定波長的功率
    
    Args:
        grid: FDTD網格
        monitor_name: Monitor名稱
        wavelength: 波長(nm)
        prop_direction: 傳播方向 ('x', 'y', 'z')
    
    Returns:
        float: 該波長的功率
    """
        # 方向索引字典
    direction_map = {
        'x': 0,
        'y': 1, 
        'z': 2
    }
    
    if prop_direction not in direction_map:
        raise ValueError(f"prop_direction 必須是 'x', 'y', 或 'z'，收到: {prop_direction}")
    
    dir_idx = direction_map[prop_direction]

    # 獲取detector數據
    monitor = getattr(grid, monitor_name)
    dataE = monitor.detector_values()['E']
    dataH = monitor.detector_values()['H']

    print(f"🔄 處理時空數據:")
    print(f"   E形狀: {dataE.shape} = (時間步, 空間點, 場分量)")
    print(f"   H形狀: {dataH.shape}")
    print(f"   轉換波長: {wavelength * 1e9:.2f} nm")
    
    # ===== 步驟1: 空間積分，時間序列 =====
    power_time_series = []

    for t in range(dataE.shape[0]):  # 遍歷時間步
        # 獲取該時間步的空間分布場
        E_t = dataE[t, :, :]  # 形狀: (100空間點, 3分量)
        H_t = dataH[t, :, :]  # 形狀: (100空間點, 3分量)
        print(f"   時間步 {t+1}/{dataE.shape[0]}: E形狀 {E_t.shape}, H形狀 {H_t.shape}")
        
        # 計算每個空間點的坡印廷矢量
        S_vec = bd.cross(E_t, bd.conj(H_t))  # 形狀: (100空間點, 3分量)
        
        # 提取傳播方向分量（假設z方向）
        S_dir = bd.real(S_vec[:, dir_idx])  # 形狀: (100空間點,)
        
        # 空間積分：對所有空間點求和
        total_power_t = bd.sum(S_dir)  # 標量
        
        power_time_series.append(total_power_t)
    
    power_array = bd.array(power_time_series)  # 形狀: (1000,) 純時間序列
    print(f"   時間序列形狀: {power_array.shape}")
    print(f"   空間積分後功率時間序列: {power_array.shape}")
    print(f"   功率範圍: {bd.min(power_array):.2e} ~ {bd.max(power_array):.2e}")
    
    # ===== 步驟2: 時間域FFT =====
    fft_result = bd.fft(power_array)  # FFT時間序列
    freqs = bd.fftfreq(len(power_array), grid.time_step)
    
    # ===== 步驟3: 找目標頻率 =====
    target_freq = bd.c0 / wavelength
    freq_idx = bd.argmin(bd.abs(freqs - target_freq))
    
    power_at_freq = bd.abs(fft_result[freq_idx])**2
    
    print(f"   目標頻率: {target_freq:.2e} Hz")
    print(f"   實際頻率: {freqs[freq_idx]:.2e} Hz")
    print(f"   該頻率功率: {power_at_freq:.2e}")
    
    return power_at_freq, power_array, fft_result, freqs


def get_result(grid, monitor_name, wavelength=1550e-9, prop_direction='z'):

    power_at_freq, power_array, fft_result, freqs = get_data(grid, monitor_name, wavelength, prop_direction)
    
    print(f"   輸入波長: {wavelength * 1e9:.2f} nm")
    # 計算源功率（如果有的話）
    P_source = None
    if hasattr(grid, 'source') and hasattr(grid.source, 'monitor_data'):
        source_data = bd.array(grid.source.monitor_data['E_field'])
        if len(source_data) > 0:
            source_fft = bd.fft(source_data)
            source_freqs = bd.fftfreq(len(source_data), grid.time_step)
            
            target_freq = bd.c0 / wavelength
            source_idx = bd.argmin(bd.abs(source_freqs - target_freq))
            P_source = bd.abs(source_fft[source_idx])**2
            
            print(f"   源功率: {P_source:.2e}")
    
    print(f"   穿透功率: {power_at_freq}")
    
    # 計算T/R比值
    if P_source is not None and P_source > 0:
        T = power_at_freq / P_source
        
        print(f"\n📈 結果:")
        print(f"   穿透率 T = {T:.4f} ({T*100:.1f}%)")
        
        return T, power_array

def get_source_power(grid, wavelength=1550e-9):
    """考慮源空間尺寸的功率計算"""
    
    source_E = bd.array(grid.source.monitor_data['E_field'])  # 中心點E
    source_H = bd.array(grid.source.monitor_data['H_field'])  # 中心點H
    print(f"🔄 處理源數據:")
    
    # 計算中心點的Poynting功率時間序列
    power_time_series = []
    for t in range(len(source_E)):
        E_t = source_E[t]  # 3分量向量
        H_t = source_H[t]  # 3分量向量
        # S_z = E_x * H_y* - E_y * H_x*
        S_z = bd.real(E_t * bd.conj(H_t[1]))  # Ex * Hy
        power_time_series.append(S_z)
    
    # 🔧 關鍵：乘以源的空間面積
    source_area = len(grid.source.x) * grid.grid_spacing  # 在2D中這是"有效長度"
    power_array = bd.array(power_time_series) * source_area
    
    print(f"   Source 有效尺寸: {len(grid.source.x)} × {len(grid.source.z)} 格點")
    print(f"   源面積: {source_area:.2f} μm²")
    
    # FFT（與檢測器相同）
    fft_result = bd.fft(power_array)
    freqs = bd.fftfreq(len(power_array), grid.time_step)
    
    target_freq = bd.c0 / wavelength
    freq_idx = bd.argmin(bd.abs(freqs - target_freq))
    
    if hasattr(freq_idx, 'item'):
        freq_idx = freq_idx.item()
    
    power_at_freq = bd.abs(fft_result[freq_idx])**2
    
    if hasattr(power_at_freq, 'item'):
        power_at_freq = power_at_freq.item()
    
    print(f"   源Poynting功率: {power_at_freq:.2e}")
    return power_at_freq