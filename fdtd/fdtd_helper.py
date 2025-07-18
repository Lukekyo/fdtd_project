# fdtd_helper.py - 完整整理版本
import numpy as np
from .backend import backend as bd

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

import numpy as np
import matplotlib.pyplot as plt

# === 簡化版診斷函數（直接可用） ===

def quick_diagnosis(grid, grid_spacing):
    """
    快速診斷您現有的FDTD設定
    """
    print("🚨 快速FDTD診斷")
    print("="*50)


    # 1. 基本資訊
    print(f"1. 基本配置:")
    print(f"   網格形狀: {grid.shape}")
    print(f"   Grid spacing: {grid_spacing*1e9:.1f} nm")
    print(f"   Time step: {grid.time_step:.2e} s")

    # 2. 場類型檢查
    E_complex = 'complex' in str(grid.E.dtype)
    H_complex = 'complex' in str(grid.H.dtype)
    print(f"\n2. 場數據類型:")
    print(f"   E場複數: {E_complex}")
    print(f"   H場複數: {H_complex}")

    # 3. 邊界條件
    print(f"\n3. 邊界條件:")
    for i, boundary in enumerate(grid.boundaries):
        btype = type(boundary).__name__
        print(f"   邊界{i}: {btype}")
        if 'Bloch' in btype:
            print(f"           k_component: {getattr(boundary, 'k_component', 'N/A')}")
            print(f"           length: {getattr(boundary, 'length', 'N/A')}")

    # 4. 源檢查
    print(f"\n4. 源配置:")
    for i, source in enumerate(grid.sources):
        stype = type(source).__name__
        print(f"   源{i}: {stype}")
        if hasattr(source, 'theta_deg'):
            print(f"         角度: {source.theta_deg}°")
        if hasattr(source, 'amplitude'):
            print(f"         振幅: {source.amplitude}")

    # 5. 檢測器
    print(f"\n5. 檢測器:")
    for i, detector in enumerate(grid.detectors):
        dtype = type(detector).__name__
        print(f"   檢測器{i}: {dtype}")
        if hasattr(detector, 'flip_sign'):
            print(f"            反向: {detector.flip_sign}")

    # 6. 場數值檢查
    E_max = float(np.max(np.abs(grid.E)))
    H_max = float(np.max(np.abs(grid.H)))
    print(f"\n6. 場強度:")
    print(f"   |E|_max: {E_max:.3e}")
    print(f"   |H|_max: {H_max:.3e}")

    if H_max > 1e-15:
        impedance = E_max / H_max
        print(f"   實際阻抗: {impedance:.1f} Ω")
        print(f"   理論阻抗: 377 Ω")
        if abs(impedance - 377) > 100:
            print("   ❌ 阻抗嚴重偏離！")
        else:
            print("   ✅ 阻抗合理")
    else:
        print("   ❌ H場為零！")

    # 7. 主要問題識別
    print(f"\n🎯 主要問題:")
    issues = []

    if not (E_complex and H_complex):
        issues.append("❌ 場不是複數類型")

    if not any('Bloch' in type(b).__name__ for b in grid.boundaries):
        issues.append("❌ 缺少BlochBoundary")

    if H_max < 1e-15:
        issues.append("❌ 磁場為零")
    elif abs(E_max/H_max - 377) > 100:
        issues.append("❌ 阻抗異常")

    if len(issues) == 0:
        print("   ✅ 基本配置正確")
    else:
        for issue in issues:
            print(f"   {issue}")

    return issues

