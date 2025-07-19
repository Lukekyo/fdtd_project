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

def analyze_results(grid, monitor_name_list, P_incident=None):
    """
    根據 monitor_name_list 回傳每個 detector 的原始數據。
    """
    results = {}
    for name in monitor_name_list:
        if not hasattr(grid, name):
            print(f"⚠️ Grid 沒有名為 {name} 的 detector")
            continue
        detector = getattr(grid, name)
        values = detector.detector_values()
        pf = values.get("power_flow", None)
        result = {"power_flow": pf, "all_values": values}
        if P_incident is not None and pf is not None:
            ratio = abs(pf) / P_incident if detector.flip_sign else pf / P_incident
            result["ratio"] = ratio
        results[name] = result
    return results
