# fdtd_helper.py

def um(val):
    """
    將微米數值轉換為公尺（m）
    Example:
        um(5) -> 5e-6
    """
    return val * 1e-6

def nm(val):
    """
    將奈米數值轉換為公尺（m）
    Example:
        nm(10) -> 10e-9
    """
    return val * 1e-9

def to_grid(val, grid_spacing):
    """
    將實際距離（以公尺為單位）轉換為 FDTD 模擬格點索引（整數）
    Example:
        to_grid(5e-6, 10e-9) -> 500
    """
    return int(val / grid_spacing + 0.5)

def from_grid(index, grid_spacing):
    """
    將格點索引轉換為實際距離（公尺）
    Example:
        from_grid(500, 10e-9) -> 5e-6
    """
    return index * grid_spacing
