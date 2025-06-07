# fdtd_helper.py

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