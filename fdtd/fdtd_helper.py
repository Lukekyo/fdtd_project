# fdtd_helper.py
import numpy as np

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

def reference_run(grid_factory, total_steps=500):
    grid, _ = grid_factory(with_structure=False)

    for _ in range(total_steps):
        grid.step()

    T_ref = np.mean(grid.T.S[-100:])
    R_ref = np.mean(grid.R.S[-100:])
    return {"T": T_ref, "R": R_ref}