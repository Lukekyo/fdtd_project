import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import time
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
from scipy.ndimage import gaussian_filter1d
import fdtd

# ==== 模擬參數 ====
fdtd.set_backend("numpy")
wavelength = fdtd.nm(1550)
grid_spacing = fdtd.nm(20)
x_span, z_span = fdtd.um(2), fdtd.um(6)
Nx = fdtd.to_grid(x_span, grid_spacing)
Nz = fdtd.to_grid(z_span, grid_spacing)
theta_deg = 0
theta = np.deg2rad(theta_deg)
k0 = 2 * np.pi / wavelength
kx = k0 * np.sin(theta)
Lx = Nx * grid_spacing
source_z = fdtd.to_grid(fdtd.um(0.5), grid_spacing)
det_z_T = fdtd.to_grid(fdtd.um(5), grid_spacing)
det_z_R = fdtd.to_grid(fdtd.um(1), grid_spacing)
start_z = fdtd.to_grid(fdtd.um(2), grid_spacing)
end_z = fdtd.to_grid(fdtd.um(3), grid_spacing)
# === 儲存資料夾 ===

# ==== make_grid 函式封裝 ====
def make_grid(with_structure=True):
    grid = fdtd.Grid(
        shape=(Nx, 1, Nz),
        grid_spacing=grid_spacing,
        permittivity=1
    )
    grid[0, :, :] = fdtd.BlochBoundary(k_component=kx, length=Lx)
    grid[:, :, :10] = fdtd.PML()
    grid[:, :, -10:] = fdtd.PML()
    grid[:, 0, source_z] = fdtd.ComplexPlaneWave(
        wavelength=wavelength,
        period=wavelength / 3e8,
        amplitude=1.0 + 0j,
        theta_deg=theta_deg,
        polarization_axis="x",
        pulse=False,
        name="source"
    )
    grid[:, 0, det_z_T] = fdtd.LineDetector(name="T")
    grid[:, 0, det_z_R] = fdtd.LineDetector(name="R", flip_sign=True)
    simfolder = None
    if with_structure:
        grid[:, 0, start_z:end_z] = fdtd.Object(n=1.5, k=0, name="n=1.5")
        simfolder = grid.save_simulation("test_bloch_xz")
    return grid, simfolder

# ==== Step 1: Reference Run ====
ref_result = fdtd.reference_run(grid_factory=make_grid, total_steps=500)
print(f"Reference Transmission: {ref_result['T']:.3e}, Reflection: {ref_result['R']:.3e}")

# ==== Step 2: 真實模擬 ====
grid, simfolder = make_grid(with_structure=True)
for t in range(500):
    grid.step()
    if t % 10 == 0:
        fig = grid.visualize(y=0, 
                             animate=True, 
                             index=t, 
                             save=True, 
                             folder=simfolder, 
                             real_field_mode=True, 
                             real_component="Ex"
                             )
        plt.title(f"t = {t}")
        ax = plt.gca()
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("z (µm)")
        plt.tight_layout()
        clear_output(wait=True)

T_struct = np.mean(grid.T.S[-100:])
R_struct = np.mean(grid.R.S[-100:])
T = T_struct / ref_result["T"]
R = R_struct / ref_result["R"]
# A = 1 - T - R

print(f"Reference Transmission: {ref_result['T']:.3e}, Reflection: {ref_result['R']:.3e}")
print(f"Transmission: {T:.3f}, Reflection: {R:.3f}")

# ==== Step 3: Compute T and R ====
# T_struct_raw = np.mean(grid.T.S[-100:])
# R_struct_raw = np.mean(grid.R.S[-100:])  # 注意：此處不再補負號

# T_struct = T_struct_raw - ref_result["T"]
# R_struct = R_struct_raw - ref_result["R"]
# P_incident = ref_result["T"] + ref_result["R"]

# T = T_struct / P_incident
# R = R_struct / P_incident

# print(f"[Raw] T_struct = {T_struct_raw:.3e}, R_struct = {R_struct_raw:.3e}")
# print(f"[Corrected] Transmission = {T:.3f}, Reflection = {R:.3f}, T+R = {T+R:.3f}")