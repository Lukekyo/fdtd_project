import os
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
import fdtd
from fdtd.fdtd_helper import um, nm, to_grid

# === 模擬參數 ===
fdtd.set_backend("numpy")
wavelength = nm(850)  # 波長 1550 nm
grid_spacing = nm(20)
x_span, y_span = um(5), um(6)
Nx = to_grid(x_span, grid_spacing)
Ny = to_grid(y_span, grid_spacing)

# === 入射角與 Bloch 相位設定 ===
theta_deg = 1
theta = np.deg2rad(theta_deg)
k0 = 2 * np.pi / wavelength
kx = k0 * np.sin(theta)
Lx = Nx * grid_spacing

# === 建立 Grid ===
grid = fdtd.Grid(
    grid_spacing=grid_spacing,
    shape=(Nx, Ny, 1),
    permittivity=1
)

# === 邊界條件 ===
grid[0, :, :] = fdtd.BlochBoundary(k_component=kx, length=Lx, name="bloch_left")
# grid[0, :, :] = fdtd.PeriodicBoundary(name="periodic")
grid[:, :10, :] = fdtd.PML(name="pml_top")
grid[:, -10:, :] = fdtd.PML(name="pml_bottom")
grid.promote_dtypes_to_complex()

# === 儲存資料夾 ===
simfolder = grid.save_simulation("test_bloch")

# === 材料：y = 3 ~ 3.5 µm，n = 1.5 ===
start_y1 = to_grid(um(3), grid_spacing)
end_y1 = to_grid(um(3.5), grid_spacing)
start_y2 = to_grid(um(2.5), grid_spacing)
end_y2 = to_grid(um(3), grid_spacing)
start_y3 = to_grid(um(2), grid_spacing)
end_y3 = to_grid(um(2.5), grid_spacing)
grid[:, start_y1:end_y1, 0] = fdtd.Object(n=1.5,k=0, name="n=1.5_1")
# grid[:, start_y2:end_y2, 0] = fdtd.Object(n=2,k=0, name="n=2")
# grid[:, start_y3:end_y3, 0] = fdtd.Object(n=1.5,k=0, name="n=1.5_2")

# === 光源：y = 4 µm，連續波 ===
source_y = to_grid(um(4), grid_spacing)

grid[:, source_y:source_y+1, 0] = fdtd.ComplexPlaneWave2D(
    wavelength=wavelength,
    period = wavelength / 3e8,
    amplitude=1.0 + 0j,
    pulse=False,
    cycle=1,
    name='plane_source'
)

# === 監視器：反射 y = 1 µm，穿透 y = 5 µm ===
det_y1 = to_grid(um(1), grid_spacing)
det_y2 = to_grid(um(5), grid_spacing)
det_x1 = to_grid(um(1), grid_spacing)
det_x2 = to_grid(um(4), grid_spacing)
grid[:, det_y1:det_y1+1, 0] = fdtd.LineDetector(name="detector_transmit")
grid[:, det_y2:det_y2+1, 0] = fdtd.LineDetector(name="detector_reflect")
grid[det_x1:det_x2, det_y1:det_y2, 0] = fdtd.BlockDetector(name="detector_field")

# === 執行模擬 ===
for t in range(400):
    grid.step()
    if t % 20 == 0:
        fig = grid.visualize(z=0, animate=True, index=t, save=True, folder=simfolder)
        plt.title(f"t = {t}")
        ax = plt.gca()
        ax.set_xticklabels([f"{x * grid_spacing * 1e6:.1f}" for x in ax.get_xticks()])
        ax.set_yticklabels([f"{y * grid_spacing * 1e6:.1f}" for y in ax.get_yticks()])
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        plt.tight_layout()
        clear_output(wait=True)

# === 儲存模擬資料 ===
grid.save_data()

# === 後處理 ===
Ez_T = np.array(grid.detector_transmit.detector_values()["E"])[..., 2]
Ez_R = np.array(grid.detector_reflect.detector_values()["E"])[..., 2]
intensity_T = np.sum(np.abs(Ez_T)**2, axis=1)
intensity_R = np.sum(np.abs(Ez_R)**2, axis=1)

# === 時間軸 ===
time_array = np.arange(len(intensity_T)) * grid.time_step * 1e15# 真實時間

# === 畫出穿透與反射 ===
plt.figure()
plt.plot(time_array, intensity_T, label="Transmitted Intensity (y=1 µm)")
plt.plot(time_array, intensity_R, label="Reflected Intensity (y=5 µm)")
plt.title("Intensity vs Time step")
plt.xlabel("Time (fs)")
plt.ylabel("Intensity (|Ez|^2 summed)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# === 計算穩態功率比（最後 100 步）===
steady_start = -100
avg_T = np.mean(intensity_T[steady_start:])
avg_R = np.mean(intensity_R[steady_start:])
total = avg_T + avg_R
print(f"Average Transmitted Power: {avg_T:.4e}")
print(f"Average Reflected Power  : {avg_R:.4e}")
print(f"Transmission Ratio       : {avg_T / total:.2%}")
print(f"Reflection Ratio         : {avg_R / total:.2%}")

# === 儲存結果 ===
df = np.load(os.path.join(simfolder, "detector_readings.npz"))
# Ez = df["detector_field (E)"][..., 2]  # 取出 E_z 分量
Ez = df["detector_field (E)"]
intensity = np.abs(Ez)**2  # 計算強度
fdtd.dB_map_2D(intensity)