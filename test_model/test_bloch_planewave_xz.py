import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import time
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
from scipy.ndimage import gaussian_filter1d
import fdtd

# === 模擬參數 ===
fdtd.set_backend("numpy")
wavelength = fdtd.nm(850)
grid_spacing = fdtd.nm(20)
x_span, z_span = fdtd.um(2), fdtd.um(6)
Nx = fdtd.to_grid(x_span, grid_spacing)
Nz = fdtd.to_grid(z_span, grid_spacing)

# === Bloch 入射角設定（沿 x 入射無用，可保留結構）===
theta_deg = 00
theta = np.deg2rad(theta_deg)
k0 = 2 * np.pi / wavelength
kx = k0 * np.sin(theta)
Lx = Nx * grid_spacing

# === 建立 Grid (XZ 模擬) ===
grid = fdtd.Grid(
    grid_spacing=grid_spacing,
    shape=(Nx, 1, Nz),  # y=1 slice，模擬 x-z 平面
    permittivity=1
)

# === 邊界條件 ===
grid[0, :, :] = fdtd.BlochBoundary(k_component=kx, length=Lx, name="bloch")
# grid[0, :, :] = fdtd.PeriodicBoundary(name="periodic")
# grid[:10, :, :10] = fdtd.PML(name="pml_right")
# grid[-10:, :, -10:] = fdtd.PML(name="pml_left")
grid[:, :, :10] = fdtd.PML(name="pml_front")
grid[:, :, -10:] = fdtd.PML(name="pml_back")
grid.promote_dtypes_to_complex()

# === 儲存資料夾 ===
simfolder = grid.save_simulation("test_bloch_xz")

# === 材料層：z = 2–3 µm，n = 1.5 ===
start_x = fdtd.to_grid(fdtd.um(0.5), grid_spacing)
end_x = fdtd.to_grid(fdtd.um(1.5), grid_spacing)
start_z = fdtd.to_grid(fdtd.um(2), grid_spacing)
end_z = fdtd.to_grid(fdtd.um(3), grid_spacing)
grid[:, 0, start_z:end_z] = fdtd.Object(n=1.5, k=0, name="n=1.5")
# grid[start_x:end_x, 0, start_z:end_z] = fdtd.Object(n=1.5, k=0, name="n=1.5")

# === 光源：z = 1 µm，橫跨 x ===
source_z = fdtd.to_grid(fdtd.um(1), grid_spacing)

grid[:, 0, source_z:source_z+1] = fdtd.ComplexPlaneWave(
    wavelength=wavelength,
    period=wavelength / 3e8,
    amplitude=1.0 + 0j,
    theta_deg=theta_deg,  # ⬅️ 入射角度
    polarization_axis="x",  # ⬅️ 根據你目前是 XZ 平面使用 Ex
    pulse=False,
    name="my_directional_source"
)

# === 偵測器：穿透 z = 5 µm，反射 z = 0.5 µm（可選）===
det_z_T = fdtd.to_grid(fdtd.um(5), grid_spacing)
det_z_R = fdtd.to_grid(fdtd.um(0.5), grid_spacing)
field_x1 = fdtd.to_grid(fdtd.um(1), grid_spacing)
field_x2 = fdtd.to_grid(fdtd.um(4), grid_spacing)
field_z1 = fdtd.to_grid(fdtd.um(1.5), grid_spacing)
field_z2 = fdtd.to_grid(fdtd.um(3.5), grid_spacing)
grid[:, 0, det_z_T:det_z_T+1] = fdtd.LineDetector(name="detector_transmit")
grid[:, 0, det_z_R:det_z_R+1] = fdtd.LineDetector(name="detector_reflect")
# grid[field_x1:field_x2, 0, field_z1:field_z2] = fdtd.BlockDetector(name="detector_field")

# === 可視化模擬域 ===
fdtd.plot_simulation_domain(grid, 
                            plane='xz', 
                            mode='n', 
                            vmin=1.0, 
                            vmax=2.0, 
                            x_range=(1, 1.25), 
                            z_range=(2, 3.2)
                            )

print("等待 2 秒自動繼續模擬...")
time.sleep(2)

# === 執行模擬 ===
for t in range(500):
    grid.step()
    if t % 10 == 0:
        # fig = grid.visualize(y=0, animate=True, index=t, save=True, folder=simfolder)
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
        # ax.set_xticklabels([f"{x * grid_spacing * 1e6:.1f}" for x in ax.get_xticks()])
        # ax.set_yticklabels([f"{z * grid_spacing * 1e6:.1f}" for z in ax.get_yticks()])
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("z (µm)")
        plt.tight_layout()
        clear_output(wait=True)

# === 儲存模擬資料 ===
grid.save_data()

# === 後處理分析 ===
Ex_T = np.array(grid.detector_transmit.detector_values()["Ex"])
Ex_R = np.array(grid.detector_reflect.detector_values()["Ex"])
intensity_T = np.sum(np.abs(Ex_T)**2, axis=1)
intensity_R = np.sum(np.abs(Ex_R)**2, axis=1)
# === 穩態功率比 ===
steady_start = -100  # 使用最後 100 個時間步長的數據計算平均值
avg_T = np.mean(intensity_T[steady_start:])
avg_R = np.mean(intensity_R[steady_start:])
total_power = avg_T + avg_R  # 穩態的總功率
# === 計算相對強度 ===
relative_intensity_T = intensity_T / total_power  # 穿透相對強度
relative_intensity_R = intensity_R / total_power  # 反射相對強度

# === 切片位置 ===
center_idx = Ex_T.shape[1] // 2  # 取中間的切片
Ex_center_real = Ex_T[:, center_idx].real  # 取 Ex_T 在中心位置的實部
smoothed = gaussian_filter1d(Ex_center_real, sigma=5)  # 平滑處理

# === 時間軸 ===
time_array = np.arange(len(intensity_T)) * grid.time_step * 1e15

# === 繪圖：穿透與反射 ===
plt.figure()
plt.plot(time_array, Ex_center_real, label="Real(Ex) (z=5 µm)")
plt.plot(time_array, smoothed, label="smoothed (z=5 µm)")
# plt.plot(time_array, relative_intensity_T, label="Transmitted Intensity (z=5 µm)")
# plt.plot(time_array, relative_intensity_R, label="Reflected Intensity (z=0.5 µm)")
plt.title("Intensity vs Time step")
plt.xlabel("Time (fs)")
plt.ylabel("Intensity (|Ez|^2 summed)")
plt.legend(loc="upper right")
plt.grid(True)
plt.tight_layout()
plt.show()

# === 穩態功率比 ===
print(f"Average Transmitted Power: {avg_T:.4e}")
print(f"Average Reflected Power  : {avg_R:.4e}")
print(f"Transmission Ratio       : {avg_T / total_power:.2%}")
print(f"Reflection Ratio         : {avg_R / total_power:.2%}")

# === 繪製 BlockDetector 圖 ===
df = np.load(os.path.join(simfolder, "detector_readings.npz"))
# Ez = df["detector_field (E)"]
# intensity = np.abs(Ez)**2
# fdtd.dB_map_2D(intensity)
