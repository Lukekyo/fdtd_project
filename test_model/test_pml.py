import os
import fdtd
import numpy as np
import matplotlib.pyplot as plt
from fdtd.fdtd_helper import um, nm, to_grid

fdtd.set_backend("numpy")

# === 模擬參數 ===
wavelength = 1550e-9  # 1550 nm
grid_spacing = 20e-9  # 20 nm
x_span = um(5)        # 5 µm
y_span = um(6)       # 10 µm

Nx = to_grid(x_span, grid_spacing)
Ny = to_grid(y_span, grid_spacing)

# === 建立 Grid（順序為 (Nx, Ny, 1)） ===
grid = fdtd.Grid(
    grid_spacing=grid_spacing,
    shape=(Nx, Ny, 1),
    permittivity=1  # 空氣背景
)

# === 邊界條件 ===
grid[:10, :, :] = fdtd.PML(name="pml_left")
grid[-10:, :, :] = fdtd.PML(name="pml_right")
grid[:, :10, :] = fdtd.PML(name="pml_top")
grid[:, -10:, :] = fdtd.PML(name="pml_bottom")

simfolder = grid.save_simulation("test")

# === 定義元件：y = 3~3.5 µm, 厚度 0.5 µm, 折射率1.5 ===
start_y = to_grid(um(3), grid_spacing)
end_y = to_grid(um(3.5), grid_spacing)
grid[:, start_y:end_y, 0] = fdtd.Object(permittivity=1.5**2, name="material")

# === 定義光源：y = 7 µm, 向 -y 發射 ===
source_y = to_grid(um(4), grid_spacing)

grid[:, source_y:source_y+1, 0] = fdtd.LineSource(
    period=wavelength / 3e8,
    name="source",
    pulse=False,
    cycle=5,
    hanning_dt=grid.time_step,
    # direction="-y"
)
# === 監視器：y = 1um, 5um ===
det_y1 = to_grid(um(1), grid_spacing)
det_y2 = to_grid(um(5), grid_spacing)
grid[:, det_y1:det_y1+1, 0] = fdtd.LineDetector(name="detector_reflect")
grid[:, det_y2:det_y2+1, 0] = fdtd.LineDetector(name="detector_transmit")

# === 執行模擬 ===
from IPython.display import clear_output
for i in range(500):
    grid.step()
    if i % 10 == 0:
        fig = grid.visualize(z=0, animate=True, index=i, save=True, folder=simfolder)
        plt.title(f"t = {i}")
        ax = plt.gca()
        ax.set_xticklabels([f"{x * grid_spacing * 1e6:.1f}" for x in ax.get_xticks()])
        ax.set_yticklabels([f"{y * grid_spacing * 1e6:.1f}" for y in ax.get_yticks()])
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        plt.tight_layout()
        clear_output(wait=True)

# === 儲存結果 ===
grid.save_data()

# === 繪圖：反射與穿透強度 ===
Ez_R = np.array(grid.detector_reflect.detector_values()["E"])[..., 2]
Ez_T = np.array(grid.detector_transmit.detector_values()["E"])[..., 2]
intensity_R = np.sum(np.abs(Ez_R)**2, axis=1)
intensity_T = np.sum(np.abs(Ez_T)**2, axis=1)

plt.figure()
plt.plot(intensity_R, label="Reflected Intensity (y=1 µm)")
plt.plot(intensity_T, label="Transmitted Intensity (y=5 µm)")
plt.title("Intensity vs Time step")
plt.xlabel("Time step")
plt.ylabel("Intensity (|Ez|^2 summed)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# === 取得 Ez 並計算穿透／反射能量 ===
detector_data = np.array(grid.detector.detector_values()["E"])
Ez_transmit = detector_data[..., 2]  # Z 向電場（穿透端）

Ez_reflect = np.array(grid.detector2.detector_values()["E"])[..., 2]  # Z 向電場（反射端）

# === 計算瞬時強度（|Ez|^2） ===
intensity_transmit = np.sum(np.abs(Ez_transmit) ** 2, axis=1)
intensity_reflect = np.sum(np.abs(Ez_reflect) ** 2, axis=1)

# === 取穩態區段，例如最後 100 個時間步 ===
steady_start = -100
avg_power_transmit = np.mean(intensity_transmit[steady_start:])
avg_power_reflect = np.mean(intensity_reflect[steady_start:])

# === 正規化為能量比例（穿透率與反射率） ===
total_power = avg_power_transmit + avg_power_reflect
transmission_ratio = avg_power_transmit / total_power
reflection_ratio = avg_power_reflect / total_power

# === 印出結果 ===
print(f"Average Transmitted Power (last 100 steps): {avg_power_transmit:.4e}")
print(f"Average Reflected Power  (last 100 steps): {avg_power_reflect:.4e}")
print(f"Transmission Ratio: {transmission_ratio:.2%}")
print(f"Reflection Ratio  : {reflection_ratio:.2%}")