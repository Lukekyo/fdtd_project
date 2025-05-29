import os
import fdtd
import numpy as np
import matplotlib.pyplot as plt
from fdtd.fdtd_helper import um, nm, to_grid, from_grid
# fdtd_helper.py

fdtd.set_backend("numpy")

# === 模擬參數 ===
grid_spacing = 20e-9  # 20 nm
x_span = um(5)        # 5 µm
y_span = um(10)       # 10 µm

Nx = to_grid(x_span, grid_spacing)
Ny = to_grid(y_span, grid_spacing)

# 建立 Grid
grid = fdtd.Grid(
    grid_spacing = grid_spacing,
    shape = (Ny, Nx, 1),
    permittivity = 1
)

# === 邊界條件 ===
grid[:10, :, :] = fdtd.PML(name="pml_top")
grid[-10:, :, :] = fdtd.PML(name="pml_bottom")
# grid[:, 0:1, 0] = fdtd.PeriodicBoundary(name="periodic_left")
grid[:, -1:, 0] = fdtd.PeriodicBoundary(name="periodic")

simfolder = grid.save_simulation("test")  # initializing environment to save simulation data
print(simfolder)

# === 定義材料層 ===
layer_A_THK = to_grid(um(3), grid_spacing)  # 1 µm
layer_B_THK = to_grid(um(0.1), grid_spacing)  # 1 µm
layer_C_THK = to_grid(um(3), grid_spacing)  # 1 µm
start = to_grid(x_span/2, grid_spacing) # 起始位置 y = 2 µm

# grid[start:start+layer_A_THK, :, 0] = fdtd.Object(permittivity=2.25, name="Layer_A")
grid[start+layer_A_THK:start+layer_A_THK+layer_B_THK, :, 0] = fdtd.Object(permittivity=4.0, name="Layer_B")
# grid[start-layer_C_THK:start, :, 0] = fdtd.Object(permittivity=6.25, name="Layer_C")

# === 放置平面波光源 ===
source_y = start-layer_A_THK - to_grid(um(0.5), grid_spacing)
source_y = start+layer_A_THK+layer_B_THK -50
grid[source_y:source_y+1, :, 0] = fdtd.LineSource(
    period=1550e-9 / 3e8,
    name="source",
    pulse=False,         # ✅ 加上這一行
    cycle=5,            # 🔁 脈衝長度（建議 5–10）
    hanning_dt=grid.time_step  # 或略大一點
)
# === 偵測器放在最底層下方記錄穿透波 ===
det_y = start+layer_A_THK+layer_B_THK +50
grid[det_y:det_y+1, :, 0] = fdtd.LineDetector(name="detector")

# === 執行模擬 ===

from IPython.display import clear_output
for i in range(300):
    grid.step()
    if i % 10 == 0:
        grid.visualize(z=0, animate=True, index=i, save=True, folder=simfolder)
        plt.title(f"t = {i}")
        clear_output(wait=True)

grid.save_data()

# === 繪製穿透強度 ===
detector_data = np.array(grid.detector.detector_values()["E"])  # 將 list 轉為 ndarray
Ez = detector_data[..., 2]  # 取出 Ez 分量
intensity = np.sum(np.abs(Ez)**2, axis=1)
plt.figure()
plt.plot(intensity)
plt.title("Transmitted Intensity vs Time")
plt.xlabel("Time step")
plt.ylabel("Intensity")
plt.grid(True)
plt.show()

# === 輸出網格與物件資訊 ===
with open(os.path.join(simfolder, "grid_info.txt"), "w") as f:
    f.write(str(grid))
    wavelength = 3e8 / grid.source.frequency
    wavelengthUnits = wavelength / grid.grid_spacing
    GD = np.array([grid.x, grid.y, grid.z])
    gridRange = [np.arange(x/grid.grid_spacing) for x in GD]
    objectRange = np.array([[gridRange[0][x.x], gridRange[1][x.y], gridRange[2][x.z]] for x in grid.objects], dtype=object).T
    f.write("\n\nGrid details (in wavelength scale):")
    f.write("\n\tGrid dimensions: ")
    f.write(str(GD / wavelength))
    f.write("\n\tSource dimensions: ")
    f.write(str(np.array([grid.source.x[-1] - grid.source.x[0] + 1, grid.source.y[-1] - grid.source.y[0] + 1, grid.source.z[-1] - grid.source.z[0] + 1]) / wavelengthUnits))
    f.write("\n\tObject dimensions: ")
    f.write(str([(max(map(max, x)) - min(map(min, x)) + 1)/wavelengthUnits for x in objectRange]))

dic = np.load(os.path.join(simfolder, "detector_readings.npz"))
import warnings; warnings.filterwarnings("ignore") # TODO: fix plot_detection to prevent warnings
fdtd.plot_detection(dic)