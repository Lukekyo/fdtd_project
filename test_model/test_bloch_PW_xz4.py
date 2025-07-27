import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
from scipy.ndimage import gaussian_filter1d
import fdtd


fdtd.set_backend("numpy")

# ==== 模擬參數 ====
wavelength = fdtd.nm(1550)
grid_spacing = fdtd.nm(20)
x_span, z_span = fdtd.um(2), fdtd.um(6)
total_steps = 500
structure_enabled = True  # 是否添加結構
Nx = fdtd.to_grid(x_span, grid_spacing)
Nz = fdtd.to_grid(z_span, grid_spacing)
theta_deg = 0
theta = np.deg2rad(theta_deg)
k0 = 2 * np.pi / wavelength
kx = k0 * np.sin(theta)
Lx = Nx * grid_spacing

pml_thickness = 10
source_z = pml_thickness + 20
det_z_R = fdtd.to_grid(fdtd.um(1.5), grid_spacing)
start_z = fdtd.to_grid(fdtd.um(2), grid_spacing)
end_z = fdtd.to_grid(fdtd.um(3), grid_spacing)
det_z_T = fdtd.to_grid(fdtd.um(4.5), grid_spacing)

print(f"📍 位置設定:")
print(f"   Source: z = {source_z} ({source_z * grid_spacing * 1e6:.2f} μm)")
print(f"   反射探測器: z = {det_z_R} ({det_z_R * grid_spacing * 1e6:.2f} μm)")
print(f"   結構: z = {start_z}-{end_z} ({start_z * grid_spacing * 1e6:.2f}-{end_z * grid_spacing * 1e6:.2f} μm)")
print(f"   穿透探測器: z = {det_z_T} ({det_z_T * grid_spacing * 1e6:.2f} μm)")

def make_grid(with_structure=structure_enabled):
    """創建FDTD網格"""
    grid = fdtd.Grid(
        shape=(Nx, 1, Nz),
        grid_spacing=grid_spacing,
        permittivity=1.0
    )

    # 邊界條件
    grid[0, :, :] = fdtd.BlochBoundary(k_component=kx, length=Lx)
    grid[:, :, :pml_thickness] = fdtd.PML(name="pml_low")
    grid[:, :, -pml_thickness:] = fdtd.PML(name="pml_high")

    # 源設定
    grid[:, 0, source_z] = fdtd.ComplexPlaneWave(
        wavelength=wavelength,
        amplitude=1.0 + 0j,
        theta_deg=theta_deg,
        polarization_axis="x",
        pulse=True,  # ===== 👇 改用連續波以便觀察週期 =====
        name="source"
    )
    # 探測器設定
    grid[:, 0, det_z_T] = fdtd.LineDetector(name="T", flip_sign=False)
    grid[:, 0, det_z_R] = fdtd.LineDetector(name="R", flip_sign=True)

    simfolder = None
    if with_structure:
        grid[:, 0, start_z:end_z] = fdtd.Object(n=1.5, k=0, name="n = 1.5")
        simfolder = grid.save_simulation("test_bloch_xz_floport")
        print(f"添加結構：n=3, 厚度={(end_z-start_z)*grid_spacing*1e6:.2f}μm")

    return grid, simfolder

def run_simulation(with_structure=True, total_steps=500, animation_interval=10):
    import time
    """運行模擬，使用 Floport 可視化"""
    print(f"\n開始模擬 ({'有結構' if with_structure else '無結構'}) - 使用 Floport 可視化")
    print("=" * 60)

    # 創建網格
    # grid, simfolder = make_grid(with_structure=with_structure)
    
    print(f"\n執行 {total_steps} 個時間步...", flush=True)
    print(f"每 {animation_interval} 步顯示一次動畫", flush=True)
    
    simulation_start = time.time()

    for t in range(total_steps):
        grid.step()
        # 每100步顯示進度，避免被清除
        # if t % 10 == 0:
        #     print(f"進度: {t}/{total_steps}", flush=True)
        # 動畫顯示 - 使用 Floport 可視化
        if USE_ANIMATION and simfolder and t % animation_interval == 0:
            try:
                fig = grid.visualize(
                    y=0, 
                    animate=True, 
                    index=t, 
                    save=True, 
                    folder=simfolder
                )
                plt.title(f"t = {t} / {total_steps}")
                ax = plt.gca()
                ax.set_xlabel("x (um)")
                ax.set_ylabel("z (um)")
                plt.tight_layout()
                clear_output(wait=True)
            except Exception as e:
                print(f"可視化錯誤: {e}")
        else:
            pass
        
    simulation_time = time.time() - simulation_start
    

    # 分析結果
    results = fdtd.analyze_results(grid, ["T", "R"])
    
    return results, grid, simulation_time

def field_visualization(grid):
    """快速場可視化，用於檢查波長"""
    print(f"\n🔍 場分佈快速檢查")
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # 1. 使用 Floport 風格總能量
    plt.sca(axes[0, 0])
    grid.visualize(y=0, show=False, cmap="Blues")
    plt.title("Total Energy Density")
    
    # 2. Ex 實部
    plt.sca(axes[0, 1])
    Ex_real = np.real(grid.E[:, 0, :, 0])
    plt.imshow(Ex_real.T, cmap="bwr", aspect="equal", origin="lower", vmin=-1, vmax=1)
    plt.title("Ex real part")
    plt.xlabel("X (grid)")
    plt.ylabel("Z (grid)")
    plt.colorbar()
    
    # 3. Z 方向剖面
    plt.sca(axes[1, 0])
    z_profile = Ex_real[Nx//2, :]  # 中間 X 位置
    z_indices = np.arange(len(z_profile))
    plt.plot(z_indices*grid_spacing, z_profile)
    plt.xlabel("Z (um)")
    plt.ylabel("Ex real part")
    plt.title("Z cross section of Ex real part")
    plt.grid(True)
    
    
    # 3. 資訊摘要
    plt.sca(axes[1, 1])
    plt.axis('off')
    info_lines = [
        f"Wavelength: {wavelength*1e9:.1f} nm",
        f"Grid spacing: {grid_spacing*1e9:.1f} nm",
        f"Total grid point: {Nx} × {Nz}",
        f"Simulation size (X(um) * Z(um)): {Nx*grid_spacing*1e6:.2f} × {Nz*grid_spacing*1e6:.2f} μm",
        "",
        f"Ex max: {np.max(np.abs(grid.E[:,:,:,0])):.3e}",
        f"Hy max: {np.max(np.abs(grid.H[:,:,:,1])):.3e}",
        f"Impedance (Ex/Hy): {np.max(np.abs(grid.E[:,:,:,0])) / np.max(np.abs(grid.H[:,:,:,1])):.1f} Ω",
        f"Simulation times (s): {simulation_time:.1f} second"
    ]
    
    for i, line in enumerate(info_lines):
        plt.text(0.1, 0.9 - i*0.08, line, fontsize=12, transform=axes[1,1].transAxes)
    
    plt.tight_layout()
    plt.show()

# ==== 主程式執行 ====
if __name__ == "__main__":
    try:
        # 選擇運行方式
        USE_ANIMATION = True  # 是否使用動畫模擬
        
        grid, simfolder = make_grid(with_structure=structure_enabled)
        grid.source.enable_monitoring()
        results, grid, simulation_time = run_simulation(
            with_structure=structure_enabled, 
            total_steps=total_steps,
            animation_interval=10  # 每n步顯示一次
        )
        # 快速場檢查
        grid.source.analyze_source_output()
        grid.source.plot_source_analysis()
        grid.source.get_source_power(grid_spacing)
        field_visualization(grid)
        
    except Exception as e:
        print(f"\n💥 程式執行錯誤: {e}")
        import traceback
        traceback.print_exc()