import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import time
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
Nx = fdtd.to_grid(x_span, grid_spacing)
Nz = fdtd.to_grid(z_span, grid_spacing)
theta_deg = 0
theta = np.deg2rad(theta_deg)
k0 = 2 * np.pi / wavelength
kx = k0 * np.sin(theta)
Lx = Nx * grid_spacing

pml_thickness = 10
source_z = pml_thickness + 5       # 在PML後面一點
det_z_R = fdtd.to_grid(fdtd.um(1.5), grid_spacing)  # 結構前面
start_z = fdtd.to_grid(fdtd.um(2), grid_spacing)     # 結構開始
end_z = fdtd.to_grid(fdtd.um(3), grid_spacing)       # 結構結束
det_z_T = fdtd.to_grid(fdtd.um(4.5), grid_spacing)  # 結構後面

print(f"📍 位置設定:")
print(f"   Source: z = {source_z} ({source_z * grid_spacing * 1e6:.2f} μm)")
print(f"   反射探測器: z = {det_z_R} ({det_z_R * grid_spacing * 1e6:.2f} μm)")
print(f"   結構: z = {start_z}-{end_z} ({start_z * grid_spacing * 1e6:.2f}-{end_z * grid_spacing * 1e6:.2f} μm)")
print(f"   穿透探測器: z = {det_z_T} ({det_z_T * grid_spacing * 1e6:.2f} μm)")

# ==== make_grid 函式封裝 ====

def make_grid(with_structure=True):
    """創建FDTD網格"""
    grid = fdtd.Grid(
    shape=(Nx, 1, Nz),
    grid_spacing=grid_spacing,
    permittivity=1.0,
    force_complex=True  # 確保複數支援
    )

    # 邊界條件
    grid[0, :, :] = fdtd.BlochBoundary(k_component=kx, length=Lx)
    grid[:, :, :pml_thickness] = fdtd.PML(name="pml_low")
    grid[:, :, -pml_thickness:] = fdtd.PML(name="pml_high")

    # Source
    grid[:, 0, source_z] = fdtd.ComplexPlaneWave(
        wavelength=wavelength,
        period=wavelength / 3e8,
        amplitude=1.0 + 0j,
        theta_deg=theta_deg,
        polarization_axis="x",
        pulse=False,
        name="source"
    )

    # 修正後的探測器設定
    grid[:, 0, det_z_T] = fdtd.LineDetector(name="T", flip_sign=False)  # 穿透
    grid[:, 0, det_z_R] = fdtd.LineDetector(name="R", flip_sign=True)   # 反射

    simfolder = None
    if with_structure:
        # 結構設定
        grid[:, 0, start_z:end_z] = fdtd.Object(n=1.5, k=0, name="structure")
        simfolder = grid.save_simulation("test_bloch_xz")
        
        print(f"添加結構：n=1.5, 厚度={(end_z-start_z)*grid_spacing*1e6:.2f}μm")

    return grid, simfolder
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


def run_simulation(with_structure=True, total_steps=500, animation_interval=10):
    """運行完整模擬"""
    print(f"\n開始模擬 ({'有結構' if with_structure else '無結構'})")
    print("=" * 60)

    # 創建網格
    grid, simfolder = make_grid(with_structure=with_structure)
    # 計算源功率
    P_incident = grid.source.get_source_power(grid_spacing=grid_spacing)
    # 運行模擬
    print(f"\n執行 {total_steps} 個時間步...")
    simulation_start = time.time()

    for t in range(total_steps):
        grid.step()
        # 動畫顯示
        if simfolder and t % animation_interval == 0:
            try:
                fig = grid.visualize(
                    y=0, 
                    animate=True, 
                    index=t, 
                    save=True, 
                    folder=simfolder, 
                    real_field_mode=True, 
                    real_component="Ex"
                )
                plt.title(f"t = {t} / {total_steps}")
                ax = plt.gca()
                ax.set_xlabel("x (μm)")
                ax.set_ylabel("z (μm)")
                plt.tight_layout()
                clear_output(wait=True)
            except Exception as e:
                print(f"可視化錯誤: {e}")
        
    simulation_time = time.time() - simulation_start
    print(f"模擬完成！耗時: {simulation_time:.1f} 秒")

    # 分析結果
    results = analyze_results(grid, ["T", "R"], P_incident)

    return results, grid


# ==== 主程式執行 ====

if __name__ == "__main__":
    try:
        results, grid = run_simulation(with_structure=True, total_steps=1200)
        P_incident = grid.source.get_source_power(grid_spacing=grid_spacing)
        print(f"\n📋 最終結果:")
        # print(f"   入射功率: P_incident = {P_incident:.3e} W")
        print(f"   穿透率: T = {results['T']['ratio']:.3f} ({results['T']['ratio']*100:.1f}%)")
        print(f"   反射率: R = {results['R']['ratio']:.3f} ({results['R']['ratio']*100:.1f}%)")
        
        # 🔧 理論比較（Fresnel公式）
        n1, n2 = 1.0, 1.5  # 真空 → 介電材料
        r_theory = ((n1 - n2) / (n1 + n2)) ** 2
        t_theory = 4*n1*n2 / (n1 + n2) ** 2
        
        print(f"\n🧮 理論比較（Fresnel公式）:")
        print(f"   理論穿透率: T = {t_theory:.3f}")
        print(f"   理論反射率: R = {r_theory:.3f}")

        plt.plot()
        plt.figure(figsize=(8,4))
        for name in ["T", "R"]:
            # 直接從 grid 取出 detector 物件
            det_obj = getattr(grid, name, None)
            if det_obj is not None and hasattr(det_obj, "S"):
                plt.plot(det_obj.S, label=f"{name} (S)")
        plt.xlabel("Time step")
        plt.ylabel("Poynting S (W)")
        plt.title("Detector S vs Time step")
        plt.legend()
        plt.tight_layout()
        plt.show()

        # det_obj = getattr(grid, "T", None)
        # if det_obj is not None and hasattr(det_obj, "E") and hasattr(det_obj, "H"):
        #     # 取每步 E, H 場 (shape=(N_line, 3))
        #     E_time_series = [E[0, :] for E in det_obj.E]  # 第0個點的 E, shape=(3,)
        #     H_time_series = [H[0, :] for H in det_obj.H]  # 第0個點的 H, shape=(3,)
        #     E_time_series = np.array(E_time_series)       # shape=(n_steps, 3)
        #     H_time_series = np.array(H_time_series)       # shape=(n_steps, 3)

        #     plt.figure(figsize=(10,6))
        #     # plt.plot(E_time_series[:, 0], label="Ex")
        #     # plt.plot(E_time_series[:, 1], label="Ey")
        #     # plt.plot(E_time_series[:, 2], label="Ez")
        #     plt.plot(H_time_series[:, 0], label="Hx")
        #     plt.plot(H_time_series[:, 1], label="Hy")
        #     plt.plot(H_time_series[:, 2], label="Hz")
        #     plt.xlabel("Time step")
        #     plt.ylabel("Field value")
        #     plt.title("Time step vs E/H (at T detector, point 0)")
        #     plt.legend()
        #     plt.tight_layout()
        #     plt.show()
            
        det_obj = getattr(grid, "T", None)
        if det_obj is not None and hasattr(det_obj, "E") and hasattr(det_obj, "H"):
            E_time_series = [E[0, :] for E in det_obj.E]
            H_time_series = [H[0, :] for H in det_obj.H]
            E_time_series = np.array(E_time_series)
            H_time_series = np.array(H_time_series)
            print("Max Ex:", np.max(np.abs(E_time_series[:, 0])))
            print("Max Ey:", np.max(np.abs(E_time_series[:, 1])))
            print("Max Ez:", np.max(np.abs(E_time_series[:, 2])))
            print("Max Hx:", np.max(np.abs(H_time_series[:, 0])))
            print("Max Hy:", np.max(np.abs(H_time_series[:, 1])))
            print("Max Hz:", np.max(np.abs(H_time_series[:, 2])))
            print("E/H (阻抗):", np.max(np.abs(E_time_series[:, 0])) / np.max(np.abs(H_time_series[:, 1])))
            print("μ₀·E/H (應為377Ω):", fdtd.backend.mu0 * np.max(np.abs(E_time_series[:, 0])) / np.max(np.abs(H_time_series[:, 1])))
            print("μ₀:", fdtd.backend.mu0)
            # ...existing code...
        print
    except Exception as e:
        print(f"\n💥 程式執行錯誤: {e}")
        import traceback
        traceback.print_exc()

    finally:
        print(f"\n✅ 程式執行完成")