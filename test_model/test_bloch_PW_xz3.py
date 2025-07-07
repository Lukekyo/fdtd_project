import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import time
import numpy as np
import fdtd

# 嘗試導入matplotlib
try:
    import matplotlib.pyplot as plt
    from IPython.display import clear_output
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False
    print("⚠️ matplotlib 或 IPython 不可用，將禁用可視化")

# ==== 模擬參數設定 ====
fdtd.set_backend("numpy")

# 基本參數
wavelength = fdtd.nm(1550)
grid_spacing = fdtd.nm(20)
x_span, z_span = fdtd.um(2), fdtd.um(6)
Nx = fdtd.to_grid(x_span, grid_spacing)
Nz = fdtd.to_grid(z_span, grid_spacing)

# 角度設定
theta_deg = 0  # 正常入射
theta = np.deg2rad(theta_deg)
k0 = 2 * np.pi / wavelength
kx = k0 * np.sin(theta)
Lx = Nx * grid_spacing

# 位置設定
pml_thickness = 5
source_z = pml_thickness + 8
structure_start_z = fdtd.to_grid(fdtd.um(2.5), grid_spacing)
structure_end_z = fdtd.to_grid(fdtd.um(3.5), grid_spacing)
det_z_R = structure_start_z - 20
det_z_T = structure_end_z + 20

# 源參數
source_amplitude = 1.0  # V/m

print("🔬" + "="*58 + "🔬")
print("        FDTD 穿透反射計算 (簡潔版)")
print("🔬" + "="*58 + "🔬")
print(f"📊 模擬參數:")
print(f"   波長: {wavelength*1e9:.0f} nm")
print(f"   網格間距: {grid_spacing*1e9:.0f} nm")
print(f"   網格大小: {Nx} × 1 × {Nz}")
print(f"   源振幅: {source_amplitude} V/m")

# def make_grid(with_structure=True):
#     """創建FDTD網格"""
    
#     print(f"\n🏗️  創建FDTD網格...")
    
#     # 創建網格
#     grid = fdtd.Grid(
#         shape=(Nx, 1, Nz),
#         grid_spacing=grid_spacing,
#         permittivity=1.0,
#         force_complex=True
#     )
    
#     # 邊界條件
#     grid[0, :, :] = fdtd.PeriodicBoundary(name="periodic_x")
#     grid[:, :, :pml_thickness] = fdtd.PML(name="pml_z_low")
#     grid[:, :, -pml_thickness:] = fdtd.PML(name="pml_z_high")
    
#     # 源設定
#     source = fdtd.ComplexPlaneWave(
#         wavelength=wavelength,
#         period=wavelength / 3e8,
#         amplitude=source_amplitude + 0j,
#         theta_deg=theta_deg,
#         polarization_axis="x",
#         pulse=False,
#         medium_n=1.0,
#         name="source"
#     )
#     grid[:, 0, source_z] = source
    
#     # 檢測器設定
#     grid[:, 0, det_z_T] = fdtd.LineDetector(name="T", flip_sign=False)
#     grid[:, 0, det_z_R] = fdtd.LineDetector(name="R", flip_sign=True)
    
#     # 結構設定
#     simfolder = None
#     if with_structure:
#         grid[:, 0, structure_start_z:structure_end_z] = fdtd.Object(n=1.5, k=0, name="structure")# 🔧 修正Object的permittivity bug
#         epsilon_correct = (1.5 + 1j * 0) ** 2  # ε = (n + ik)² = 2.25
#         inverse_epsilon_correct = 1.0 / epsilon_correct  # 1/ε = 0.444
#         grid.inverse_permittivity[:, 0, structure_start_z:structure_end_z, :] = inverse_epsilon_correct

#         print(f"✅ 修正後 inverse_permittivity: {inverse_epsilon_correct}")
#         print(f"   對應 permittivity: {1/inverse_epsilon_correct}")
#         print(f"   對應 n: {abs((1/inverse_epsilon_correct)**0.5)}")
#         if VISUALIZATION_AVAILABLE:
#             simfolder = grid.save_simulation("2D_transmission_reflection")
#         structure_thickness = (structure_end_z - structure_start_z) * grid_spacing * 1e6
#         print(f"   ✅ 添加結構：n=1.5, 厚度={structure_thickness:.2f}μm")
    
#     print(f"   ✅ 網格創建完成")
    
#     return grid, simfolder

def make_grid(with_structure=True):
    """創建FDTD網格（臨時解決方案）"""
    
    # 創建網格
    grid = fdtd.Grid(
        shape=(Nx, 1, Nz),
        grid_spacing=grid_spacing,
        permittivity=1.0,
        force_complex=True
    )
    
    # 邊界條件（保持不變）
    grid[0, :, :] = fdtd.PeriodicBoundary(name="periodic_x")
    grid[:, :, :pml_thickness] = fdtd.PML(name="pml_z_low")
    grid[:, :, -pml_thickness:] = fdtd.PML(name="pml_z_high")
    
    # 源設定（保持不變）
    source = fdtd.ComplexPlaneWave(
        wavelength=wavelength,
        period=wavelength / 3e8,
        amplitude=source_amplitude + 0j,
        theta_deg=theta_deg,
        polarization_axis="x",
        pulse=False,
        medium_n=1.0,
        name="source"
    )
    grid[:, 0, source_z] = source
    
    # 檢測器設定（保持不變）
    grid[:, 0, det_z_T] = fdtd.LineDetector(name="T", flip_sign=False)
    grid[:, 0, det_z_R] = fdtd.LineDetector(name="R", flip_sign=True)
    
    # 🔧 手動設定結構（完全避免 fdtd.Object）
    simfolder = None
    if with_structure:
        # 不使用 fdtd.Object，手動設定材料
        n, k = 1.5, 0.0
        epsilon = (n + 1j * k) ** 2  # ε = 2.25
        inverse_epsilon = 1.0 / epsilon  # 1/ε = 0.444
        
        # 直接設定 inverse_permittivity
        grid.inverse_permittivity[:, 0, structure_start_z:structure_end_z, :] = inverse_epsilon
        
        # 驗證設定
        sample = grid.inverse_permittivity[50, 0, (structure_start_z + structure_end_z)//2, 0]
        print(f"✅ 手動設定結構：n={n}, ε={epsilon:.3f}, 1/ε={inverse_epsilon:.3f}")
        print(f"   驗證樣本值: {sample:.6f}")
        
        structure_thickness = (structure_end_z - structure_start_z) * grid_spacing * 1e6
        print(f"   結構厚度={structure_thickness:.2f}μm")
        
        if VISUALIZATION_AVAILABLE:
            simfolder = grid.save_simulation("2D_transmission_reflection")
    
    print(f"   ✅ 網格創建完成")
    
    return grid, simfolder

def main():
    """主程式執行流程"""
    
    start_time = time.time()
    
    try:
        # ==== 步驟1: 創建網格 ====
        print(f"\n🔧 步驟1：網格設置")
        print("=" * 50)
        
        grid, simfolder = make_grid(with_structure=True)
        
        # ==== 步驟2: 真實模擬（保留原本動態） ====
        print(f"\n⚡ 步驟2：FDTD模擬執行")
        print("=" * 50)
        
        total_steps = 1200
        print(f"模擬步數: {total_steps}")
        
        simulation_start = time.time()
        
        for t in range(total_steps):
            grid.step()
            
            # 文字進度條
            if t % 100 == 0:
                elapsed = time.time() - simulation_start
                if t > 0:
                    eta = elapsed * (total_steps - t) / t
                    print(f"   進度: {t:3d}/{total_steps} ({t/total_steps*100:5.1f}%) - 已用{elapsed:.1f}s, 剩餘{eta:.1f}s")
                else:
                    print(f"   進度: {t:3d}/{total_steps} ({t/total_steps*100:5.1f}%) - 已用{elapsed:.1f}s")
        
        simulation_time = time.time() - simulation_start
        print(f"✅ 模擬完成！耗時: {simulation_time:.1f} 秒")
        
        # 🚨 在這裡加入診斷代碼！
        # ==== 步驟3: 結果計算 ====
        print(f"\n📊 步驟3：結果計算")
        print("=" * 50)
        # 1. 磁場問題診斷
        fdtd.comprehensive_magnetic_field_diagnosis_fixed(grid, grid_spacing)
        
        # 2. 參數診斷
        params, units_ok, fix_success = fdtd.comprehensive_parameter_diagnosis_fixed(grid, grid_spacing)
        fdtd.comprehensive_impedance_error_analysis(grid, grid_spacing)
        
        # 3. 根據診斷結果決定後續處理
        if fix_success:
            print(f"\n🎉 診斷成功！發現H場數值錯誤可以修正")
            print(f"建議：檢查FDTD庫的時間步長或更新係數設定")
        elif not units_ok:
            print(f"\n⚠️ 發現阻抗問題，H場數值異常")
            print(f"更新係數比例: {params['coeff_ratio']:.2e}")
        else:
            print(f"\n❓ 診斷結果不明確，可能需要更深入調查")
        
        # ==== 步驟3: 原有的結果計算 ====
        print(f"\n📊 步驟3：結果計算")
        print("=" * 50)
        
        # 原有的計算
        results = fdtd.calculate_simple_TR_corrected(grid, grid_spacing)
        
        if results is None:
            print("❌ 計算失敗")
            return None
        
        # 顯示原始結果
        print(f"\n📋 原始計算結果:")
        print(f"   T = {results['T']:.3f} ({results['T']*100:.1f}%)")
        print(f"   R = {results['R']:.3f} ({results['R']*100:.1f}%)")
        print(f"   能量守恆: {results['energy_conservation']:.3f}")
        
        # 如果診斷發現問題，提供修正建議
        if not units_ok and params:
            print(f"\n💡 基於診斷的修正建議:")
            correction_factor = 377.0 / (np.max(np.abs(grid.E)) / np.max(np.abs(grid.H)))
            print(f"   H場應該除以: {correction_factor:.2f}")
            print(f"   這會大幅改善T/R計算結果")
        
        # ... 其餘原有代碼 ...
        
    except Exception as e:
        print(f"\n💥 程式執行錯誤: {e}")
        import traceback
        traceback.print_exc()
        return None


# ==== 程式入口 ====
if __name__ == "__main__":
    print("🚀 啟動FDTD穿透反射計算程式（調試版）...")
    
    final_results = main()
    
    if final_results:
        print("\n✅ 程式執行完成！")
        
        # 一行摘要
        print(f"\n📋 一行摘要:")
        print(f"T={final_results['T']:.3f} ({final_results['T']*100:.1f}%), " + 
              f"R={final_results['R']:.3f} ({final_results['R']*100:.1f}%), " +
              f"A={final_results['A']:.3f} ({final_results['A']*100:.1f}%) " +
              f"[方法: {final_results['method_used']}]")
        
    else:
        print("\n❌ 程式執行失敗")