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

# 將這些函數加入到您的 test_bloch_PW_xz3.py 文件中

def check_actual_source_power(grid, grid_spacing):
    """檢查源實際注入的功率"""
    
    print("\n🔬 檢查實際源功率")
    print("=" * 40)
    
    # 獲取源位置的實際E場值
    source = grid.sources[0]
    source_z = source.z[0] if hasattr(source, 'z') else 0
    
    print(f"源類型: {source.__class__.__name__}")
    print(f"源amplitude設定: {getattr(source, 'amplitude', 'N/A')}")
    print(f"源位置: z={source_z}")
    print(f"源x範圍: {source.x[0]} 到 {source.x[-1]} (共{len(source.x)}個點)")
    
    # 檢查源位置的實際E場值
    total_E_squared = 0
    E_values = []
    
    for i, x_idx in enumerate(source.x[:10]):  # 只檢查前10個點
        E_val = grid.E[x_idx, 0, source_z, 0]  # Ex分量
        E_magnitude = abs(E_val)
        E_values.append(E_magnitude)
        total_E_squared += E_magnitude**2
        
        if i < 5:  # 只打印前5個
            print(f"   點{i}: E[{x_idx},0,{source_z}] = {E_magnitude:.6e} V/m")
    
    # 計算所有點的總和
    total_E_squared_all = 0
    for x_idx in source.x:
        E_val = grid.E[x_idx, 0, source_z, 0]
        total_E_squared_all += abs(E_val)**2
    
    # 實際源功率計算
    Z0 = 377.0
    actual_source_power = 0.5 * total_E_squared_all * grid_spacing / Z0
    
    print(f"\n💡 基於實際E場的源功率:")
    print(f"   源點數: {len(source.x)}")
    print(f"   平均|E|: {(total_E_squared_all/len(source.x))**0.5:.6e} V/m")
    print(f"   Σ|E|²: {total_E_squared_all:.6e}")
    print(f"   實際源功率: {actual_source_power:.6e} W/m")
    
    # 與原始amplitude比較
    if hasattr(source, 'amplitude'):
        expected_total_E_squared = len(source.x) * abs(source.amplitude)**2
        print(f"\n🔍 與設定值比較:")
        print(f"   預期Σ|E|² (如果每點都是amplitude): {expected_total_E_squared:.6e}")
        print(f"   實際/預期比例: {total_E_squared_all/expected_total_E_squared:.2f}")
    
    return actual_source_power

def debug_detector_power_detailed(grid, grid_spacing):
    """詳細調試檢測器功率計算"""
    
    print("\n🔬 詳細檢測器功率調試")
    print("=" * 40)
    
    results = {}
    
    # 檢查T檢測器
    if hasattr(grid, 'T'):
        detector = grid.T
        
        print(f"T檢測器:")
        print(f"   位置: z={detector.z[0] if hasattr(detector, 'z') else 'N/A'}")
        print(f"   點數: {len(detector.x) if hasattr(detector, 'x') else 'N/A'}")
        
        if len(detector.S) > 0:
            recent_power = detector.S[-1]
            print(f"   最新功率: {recent_power:.6e} W/m")
            results['T'] = recent_power
            
            # 檢查原始E, H場值
            if len(detector.E) > 0 and len(detector.H) > 0:
                E_latest = np.array(detector.E[-1])  # 最新的E場
                H_latest = np.array(detector.H[-1])  # 最新的H場
                
                print(f"   E場形狀: {E_latest.shape}")
                print(f"   H場形狀: {H_latest.shape}")
                
                if len(E_latest.shape) == 2:  # (N_points, 3)
                    E_max = np.max(np.abs(E_latest))
                    H_max = np.max(np.abs(H_latest))
                    
                    print(f"   E場最大值: {E_max:.6e} V/m")
                    print(f"   H場最大值: {H_max:.6e} A/m")
                    print(f"   阻抗比: {E_max/H_max:.1f} Ω" if H_max > 0 else "   H場為零")
                    
                    # 手動計算Poynting向量檢驗
                    mu0 = 4e-7 * np.pi
                    S_z_manual = np.real(E_latest[:, 0] * np.conj(H_latest[:, 1]) - 
                                       E_latest[:, 1] * np.conj(H_latest[:, 0])) / mu0
                    S_total_manual = np.sum(S_z_manual) * grid_spacing
                    
                    print(f"   S_z範圍: [{np.min(S_z_manual):.2e}, {np.max(S_z_manual):.2e}] W/m²")
                    print(f"   手動計算總功率: {S_total_manual:.6e} W/m")
                    
                    results['T_manual'] = S_total_manual
                    
                    # 檢查前幾個點的詳細計算
                    print(f"\n   前5個點的詳細計算:")
                    for i in range(min(5, len(S_z_manual))):
                        Ex = E_latest[i, 0]
                        Ey = E_latest[i, 1] 
                        Hx = H_latest[i, 0]
                        Hy = H_latest[i, 1]
                        S_point = (Ex * np.conj(Hy) - Ey * np.conj(Hx)) / mu0
                        print(f"     點{i}: Ex={abs(Ex):.2e}, Hy={abs(Hy):.2e}, S_z={np.real(S_point):.2e}")
    
    # 同樣檢查R檢測器
    if hasattr(grid, 'R'):
        print(f"\nR檢測器:")
        detector = grid.R
        if len(detector.S) > 0:
            recent_power = detector.S[-1]
            print(f"   最新功率: {recent_power:.6e} W/m")
            results['R'] = recent_power
    
    return results

def comprehensive_power_analysis(grid, grid_spacing):
    """綜合功率分析"""
    
    print("\n🚨 綜合功率分析")
    print("=" * 80)
    
    # 1. 檢查實際源功率
    actual_source_power = check_actual_source_power(grid, grid_spacing)
    
    # 2. 詳細檢測器調試
    detector_results = debug_detector_power_detailed(grid, grid_spacing)
    
    # 3. 比較和分析
    print(f"\n📊 功率比較分析:")
    print(f"   實際源功率: {actual_source_power:.6e} W/m")
    
    if 'T' in detector_results:
        T_power = detector_results['T']
        ratio = abs(T_power) / actual_source_power if actual_source_power > 0 else float('inf')
        
        print(f"   T檢測器功率: {T_power:.6e} W/m")
        print(f"   T/源功率比例: {ratio:.2e}")
        
        if ratio > 10:
            print(f"   ❌ T檢測器功率異常大！")
            print(f"      可能原因:")
            print(f"      - 源功率計算錯誤（amplitude定義問題）")
            print(f"      - H場數值錯誤（阻抗問題）")
            print(f"      - Poynting向量計算錯誤")
        elif ratio < 0.1:
            print(f"   ❌ T檢測器功率異常小！")
            print(f"      可能原因:")
            print(f"      - H場太小或為零") 
            print(f"      - 檢測器位置不對")
        else:
            print(f"   ✅ T檢測器功率合理")
    
    if 'R' in detector_results:
        R_power = detector_results['R']
        print(f"   R檢測器功率: {R_power:.6e} W/m")
    
    # 4. 理論期望檢查
    print(f"\n🎯 理論期望檢查:")
    
    # Fresnel反射理論 (n1=1 -> n2=1.5)
    r_theory = (1.0 - 1.5) / (1.0 + 1.5)  # r = -0.2
    R_theory = abs(r_theory)**2  # R = 0.04 = 4%
    T_theory = 1 - R_theory      # T = 0.96 = 96%
    
    print(f"   理論穿透率: T = {T_theory:.1%}")
    print(f"   理論反射率: R = {R_theory:.1%}")
    
    if actual_source_power > 0:
        expected_T_power = T_theory * actual_source_power
        expected_R_power = R_theory * actual_source_power
        
        print(f"   期望T功率: {expected_T_power:.6e} W/m")
        print(f"   期望R功率: {expected_R_power:.6e} W/m")
        
        if 'T' in detector_results:
            T_error = abs(detector_results['T'] - expected_T_power) / expected_T_power * 100
            print(f"   T功率誤差: {T_error:.1f}%")
    
    return {
        'actual_source_power': actual_source_power,
        'detector_results': detector_results
    }

# 修正您的主程式
def main_with_power_debug():
    """修正版主程式（包含功率調試）"""
    
    start_time = time.time()
    
    try:
        # ==== 步驟1: 創建網格 ====
        print(f"\n🔧 步驟1：網格設置")
        print("=" * 50)
        
        grid, simfolder = make_grid(with_structure=True)
        
        # ==== 步驟2: FDTD模擬執行 ====
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
        
        # ==== 🔬 新增：功率調試分析 ====
        power_results = comprehensive_power_analysis(grid, grid_spacing)
        
        # ==== 步驟3: 原有的結果計算 ====
        print(f"\n📊 步驟3：傳統結果計算")
        print("=" * 50)
        
        # 原有的計算方法
        results = fdtd.calculate_simple_TR_corrected(grid, grid_spacing)
        
        if results is None:
            print("❌ 計算失敗")
            return None
        
        # 顯示比較結果
        print(f"\n📋 計算方法比較:")
        print(f"   傳統方法 - T: {results['T']:.3f}, R: {results['R']:.3f}")
        
        # 基於實際源功率的修正計算
        if power_results['actual_source_power'] > 0:
            actual_source = power_results['actual_source_power']
            if 'T' in power_results['detector_results'] and 'R' in power_results['detector_results']:
                T_corrected = abs(power_results['detector_results']['T']) / actual_source
                R_corrected = abs(power_results['detector_results']['R']) / actual_source
                
                print(f"   修正方法 - T: {T_corrected:.3f}, R: {R_corrected:.3f}")
                print(f"   能量守恆: {T_corrected + R_corrected:.3f}")
                
                # 與理論對比
                r_theory = (1.0 - 1.5) / (1.0 + 1.5)
                R_theory = abs(r_theory)**2
                T_theory = 1 - R_theory
                
                print(f"   理論值 - T: {T_theory:.3f}, R: {R_theory:.3f}")
                
                T_error = abs(T_corrected - T_theory) / T_theory * 100
                R_error = abs(R_corrected - R_theory) / R_theory * 100
                
                print(f"   誤差 - T: ±{T_error:.1f}%, R: ±{R_error:.1f}%")
                
                return {
                    'T': T_corrected,
                    'R': R_corrected,
                    'T_error': T_error,
                    'R_error': R_error,
                    'method': 'corrected_source_power'
                }
        
        return results
        
    except Exception as e:
        print(f"\n💥 程式執行錯誤: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    print("🚀 啟動FDTD穿透反射計算程式（功率調試版）...")
    
    final_results = main_with_power_debug()  # 使用新的調試版本
    
    if final_results:
        print("\n✅ 程式執行完成！")
        
        # 顯示最終結果
        if 'method' in final_results and final_results['method'] == 'corrected_source_power':
            print(f"\n📋 修正後結果:")
            print(f"T={final_results['T']:.3f} ({final_results['T']*100:.1f}%), " + 
                  f"R={final_results['R']:.3f} ({final_results['R']*100:.1f}%), " +
                  f"誤差: T±{final_results['T_error']:.1f}%, R±{final_results['R_error']:.1f}%")
        else:
            print(f"\n📋 傳統方法結果:")
            print(f"T={final_results['T']:.3f} ({final_results['T']*100:.1f}%), " + 
                  f"R={final_results['R']:.3f} ({final_results['R']*100:.1f}%)")
    else:
        print("\n❌ 程式執行失敗")