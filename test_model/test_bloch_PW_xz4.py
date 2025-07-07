"""
測試 ComplexPlaneWave + Bloch Boundary 組合
避免使用有問題的 PML，專注測試您的實際使用場景
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# 導入您修復後的 fdtd 庫
sys.path.append('.')
import fdtd

def test_complex_planewave_stability():
    """測試 ComplexPlaneWave 的數值穩定性"""
    
    print("🌊 測試 ComplexPlaneWave 穩定性")
    print("="*50)
    
    # 基本參數
    wavelength = 1550e-9
    grid_spacing = wavelength / 15  # λ/15 解析度
    
    # 創建較小的測試網格
    grid = fdtd.Grid(
        shape=(60, 1, 120),  # 較小的網格
        grid_spacing=grid_spacing,
        courant_number=0.3,  # 保守的 Courant 數
        force_complex=True
    )
    
    print(f"📏 測試參數:")
    print(f"   波長: {wavelength*1e9:.0f} nm")
    print(f"   網格間距: {grid_spacing*1e9:.1f} nm")
    print(f"   網格形狀: {grid.shape}")
    print(f"   時間步: {grid.time_step:.3e} s")
    print(f"   Courant數: {grid.courant_number:.3f}")
    
    # 添加 ComplexPlaneWave
    source_z = 20
    try:
        source = fdtd.ComplexPlaneWave(
            wavelength=wavelength,
            period=wavelength / 3e8,
            amplitude=0.05 + 0j,  # 小振幅
            theta_deg=0,
            polarization_axis="x",
            medium_n=1.0,
            name="plane_wave"
        )
        grid[:, 0, source_z] = source
        print(f"   ✅ ComplexPlaneWave 源創建成功")
        print(f"   源位置: z = {source_z}")
        print(f"   源振幅: {source.amplitude}")
    except Exception as e:
        print(f"   ❌ ComplexPlaneWave 創建失敗: {e}")
        return False, None
    
    # 添加檢測器 (避免使用 PML 區域)
    det_z1 = 40
    det_z2 = 80
    
    grid[30, 0, det_z1] = fdtd.LineDetector(name="detector1")
    grid[30, 0, det_z2] = fdtd.LineDetector(name="detector2")
    
    print(f"   檢測器位置: z = {det_z1}, {det_z2}")
    
    # 運行穩定性測試
    print(f"\n⚡ 運行穩定性測試...")
    
    max_steps = 150
    field_history = []
    
    for step in range(max_steps):
        grid.step()
        
        # 監控場值
        E_max = float(np.max(np.abs(grid.E)))
        H_max = float(np.max(np.abs(grid.H)))
        
        field_history.append({
            'step': step,
            'E_max': E_max,
            'H_max': H_max
        })
        
        # 檢查數值穩定性
        if np.isnan(E_max) or np.isnan(H_max):
            print(f"❌ 步驟 {step}: 出現 NaN！")
            return False, field_history
            
        if np.isinf(E_max) or np.isinf(H_max):
            print(f"❌ 步驟 {step}: 出現 Inf！")
            return False, field_history
            
        # 更寬鬆的閾值檢查
        if E_max > 100 or H_max > 1:
            print(f"❌ 步驟 {step}: 場值異常大！")
            print(f"   E_max = {E_max:.3e}")
            print(f"   H_max = {H_max:.3e}")
            return False, field_history
        
        # 報告進度
        if step % 30 == 0:
            print(f"   步驟 {step:3d}: E_max={E_max:.3e}, H_max={H_max:.3e}")
    
    print(f"✅ ComplexPlaneWave 穩定性測試通過！")
    print(f"   最終 E_max: {E_max:.3e}")
    print(f"   最終 H_max: {H_max:.3e}")
    
    return True, field_history

def test_bloch_boundary_stability():
    """測試 Bloch Boundary 的穩定性"""
    
    print("\n🔄 測試 Bloch Boundary 穩定性")
    print("="*50)
    
    # 創建測試網格
    wavelength = 1550e-9
    grid_spacing = wavelength / 12
    
    grid = fdtd.Grid(
        shape=(40, 1, 80),
        grid_spacing=grid_spacing,
        courant_number=0.25,  # 更保守
        force_complex=True
    )
    
    print(f"📏 Bloch 測試參數:")
    print(f"   網格形狀: {grid.shape}")
    print(f"   Courant數: {grid.courant_number:.3f}")
    
    # 設定 Bloch 邊界條件
    try:
        # x 方向 Bloch 邊界
        k_component = 0.1  # 小的波數分量
        length = grid.Nx * grid_spacing
        
        grid[0, :, :] = fdtd.BlochBoundary(
            k_component=k_component,
            length=length,
            name="bloch_x"
        )
        
        print(f"   ✅ Bloch 邊界創建成功")
        print(f"   k_component: {k_component}")
        print(f"   length: {length*1e6:.2f} μm")
        
    except Exception as e:
        print(f"   ❌ Bloch 邊界創建失敗: {e}")
        return False, None
    
    # 添加簡單源
    grid[20, 0, 20] = fdtd.PointSource(
        period=50,
        amplitude=0.02,  # 很小的振幅
        name="test_source"
    )
    
    # 運行測試
    print(f"\n⚡ 運行 Bloch 邊界測試...")
    
    max_steps = 100
    
    for step in range(max_steps):
        grid.step()
        
        E_max = float(np.max(np.abs(grid.E)))
        H_max = float(np.max(np.abs(grid.H)))
        
        # 檢查穩定性
        if np.isnan(E_max) or np.isnan(H_max) or E_max > 10 or H_max > 0.1:
            print(f"❌ 步驟 {step}: Bloch 邊界不穩定！")
            print(f"   E_max = {E_max:.3e}, H_max = {H_max:.3e}")
            return False, None
        
        if step % 25 == 0:
            print(f"   步驟 {step:3d}: E_max={E_max:.3e}, H_max={H_max:.3e}")
    
    print(f"✅ Bloch Boundary 穩定性測試通過！")
    return True, None

def test_combined_complex_planewave_bloch():
    """測試 ComplexPlaneWave + Bloch Boundary 組合"""
    
    print("\n🎯 測試 ComplexPlaneWave + Bloch Boundary 組合")
    print("="*60)
    
    # 參數設定
    wavelength = 1550e-9
    grid_spacing = wavelength / 20
    
    # 創建網格
    grid = fdtd.Grid(
        shape=(80, 1, 160),
        grid_spacing=grid_spacing,
        courant_number=0.3,
        force_complex=True
    )
    
    print(f"📏 組合測試參數:")
    print(f"   波長: {wavelength*1e9:.0f} nm")
    print(f"   網格: {grid.shape}")
    print(f"   解析度: λ/{wavelength/grid_spacing:.0f}")
    
    # 添加 Bloch 邊界 (x 方向)
    k_x = 2 * np.pi / wavelength * 0.1  # 小的橫向波數
    length_x = grid.Nx * grid_spacing
    
    grid[0, :, :] = fdtd.BlochBoundary(
        k_component=k_x,
        length=length_x,
        name="bloch_x"
    )
    
    # 添加 ComplexPlaneWave
    source_z = 30
    source = fdtd.ComplexPlaneWave(
        wavelength=wavelength,
        period=wavelength / 3e8,
        amplitude=0.03 + 0j,
        theta_deg=0,
        polarization_axis="x",
        name="plane_wave"
    )
    grid[:, 0, source_z] = source
    
    # 添加介電結構測試
    structure_z_start = 80
    structure_z_end = 100
    structure_thickness = (structure_z_end - structure_z_start) * grid_spacing * 1e6
    
    # n=1.3 的小折射率差
    grid.inverse_permittivity[:, 0, structure_z_start:structure_z_end, :] = 1.0 / (1.3**2)
    
    print(f"   介電結構: z={structure_z_start}:{structure_z_end} (厚度 {structure_thickness:.2f}μm)")
    print(f"   折射率: n=1.3")
    
    # 添加檢測器
    refl_z = source_z + 15
    trans_z = structure_z_end + 20
    
    grid[40, 0, refl_z] = fdtd.LineDetector(name="R", flip_sign=True)
    grid[40, 0, trans_z] = fdtd.LineDetector(name="T", flip_sign=False)
    
    print(f"   反射檢測器: z={refl_z}")
    print(f"   穿透檢測器: z={trans_z}")
    
    # 運行組合測試
    print(f"\n⚡ 運行組合模擬...")
    
    total_steps = 500
    check_interval = 100
    
    for step in range(total_steps):
        grid.step()
        
        if step % check_interval == 0:
            E_max = float(np.max(np.abs(grid.E)))
            H_max = float(np.max(np.abs(grid.H)))
            
            print(f"   步驟 {step:3d}: E_max={E_max:.3e}, H_max={H_max:.3e}")
            
            # 檢查穩定性
            if np.isnan(E_max) or np.isnan(H_max):
                print(f"❌ 組合測試失敗：NaN 出現於步驟 {step}")
                return False, None
                
            if E_max > 50 or H_max > 0.5:
                print(f"❌ 組合測試失敗：場值過大於步驟 {step}")
                return False, None
    
    # 分析結果
    print(f"\n📊 組合測試結果分析...")
    
    if len(grid.T.S) > 50 and len(grid.R.S) > 50:
        # 取穩態數據
        steady_steps = 50
        T_data = np.array(grid.T.S[-steady_steps:])
        R_data = np.array(grid.R.S[-steady_steps:])
        
        P_T = np.mean(np.real(T_data))
        P_R = np.mean(np.abs(R_data))
        
        print(f"   穿透功率: {P_T:.3e}")
        print(f"   反射功率: {P_R:.3e}")
        
        # 簡單的能量分析
        if abs(P_T) > 1e-20 and abs(P_R) > 1e-20:
            ratio = P_T / (P_T + P_R)
            print(f"   穿透比例: {ratio:.3f}")
            
            if 0.1 < ratio < 0.9:
                print(f"   ✅ 能量分布合理")
            else:
                print(f"   ⚠️ 能量分布可能有問題")
        else:
            print(f"   ⚠️ 檢測器功率過小，可能需要調整")
    
    print(f"✅ ComplexPlaneWave + Bloch Boundary 組合測試完成！")
    return True, {
        'final_E_max': E_max,
        'final_H_max': H_max,
        'P_T': P_T if 'P_T' in locals() else 0,
        'P_R': P_R if 'P_R' in locals() else 0
    }

def plot_stability_comparison(complex_history, title="場穩定性測試"):
    """繪製穩定性對比圖"""
    
    if not complex_history:
        print("⚠️ 沒有歷史數據可繪製")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    steps = [h['step'] for h in complex_history]
    E_max = [h['E_max'] for h in complex_history]
    H_max = [h['H_max'] for h in complex_history]
    
    # E 場圖
    ax1.semilogy(steps, E_max, 'b-', linewidth=2, label='E_max')
    ax1.set_title('電場最大值')
    ax1.set_xlabel('時間步')
    ax1.set_ylabel('|E|_max (V/m)')
    ax1.grid(True)
    ax1.legend()
    
    # H 場圖
    ax2.semilogy(steps, H_max, 'r-', linewidth=2, label='H_max')
    ax2.set_title('磁場最大值')
    ax2.set_xlabel('時間步')
    ax2.set_ylabel('|H|_max (A/m)')
    ax2.grid(True)
    ax2.legend()
    
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig('complex_planewave_stability.png', dpi=150)
    plt.show()
    
    print(f"📈 穩定性圖已保存為 'complex_planewave_stability.png'")

def main():
    """主測試程序"""
    
    print("🚀 ComplexPlaneWave + Bloch Boundary 完整測試")
    print("="*70)
    
    success_count = 0
    total_tests = 3
    
    try:
        # 測試 1: ComplexPlaneWave 穩定性
        print("\n🧪 測試 1/3: ComplexPlaneWave 穩定性")
        test1_ok, complex_history = test_complex_planewave_stability()
        if test1_ok:
            success_count += 1
            print("✅ 測試 1 通過")
        else:
            print("❌ 測試 1 失敗")
        
        # 測試 2: Bloch Boundary 穩定性
        print("\n🧪 測試 2/3: Bloch Boundary 穩定性")
        test2_ok, _ = test_bloch_boundary_stability()
        if test2_ok:
            success_count += 1
            print("✅ 測試 2 通過")
        else:
            print("❌ 測試 2 失敗")
        
        # 測試 3: 組合測試
        if test1_ok and test2_ok:
            print("\n🧪 測試 3/3: 組合功能測試")
            test3_ok, results = test_combined_complex_planewave_bloch()
            if test3_ok:
                success_count += 1
                print("✅ 測試 3 通過")
            else:
                print("❌ 測試 3 失敗")
        else:
            print("\n⏭️ 跳過測試 3（前面的測試失敗）")
        
        # 繪製結果
        if complex_history and test1_ok:
            try:
                plot_stability_comparison(complex_history)
            except Exception as e:
                print(f"⚠️ 繪圖失敗: {e}")
        
    except Exception as e:
        print(f"\n💥 測試過程出錯: {e}")
        import traceback
        traceback.print_exc()
    
    # 總結
    print(f"\n📋 測試總結")
    print("="*40)
    print(f"通過測試: {success_count}/{total_tests}")
    
    if success_count == total_tests:
        print("🎉 所有測試通過！您的修復非常成功！")
        print("\n💡 下一步建議:")
        print("   1. 現在可以安全使用 ComplexPlaneWave + Bloch")
        print("   2. 避免使用 PML (可能仍有問題)")
        print("   3. 可以開始真實的光子晶體模擬")
        
    elif success_count >= 2:
        print("👍 大部分測試通過，基本功能正常")
        print("⚠️ 建議先解決失敗的測試再進行複雜模擬")
        
    else:
        print("😞 多數測試失敗，需要進一步調試")
        print("🔧 建議檢查:")
        print("   1. 確認 grid.py 修復正確")
        print("   2. 檢查 ComplexPlaneWave 實現")
        print("   3. 檢查 Bloch boundary 實現")
    
    return success_count == total_tests

if __name__ == "__main__":
    success = main()
    
    if success:
        print(f"\n✅ 完整測試成功 - 可以開始使用了！")
    else:
        print(f"\n❌ 測試未完全通過 - 需要進一步調試")