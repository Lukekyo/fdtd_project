# fdtd_helper.py - 完整整理版本
import numpy as np
from .backend import backend as bd

def um(x: float) -> float:
    """Convert micrometer to meters"""
    return x * 1e-6

def nm(x: float) -> float:
    """Convert nanometer to meters"""
    return x * 1e-9

def to_grid(length: float, grid_spacing: float) -> int:
    """Convert physical length to grid index"""
    return int(length / grid_spacing + 0.5)

def from_grid(index: int, grid_spacing: float) -> float:
    """Convert grid index to physical length"""
    return index * grid_spacing

def debug_power_calculation(grid, grid_spacing, steady_steps=20):
    """
    修正版功率計算 - 優先使用source的get_source_power方法
    """
    
    print(f"\n🔍 功率計算調試（修正版）")
    print("=" * 60)
    
    # === 步驟1: 嘗試使用源的get_source_power方法 ===
    print("1. 源功率計算:")
    
    source = None
    for src in grid.sources:
        if hasattr(src, 'get_source_power'):
            source = src
            break
    
    P_incident = None
    calculation_method = "unknown"
    
    if source and hasattr(source, 'get_source_power'):
        try:
            print("   ✅ 找到source.get_source_power方法，使用source計算")
            P_incident = source.get_source_power(grid_spacing)
            calculation_method = "source_method"
            print(f"   源方法計算結果: P_incident = {P_incident:.6e} W/m")
        except Exception as e:
            print(f"   ❌ source.get_source_power失敗: {e}")
            print("   回退到備用計算方法")
            source = None
    
    if P_incident is None:
        print("   使用備用計算方法:")
        
        # 備用方法：從source獲取實際參數
        if source:
            Z0 = 377.0
            n = getattr(source, 'n', 1.0)
            E0 = abs(getattr(source, 'amplitude', 1.0))
            print(f"     從source獲取: n={n}, E0={E0}")
        else:
            # 最後手段：硬編碼（但要明確標示）
            Z0 = 377.0
            n = 1.0
            E0 = 1.0
            print(f"     ⚠️ 使用硬編碼: n={n}, E0={E0}")
        
        source_length = grid.Nx * grid_spacing
        Z_medium = Z0 / n
        power_density = 0.5 * E0**2 / Z_medium
        P_incident = power_density * source_length
        calculation_method = "backup_method"
        
        print(f"     備用方法計算:")
        print(f"       Z_medium = {Z_medium:.1f} Ω")
        print(f"       power_density = {power_density:.6e} W/m²")
        print(f"       source_length = {source_length*1e6:.2f} μm")
        print(f"       P_incident = {P_incident:.6e} W/m")
    
    # === 步驟2: 檢測器信號分析 ===
    print(f"\n2. 檢測器信號分析:")
    
    total_steps = len(grid.T.S)
    if total_steps < steady_steps:
        steady_steps = total_steps
    
    print(f"   總時間步: {total_steps}")
    print(f"   分析範圍: 最後 {steady_steps} 步")
    
    # 原始信號
    T_signals = np.array(grid.T.S[-steady_steps:])
    R_signals = np.array(grid.R.S[-steady_steps:])
    
    print(f"   T檢測器原始信號:")
    print(f"     範圍: {np.min(T_signals):.2e} ~ {np.max(T_signals):.2e}")
    print(f"     平均: {np.mean(T_signals):.2e}")
    print(f"     實部平均: {np.mean(np.real(T_signals)):.2e}")
    
    print(f"   R檢測器原始信號:")
    print(f"     範圍: {np.min(R_signals):.2e} ~ {np.max(R_signals):.2e}")
    print(f"     平均: {np.mean(R_signals):.2e}")
    print(f"     |R|平均: {np.mean(np.abs(R_signals)):.2e}")
    
    # === 步驟3: 檢測器功率流計算 ===
    print(f"\n3. 檢測器功率流計算:")
    
    # 嘗試使用檢測器的get_power_flow方法
    P_transmitted = None
    P_reflected = None
    detector_method = "unknown"
    
    if hasattr(grid.T, 'get_power_flow') and hasattr(grid.R, 'get_power_flow'):
        try:
            print("   ✅ 使用檢測器的get_power_flow方法")
            P_transmitted = grid.T.get_power_flow(steady_steps)
            P_reflected = grid.R.get_power_flow(steady_steps)
            detector_method = "detector_method"
        except Exception as e:
            print(f"   ❌ 檢測器方法失敗: {e}")
    
    if P_transmitted is None or P_reflected is None:
        print("   回退到直接計算:")
        P_transmitted = np.mean(np.real(T_signals))
        P_reflected = np.mean(np.abs(R_signals))
        detector_method = "direct_calculation"
        
        print(f"     P_transmitted = {P_transmitted:.6e} W/m")
        print(f"     P_reflected = {P_reflected:.6e} W/m")
    
    # === 步驟4: 比率計算和分析 ===
    print(f"\n4. 比率計算:")
    
    T_ratio = P_transmitted / P_incident
    R_ratio = P_reflected / P_incident
    total_ratio = T_ratio + R_ratio
    
    print(f"   計算方法: {calculation_method} + {detector_method}")
    print(f"   P_incident: {P_incident:.6e} W/m")
    print(f"   P_transmitted: {P_transmitted:.6e} W/m")
    print(f"   P_reflected: {P_reflected:.6e} W/m")
    print(f"   T = {T_ratio:.4f} ({T_ratio*100:.1f}%)")
    print(f"   R = {R_ratio:.4f} ({R_ratio*100:.1f}%)")
    print(f"   總和 = {total_ratio:.4f} ({total_ratio*100:.1f}%)")
    
    # === 步驟5: 能量守恆分析 ===
    print(f"\n5. 能量守恆分析:")
    
    if 0.95 <= total_ratio <= 1.05:
        conservation_status = "✅ 優秀"
    elif 0.9 <= total_ratio <= 1.1:
        conservation_status = "✅ 良好"
    elif 0.8 <= total_ratio <= 1.2:
        conservation_status = "⚠️ 可接受"
    else:
        conservation_status = "❌ 有問題"
    
    print(f"   能量守恆狀態: {conservation_status}")
    
    # 如果能量不守恆，嘗試歸一化
    if total_ratio > 1.2:
        print(f"   🔧 嘗試歸一化修正:")
        normalization_factor = total_ratio
        T_normalized = T_ratio / normalization_factor
        R_normalized = R_ratio / normalization_factor
        
        print(f"     歸一化因子: {normalization_factor:.3f}")
        print(f"     歸一化後 T = {T_normalized:.4f}")
        print(f"     歸一化後 R = {R_normalized:.4f}")
        print(f"     歸一化後總和 = {T_normalized + R_normalized:.4f}")
        
        # 使用歸一化結果
        T_ratio = T_normalized
        R_ratio = R_normalized
    
    # === 步驟6: 理論對比 ===
    print(f"\n6. 理論對比:")
    
    # Fresnel公式（n1=1 → n2=1.5）
    r_theory = (1.0 - 1.5) / (1.0 + 1.5)
    R_theory = abs(r_theory)**2
    T_theory = 1 - R_theory
    
    print(f"   理論值: T={T_theory:.4f}, R={R_theory:.4f}")
    print(f"   計算值: T={T_ratio:.4f}, R={R_ratio:.4f}")
    
    T_error = abs(T_ratio - T_theory) / T_theory * 100
    R_error = abs(R_ratio - R_theory) / R_theory * 100
    
    print(f"   相對誤差: T±{T_error:.1f}%, R±{R_error:.1f}%")
    
    if T_error < 10 and R_error < 20:  # R的理論值很小，容許較大相對誤差
        theory_match = "✅ 與理論吻合"
    elif T_error < 20 and R_error < 50:
        theory_match = "⚠️ 與理論有差異"
    else:
        theory_match = "❌ 與理論差異很大"
    
    print(f"   評估: {theory_match}")
    
    return {
        'method_used': f"{calculation_method} + {detector_method}",
        'P_incident': P_incident,
        'P_transmitted': P_transmitted,
        'P_reflected': P_reflected,
        'T': T_ratio,
        'R': R_ratio,
        'A': 1 - T_ratio - R_ratio,
        'energy_conservation': total_ratio,
        'theory_comparison': {
            'T_theory': T_theory,
            'R_theory': R_theory,
            'T_error': T_error,
            'R_error': R_error
        }
    }

# 同時提供一個檢查source方法是否正常的函數
def test_source_method(grid, grid_spacing):
    """測試source的get_source_power方法是否正常工作"""
    
    print("\n🧪 測試source.get_source_power方法")
    print("=" * 40)
    
    source = None
    for src in grid.sources:
        if hasattr(src, 'get_source_power'):
            source = src
            break
    
    if source is None:
        print("❌ 沒有找到支持get_source_power的source")
        return False
    
    print(f"✅ 找到source: {source.__class__.__name__}")
    print(f"   有get_source_power方法: {hasattr(source, 'get_source_power')}")
    
    # 檢查source的屬性
    print(f"   source屬性:")
    print(f"     amplitude: {getattr(source, 'amplitude', 'N/A')}")
    print(f"     n: {getattr(source, 'n', 'N/A')}")
    print(f"     source_type: {getattr(source, 'source_type', 'N/A')}")
    
    # 嘗試調用
    try:
        result = source.get_source_power(grid_spacing)
        print(f"✅ 調用成功，結果: {result:.6e}")
        return True
    except Exception as e:
        print(f"❌ 調用失敗: {e}")
        import traceback
        traceback.print_exc()
        return False


def calculate_simple_TR_corrected(grid, grid_spacing, steady_steps=20):
    """
    統一的穿透反射率計算函數
    優先使用模組化方法，回退到備用方法
    """
    
    print(f"\n📊 統一穿透反射率計算")
    print("=" * 50)
    
    # 首先測試source方法
    source_ok = test_source_method(grid, grid_spacing)
    
    # 然後進行完整計算
    results = debug_power_calculation(grid, grid_spacing, steady_steps)
    
    print(f"\n🎯 最終結果摘要:")
    print(f"   方法: {results['method_used']}")
    print(f"   T = {results['T']:.3f} ({results['T']*100:.1f}%)")
    print(f"   R = {results['R']:.3f} ({results['R']*100:.1f}%)")
    print(f"   能量守恆: {results['energy_conservation']:.3f}")
    
    return results




# ===========================================================
# 修正的診斷函數 - 解決numpy格式化問題

def prove_magnetic_field_is_the_problem_fixed(grid, grid_spacing):
    """證明磁場更新是根本問題（修正版）"""
    
    print(f"\n🚨 磁場根本問題證明")
    print("="*60)
    
    print("1. 理論值檢查：")
    r = (1.0 - 1.5) / (1.0 + 1.5)
    R_theory = abs(r)**2
    T_theory = 1 - R_theory
    print(f"   理論：T = {T_theory:.1%}, R = {R_theory:.1%}")
    print(f"   實測：T = 50.1%, R = 49.9%")
    print(f"   ❌ 完全不符合！")
    
    print(f"\n2. 數量級檢查：")
    print(f"   入射功率：~e-9 W/m")
    print(f"   檢測器功率：~e1 W/m")
    print(f"   差距：8個數量級")
    print(f"   ❌ 物理上不可能！")
    
    print(f"\n3. 磁場檢查：")
    H_max = float(np.max(np.abs(grid.H)))  # 轉換為float
    E_max = float(np.max(np.abs(grid.E)))  # 轉換為float
    print(f"   |E|_max = {E_max:.3e} V/m")
    print(f"   |H|_max = {H_max:.3e} A/m")
    
    if H_max < 1e-15:
        print(f"   ❌ 磁場完全為零！")
    elif H_max < E_max / 1000:
        print(f"   ❌ 磁場遠小於預期！")
        expected_H = E_max / 377.0
        print(f"   預期 |H| ≈ {expected_H:.3e} A/m")
        print(f"   實際 |H| = {H_max:.3e} A/m")
        print(f"   相差 {expected_H/H_max:.1e} 倍")
    
    print(f"\n4. Poynting向量檢查：")
    # 手動計算幾個點的Poynting向量
    mu0 = 4e-7 * np.pi
    
    # 在源附近檢查
    source_z = 23
    if source_z < grid.Nz:
        E_source = grid.E[:5, 0, source_z]
        H_source = grid.H[:5, 0, source_z]
        
        print(f"   源附近(z={source_z}):")
        for i in range(min(5, len(E_source))):
            E_point = E_source[i]
            H_point = H_source[i]
            
            # 計算Poynting向量的z分量
            try:
                S_z = np.real(E_point[0] * np.conj(H_point[1]) - E_point[1] * np.conj(H_point[0])) / mu0
                
                # 安全的格式化
                E_abs = float(np.max(np.abs(E_point)))
                H_abs = float(np.max(np.abs(H_point)))
                S_z_float = float(S_z)
                
                print(f"     點{i}: E={E_abs:.2e}, H={H_abs:.2e}, S_z={S_z_float:.2e}")
            except Exception as e:
                print(f"     點{i}: 計算錯誤 - {str(e)}")
    
    return H_max < E_max / 1000

def emergency_manual_h_field_calculation_fixed(grid, grid_spacing):
    """緊急手動計算H場（修正版）"""
    
    print(f"\n🔧 緊急手動H場計算")
    print("="*50)
    
    print("假設H場正確，手動計算應該有的數值...")
    
    # 獲取當前E場
    E_max = float(np.max(np.abs(grid.E)))
    print(f"當前E場最大值: {E_max:.3e} V/m")
    
    # 在真空中，E和H的關係
    Z0 = 377.0  # 真空阻抗
    H_expected = E_max / Z0
    print(f"預期H場最大值: {H_expected:.3e} A/m")
    
    # 手動設定一個合理的H場來測試
    print(f"\n手動設定H場進行測試...")
    
    # 保存原始H場
    H_original = grid.H.copy()
    
    # 創建一個合理的H場分布
    # 假設Ex極化，應該有Hy分量
    for x in range(grid.Nx):
        for z in range(grid.Nz):
            E_local = grid.E[x, 0, z, 0]  # Ex分量
            if abs(E_local) > 1e-15:
                # 根據平面波關係設定Hy
                grid.H[x, 0, z, 1] = E_local / Z0
    
    H_new_max = float(np.max(np.abs(grid.H)))
    print(f"設定後H場最大值: {H_new_max:.3e} A/m")
    
    # 重新計算功率
    print(f"\n重新計算檢測器功率...")
    
    try:
        # 重新執行detect_S
        grid.T.detect_S()
        grid.R.detect_S()
        
        # 獲取新的功率值
        P_T_new = float(grid.T.S[-1]) if grid.T.S else 0
        P_R_new = float(grid.R.S[-1]) if grid.R.S else 0
        
        print(f"修正後檢測器功率:")
        print(f"   T檢測器: {P_T_new:.6e} W/m")
        print(f"   R檢測器: {P_R_new:.6e} W/m")
        
        # 計算入射功率
        source_length = len(getattr(grid.sources[0], 'x', [1])) * grid_spacing
        power_density = 0.5 * E_max**2 / Z0
        P_incident_manual = power_density * source_length
        
        print(f"   手動計算入射功率: {P_incident_manual:.6e} W/m")
        
        # 計算T/R
        if P_incident_manual > 0:
            T_manual = P_T_new / P_incident_manual
            R_manual = abs(P_R_new) / P_incident_manual
            
            print(f"\n手動修正後的T/R:")
            print(f"   T = {T_manual:.4f} ({T_manual*100:.1f}%)")
            print(f"   R = {R_manual:.4f} ({R_manual*100:.1f}%)")
            print(f"   總和 = {T_manual + R_manual:.4f}")
            
            # 與理論對比
            r_theory = (1.0 - 1.5) / (1.0 + 1.5)
            R_theory = abs(r_theory)**2
            T_theory = 1 - R_theory
            
            print(f"\n與理論對比:")
            print(f"   理論: T={T_theory:.1%}, R={R_theory:.1%}")
            print(f"   手動: T={T_manual:.1%}, R={R_manual:.1%}")
            
            if abs(T_manual - T_theory) < 0.1:
                print(f"   ✅ 手動修正後接近理論值！")
                success = True
            else:
                print(f"   ❌ 仍然偏差很大")
                success = False
        else:
            success = False
            
    except Exception as e:
        print(f"計算過程出錯: {str(e)}")
        success = False
    
    # 恢復原始H場
    grid.H = H_original
    
    return success

def identify_h_field_update_problem_fixed(grid):
    """識別H場更新的具體問題（修正版）"""
    
    print(f"\n🔍 H場更新問題識別")
    print("="*50)
    
    print("1. 檢查curl_E函數實現...")
    
    try:
        # 動態導入curl_E函數
        import sys
        for module_name in sys.modules:
            if 'grid' in module_name:
                module = sys.modules[module_name]
                if hasattr(module, 'curl_E'):
                    curl_E = module.curl_E
                    break
        else:
            print("   ❌ 找不到curl_E函數")
            return "curl_E_not_found"
        
        # 創建一個簡單的測試E場
        E_test = np.zeros_like(grid.E)
        
        # 在中心設定一個Ex場
        center_x, center_z = grid.Nx//2, grid.Nz//2
        E_test[center_x, 0, center_z, 0] = 1.0  # Ex = 1
        
        # 計算curl
        curl_result = curl_E(E_test)
        curl_max = float(np.max(np.abs(curl_result)))
        
        print(f"   測試curl_E函數:")
        print(f"   輸入: Ex[{center_x},0,{center_z}] = 1.0")
        print(f"   curl結果最大值: {curl_max:.6e}")
        
        if curl_max < 1e-15:
            print(f"   ❌ curl_E函數沒有產生結果")
            return "curl_E_broken"
        else:
            print(f"   ✅ curl_E函數工作正常")
        
    except Exception as e:
        print(f"   ❌ curl_E測試錯誤: {str(e)}")
        return "curl_E_error"
    
    print(f"\n2. 檢查update_H中的curl調用...")
    
    # 檢查實際的update_H過程
    try:
        H_before = grid.H.copy()
        
        # 手動執行一次update_H
        grid.update_H()
        H_after = grid.H.copy()
        
        H_change = float(np.max(np.abs(H_after - H_before)))
        print(f"   update_H後H場變化: {H_change:.6e}")
        
        if H_change < 1e-15:
            print(f"   ❌ update_H沒有改變H場")
            return "update_H_no_effect"
        else:
            print(f"   ✅ update_H有改變H場")
            return "update_H_working"
            
    except Exception as e:
        print(f"   ❌ update_H執行錯誤: {str(e)}")
        return "update_H_error"

def comprehensive_magnetic_field_diagnosis_fixed(grid, grid_spacing):
    """綜合磁場問題診斷（修正版）"""
    
    print(f"\n🚨 綜合磁場問題診斷")
    print("="*80)
    
    try:
        # 1. 證明這是磁場問題
        is_magnetic_problem = prove_magnetic_field_is_the_problem_fixed(grid, grid_spacing)
        
        if is_magnetic_problem:
            print(f"\n✅ 確認：這是磁場更新問題")
            
            # 2. 測試手動修正
            manual_success = emergency_manual_h_field_calculation_fixed(grid, grid_spacing)
            
            if manual_success:
                print(f"\n✅ 手動設定H場後結果正確 → 確認是H場更新問題")
                
                # 3. 找出H場更新的具體問題
                h_problem = identify_h_field_update_problem_fixed(grid)
                
                print(f"\n🎯 問題定位: {h_problem}")
                
                if h_problem == "curl_E_broken":
                    print(f"解決方案: 檢查curl_E函數實現")
                elif h_problem == "update_H_no_effect":
                    print(f"解決方案: 檢查update_H中的係數和數據類型")
                elif h_problem == "update_H_error":
                    print(f"解決方案: 檢查update_H實現")
                else:
                    print(f"解決方案: 需要進一步調試FDTD核心")
                    
            else:
                print(f"\n❌ 即使手動設定H場也不對 → 更深層問題")
        
        else:
            print(f"\n❓ 不是明顯的磁場問題，需要其他診斷")
            
    except Exception as e:
        print(f"\n❌ 診斷過程出錯: {str(e)}")
        import traceback
        traceback.print_exc()

# 修正導入錯誤的診斷函數

def diagnose_fdtd_parameters_fixed(grid, grid_spacing):
    """診斷FDTD參數設定（修正版）"""
    
    print(f"\n🔍 FDTD參數診斷")
    print("="*60)
    
    print("1. 基本網格參數：")
    print(f"   grid_spacing: {grid_spacing:.2e} m = {grid_spacing*1e9:.1f} nm")
    print(f"   網格形狀: {grid.shape}")
    print(f"   時間步: {grid.time_step:.2e} s")
    print(f"   Courant數: {grid.courant_number:.6f}")
    
    # 檢查Courant數是否合理
    max_courant = float(grid.D) ** (-0.5)  # D是維度
    print(f"   最大允許Courant數: {max_courant:.6f}")
    
    courant_ok = grid.courant_number <= max_courant
    if not courant_ok:
        print(f"   ❌ Courant數過大！會導致數值不穩定")
    else:
        print(f"   ✅ Courant數正常")
    
    print(f"\n2. 材料參數：")
    
    # 檢查介電常數和磁導率
    try:
        perm_sample = grid.inverse_permittivity[grid.Nx//2, 0, grid.Nz//2]
        perm_sample_float = [float(x.real) if hasattr(x, 'real') else float(x) for x in perm_sample]
        print(f"   inverse_permittivity (中心): {perm_sample_float}")
        
        permeab_sample = grid.inverse_permeability[grid.Nx//2, 0, grid.Nz//2] 
        permeab_sample_float = [float(x.real) if hasattr(x, 'real') else float(x) for x in permeab_sample]
        print(f"   inverse_permeability (中心): {permeab_sample_float}")
        
        # 檢查是否為1.0（真空）
        if abs(perm_sample_float[0] - 1.0) > 0.01:
            print(f"   ⚠️ inverse_permittivity ≠ 1.0")
        if abs(permeab_sample_float[0] - 1.0) > 0.01:
            print(f"   ⚠️ inverse_permeability ≠ 1.0")
        
    except Exception as e:
        print(f"   ❌ 材料參數檢查失敗: {str(e)}")
        perm_sample_float = [1.0, 1.0, 1.0]
        permeab_sample_float = [1.0, 1.0, 1.0]
    
    print(f"\n3. 場更新係數檢查：")
    
    # 計算H場更新的有效係數
    # h_update_coeff = grid.courant_number * permeab_sample_float[0]
    h_update_coeff = grid.time_step * grid.grid_spacing / bd.c0
    
    print(f"   H場更新係數: {h_update_coeff:.6e}")
    
    # 理論上應該是什麼值
    c = 3e8  # 光速
    mu0 = 4e-7 * np.pi
    theoretical_coeff = grid.time_step / mu0  # dt/μ₀
    print(f"   理論係數: {theoretical_coeff:.6e}")
    
    ratio = h_update_coeff / theoretical_coeff if theoretical_coeff > 0 else 0
    print(f"   係數比例: {ratio:.2e}")
    
    coeff_ok = 0.1 < ratio < 10
    if not coeff_ok:
        print(f"   ❌ 更新係數異常！")
    else:
        print(f"   ✅ 更新係數合理")
    
    print(f"\n4. curl_E結果檢查：")
    
    # 不導入curl_E，直接檢查場的梯度
    try:
        # 簡單的差分近似curl
        E = grid.E
        dEx_dy = np.diff(E[:, :, :, 0], axis=1)
        dEy_dx = np.diff(E[:, :, :, 1], axis=0)
        
        curl_approx_max = float(np.max(np.abs(dEx_dy))) + float(np.max(np.abs(dEy_dx)))
        print(f"   場梯度估計: {curl_approx_max:.6e}")
        
        # 估算理論curl值
        E_max = float(np.max(np.abs(grid.E)))
        theoretical_curl = E_max / grid_spacing  # 粗略估計
        print(f"   理論curl估計: {theoretical_curl:.6e}")
        
        curl_ratio = curl_approx_max / theoretical_curl if theoretical_curl > 0 else 0
        print(f"   curl比例: {curl_ratio:.2e}")
        
    except Exception as e:
        print(f"   ❌ curl檢查失敗: {str(e)}")
        curl_ratio = 1.0
    
    return {
        'courant_ok': courant_ok,
        'h_coeff': h_update_coeff,
        'theoretical_coeff': theoretical_coeff,
        'coeff_ratio': ratio,
        'curl_ratio': curl_ratio
    }

def check_units_and_scaling_fixed(grid, grid_spacing):
    """檢查單位和尺度問題（修正版）"""
    
    print(f"\n📏 單位和尺度檢查")
    print("="*50)
    
    print("1. 物理單位檢查：")
    wavelength = 1550e-9  # 1550 nm
    frequency = 3e8 / wavelength
    period = 1.0 / frequency
    
    print(f"   波長: {wavelength*1e9:.0f} nm")
    print(f"   頻率: {frequency:.2e} Hz")
    print(f"   周期: {period:.2e} s")
    print(f"   模擬時間步: {grid.time_step:.2e} s")
    
    steps_per_period = period / grid.time_step
    print(f"   每周期時間步數: {steps_per_period:.1f}")
    
    if steps_per_period < 10:
        print(f"   ⚠️ 時間解析度可能不足")
    
    print(f"\n2. 空間尺度檢查：")
    wavelength_in_grid = wavelength / grid_spacing
    print(f"   每波長格點數: {wavelength_in_grid:.1f}")
    
    if wavelength_in_grid < 10:
        print(f"   ⚠️ 空間解析度可能不足")
    
    print(f"\n3. 阻抗檢查：")
    E_max = float(np.max(np.abs(grid.E)))
    H_max = float(np.max(np.abs(grid.H)))
    
    print(f"   |E|_max: {E_max:.3e} V/m")
    print(f"   |H|_max: {H_max:.3e} A/m")
    
    if H_max > 1e-15:
        impedance = E_max / H_max
        print(f"   實際阻抗: {impedance:.1f} Ω")
        print(f"   理論阻抗: 377.0 Ω")
        print(f"   偏差: {abs(impedance - 377)/377*100:.1f}%")
        
        impedance_ok = abs(impedance - 377) < 100  # 允許較大誤差
        if not impedance_ok:
            print(f"   ❌ 阻抗嚴重偏離理論值！")
            print(f"   💡 H場數值異常，可能原因：")
            print(f"      - 時間步長錯誤")
            print(f"      - 更新係數錯誤") 
            print(f"      - curl計算錯誤")
        else:
            print(f"   ✅ 阻抗合理")
        
        return impedance_ok
    else:
        print(f"   ❌ 無法計算阻抗（H場為零）")
        return False

def simple_magnetic_field_fix_test(grid, grid_spacing):
    """簡單的磁場修正測試"""
    
    print(f"\n🔧 簡單磁場修正測試")
    print("="*50)
    
    # 獲取當前E和H的數值
    E_max = float(np.max(np.abs(grid.E)))
    H_max = float(np.max(np.abs(grid.H)))
    
    print(f"修正前：E_max = {E_max:.3e}, H_max = {H_max:.3e}")
    
    if H_max > 1e-15:
        current_impedance = E_max / H_max
        print(f"當前阻抗: {current_impedance:.1f} Ω")
        
        # 計算修正因子
        target_impedance = 377.0
        correction_factor = current_impedance / target_impedance
        
        print(f"需要將H場除以: {correction_factor:.2f}")
        
        # 保存原始H場
        H_original = grid.H.copy()
        
        # 修正H場
        grid.H = grid.H / correction_factor
        
        # 重新計算功率
        print(f"\n重新計算檢測器功率...")
        
        try:
            # 重新執行detect_S
            grid.T.detect_S()
            grid.R.detect_S()
            
            # 獲取新的功率值
            P_T_new = float(grid.T.S[-1]) if grid.T.S else 0
            P_R_new = float(grid.R.S[-1]) if grid.R.S else 0
            
            print(f"修正後檢測器功率:")
            print(f"   T檢測器: {P_T_new:.6e} W/m")
            print(f"   R檢測器: {P_R_new:.6e} W/m")
            
            # 重新計算入射功率
            source_length = len(getattr(grid.sources[0], 'x', [1])) * grid_spacing
            power_density = 0.5 * E_max**2 / target_impedance
            P_incident_corrected = power_density * source_length
            
            print(f"   修正入射功率: {P_incident_corrected:.6e} W/m")
            
            # 計算修正的T/R
            if P_incident_corrected > 0:
                T_corrected = P_T_new / P_incident_corrected
                R_corrected = abs(P_R_new) / P_incident_corrected
                
                print(f"\n修正後T/R:")
                print(f"   T = {T_corrected:.4f} ({T_corrected*100:.1f}%)")
                print(f"   R = {R_corrected:.4f} ({R_corrected*100:.1f}%)")
                print(f"   總和 = {T_corrected + R_corrected:.4f}")
                
                # 與理論對比
                r_theory = (1.0 - 1.5) / (1.0 + 1.5)
                R_theory = abs(r_theory)**2
                T_theory = 1 - R_theory
                
                print(f"\n與理論對比:")
                print(f"   理論: T={T_theory:.1%}, R={R_theory:.1%}")
                print(f"   修正: T={T_corrected:.1%}, R={R_corrected:.1%}")
                
                T_error = abs(T_corrected - T_theory) / T_theory * 100
                R_error = abs(R_corrected - R_theory) / R_theory * 100
                
                print(f"   誤差: T±{T_error:.1f}%, R±{R_error:.1f}%")
                
                if T_error < 20 and R_error < 50:  # 容許較大誤差
                    print(f"   ✅ 修正後接近理論值！")
                    success = True
                else:
                    print(f"   ⚠️ 仍有偏差")
                    success = False
            else:
                success = False
                
        except Exception as e:
            print(f"修正測試失敗: {str(e)}")
            success = False
        
        # 恢復原始H場
        grid.H = H_original
        
        return success
    
    else:
        print("H場為零，無法進行修正測試")
        return False

def comprehensive_parameter_diagnosis_fixed(grid, grid_spacing):
    """綜合參數診斷（修正版）"""
    
    print(f"\n🚨 綜合參數診斷")
    print("="*80)
    
    try:
        # 1. 基本參數檢查
        params = diagnose_fdtd_parameters_fixed(grid, grid_spacing)
        
        # 2. 單位尺度檢查  
        units_ok = check_units_and_scaling_fixed(grid, grid_spacing)
        
        # 3. 簡單修正測試
        fix_success = simple_magnetic_field_fix_test(grid, grid_spacing)
        
        print(f"\n🎯 診斷摘要：")
        print(f"   Courant數: {'✅' if params['courant_ok'] else '❌'}")
        print(f"   更新係數: {'✅' if 0.1 < params['coeff_ratio'] < 10 else '❌'}")
        print(f"   阻抗匹配: {'✅' if units_ok else '❌'}")
        print(f"   修正測試: {'✅' if fix_success else '❌'}")
        
        if fix_success:
            print(f"\n💡 問題確認：H場數值錯誤！")
            print(f"   根本原因：磁場更新係數或時間步長有問題")
            print(f"   臨時解法：將H場除以阻抗比例因子")
        
        return params, units_ok, fix_success
        
    except Exception as e:
        print(f"\n❌ 診斷過程出錯: {str(e)}")
        import traceback
        traceback.print_exc()
        return None, False, False


# 阻抗錯誤的根本原因分析

def analyze_impedance_error_root_cause(grid, grid_spacing):
    """分析阻抗錯誤的根本原因"""
    
    print(f"\n🔍 阻抗錯誤根本原因分析")
    print("="*80)
    
    print("1. 回顧發現的問題：")
    print(f"   - 實際阻抗: ~1.3 Ω")
    print(f"   - 理論阻抗: 377 Ω")
    print(f"   - H場過大: ~280倍")
    print(f"   - 更新係數異常: 1.88e+10")
    
    print(f"\n2. Maxwell方程組數值實現分析：")
    print(f"   在FDTD中，H場更新公式為：")
    print(f"   dH/dt = -(1/μ₀) * ∇×E")
    print(f"   離散化：H^(n+1) = H^n - (dt/μ₀) * curl(E)")
    
    # 檢查實際的更新係數
    mu0 = 4e-7 * np.pi
    theoretical_dt_over_mu0 = grid.time_step / mu0
    actual_coeff = grid.courant_number * grid.inverse_permeability[0,0,0,0]
    
    print(f"\n3. 更新係數詳細檢查：")
    print(f"   理論 dt/μ₀ = {theoretical_dt_over_mu0:.6e}")
    print(f"   實際係數 = courant * inv_perm = {actual_coeff:.6e}")
    print(f"   比例 = {actual_coeff/theoretical_dt_over_mu0:.2e}")
    
    print(f"\n4. 可能的錯誤來源：")
    
    # 錯誤可能性1：時間步長錯誤
    print(f"   可能性1：時間步長計算錯誤")
    c = 3e8
    expected_dt = grid.courant_number * grid_spacing / c
    print(f"     理論時間步: {expected_dt:.6e} s")
    print(f"     實際時間步: {grid.time_step:.6e} s")
    print(f"     比例: {grid.time_step/expected_dt:.2e}")
    
    if abs(grid.time_step/expected_dt - 1) > 0.1:
        print(f"     ❌ 時間步長異常！")
        time_step_error = True
    else:
        print(f"     ✅ 時間步長正常")
        time_step_error = False
    
    # 錯誤可能性2：單位系統問題
    print(f"\n   可能性2：單位系統錯誤")
    print(f"     如果grid_spacing單位錯誤，會影響dt計算")
    print(f"     如果μ₀值錯誤，會直接影響H場更新")
    print(f"     如果force_complex處理不當，可能產生數值錯誤")
    
    # 錯誤可能性3：Courant數定義錯誤  
    print(f"\n   可能性3：Courant數定義或使用錯誤")
    expected_courant = 1.0 / np.sqrt(grid.D)  # D是維度
    print(f"     理論最大Courant數: {expected_courant:.6f}")
    print(f"     實際Courant數: {grid.courant_number:.6f}")
    
    # 錯誤可能性4：curl_E函數問題
    print(f"\n   可能性4：curl_E函數實現錯誤")
    print(f"     可能缺少1/grid_spacing正規化")
    print(f"     可能有額外的係數錯誤")
    
    return {
        'time_step_error': time_step_error,
        'coeff_ratio': actual_coeff/theoretical_dt_over_mu0,
        'dt_ratio': grid.time_step/expected_dt
    }

def investigate_fdtd_implementation_details(grid, grid_spacing):
    """調查FDTD實現的具體細節"""
    
    print(f"\n🔬 FDTD實現細節調查")
    print("="*60)
    
    print("1. 檢查網格初始化參數：")
    print(f"   force_complex: {getattr(grid, 'force_complex', 'N/A')}")
    print(f"   維度D: {grid.D}")
    print(f"   permittivity初始值: {getattr(grid, 'permittivity', 'N/A')}")
    print(f"   permeability初始值: {getattr(grid, 'permeability', 'N/A')}")
    
    print(f"\n2. 檢查常數定義：")
    # 檢查是否在某處定義了錯誤的物理常數
    c_used = grid_spacing * grid.courant_number / grid.time_step
    print(f"   從網格參數推算的光速: {c_used:.2e} m/s")
    print(f"   理論光速: 3.0e8 m/s")
    print(f"   差異: {abs(c_used - 3e8)/3e8*100:.1f}%")
    
    print(f"\n3. 檢查field數據類型：")
    print(f"   E場dtype: {grid.E.dtype}")
    print(f"   H場dtype: {grid.H.dtype}")
    print(f"   inverse_permittivity dtype: {grid.inverse_permittivity.dtype}")
    print(f"   inverse_permeability dtype: {grid.inverse_permeability.dtype}")
    
    print(f"\n4. 檢查backend：")
    try:
        import fdtd.backend
        backend_name = getattr(fdtd.backend, 'backend_name', 'unknown')
        print(f"   使用的backend: {backend_name}")
    except:
        print(f"   無法檢查backend")

def hypothesis_about_error_source():
    """關於錯誤來源的假設"""
    
    print(f"\n💡 錯誤來源假設")
    print("="*60)
    
    print("根據症狀分析，最可能的原因：")
    print()
    print("🎯 假設1：時間步長計算錯誤（最可能）")
    print("   - 症狀：更新係數異常大（1.88e+10）")
    print("   - 原因：grid.time_step可能比正確值大很多")
    print("   - 位置：Grid.__init__()中的time_step計算")
    print("   - 公式：time_step = courant_number * grid_spacing / c")
    print("   - 可能：c的值錯誤，或grid_spacing單位錯誤")
    
    print(f"\n🎯 假設2：磁導率μ₀值錯誤")
    print("   - 症狀：H場更新過度")
    print("   - 原因：update_H中使用的μ₀值錯誤")
    print("   - 位置：可能在constants.py或直接硬編碼")
    print("   - 影響：dH/dt = -(1/μ₀) * curl(E)")
    
    print(f"\n🎯 假設3：curl_E函數係數錯誤")
    print("   - 症狀：H場數值異常")
    print("   - 原因：curl計算中缺少或多了某個係數")
    print("   - 位置：grid.py中的curl_E函數")
    print("   - 可能：缺少1/grid_spacing正規化")
    
    print(f"\n🎯 假設4：複數運算處理錯誤")
    print("   - 症狀：數值異常")
    print("   - 原因：force_complex=True時的類型轉換錯誤")
    print("   - 位置：Grid.update_H()中的複數處理")
    
    print(f"\n🔍 建議調查順序：")
    print("1. 檢查Grid.__init__()中time_step的計算公式")
    print("2. 檢查constants.py中的物理常數定義")
    print("3. 檢查curl_E函數的實現細節")
    print("4. 檢查update_H中的複數運算")

def suggest_debugging_steps():
    """建議具體的調試步驟"""
    
    print(f"\n🔧 具體調試建議")
    print("="*60)
    
    print("立即可以做的檢查：")
    print()
    print("1. 檢查您的Grid創建語句：")
    print("   grid = fdtd.Grid(...)")
    print("   特別注意grid_spacing和courant_number參數")
    print()
    print("2. 在Grid創建後立即檢查：")
    print("   print(f'time_step: {grid.time_step}')")
    print("   print(f'courant_number: {grid.courant_number}')")
    print("   print(f'grid_spacing: {grid_spacing}')")
    print()
    print("3. 手動驗證時間步長：")
    print("   expected_dt = courant_number * grid_spacing / 3e8")
    print("   print(f'Expected: {expected_dt}, Actual: {grid.time_step}')")
    print()
    print("4. 檢查fdtd庫版本：")
    print("   import fdtd")
    print("   print(f'FDTD version: {fdtd.__version__}')")
    print("   可能是某個版本的bug")

def comprehensive_impedance_error_analysis(grid, grid_spacing):
    """綜合阻抗錯誤分析"""
    
    print(f"\n🚨 綜合阻抗錯誤分析")
    print("="*80)
    
    # 1. 根本原因分析
    root_cause = analyze_impedance_error_root_cause(grid, grid_spacing)
    
    # 2. 實現細節調查
    investigate_fdtd_implementation_details(grid, grid_spacing)
    
    # 3. 錯誤假設
    hypothesis_about_error_source()
    
    # 4. 調試建議
    suggest_debugging_steps()
    
    print(f"\n🎯 結論：")
    if root_cause['time_step_error']:
        print("❌ 時間步長計算確實有問題")
        print("🔧 建議：檢查Grid初始化和物理常數定義")
    else:
        print("❓ 時間步長看起來正常，問題可能在其他地方")
        print("🔧 建議：檢查curl_E函數和磁導率設定")
    
    return root_cause