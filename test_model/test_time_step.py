# 物理頻率 vs 模擬頻率完整解釋

import numpy as np
import matplotlib.pyplot as plt

def explain_why_period_exists():
    """
    解釋為什麼FDTD需要period這個概念
    """
    
    print("🤔 為什麼還要period這東西？")
    print("=" * 60)
    
    print("FDTD是「離散」的數值方法，需要把「連續」的物理世界數位化：")
    print()
    print("📍 現實世界（連續）:")
    print("   - 光波頻率：f = 1.94×10¹⁴ Hz （1550nm光）")
    print("   - 光波週期：T = 5.17×10⁻¹⁵ 秒")
    print("   - 波形：E(t) = E₀ cos(2πft)")
    print()
    print("📍 FDTD世界（離散）:")
    print("   - 時間步長：Δt ≈ 10⁻¹⁶ 秒（由Courant條件決定）")
    print("   - 時間變成：t₀, t₁, t₂, t₃, ... = 0×Δt, 1×Δt, 2×Δt, 3×Δt, ...")
    print("   - 波形變成：E[0], E[1], E[2], E[3], ... （數字序列）")
    print()
    print("❓ 問題：如何把連續的5.17×10⁻¹⁵秒週期，轉換成離散的「第幾步」？")
    print()
    print("💡 答案：period = T / Δt = 一個物理週期需要多少個時間步")
    print()
    
    # 實際計算示例
    wavelength = 1550e-9  # 1550 nm
    c = 3e8
    f_physical = c / wavelength  # Hz
    T_physical = 1.0 / f_physical  # 秒
    
    # FDTD參數
    grid_spacing = 155e-9  # 155 nm
    courant = 0.5
    dt = courant * grid_spacing / c  # 時間步長
    
    period_timesteps = T_physical / dt  # 這就是period！
    
    print(f"📊 實際數字:")
    print(f"   物理週期 T = {T_physical:.2e} 秒")
    print(f"   時間步長 Δt = {dt:.2e} 秒")
    print(f"   period = T/Δt = {period_timesteps:.1f} 時間步")
    print()
    print(f"🎯 意義：每{period_timesteps:.1f}個時間步，光波完成一個完整週期")
    
    return {
        'T_physical': T_physical,
        'dt': dt,
        'period_timesteps': period_timesteps,
        'f_physical': f_physical
    }

def explain_frequency_types():
    """
    詳細解釋物理頻率和模擬頻率
    """
    
    print(f"\n🔄 物理頻率 vs 模擬頻率")
    print("=" * 60)
    
    print("有兩套完全不同的「頻率」概念：")
    print()
    
    print("📍 物理頻率（Physical Frequency）：")
    print("   - 定義：光波在現實世界中每秒振盪幾次")
    print("   - 公式：f = c / λ")
    print("   - 單位：Hz = 1/秒")
    print("   - 例子：1550nm光的頻率 = 1.94×10¹⁴ Hz")
    print("   - 用途：計算物理公式，如 ω = 2πf, E = hf")
    print()
    
    print("📍 模擬頻率（Simulation Frequency）：")
    print("   - 定義：在FDTD模擬中，每個時間步振盪多少")
    print("   - 公式：f_sim = 1 / period")
    print("   - 單位：1/timestep（無量綱）")
    print("   - 例子：period=50時，f_sim = 1/50 = 0.02 /timestep")
    print("   - 用途：FDTD數值計算，控制波形重複")
    print()
    
    print("🔗 關係：")
    print("   物理頻率 × 時間步長 = 模擬頻率 × 物理週期")
    print("   f_physical × dt = f_sim × T_physical")
    print("   或者：f_sim = f_physical × dt")
    
    return None

def demonstrate_frequency_conversion():
    """
    演示頻率轉換的具體過程
    """
    
    print(f"\n🔄 頻率轉換具體過程")
    print("=" * 60)
    
    # 設定參數
    wavelength = 1550e-9  # nm
    grid_spacing = 155e-9  # nm 
    courant = 0.5
    c = 3e8
    
    # 步驟1：計算物理參數
    print("步驟1️⃣：從波長計算物理參數")
    f_physical = c / wavelength
    T_physical = 1.0 / f_physical
    omega_physical = 2 * np.pi * f_physical
    
    print(f"   波長 λ = {wavelength*1e9:.0f} nm")
    print(f"   物理頻率 f = c/λ = {f_physical:.2e} Hz")
    print(f"   物理週期 T = 1/f = {T_physical:.2e} 秒")
    print(f"   物理角頻率 ω = 2πf = {omega_physical:.2e} rad/s")
    
    # 步驟2：計算FDTD參數
    print(f"\n步驟2️⃣：計算FDTD時間步長")
    dt = courant * grid_spacing / c
    print(f"   網格間距 Δx = {grid_spacing*1e9:.0f} nm")
    print(f"   Courant數 = {courant}")
    print(f"   時間步長 Δt = {dt:.2e} 秒")
    
    # 步驟3：轉換到模擬單位
    print(f"\n步驟3️⃣：轉換到模擬單位")
    period_timesteps = T_physical / dt
    f_simulation = 1.0 / period_timesteps
    omega_simulation = 2 * np.pi / period_timesteps
    
    print(f"   period = T/Δt = {period_timesteps:.1f} timesteps")
    print(f"   模擬頻率 f_sim = 1/period = {f_simulation:.6f} /timestep")
    print(f"   模擬角頻率 ω_sim = 2π/period = {omega_simulation:.6f} rad/timestep")
    
    # 步驟4：驗證轉換
    print(f"\n步驟4️⃣：驗證轉換正確性")
    check1 = f_physical * dt
    check2 = f_simulation
    print(f"   f_physical × Δt = {check1:.6f}")
    print(f"   f_simulation = {check2:.6f}")
    print(f"   差異 = {abs(check1 - check2):.2e} {'✅ 正確' if abs(check1-check2)<1e-10 else '❌ 錯誤'}")
    
    # 步驟5：時間演化比較
    print(f"\n步驟5️⃣：時間演化公式比較")
    print(f"   物理公式：E(t) = E₀ exp(i × {omega_physical:.2e} × t)")
    print(f"   模擬公式：E[n] = E₀ exp(i × {omega_simulation:.6f} × n)")
    print(f"   其中：t = n × Δt")
    
    return {
        'f_physical': f_physical,
        'f_simulation': f_simulation,
        'period_timesteps': period_timesteps,
        'dt': dt,
        'omega_physical': omega_physical,
        'omega_simulation': omega_simulation
    }

def visualize_frequency_concepts():
    """
    視覺化展示兩種頻率概念
    """
    
    print(f"\n📊 視覺化展示")
    print("=" * 60)
    
    # 參數設定
    wavelength = 1550e-9
    c = 3e8
    grid_spacing = 155e-9
    courant = 0.5
    
    # 計算參數
    f_physical = c / wavelength
    dt = courant * grid_spacing / c
    period_timesteps = (1.0 / f_physical) / dt
    omega_simulation = 2 * np.pi / period_timesteps
    
    # 建立時間軸
    max_timesteps = int(2 * period_timesteps)  # 兩個週期
    timesteps = np.arange(max_timesteps)
    times = timesteps * dt
    
    # 計算波形
    E_continuous = np.cos(2 * np.pi * f_physical * times)  # 連續時間公式
    E_discrete = np.cos(omega_simulation * timesteps)      # 離散時間公式
    
    print(f"生成了{max_timesteps}個時間步的波形數據")
    print(f"物理時間範圍：0 到 {times[-1]:.2e} 秒")
    print(f"包含 {max_timesteps/period_timesteps:.1f} 個完整週期")
    
    # 檢查兩種計算的一致性
    difference = np.max(np.abs(E_continuous - E_discrete))
    print(f"兩種計算方法的最大差異：{difference:.2e}")
    print(f"一致性檢查：{'✅ 通過' if difference < 1e-10 else '❌ 失敗'}")
    
    return {
        'timesteps': timesteps,
        'times': times,
        'E_continuous': E_continuous,
        'E_discrete': E_discrete,
        'period_timesteps': period_timesteps
    }

def explain_practical_usage():
    """
    解釋實際使用中的考量
    """
    
    print(f"\n🛠️ 實際使用考量")
    print("=" * 60)
    
    print("在ComplexPlaneWave中，你需要兩種頻率：")
    print()
    
    print("1️⃣ 物理頻率 - 用於物理公式：")
    print("   - 計算 ω = 2πf （角頻率）")
    print("   - 計算 k = ω/c = 2π/λ （波數）")
    print("   - 計算 hanning_dt = 0.5/f （脈衝寬度）")
    print("   - 與wavelength保持一致性")
    print()
    
    print("2️⃣ 模擬頻率 - 用於數值計算：")
    print("   - 計算 time_phase = 2π × timestep / period")
    print("   - 控制波形在時間軸上的重複")
    print("   - 確保數值穩定性")
    print()
    
    print("🎯 完整的ComplexPlaneWave設定流程：")
    
    code_example = '''
# 第一步：從wavelength開始
wavelength = 1550e-9  # 物理波長

# 第二步：計算物理頻率
f_physical = 3e8 / wavelength  # Hz
omega_physical = 2 * np.pi * f_physical  # rad/s

# 第三步：等Grid創建後，轉換到模擬單位
# （在_register_grid中）
dt = grid.time_step
period_timesteps = (1.0 / f_physical) / dt  # timesteps
f_simulation = 1.0 / period_timesteps  # 1/timestep

# 第四步：在update_E中使用
def update_E(self):
    q = self.grid.time_steps_passed
    time_phase = 2 * np.pi * q / self.period  # 使用period_timesteps
    # 或者 time_phase = self.omega * q * self.grid.time_step
    '''
    
    print(code_example)
    
    print("❗ 關鍵理解：")
    print("   - period不是多餘的，它是物理世界到數位世界的橋樑")
    print("   - 沒有period，你無法知道「第100個timestep對應多少物理時間」")
    print("   - period = 物理週期 ÷ 時間步長，單位是timesteps")

def main_explanation():
    """
    主要解釋函數
    """
    
    print("🌊 完整解釋：為什麼需要period？物理頻率vs模擬頻率")
    print("=" * 80)
    
    # 1. 解釋period存在的必要性
    params1 = explain_why_period_exists()
    
    # 2. 解釋兩種頻率概念
    explain_frequency_types()
    
    # 3. 演示轉換過程
    params2 = demonstrate_frequency_conversion()
    
    # 4. 視覺化
    vis_data = visualize_frequency_concepts()
    
    # 5. 實際使用
    explain_practical_usage()
    
    print(f"\n" + "="*80)
    print("🎯 總結")
    print("="*80)
    print("1. period是必須的：它是物理時間到模擬時間的轉換因子")
    print("2. 物理頻率：用於物理公式，單位Hz")
    print("3. 模擬頻率：用於FDTD計算，單位1/timestep") 
    print("4. 轉換關係：f_sim = f_physical × dt")
    print("5. period = 物理週期 / 時間步長 = 1/(f_sim)")
    print("6. 兩套頻率都需要，各有用途，缺一不可！")

if __name__ == "__main__":
    main_explanation()