# Floport的period轉換機制解析

import numpy as np

def explain_floport_period_conversion():
    """
    解釋Floport在_register_grid中調整period的原因和機制
    """
    
    print("🔄 Floport的period轉換機制")
    print("=" * 60)
    
    print("📍 Floport的巧妙設計：支援兩種period輸入方式")
    print()
    print("方式1️⃣：直接輸入timesteps（整數）")
    print("   source = LineSource(period=15)  # 15個時間步")
    print("   → _handle_time(15) → 15 (不變)")
    print()
    print("方式2️⃣：輸入物理時間（浮點數）")
    print("   source = LineSource(period=5.17e-15)  # 5.17飛秒")
    print("   → _handle_time(5.17e-15) → 轉換成對應的timesteps")
    print()
    
    print("🔧 `grid._handle_time()` 函數邏輯：")
    code = '''
def _handle_time(self, time: Number) -> int:
    if not isinstance(time, int):
        # 物理時間(秒) → timesteps
        return int(float(time) / self.time_step + 0.5)
    else:
        # 已經是timesteps，直接返回
        return time
    '''
    print(code)
    
    return None

def demonstrate_period_conversion():
    """
    實際演示period轉換過程
    """
    
    print(f"\n📊 period轉換實例")
    print("=" * 60)
    
    # 模擬FDTD參數
    grid_spacing = 155e-9  # 155 nm
    courant = 0.5
    time_step = courant * grid_spacing / 3e8
    
    print(f"FDTD設定：")
    print(f"   grid_spacing = {grid_spacing*1e9:.0f} nm")
    print(f"   time_step = {time_step:.2e} s")
    print()
    
    # 測試不同的period輸入
    test_cases = [
        ("整數timesteps", 15, "timesteps"),
        ("浮點timesteps", 15.0, "timesteps"), 
        ("物理時間(秒)", 5.17e-15, "秒"),
        ("wavelength/c", 1550e-9/3e8, "秒")
    ]
    
    def handle_time_simulation(time_value, time_step):
        """模擬_handle_time函數"""
        if not isinstance(time_value, int):
            return int(float(time_value) / time_step + 0.5)
        return time_value
    
    print(f"{'輸入類型':<15} {'輸入值':<15} {'轉換後period':<15} {'最終frequency':<15}")
    print("-" * 70)
    
    for case_name, input_period, unit in test_cases:
        # 模擬_register_grid中的轉換
        converted_period = handle_time_simulation(input_period, time_step)
        final_frequency = 1.0 / converted_period
        
        if unit == "秒":
            input_str = f"{input_period:.2e} s"
        else:
            input_str = f"{input_period}"
            
        print(f"{case_name:<15} {input_str:<15} {converted_period:<15} {final_frequency:<15.6f}")
    
    return None

def explain_why_reconvert_frequency():
    """
    解釋為什麼要重新計算frequency
    """
    
    print(f"\n❓ 為什麼要重新計算frequency？")
    print("=" * 60)
    
    print("原因1️⃣：單位統一")
    print("   初始化時：frequency = 1.0 / period")
    print("   如果period是秒 → frequency單位是Hz")
    print("   如果period是timesteps → frequency單位是1/timestep")
    print("   _register_grid後：period統一成timesteps")
    print("   所以frequency也要重新計算成1/timestep")
    print()
    
    print("原因2️⃣：數值一致性")
    print("   轉換前：period=5.17e-15秒，frequency=1.94e14 Hz")
    print("   轉換後：period=33.4 timesteps，frequency=0.03 /timestep")
    print("   保證在update_E中使用時單位正確")
    print()
    
    print("原因3️⃣：避免重複計算")
    print("   如果不重新計算，每次update_E都要轉換")
    print("   預先轉換一次，後續計算更高效")
    
    return None

def compare_with_complexplanewave():
    """
    與ComplexPlaneWave比較
    """
    
    print(f"\n🆚 與你的ComplexPlaneWave比較")
    print("=" * 60)
    
    print("Floport LineSource的做法：")
    print("✅ 支援兩種period輸入方式")
    print("✅ 自動單位轉換")  
    print("✅ 統一的內部表示")
    print("❌ 沒有wavelength概念")
    print("❌ 無法直接控制光學參數")
    print()
    
    print("你的ComplexPlaneWave應該採用：")
    print("💡 方案1：學習Floport的轉換機制")
    code1 = '''
def _register_grid(self, grid, x, y, z):
    # ... 其他註冊代碼 ...
    
    # 統一period單位（學習Floport）
    self.period = grid._handle_time(self.period)
    
    # 重新計算模擬頻率
    self.simulation_frequency = 1.0 / self.period
    
    # 但保持物理頻率不變
    self.physical_frequency = 3e8 / self.wavelength
    '''
    print(code1)
    
    print("💡 方案2：明確區分兩種頻率")
    code2 = '''
def __init__(self, wavelength, period=None, ...):
    self.wavelength = wavelength
    self.physical_frequency = 3e8 / wavelength
    
    if period is None:
        # 自動計算period
        physical_period = 1.0 / self.physical_frequency
        self.period = physical_period  # 先設為秒，等註冊時轉換
    else:
        self.period = period
        
def _register_grid(self, grid, x, y, z):
    # 轉換period（如果需要）
    self.period = grid._handle_time(self.period)
    self.simulation_frequency = 1.0 / self.period
    '''
    print(code2)

def main_analysis():
    """
    主要分析函數
    """
    
    print("🔍 Floport的period調整機制完整分析")
    print("=" * 80)
    
    # 1. 解釋轉換機制
    explain_floport_period_conversion()
    
    # 2. 演示轉換過程
    demonstrate_period_conversion()
    
    # 3. 解釋重新計算的原因
    explain_why_reconvert_frequency()
    
    # 4. 與ComplexPlaneWave比較
    compare_with_complexplanewave()
    
    print(f"\n" + "="*80)
    print("🎯 你的發現總結")
    print("="*80)
    print("1. Floport設計得很巧妙：支援兩種period輸入方式")
    print("2. _register_grid是單位統一的關鍵步驟")
    print("3. period轉換後，frequency也必須重新計算")
    print("4. 這確保了update_E中的計算單位正確")
    print("5. 你的ComplexPlaneWave也應該採用類似機制")
    print("6. 關鍵是要區分物理頻率和模擬頻率！")

if __name__ == "__main__":
    main_analysis()