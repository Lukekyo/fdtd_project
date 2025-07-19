import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos

# 解決中文顯示問題
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

# 嘗試設置中文字體，順序：最常見 -> 最兼容
def setup_chinese_font():
    """設置中文字體，自動檢測可用字體"""
    font_candidates = [
        'SimHei',           # Windows 黑體
        'Microsoft YaHei',  # Windows 微軟雅黑
        'Arial Unicode MS', # macOS
        'PingFang SC',      # macOS 蘋方
        'WenQuanYi Micro Hei', # Linux 文泉驛微米黑
        'Noto Sans CJK SC', # Google Noto 字體
        'DejaVu Sans'       # 後備字體
    ]
    
    # 檢查可用字體
    from matplotlib.font_manager import FontProperties
    import matplotlib.font_manager as fm
    
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    
    for font in font_candidates:
        if font in available_fonts:
            plt.rcParams['font.sans-serif'] = [font]
            plt.rcParams['axes.unicode_minus'] = False
            return True, font
    
    # 如果沒有中文字體，使用默認字體
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    return False, 'DejaVu Sans'

# 設置字體
try:
    USE_CHINESE, current_font = setup_chinese_font()
    if not USE_CHINESE:
        print(f"⚠️  未檢測到中文字體，使用 {current_font}")
except:
    USE_CHINESE = False
    current_font = 'default'

# 標籤字典
LABELS = {
    'sudden_signal': '突然開關信號' if USE_CHINESE else 'Sudden Switch Signal',
    'hanning_signal': '漢寧窗調制信號' if USE_CHINESE else 'Hanning Windowed Signal', 
    'envelope': '漢寧窗包絡' if USE_CHINESE else 'Hanning Window Envelope',
    'carrier': '載波信號' if USE_CHINESE else 'Carrier Signal',
    'time_steps': '時間步數 (q)' if USE_CHINESE else 'Time Steps (q)',
    'amplitude': '振幅' if USE_CHINESE else 'Amplitude',
    'signal_comparison': '信號對比：突變 vs 漢寧窗調制' if USE_CHINESE else 'Signal Comparison: Sudden vs Hanning Windowed',
    'boundary_zoom': '邊界區域放大：清楚看出平滑效果' if USE_CHINESE else 'Boundary Region Zoom: Clear Smoothing Effect',
    'cw': '連續波 (CW)' if USE_CHINESE else 'Continuous Wave (CW)',
    'pulse': '脈衝 (5週期漢寧窗)' if USE_CHINESE else 'Pulse (5-cycle Hanning)',
    'fdtd_application': 'FDTD中的實際應用：連續波 vs 脈衝' if USE_CHINESE else 'FDTD Application: CW vs Pulse',
    'frequency': '頻率' if USE_CHINESE else 'Frequency',
    'power_density': '功率密度' if USE_CHINESE else 'Power Density',
    'freq_comparison': '頻域對比：漢寧窗減少了不需要的頻率成分' if USE_CHINESE else 'Frequency Domain: Hanning Reduces Unwanted Components',
    'sudden_spectrum': '突變信號頻譜' if USE_CHINESE else 'Sudden Signal Spectrum',
    'hanning_spectrum': '漢寧窗信號頻譜' if USE_CHINESE else 'Hanning Signal Spectrum'
}

# ========================
# 核心漢寧窗函數
# ========================

def hanning_window(t, duration):
    """
    計算漢寧窗函數值
    
    Args:
        t: 當前時間
        duration: 總持續時間
    
    Returns:
        漢寧窗函數值 (0到1之間)
    """
    if t < 0 or t > duration:
        return 0.0
    
    # 漢寧窗公式：0.5 * (1 - cos(2π * t / T))
    normalized_t = t / duration
    return 0.5 * (1 - cos(2 * pi * normalized_t))

def hanning_modulated_signal(frequency, t, cycles=5, amplitude=1.0, phase=0.0):
    """
    生成漢寧窗調制的信號
    
    Args:
        frequency: 載波頻率
        t: 時間點
        cycles: 總週期數
        amplitude: 振幅
        phase: 相位偏移
        
    Returns:
        調制後的信號值
    """
    period = 1.0 / frequency
    total_duration = cycles * period
    
    if t < 0 or t > total_duration:
        return 0.0
    
    # 載波信號
    carrier = amplitude * sin(2 * pi * frequency * t + phase)
    
    # 漢寧窗包絡
    envelope = hanning_window(t, total_duration)
    
    return carrier * envelope

# ========================
# FDTD 源類實現
# ========================

class PointSource:
    """FDTD中的點源實現"""
    
    def __init__(self, period=15, amplitude=1.0, phase_shift=0.0, 
                 pulse=False, cycle=5, hanning_dt=None):
        """
        初始化點源
        
        Args:
            period: 週期（時間步數或秒）
            amplitude: 振幅
            phase_shift: 相位偏移
            pulse: 是否使用脈衝（漢寧窗）
            cycle: 脈衝週期數
            hanning_dt: 漢寧窗時間步長
        """
        self.period = period
        self.amplitude = amplitude
        self.phase_shift = phase_shift
        self.pulse = pulse
        self.cycle = cycle
        self.frequency = 1.0 / period
        self.hanning_dt = hanning_dt if hanning_dt is not None else 0.5 / self.frequency
        
        # 模擬網格時間
        self.time_steps_passed = 0
    
    def get_signal_value(self, q=None):
        """
        獲取當前時間步的信號值
        
        Args:
            q: 時間步數（可選，默認使用內部計數器）
        
        Returns:
            信號值
        """
        if q is None:
            q = self.time_steps_passed
        
        if self.pulse:
            # 脈衝模式：使用漢寧窗
            t1 = int(2 * pi / (self.frequency * self.hanning_dt / self.cycle))
            if q < t1:
                # 使用漢寧窗調制
                time_normalized = q * self.hanning_dt
                signal = hanning_modulated_signal(
                    self.frequency, time_normalized, self.cycle, self.amplitude
                )
                return signal
            else:
                return 0.0
        else:
            # 連續波模式
            return self.amplitude * sin(2 * pi * q / self.period + self.phase_shift)
    
    def step(self):
        """推進一個時間步"""
        self.time_steps_passed += 1
        return self.get_signal_value()

# ========================
# 信號比較和可視化函數
# ========================

def generate_signal_comparison():
    """生成信號對比數據"""
    
    # 參數設置
    total_points = 300
    pulse_start = 50
    pulse_end = 250
    frequency = 0.05  # 載波頻率
    
    time_points = np.arange(total_points + 1)
    
    # 初始化信號數組
    sudden_signal = np.zeros(len(time_points))
    hanning_signal = np.zeros(len(time_points))
    envelope = np.zeros(len(time_points))
    carrier = np.sin(2 * np.pi * frequency * time_points)
    
    # 生成信號
    for i, t in enumerate(time_points):
        # 突然開關信號
        if pulse_start <= t <= pulse_end:
            sudden_signal[i] = carrier[i]
        
        # 漢寧窗包絡
        if pulse_start <= t <= pulse_end:
            normalized = (t - pulse_start) / (pulse_end - pulse_start)
            envelope[i] = 0.5 * (1 - np.cos(2 * np.pi * normalized))
        
        # 漢寧窗調制信號
        hanning_signal[i] = sudden_signal[i] * envelope[i]
    
    return {
        'time': time_points,
        'sudden': sudden_signal,
        'hanning': hanning_signal,
        'envelope': envelope,
        'carrier': carrier
    }

def plot_signal_comparison(data):
    """繪製信號對比圖"""
    
    plt.figure(figsize=(12, 8))
    
    # 主要對比圖
    plt.subplot(2, 1, 1)
    plt.plot(data['time'], data['sudden'], 'r-', linewidth=2, label=LABELS['sudden_signal'])
    plt.plot(data['time'], data['hanning'], 'b-', linewidth=2, label=LABELS['hanning_signal'])
    plt.plot(data['time'], data['envelope'], 'g-', linewidth=2, label=LABELS['envelope'])
    plt.plot(data['time'], data['carrier'], 'orange', linewidth=1, alpha=0.7, label=LABELS['carrier'])
    
    plt.xlabel(LABELS['time_steps'])
    plt.ylabel(LABELS['amplitude'])
    plt.title(LABELS['signal_comparison'])
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 放大邊界區域
    plt.subplot(2, 1, 2)
    start_zoom = 40
    end_zoom = 70
    mask = (data['time'] >= start_zoom) & (data['time'] <= end_zoom)
    
    plt.plot(data['time'][mask], data['sudden'][mask], 'r-', linewidth=2, label=LABELS['sudden_signal'])
    plt.plot(data['time'][mask], data['hanning'][mask], 'b-', linewidth=2, label=LABELS['hanning_signal'])
    plt.plot(data['time'][mask], data['envelope'][mask], 'g-', linewidth=2, label=LABELS['envelope'])
    
    plt.xlabel(LABELS['time_steps'])
    plt.ylabel(LABELS['amplitude'])
    plt.title(LABELS['boundary_zoom'])
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return plt

def simulate_fdtd_sources():
    """模擬FDTD中的不同源類型"""
    
    # 創建不同類型的源
    cw_source = PointSource(period=25, amplitude=1.0, pulse=False)
    pulse_source = PointSource(period=25, amplitude=1.0, pulse=True, cycle=5)
    
    # 模擬200個時間步
    time_steps = 200
    cw_values = []
    pulse_values = []
    
    for step in range(time_steps):
        cw_values.append(cw_source.step())
        pulse_values.append(pulse_source.step())
    
    return {
        'time': np.arange(time_steps),
        'cw': np.array(cw_values),
        'pulse': np.array(pulse_values)
    }

def plot_fdtd_sources(data):
    """繪製FDTD源對比"""
    
    plt.figure(figsize=(12, 6))
    
    plt.plot(data['time'], data['cw'], 'g-', linewidth=1, label=LABELS['cw'])
    plt.plot(data['time'], data['pulse'], 'purple', linewidth=2, label=LABELS['pulse'])
    
    plt.xlabel(LABELS['time_steps'])
    plt.ylabel(LABELS['amplitude'])
    plt.title(LABELS['fdtd_application'])
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    return plt

# ========================
# 頻域分析（模擬）
# ========================

def simulate_frequency_domain():
    """模擬頻域效果"""
    
    frequencies = np.arange(0, 101)
    center_freq = 20
    
    # 模擬突變信號頻譜（有旁瓣）
    sudden_spectrum = np.exp(-0.1 * np.abs(frequencies - center_freq))
    
    # 添加旁瓣效果
    for i, f in enumerate(frequencies):
        if f != center_freq:
            sudden_spectrum[i] += 0.3 * np.exp(-0.05 * np.abs(f - center_freq)) * \
                                 np.abs(np.sin(0.3 * (f - center_freq)))
    
    # 模擬漢寧窗信號頻譜（更集中）
    hanning_spectrum = np.exp(-0.2 * np.abs(frequencies - center_freq))
    
    # 較小的旁瓣
    for i, f in enumerate(frequencies):
        if f != center_freq:
            hanning_spectrum[i] += 0.1 * np.exp(-0.1 * np.abs(f - center_freq)) * \
                                  np.abs(np.sin(0.5 * (f - center_freq)))
    
    return {
        'frequencies': frequencies,
        'sudden': sudden_spectrum,
        'hanning': hanning_spectrum
    }

def plot_frequency_comparison(data):
    """繪製頻域對比"""
    
    plt.figure(figsize=(10, 6))
    
    plt.fill_between(data['frequencies'], data['sudden'], alpha=0.3, color='red', label=LABELS['sudden_spectrum'])
    plt.fill_between(data['frequencies'], data['hanning'], alpha=0.3, color='blue', label=LABELS['hanning_spectrum'])
    plt.plot(data['frequencies'], data['sudden'], 'r-', linewidth=2)
    plt.plot(data['frequencies'], data['hanning'], 'b-', linewidth=2)
    
    plt.xlabel(LABELS['frequency'])
    plt.ylabel(LABELS['power_density'])
    plt.title(LABELS['freq_comparison'])
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    return plt

# ========================
# 主函數和使用示例
# ========================

def main():
    """主函數：演示所有功能"""
    
    print("🔧 FDTD Hanning Window Signal Python Implementation")
    print("="*60)
    
    # 檢查中文支持
    print(f"📝 當前字體: {current_font}")
    if USE_CHINESE:
        print("✅ 中文顯示支持已啟用")
    else:
        print("⚠️  中文顯示不支持，使用英文標籤")
        print("💡 解決方案:")
        print("   1. 安裝中文字體包")
        print("   2. Windows: 通常已有中文字體")
        print("   3. Linux: sudo apt install fonts-wqy-microhei")
        print("   4. macOS: 通常已有中文字體")
        print("   5. 或在 Jupyter 中運行: !pip install matplotlib --upgrade")
    
    # 1. 基本漢寧窗測試
    print("\n1. 測試漢寧窗函數:")
    for t in [0, 0.25, 0.5, 0.75, 1.0]:
        value = hanning_window(t, 1.0)
        print(f"   t={t:.2f}: hanning={value:.4f}")
    
    # 2. 信號對比
    print("\n2. 生成信號對比圖...")
    signal_data = generate_signal_comparison()
    plot_signal_comparison(signal_data)
    plt.show()
    
    # 3. FDTD源模擬
    print("\n3. 模擬FDTD源...")
    fdtd_data = simulate_fdtd_sources()
    plot_fdtd_sources(fdtd_data)
    plt.show()
    
    # 4. 頻域分析
    print("\n4. 頻域分析...")
    freq_data = simulate_frequency_domain()
    plot_frequency_comparison(freq_data)
    plt.show()
    
    # 5. 實際使用示例
    print("\n5. 實際使用示例:")
    print("\n   # 創建連續波源")
    print("   cw_source = PointSource(period=25, pulse=False)")
    print("   signal_value = cw_source.step()")
    
    print("\n   # 創建脈衝源")
    print("   pulse_source = PointSource(period=25, pulse=True, cycle=5)")
    print("   signal_value = pulse_source.step()")
    
    cw_source = PointSource(period=25, pulse=False)
    pulse_source = PointSource(period=25, pulse=True, cycle=5)
    
    print(f"\n   連續波第10步信號值: {cw_source.get_signal_value(10):.4f}")
    print(f"   脈衝第10步信號值: {pulse_source.get_signal_value(10):.4f}")
    
    print("\n🎯 關鍵理解要點:")
    print("   ✅ 漢寧窗提供平滑的開關效果")
    print("   ✅ 減少頻域中的旁瓣")
    print("   ✅ 提高FDTD數值穩定性")
    print("   ✅ 更符合物理實際")

# 額外的字體安裝函數
def install_chinese_fonts():
    """提供安裝中文字體的指導"""
    print("\n🔧 中文字體安裝指南:")
    print("="*40)
    print("🖥️  Windows:")
    print("   通常已預裝中文字體，如果仍有問題:")
    print("   1. 控制面板 -> 字體 -> 檢查是否有 SimHei 或 Microsoft YaHei")
    print("   2. 重啟 Python kernel")
    
    print("\n🍎 macOS:")
    print("   1. 已預裝 Arial Unicode MS")
    print("   2. 如有問題，可安裝 PingFang SC")
    
    print("\n🐧 Linux (Ubuntu/Debian):")
    print("   sudo apt update")
    print("   sudo apt install fonts-wqy-microhei fonts-noto-cjk")
    print("   fc-cache -fv  # 刷新字體緩存")
    
    print("\n📦 通用解決方案 (pip):")
    print("   pip install matplotlib --upgrade --force-reinstall")
    
    print("\n🔍 調試步驟:")
    print("   1. 檢查可用字體:")
    print("      import matplotlib.font_manager as fm")
    print("      fonts = [f.name for f in fm.fontManager.ttflist]")
    print("      chinese_fonts = [f for f in fonts if any(c in f for c in ['黑體', 'YaHei', 'SimHei'])]")
    print("      print(chinese_fonts)")
    
    print("\n   2. 手動設置字體:")
    print("      plt.rcParams['font.family'] = ['你的字體名稱']")
    
    cw_source = PointSource(period=25, pulse=False)
    pulse_source = PointSource(period=25, pulse=True, cycle=5)
    
    print(f"\n   連續波第10步信號值: {cw_source.get_signal_value(10):.4f}")
    print(f"   脈衝第10步信號值: {pulse_source.get_signal_value(10):.4f}")
    
    print("\n🎯 關鍵理解要點:")
    print("   ✅ 漢寧窗提供平滑的開關效果")
    print("   ✅ 減少頻域中的旁瓣")
    print("   ✅ 提高FDTD數值穩定性")
    print("   ✅ 更符合物理實際")

if __name__ == "__main__":
    main()