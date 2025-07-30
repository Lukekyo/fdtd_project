import numpy as np
import matplotlib.pyplot as plt

def explain_fft_power_calculation():
    """
    解釋為什麼 FFT 結果要做絕對值平方
    """
    
    print("🔍 FFT 功率計算解釋")
    print("="*50)
    
    # 創建測試信號
    dt = 0.01
    N = 1000
    time = np.arange(N) * dt
    
    # 信號: 1Hz 振幅=2, 3Hz 振幅=1
    signal = 2.0 * np.sin(2*np.pi*1*time) + 1.0 * np.sin(2*np.pi*3*time)
    
    # FFT
    fft_result = np.fft.fft(signal)
    freqs = np.fft.fftfreq(N, dt)
    
    # 找到 1Hz 和 3Hz 的位置
    idx_1hz = np.argmin(np.abs(freqs - 1.0))
    idx_3hz = np.argmin(np.abs(freqs - 3.0))
    
    print(f"📊 測試信號分析:")
    print(f"   1Hz 分量 - 理論振幅: 2.0")
    print(f"   3Hz 分量 - 理論振幅: 1.0")
    
    print(f"\n🔢 FFT 結果 (複數):")
    fft_1hz = fft_result[idx_1hz]
    fft_3hz = fft_result[idx_3hz]
    
    print(f"   1Hz FFT係數: {fft_1hz:.1f}")
    print(f"   3Hz FFT係數: {fft_3hz:.1f}")
    
    print(f"\n📏 不同的「大小」計算:")
    
    # 1. 只取實部 (錯誤!)
    real_1hz = np.real(fft_1hz)
    real_3hz = np.real(fft_3hz)
    print(f"   只取實部:")
    print(f"     1Hz: {real_1hz:.1f}")
    print(f"     3Hz: {real_3hz:.1f}")
    print(f"     ❌ 不對！會遺失虛部信息")
    
    # 2. 絕對值 (振幅)
    abs_1hz = np.abs(fft_1hz)
    abs_3hz = np.abs(fft_3hz)
    print(f"   絕對值 |FFT|:")
    print(f"     1Hz: {abs_1hz:.1f}")
    print(f"     3Hz: {abs_3hz:.1f}")
    print(f"     ✅ 這是振幅！")
    
    # 3. 絕對值平方 (功率)
    power_1hz = np.abs(fft_1hz)**2
    power_3hz = np.abs(fft_3hz)**2
    print(f"   絕對值平方 |FFT|²:")
    print(f"     1Hz: {power_1hz:.0f}")
    print(f"     3Hz: {power_3hz:.0f}")
    print(f"     ✅ 這是功率！")
    
    print(f"\n⚡ 物理意義:")
    print(f"   - 複數 FFT 係數 = 振幅 × e^(iφ)")
    print(f"   - |FFT| = 振幅")
    print(f"   - |FFT|² = 功率 ∝ 振幅²")
    
    print(f"\n🎯 為什麼用功率？")
    print(f"   1. 功率是物理可測量的量")
    print(f"   2. 功率與能量成正比")
    print(f"   3. 相位信息在功率分析中不重要")
    print(f"   4. 符合 Parseval 定理")


def complex_number_magnitude_demo():
    """
    示範複數的絕對值計算
    """
    print(f"\n📐 複數絕對值計算示範:")
    print("="*40)
    
    # 幾個複數例子
    examples = [
        3 + 4j,      # |z| = 5
        1 + 0j,      # |z| = 1  
        0 + 2j,      # |z| = 2
        -1 - 1j      # |z| = √2
    ]
    
    for z in examples:
        magnitude = np.abs(z)
        power = np.abs(z)**2
        
        print(f"   複數: {z}")
        print(f"     |z| = √({z.real}² + {z.imag}²) = {magnitude:.2f}")
        print(f"     |z|² = {power:.2f}")
        print()


def why_power_not_amplitude():
    """
    解釋為什麼用功率而不是振幅
    """
    print(f"\n🤔 為什麼用 |FFT|² 而不是 |FFT|？")
    print("="*50)
    
    print(f"📚 理論原因:")
    print(f"   1. 功率守恆: Parseval定理")
    print(f"      ∑|x[n]|² = (1/N)∑|X[k]|²")
    print(f"   2. 能量密度正比於振幅平方")
    print(f"   3. 與實驗測量一致")
    
    print(f"\n🔬 實際原因:")
    print(f"   - Detector 測量的是功率，不是振幅")
    print(f"   - 穿透率 T = P_transmitted / P_incident")
    print(f"   - 功率譜密度 = |FFT|²")
    
    print(f"\n⚖️ 對比:")
    print(f"   振幅: 告訴你「場有多強」")
    print(f"   功率: 告訴你「能量流有多大」")
    print(f"   ✅ 我們要的是能量流！")


# 執行說明
if __name__ == "__main__":
    explain_fft_power_calculation()
    complex_number_magnitude_demo()
    why_power_not_amplitude()