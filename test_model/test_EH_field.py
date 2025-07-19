import numpy as np
import matplotlib.pyplot as plt

# 基本參數設定
c = 3e8  # 光速 (m/s)
mu0 = 4 * np.pi * 1e-7  # 真空磁導率 (H/m)
eps0 = 1 / (mu0 * c**2)  # 真空介電常數 (F/m)
Z0 = np.sqrt(mu0 / eps0)  # 真空阻抗 (Ω) ≈ 377

# 網格設定
nx = 200  # 網格點數
grid_spacing = 1e-3  # 網格間距 (m)
courant_number = 0.5  # Courant 數
time_step = courant_number * grid_spacing / c  # 時間步長

print(f"網格參數:")
print(f"網格間距: {grid_spacing*1e3:.1f} mm")
print(f"Courant 數: {courant_number}")
print(f"時間步長: {time_step*1e12:.2f} ps")
print(f"真空阻抗: {Z0:.1f} Ω")
print()

# 初始化電磁場
Ex = np.zeros(nx)
Hy = np.zeros(nx)

# 儲存結果
time_steps = 100
results = []

# 時間迴圈
for n in range(time_steps):
    # 在中心點加入高斯脈衝源
    source_position = nx // 2
    if n < 40:  # 前40個時間步加入源
        pulse = np.exp(-((n - 20) / 10)**2)
        Ex[source_position] = pulse
    
    # 更新 Ex (除了最後一點)
    for i in range(nx - 1):
        Ex[i] = Ex[i] + (time_step / eps0) * (Hy[i] - Hy[i + 1]) / grid_spacing
    
    # 更新 Hy (除了第一點)
    for i in range(1, nx):
        Hy[i] = Hy[i] + (time_step / mu0) * (Ex[i - 1] - Ex[i]) / grid_spacing
    
    # 記錄幾個特定時間步的結果
    if n % 10 == 0:
        # 找到場強最大的位置
        max_E_idx = np.argmax(np.abs(Ex))
        max_H_idx = np.argmax(np.abs(Hy))
        
        # 計算阻抗 (避免除零)
        E_val = Ex[max_E_idx]
        H_val = Hy[max_E_idx] if abs(Hy[max_E_idx]) > 1e-12 else 1e-12
        impedance = E_val / H_val if abs(H_val) > 1e-12 else 0
        
        results.append({
            'time_step': n,
            'max_E': E_val,
            'max_H': H_val,
            'impedance': impedance,
            'E_position': max_E_idx,
            'H_position': max_H_idx
        })
        
        print(f"時間步 {n:3d}: E_max = {E_val:8.3e} V/m, H_max = {H_val:8.3e} A/m, E/H = {impedance:8.1f} Ω")

print(f"\n理論真空阻抗: {Z0:.1f} Ω")
print(f"預期的 E/H 比值應該接近 {Z0:.1f} Ω")

# 繪製最後時間步的場分佈
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(Ex, 'b-', label='Ex (V/m)')
plt.ylabel('Ex (V/m)')
plt.title('電磁場分佈')
plt.grid(True)
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(Hy, 'r-', label='Hy (A/m)')
plt.ylabel('Hy (A/m)')
plt.grid(True)
plt.legend()

plt.subplot(3, 1, 3)
# 計算阻抗分佈 (避免除零)
impedance_dist = np.where(np.abs(Hy) > 1e-12, Ex / Hy, 0)
plt.plot(impedance_dist, 'g-', label='E/H (Ω)')
plt.axhline(y=Z0, color='k', linestyle='--', label=f'理論值 {Z0:.0f} Ω')
plt.ylabel('阻抗 (Ω)')
plt.xlabel('網格位置')
plt.grid(True)
plt.legend()
plt.ylim(-500, 500)  # 限制 y 軸範圍避免奇異值

plt.tight_layout()
plt.show()

# 總結分析
print(f"\n總結分析:")
print(f"1. 網格設定: dx = {grid_spacing*1e3:.1f} mm, dt = {time_step*1e12:.2f} ps")
print(f"2. Courant 數: {courant_number} (< 1, 穩定)")
print(f"3. 理論阻抗: {Z0:.1f} Ω")
print(f"4. 如果數值正確，E/H 比值應該接近 {Z0:.1f} Ω")