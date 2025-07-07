import fdtd

# 創建一個源
source = fdtd.ComplexPlaneWave(wavelength=1550e-9, period=1e-15, amplitude=1.0)

# 檢查是否有方法
print(hasattr(source, 'get_source_power'))