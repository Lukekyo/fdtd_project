""" Sources are objects that inject the fields into the grid.

Available sources:

- PointSource
- LineSource

"""
## Imports

# other
from math import pi, sin

# typing
from .typing_ import Tuple, Number, ListOrSlice, List
from numpy import ndarray

# relatvie
from .grid import Grid
from .backend import backend as bd
from .waveforms import *
from .detectors import CurrentDetector

## PointSource class
class PointSource:
    """A source placed at a single point (grid cell) in the grid"""

    def __init__(
        self,
        period: Number = 15,
        amplitude: float = 1.0,
        phase_shift: float = 0.0,
        name: str = None,
        pulse: bool = False,
        cycle: int = 5,
        hanning_dt: float = 10.0,
    ):
        """Create a LineSource with a gaussian profile

        Args:
            period: The period of the source. The period can be specified
                as integer [timesteps] or as float [seconds]
            amplitude: The electric field amplitude in simulation units
            phase_shift: The phase offset of the source.
            name: name of the source.
            pulse: Set True to use a Hanning window pulse instead of continuous wavefunction.
            cycle: cycles for Hanning window pulse.
            hanning_dt: timestep used for Hanning window pulse width (optional).

        """
        self.grid = None
        self.period = period
        self.amplitude = amplitude
        self.phase_shift = phase_shift
        self.name = name
        self.pulse = pulse
        self.cycle = cycle
        self.frequency = 1.0 / period
        self.hanning_dt = hanning_dt if hanning_dt is not None else 0.5 / self.frequency

    def _register_grid(self, grid: Grid, x: Number, y: Number, z: Number):
        """Register a grid for the source.

        Args:
            grid: the grid to place the source into.
            x: The x-location of the source in the grid
            y: The y-location of the source in the grid
            z: The z-location of the source in the grid

        Note:
            As its name suggests, this source is a POINT source.
            Hence it should be placed at a single coordinate tuple
            int the grid.
        """
        self.grid = grid
        self.grid.sources.append(self)
        if self.name is not None:
            if not hasattr(grid, self.name):
                setattr(grid, self.name, self)
            else:
                raise ValueError(
                    f"The grid already has an attribute with name {self.name}"
                )

        try:
            (x,), (y,), (z,) = x, y, z
        except (TypeError, ValueError):
            raise ValueError("a point source should be placed on a single grid cell.")
        self.x, self.y, self.z = grid._handle_tuple((x, y, z))
        self.period = grid._handle_time(self.period)
        self.frequency = 1.0 / self.period

    def update_E(self):
        """Add the source to the electric field"""
        q = self.grid.time_steps_passed
        # if pulse
        if self.pulse:
            t1 = int(2 * pi / (self.frequency * self.hanning_dt / self.cycle))
            if q < t1:
                src = self.amplitude * hanning(
                    self.frequency, q * self.hanning_dt, self.cycle
                )
            else:
                # src = - self.grid.E[self.x, self.y, self.z, 2]
                src = 0
        # if not pulse
        else:
            src = self.amplitude * sin(2 * pi * q / self.period + self.phase_shift)
        self.grid.E[self.x, self.y, self.z, 2] += src

    def update_H(self):
        """Add the source to the magnetic field"""

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(period={self.period}, "
            f"amplitude={self.amplitude}, phase_shift={self.phase_shift}, "
            f"name={repr(self.name)})"
        )

    def __str__(self):
        s = "    " + repr(self) + "\n"
        x = f"{self.x}"
        y = f"{self.y}"
        z = f"{self.z}"
        s += f"        @ x={x}, y={y}, z={z}\n"
        return s


## LineSource class
class LineSource:
    """A source along a line in the FDTD grid"""

    def __init__(
        self,
        period: Number = 15,
        amplitude: float = 1.0,
        phase_shift: float = 0.0,
        name: str = None,
        pulse: bool = False,
        cycle: int = 5,
        hanning_dt: float = 10.0,
    ):
        """Create a LineSource with a gaussian profile

        Args:
            period: The period of the source. The period can be specified
                as integer [timesteps] or as float [seconds]
            amplitude: The amplitude of the source in simulation units
            phase_shift: The phase offset of the source.
            pulse: Set True to use a Hanning window pulse instead of continuous wavefunction.
            cycle: cycles for Hanning window pulse.
            hanning_dt: timestep used for Hanning window pulse width (optional).

        """
        self.grid = None
        self.period = period
        self.amplitude = amplitude
        self.phase_shift = phase_shift
        self.name = name
        self.pulse = pulse
        self.cycle = cycle
        self.frequency = 1.0 / period
        self.hanning_dt = hanning_dt if hanning_dt is not None else 0.5 / self.frequency

    def _register_grid(
        self, grid: Grid, x: ListOrSlice, y: ListOrSlice, z: ListOrSlice
    ):
        """Register a grid for the source.

        Args:
            grid: the grid to place the source into.
            x: The x-location of the source in the grid
            y: The y-location of the source in the grid
            z: The z-location of the source in the grid

        Note:
            As its name suggests, this source is a LINE source.
            Hence the source spans the diagonal of the cube
            defined by the slices in the grid.
        """
        self.grid = grid
        self.grid.sources.append(self)
        if self.name is not None:
            if not hasattr(grid, self.name):
                setattr(grid, self.name, self)
            else:
                raise ValueError(
                    f"The grid already has an attribute with name {self.name}"
                )

        self.x, self.y, self.z = self._handle_slices(x, y, z) # convert slices to lists

        self.period = grid._handle_time(self.period) # convert to time
        self.frequency = 1.0 / self.period # convert to frequency

        L = len(self.x) # length of the line
        vect = bd.array(
            (bd.array(self.x) - self.x[L // 2]) ** 2
            + (bd.array(self.y) - self.y[L // 2]) ** 2
            + (bd.array(self.z) - self.z[L // 2]) ** 2,
            bd.float,
        ) # distance from the center of the line, only work for uniform grid

        self.profile = bd.exp(-(vect ** 2) / (2 * (0.5 * vect.max()) ** 2)) # gaussian profile
        self.profile /= self.profile.sum() # normalize the profile
        self.profile *= self.amplitude # scale the profile to the amplitude

    def _handle_slices(
        self, x: ListOrSlice, y: ListOrSlice, z: ListOrSlice
    ) -> Tuple[List, List, List]:
        """Convert slices in the grid to lists

        This is necessary to make the source span the diagonal of the volume
        defined by the slices.

        Args:
            x: The x-location of the volume in the grid
            y: The y-location of the volume in the grid
            z: The z-location of the volume in the grid

        Returns:
            x, y, z: the x, y and z coordinates of the source as lists

        """

        # if list-indices were chosen:
        if isinstance(x, list) and isinstance(y, list) and isinstance(z, list):
            if len(x) != len(y) or len(y) != len(z) or len(z) != len(x):
                raise IndexError(
                    "sources require grid to be indexed with slices or equal length list-indices"
                )
            x = [self.grid._handle_distance(_x) for _x in x]
            y = [self.grid._handle_distance(_y) for _y in y]
            z = [self.grid._handle_distance(_z) for _z in z]
            return x, y, z

        # if a combination of list-indices and slices were chosen,
        # convert the list-indices to slices.
        # TODO: maybe issue a warning here?
        if isinstance(x, list):
            x = slice(
                self.grid._handle_distance(x[0]),
                self.grid._handle_distance(x[-1]),
                None,
            )
        if isinstance(y, list):
            y = slice(
                self.grid._handle_distance(y[0]),
                self.grid._handle_distance(y[-1]),
                None,
            )
        if isinstance(z, list):
            z = slice(
                self.grid._handle_distance(z[0]),
                self.grid._handle_distance(z[-1]),
                None,
            )

        # if we get here, we can assume slices:
        x0 = self.grid._handle_distance(x.start if x.start is not None else 0)
        y0 = self.grid._handle_distance(y.start if y.start is not None else 0)
        z0 = self.grid._handle_distance(z.start if z.start is not None else 0)
        x1 = self.grid._handle_distance(x.stop if x.stop is not None else self.grid.Nx)
        y1 = self.grid._handle_distance(y.stop if y.stop is not None else self.grid.Ny)
        z1 = self.grid._handle_distance(z.stop if z.stop is not None else self.grid.Nz)

        # we can now convert these coordinates into index lists
        m = max(abs(x1 - x0), abs(y1 - y0), abs(z1 - z0))
        if m < 2:
            raise ValueError("a LineSource should consist of at least two gridpoints")
        x = [v.item() for v in bd.array(bd.linspace(x0, x1, m, endpoint=False), bd.int)]
        y = [v.item() for v in bd.array(bd.linspace(y0, y1, m, endpoint=False), bd.int)]
        z = [v.item() for v in bd.array(bd.linspace(z0, z1, m, endpoint=False), bd.int)]

        return x, y, z

    def update_E(self):
        """Add the source to the electric field"""
        q = self.grid.time_steps_passed
        # if pulse
        if self.pulse:
            t1 = int(2 * pi / (self.frequency * self.hanning_dt / self.cycle))
            if q < t1:
                vect = self.profile * hanning(
                    self.frequency, q * self.hanning_dt, self.cycle
                )
            else:
                # src = - self.grid.E[self.x, self.y, self.z, 2]
                vect = self.profile * 0
        # if not pulse
        else:
            vect = self.profile * sin(2 * pi * q / self.period + self.phase_shift)
        # do not use list indexing here, as this is much slower especially for torch backend
        # DISABLED: self.grid.E[self.x, self.y, self.z, 2] = vect
        for x, y, z, value in zip(self.x, self.y, self.z, vect):
            self.grid.E[x, y, z, 2] += value

    def update_H(self):
        """Add the source to the magnetic field"""

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(period={self.period}, "
            f"amplitude={self.amplitude}, phase_shift={self.phase_shift}, "
            f"name={repr(self.name)})"
        )

    def __str__(self):
        s = "    " + repr(self) + "\n"
        x = f"[{self.x[0]}, ... , {self.x[-1]}]"
        y = f"[{self.y[0]}, ... , {self.y[-1]}]"
        z = f"[{self.z[0]}, ... , {self.z[-1]}]"
        s += f"        @ x={x}, y={y}, z={z}\n"
        return s

## PlaneSource class
class PlaneSource:
    """A source along a plane in the FDTD grid"""

    def __init__(
        self,
        period: Number = 15,
        amplitude: float = 1.0,
        phase_shift: float = 0.0,
        name: str = None,
        polarization: str = 'z',
    ):
        """Create a PlaneSource.

        Args:
            period: The period of the source. The period can be specified
                as integer [timesteps] or as float [seconds]
            amplitude: The amplitude of the source in simulation units
            phase_shift: The phase offset of the source.
            polarization: Axis of E-field polarization ('x','y',or 'z')
        """
        self.grid = None
        self.period = period
        self.amplitude = amplitude
        self.phase_shift = phase_shift
        self.name = name
        self.polarization = polarization

    def _register_grid(
        self, grid: Grid, x: ListOrSlice, y: ListOrSlice, z: ListOrSlice
    ):
        """Register a grid for the source.

        Args:
            grid: the grid to place the source into.
            x: The x-location of the source in the grid
            y: The y-location of the source in the grid
            z: The z-location of the source in the grid

        Note:
        As its name suggests, this source is a LINE source.
            Hence the source spans the diagonal of the cube
            defined by the slices in the grid.
        """
        self.grid = grid
        self.grid.sources.append(self)
        if self.name is not None:
            if not hasattr(grid, self.name):
                setattr(grid, self.name, self)
            else:
                raise ValueError(
                    f"The grid already has an attribute with name {self.name}"
                )

        self.x, self.y, self.z = self._handle_slices(x, y, z)

        self.period = grid._handle_time(self.period)
        self.frequency = 1.0 / self.period

        x = bd.arange(self.x.start, self.x.stop, 1) - (self.x.start + self.x.stop) // 2
        y = bd.arange(self.y.start, self.y.stop, 1) - (self.y.start + self.y.stop) // 2
        z = bd.arange(self.z.start, self.z.stop, 1) - (self.z.start + self.z.stop) // 2
        xvec, yvec, zvec = bd.broadcast_arrays( # broadcast the x, y and z arrays
            x[:, None, None], y[None, :, None], z[None, None, :]
        )
        _xvec = bd.array(xvec, float)
        _yvec = bd.array(yvec, float)
        _zvec = bd.array(zvec, float)

        profile = bd.ones(_xvec.shape)
        self.profile = self.amplitude * profile

    def _handle_slices(
        self, x: ListOrSlice, y: ListOrSlice, z: ListOrSlice
    ) -> Tuple[List, List, List]:
        """Validate slices and calculate center of plane

        Args:
            x: The x-location of the volume in the grid
            y: The y-location of the volume in the grid
            z: The z-location of the volume in the grid

        Returns:
            x, y, z: the x, y and z coordinates of the source as slices

        """
        # ensure all slices
        if not isinstance(x, slice):
            if isinstance(x, list):
                (x,) = x
            x = slice(
                self.grid._handle_distance(x), self.grid._handle_distance(x) + 1, None
            )
        if not isinstance(y, slice):
            if isinstance(y, list):
                (y,) = y
            y = slice(
                self.grid._handle_distance(y), self.grid._handle_distance(y) + 1, None
            )
        if not isinstance(z, slice):
            if isinstance(z, list):
                (z,) = z
            z = slice(
                self.grid._handle_distance(z), self.grid._handle_distance(z) + 1, None
            )

        # if we get here, we can assume slices:
        x0 = self.grid._handle_distance(x.start if x.start is not None else 0)
        y0 = self.grid._handle_distance(y.start if y.start is not None else 0)
        z0 = self.grid._handle_distance(z.start if z.start is not None else 0)
        x1 = self.grid._handle_distance(x.stop if x.stop is not None else self.grid.Nx)
        y1 = self.grid._handle_distance(y.stop if y.stop is not None else self.grid.Ny)
        z1 = self.grid._handle_distance(z.stop if z.stop is not None else self.grid.Nz)

        # make sure all slices have a start, stop and no step:
        x = (
            slice(x0, x1)
            if x0 < x1
            else (slice(x1, x0) if x0 > x1 else slice(x0, x0 + 1))
        )
        y = (
            slice(y0, y1)
            if y0 < y1
            else (slice(y1, y0) if y0 > y1 else slice(y0, y0 + 1))
        )
        z = (
            slice(z0, z1)
            if z0 < z1
            else (slice(z1, z0) if z0 > z1 else slice(z0, z0 + 1))
        )

        if [x.stop - x.start, y.stop - y.start, z.stop - z.start].count(0) > 0:
            raise ValueError(
                "Given location for PlaneSource results in slices of length 0!"
            )
        if [x.stop - x.start, y.stop - y.start, z.stop - z.start].count(1) == 0:
            raise ValueError("Given location for PlaneSource is not a 2D plane!")
        if [x.stop - x.start, y.stop - y.start, z.stop - z.start].count(1) > 1:
            raise ValueError(
                "Given location for PlaneSource should have no more than one dimension in which it's flat.\n"
                "Use a LineSource for lower dimensional sources."
            )

        self._Epol = 'xyz'.index(self.polarization)
        if (x.stop - x.start == 1 and self.polarization == 'x') or \
           (y.stop - y.start == 1 and self.polarization == 'y') or \
           (z.stop - z.start == 1 and self.polarization == 'z'):
            raise ValueError(
                "PlaneSource cannot be polarized perpendicular to the orientation of the plane."
            )
        _Hpols = [(z,1,2), (z,0,2), (y,0,1)][self._Epol]
        if _Hpols[0].stop - _Hpols[0].start == 1:
            self._Hpol = _Hpols[1]
        else:
            self._Hpol = _Hpols[2]

        return x, y, z

    def update_E(self):
        """Add the source to the electric field"""
        q = self.grid.time_steps_passed
        vect = self.profile * sin(2 * pi * q / self.period + self.phase_shift)
        self.grid.E[self.x, self.y, self.z, self._Epol] = vect

    def update_H(self):
        """Add the source to the magnetic field"""
        q = self.grid.time_steps_passed
        vect = self.profile * sin(2 * pi * q / self.period + self.phase_shift)
        self.grid.H[self.x, self.y, self.z, self._Hpol] = vect

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(period={self.period}, "
            f"amplitude={self.amplitude}, phase_shift={self.phase_shift}, "
            f"name={repr(self.name)}, polarization={repr(self.polarization)})"
        )

    def __str__(self):
        s = "    " + repr(self) + "\n"
        x = f"[{self.x.start}, ... , {self.x.stop}]"
        y = f"[{self.y.start}, ... , {self.y.stop}]"
        z = f"[{self.z.start}, ... , {self.z.stop}]"
        s += f"        @ x={x}, y={y}, z={z}\n"
        return s

class ComplexPlaneWave:
    """A 2D plane wave source (Ez polarization)"""

    def __init__(
        self,
        # 支援兩種介面方式
        wavelength: float = None,           # 原有的單波長介面
        wavelength_min: float = None,       # 新增：Lumerical 風格介面
        wavelength_max: float = None,       # 新增：Lumerical 風格介面
        
        # 其他參數保持不變
        amplitude: complex = 1.0 + 0.0j,
        phase_shift: float = 0.0,
        theta_deg: float = 0.0,
        polarization_axis: str = "x",
        name: str = None,
    
        # 脈衝參數 - 改為自動計算
        pulse: bool = None,                 # None = 自動決定
        cycle: int = None,                  # None = 自動計算
        optimize_for_short_pulse: bool = True,  # 新增：模仿 Lumerical 選項
        
        hanning_dt: float = None,
        medium_n: float = None
    ):
        """
        創建ComplexPlaneWave，支援 Lumerical 風格的波長設定
        
        使用方式：
        1. 原有方式：wavelength=1550e-9  
        2. Lumerical 風格：wavelength_min=1500e-9, wavelength_max=1600e-9
        3. 單頻 Lumerical 風格：wavelength_min=1550e-9, wavelength_max=1550e-9
        """
    
        # 解析波長設定
        self._parse_wavelength_input(wavelength, wavelength_min, wavelength_max)
    
        # 基本參數
        self.amplitude = amplitude
        self.phase_shift = phase_shift
        self.theta_deg = theta_deg
        self.theta = bd.deg2rad(theta_deg)
        self.polarization_axis = polarization_axis.lower()
        self.name = name
        self.n = medium_n if medium_n is not None else 1.0
        self.optimize_for_short_pulse = optimize_for_short_pulse

        # 計算頻率參數（使用中心波長）
        self.frequency = bd.c0 / self.n / self.wavelength  # Hz
        self.period = 1.0 / self.frequency  # 秒
        self.omega = 2 * bd.pi * self.frequency  # rad/s
        self.k = 2 * np.pi / self.wavelength  # 1/m
        
        # Lumerical 風格：總是使用脈衝，自動計算參數
        self._setup_lumerical_pulse_params(pulse, cycle)
            
        # 初始化其他屬性（保持原有結構）
        self.grid = None
        self.period_sim = None
        self.frequency_sim = None
        self.omega_sim = None
        self.hanning_dt_sim = None
        self.hanning_dt = hanning_dt
        
        # 初始化其他屬性
        self.grid = None
        self.period_sim = None      # 將在_register_grid中計算
        self.frequency_sim = None   # 將在_register_grid中計算
        self.omega_sim = None       # 將在_register_grid中計算
        self.hanning_dt_sim = None       # 將在_register_grid中計算

        # monitor source properties（保持原有）
        self.monitoring_enabled = False
        self.monitor_data = {
            'timesteps': [],
            'E_field': [],
            'envelope': [],
            'phase': []
        }
        self.monitor_center_x = None
        self.monitor_center_z = None
    
    def _parse_wavelength_input(self, wavelength, wavelength_min, wavelength_max):
        """解析波長輸入參數"""
        
        # 計算設定了多少種方式
        old_style = wavelength is not None
        new_style = wavelength_min is not None or wavelength_max is not None
        
        if old_style and new_style:
            raise ValueError("不能同時使用舊式 (wavelength) 和新式 (wavelength_min/max) 參數")
        
        if old_style:
            # 舊式：單一波長
            self.wavelength = wavelength
            self.wavelength_min = wavelength
            self.wavelength_max = wavelength
            self.input_style = "legacy"
            
        elif new_style:
            # 新式：Lumerical 風格
            if wavelength_min is None or wavelength_max is None:
                raise ValueError("使用 Lumerical 風格時，必須同時設定 wavelength_min 和 wavelength_max")
                
            if wavelength_min <= 0 or wavelength_max <= 0:
                raise ValueError("波長必須大於 0")
            if wavelength_min > wavelength_max:
                raise ValueError("wavelength_min 不能大於 wavelength_max")
                
            self.wavelength_min = wavelength_min
            self.wavelength_max = wavelength_max
            self.wavelength = (wavelength_min + wavelength_max) / 2  # 中心波長
            self.input_style = "lumerical"
            
        else:
            raise ValueError("必須設定 wavelength 或 (wavelength_min, wavelength_max)")
        
        # 計算頻譜特性
        self.wavelength_span = self.wavelength_max - self.wavelength_min
        self.frequency_center = bd.c0 / self.n / self.wavelength
        self.frequency_min = bd.c0 / self.n / self.wavelength_max  # 注意反向
        self.frequency_max = bd.c0 / self.n / self.wavelength_min
        self.frequency_bandwidth = self.frequency_max - self.frequency_min
        self.relative_bandwidth = self.wavelength_span / self.wavelength if self.wavelength > 0 else 0
        
        # 判斷是否為單頻
        tolerance = 1e-12
        self.is_single_frequency = abs(self.wavelength_max - self.wavelength_min) < tolerance

    def _setup_lumerical_pulse_params(self, pulse_override, cycle_override):
        """設定 Lumerical 風格的脈衝參數"""
        
        # Lumerical 風格：總是使用脈衝
        self.pulse = True  # 強制使用脈衝
        
        if cycle_override is not None:
            # 用戶手動指定 cycle
            self.cycle = cycle_override
            self.pulse_mode = "manual"
        else:
            # 自動計算 cycle（模仿 Lumerical 邏輯）
            if self.is_single_frequency:
                # 單頻：較多週期，產生更窄的頻譜
                if self.optimize_for_short_pulse:
                    self.cycle = 5   # 對應 Lumerical 的 "standard" 脈衝
                else:
                    self.cycle = 10  # 更長脈衝，更精確
                self.pulse_mode = "single_freq_auto"
            else:
                # 寬頻：根據相對頻寬調整
                if self.optimize_for_short_pulse:
                    # 短脈衝，對應 Lumerical 的 "broadband" 脈衝
                    self.cycle = max(2, min(4, int(1.0 / max(self.relative_bandwidth, 0.1))))
                else:
                    # 稍長脈衝，減少頻率範圍外的功率
                    self.cycle = max(3, min(8, int(1.5 / max(self.relative_bandwidth, 0.05))))
                self.pulse_mode = "broadband_auto"
        
        # 設定脈衝類型標籤（用於顯示）
        if self.is_single_frequency:
            self.lumerical_pulse_type = "standard"
        else:
            self.lumerical_pulse_type = "broadband"

    def get_lumerical_info(self):
        """返回 Lumerical 風格的資訊"""
        return {
            'input_style': self.input_style,
            'is_single_frequency': self.is_single_frequency,
            'lumerical_pulse_type': self.lumerical_pulse_type,
            'pulse_mode': self.pulse_mode,
            'optimize_for_short_pulse': self.optimize_for_short_pulse,
            
            # 波長資訊
            'wavelength_min_nm': self.wavelength_min * 1e9,
            'wavelength_max_nm': self.wavelength_max * 1e9,
            'wavelength_center_nm': self.wavelength * 1e9,
            'wavelength_span_nm': self.wavelength_span * 1e9,
            
            # 頻率資訊
            'frequency_center_THz': self.frequency / 1e12,
            'frequency_min_THz': self.frequency_min / 1e12,
            'frequency_max_THz': self.frequency_max / 1e12,
            'frequency_bandwidth_THz': self.frequency_bandwidth / 1e12,
            'relative_bandwidth_percent': self.relative_bandwidth * 100,
            
            # 脈衝資訊
            'pulse_cycles': self.cycle,
            'pulse_duration_periods': self.cycle,
            'pulse_duration_fs': self.cycle * self.period * 1e15,
        }

    def print_lumerical_info(self):
        """印出 Lumerical 風格的資訊"""
        info = self.get_lumerical_info()
        
        print(f"📊 ComplexPlaneWave: {self.name or 'unnamed'}")
        print(f"   Input Style: {info['input_style']}")
        print(f"   Source Type: {'Single Frequency' if info['is_single_frequency'] else 'Broadband'}")
        print(f"   Pulse Type: {info['lumerical_pulse_type']}")
        print(f"   Pulse Mode: {info['pulse_mode']}")
        print()
        
        if info['is_single_frequency']:
            print(f"🎯 Single Frequency Setting:")
            print(f"   Wavelength: {info['wavelength_center_nm']:.1f} nm")
            print(f"   Frequency: {info['frequency_center_THz']:.2f} THz")
        else:
            print(f"🌈 Broadband Setting:")
            print(f"   Wavelength Range: {info['wavelength_min_nm']:.1f} - {info['wavelength_max_nm']:.1f} nm")
            print(f"   Center Wavelength: {info['wavelength_center_nm']:.1f} nm")
            print(f"   Span: {info['wavelength_span_nm']:.1f} nm ({info['relative_bandwidth_percent']:.1f}%)")
            print(f"   Frequency Range: {info['frequency_min_THz']:.2f} - {info['frequency_max_THz']:.2f} THz")
        print()
        
        print(f"⚡ Pulse Configuration:")
        print(f"   Always Pulse: {self.pulse}")
        print(f"   Cycles: {info['pulse_cycles']}")
        print(f"   Duration: {info['pulse_duration_fs']:.1f} fs")
        print(f"   Optimize for Short Pulse: {info['optimize_for_short_pulse']}")
        print()
        
        print(f"💡 How it works (Lumerical style):")
        if info['is_single_frequency']:
            print(f"   • Longer pulse ({info['pulse_cycles']} cycles) creates narrow spectrum")
            print(f"   • FFT extracts steady-state response at {info['frequency_center_THz']:.2f} THz")
            print(f"   • Result: Perfect CW response at target frequency")
        else:
            print(f"   • Shorter pulse ({info['pulse_cycles']} cycles) creates broad spectrum") 
            print(f"   • FFT extracts responses across {info['frequency_min_THz']:.2f}-{info['frequency_max_THz']:.2f} THz")
            print(f"   • Result: Broadband frequency response from single simulation")

    def _register_grid(self, grid: Grid, x: ListOrSlice, y: ListOrSlice, z: ListOrSlice):
        self.grid = grid
        self.grid.sources.append(self)
        
        if self.name is not None:
            if hasattr(grid, self.name):
                raise ValueError(f"The grid already has an attribute with name {self.name}")
            setattr(grid, self.name, self)

        self.x, self.y, self.z = self._handle_slices(x, y, z) # convert slices to lists

        # period, frequency and omega in simulation time steps
        self.period_sim = self.period / grid.time_step  # timesteps
        self.frequency_sim = 1.0 / self.period_sim         # 1/timestep  
        self.omega_sim = 2 * bd.pi / self.period_sim       # rad/timestep
        # 將 hanning_dt 轉成模擬步數
        self.hanning_dt_sim = 0.5 / self.frequency_sim
        # === 空間相位設定 計算 kx, kz ===
        kx = self.k * bd.sin(self.theta)
        kz = self.k * bd.cos(self.theta)
        # 計算每個格點的空間相位
        self.spatial_phase = {
            (xi, zi): (kx * xi * grid.grid_spacing + kz * zi * grid.grid_spacing)
            for xi, zi in zip(self.x, self.z)
        }

        self.pol_index = {"x": 0, "y": 1, "z": 2}[self.polarization_axis]
        
    def _handle_slices(
        self, x: ListOrSlice, y: ListOrSlice, z: ListOrSlice
    ) -> Tuple[List, List, List]:
        """Convert slices in the grid to lists

        This is necessary to make the source span the diagonal of the volume
        defined by the slices.

        Args:
            x: The x-location of the volume in the grid
            y: The y-location of the volume in the grid
            z: The z-location of the volume in the grid

        Returns:
            x, y, z: the x, y and z coordinates of the source as lists

        """

        # if list-indices were chosen:
        if isinstance(x, list) and isinstance(y, list) and isinstance(z, list):
            if len(x) != len(y) or len(y) != len(z) or len(z) != len(x):
                raise IndexError(
                    "sources require grid to be indexed with slices or equal length list-indices"
                )
            x = [self.grid._handle_distance(_x) for _x in x]
            y = [self.grid._handle_distance(_y) for _y in y]
            z = [self.grid._handle_distance(_z) for _z in z]
            return x, y, z

        # if a combination of list-indices and slices were chosen,
        # convert the list-indices to slices.
        # TODO: maybe issue a warning here?
        if isinstance(x, list):
            x = slice(
                self.grid._handle_distance(x[0]),
                self.grid._handle_distance(x[-1]),
                None,
            )
        if isinstance(y, list):
            y = slice(
                self.grid._handle_distance(y[0]),
                self.grid._handle_distance(y[-1]),
                None,
            )
        if isinstance(z, list):
            z = slice(
                self.grid._handle_distance(z[0]),
                self.grid._handle_distance(z[-1]),
                None,
            )

        # if we get here, we can assume slices:
        x0 = self.grid._handle_distance(x.start if x.start is not None else 0)
        y0 = self.grid._handle_distance(y.start if y.start is not None else 0)
        z0 = self.grid._handle_distance(z.start if z.start is not None else 0)
        x1 = self.grid._handle_distance(x.stop if x.stop is not None else self.grid.Nx)
        y1 = self.grid._handle_distance(y.stop if y.stop is not None else self.grid.Ny)
        z1 = self.grid._handle_distance(z.stop if z.stop is not None else self.grid.Nz)

        # we can now convert these coordinates into index lists
        m = max(abs(x1 - x0), abs(y1 - y0), abs(z1 - z0))
        if m < 2:
            raise ValueError("a LineSource should consist of at least two gridpoints")
        x = [v.item() for v in bd.array(bd.linspace(x0, x1, m, endpoint=False), bd.int)]
        y = [v.item() for v in bd.array(bd.linspace(y0, y1, m, endpoint=False), bd.int)]
        z = [v.item() for v in bd.array(bd.linspace(z0, z1, m, endpoint=False), bd.int)]

        return x, y, z

    def update_E(self):

        q = self.grid.time_steps_passed
                
        # 載波相位
        time_phase = self.omega_sim * q + self.phase_shift

        if self.pulse:
            pulse_duration = int(self.cycle * self.period_sim)
            
            if q < pulse_duration:
                # 即時計算hanning窗值
                envelope = 0.5 * (1 - bd.cos(2 * bd.pi * q / pulse_duration))
            else:
                envelope = 0
        else:
            envelope = 1
        
        # 更新每個格點的電場
        for xi, zi in zip(self.x, self.z):
            spatial_phase = self.spatial_phase[(xi, zi)]
            total_phase = time_phase + spatial_phase
            
            field_value = self.amplitude * bd.exp(1j * total_phase) * envelope
            self.grid.E[xi, 0, zi, self.pol_index] += field_value

            # 如果是中心點，記錄數據
            if self.monitoring_enabled: 
                if xi == self.monitor_center_x and zi == self.monitor_center_z:
                    self._record_source_output(q, envelope, field_value, time_phase)

    def update_H(self):
        pass

    # 新增便捷創建方法
    @classmethod
    def single_frequency(cls, wavelength: float, **kwargs):
        """便捷方法：創建單頻光源（Lumerical 風格）"""
        return cls(wavelength_min=wavelength, wavelength_max=wavelength, **kwargs)

    @classmethod
    def broadband(cls, wavelength_min: float, wavelength_max: float, **kwargs):
        """便捷方法：創建寶頻光源（Lumerical 風格）"""
        return cls(wavelength_min=wavelength_min, wavelength_max=wavelength_max, **kwargs)

    @classmethod
    def broadband_center_span(cls, wavelength_center: float, wavelength_span: float, **kwargs):
        """便捷方法：創建寬頻光源（中心+跨度方式）"""
        wavelength_min = wavelength_center - wavelength_span / 2
        wavelength_max = wavelength_center + wavelength_span / 2
        return cls(wavelength_min=wavelength_min, wavelength_max=wavelength_max, **kwargs)

    def get_source_power(self, grid_spacing):
        """功率計算"""
        Z_medium = bd.eta0 / self.n
        E0 = abs(self.amplitude)
        print(f"eta0 = {bd.eta0:.2f} Ω, n = {self.n}, Z_medium = {Z_medium:.2f} Ω")
        
        source_length = len(getattr(self, 'x', [1])) * grid_spacing
        power_density = 0.5 * E0**2 / Z_medium
        P_incident = power_density * source_length
        
        print(f"🔋 源功率: {P_incident:.6e} W/m (長度: {source_length*1e6:.2f}μm)")
        return P_incident
    
    def enable_monitoring(self):
        """啟用監測功能"""
        self.monitoring_enabled = True
        # 將座標轉換為排序的列表
        x_list = sorted(list(self.x))
        z_list = sorted(list(self.z))
        # 取中心點
        self.monitor_center_x = x_list[len(x_list) // 2]
        self.monitor_center_z = z_list[len(z_list) // 2]
        # print(f" 監測已啟用")。

    def _record_source_output(self, q, envelope, field_value, time_phase):
        """紀錄source輸出數據"""
        if not self.monitoring_enabled:
            return
        
        # 記錄當前狀態
        self.monitor_data['timesteps'].append(q)
        self.monitor_data['E_field'].append(field_value)
        self.monitor_data['envelope'].append(envelope)
        # 記錄載波相位
        self.monitor_data['carrier_phase'].append(time_phase)
    
    def get_monitor_data(self):
        """獲取監測數據"""
        if not self.monitoring_enabled:
            print("監測未啟用，請先調用 enable_monitoring()")
            return None
        return self.monitor_data.copy()
    
    def plot_source_analysis(self):
        """分析並繪製source輸出"""
        if not self.monitoring_enabled or len(self.monitor_data['timesteps']) == 0:
            print("沒有監測數據")
            return
        
        import matplotlib.pyplot as plt

        timesteps = bd.array(self.monitor_data['timesteps'])
        E_field = bd.array(self.monitor_data['E_field'])

        # 繪製
        plt.figure(figsize=(12, 6))
        plt.plot(timesteps, bd.real(E_field), 'b-', label='Real part', linewidth=1.5)
        plt.plot(timesteps, bd.imag(E_field), 'r--', label='Imag part', linewidth=1.5)
        plt.plot(timesteps, bd.abs(E_field), 'g:', label='|E|', linewidth=2)
        
        plt.title(f'Source E_field vs timesteps ({"Pulse" if self.pulse else "CW"})')
        plt.xlabel('Timesteps')
        plt.ylabel('Electric Field')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    def analyze_source_output(self):
        """簡單的數值統計"""
        if not self.monitoring_enabled or len(self.monitor_data['timesteps']) == 0:
            print("沒有監測數據")
            return
        
        E_field = bd.array(self.monitor_data['E_field'])
        
        print(f"Source Output Analysis:")
        print(f"  數據點數: {len(self.monitor_data['timesteps'])}")
        print(f"  最大振幅: {bd.max(bd.abs(E_field)):.3f}")
        print(f"  最小振幅: {bd.min(bd.abs(E_field)):.3f}")

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(period={self.period}, "
            f"amplitude={self.amplitude}, phase_shift={self.phase_shift}, "
            f"name={repr(self.name)})"
        )

    def __str__(self):
        s = "    " + repr(self) + "\n"
        if hasattr(self, 'x') and hasattr(self, 'y') and hasattr(self, 'z'):
            x = f"[{self.x[0]}, ... , {self.x[-1]}]" if len(self.x) > 1 else f"[{self.x[0]}]"
            y = f"[{self.y[0]}, ... , {self.y[-1]}]" if len(self.y) > 1 else f"[{self.y[0]}]"
            z = f"[{self.z[0]}, ... , {self.z[-1]}]" if len(self.z) > 1 else f"[{self.z[0]}]"
            s += f"        @ x={x}, y={y}, z={z}\n"
            if hasattr(self, 'source_type'):
                s += f"        type: {self.source_type}\n"
        return s

class SoftArbitraryPointSource:
    r"""

    A source placed at a single point (grid cell) in the grid.
    This source is special: it's both a source and a detector.

    Unlike the other sources, the input is a voltage, not an electric field.
    (really? why? should we convert back and forth?)

    For electrical measurements I've only needed a single-index source,
    so I don't know how the volume/line sources above work.
    We want the FFT function to operate over any detector.
    Maybe all sources should take an arbitary waveform argument?

    Each index in the *waveform* array represents 1 value at a timestep.

    There are many different *geometries* of "equivalent sources".
    The detector/source paradigm used in /fdtd might perhaps not correspond to this in an ideal
    fashion.

    It's not intuitively clear to me what a "soft" source would imply in the optical case, or what
    impedance even means for a laser.

    /fdtd/ seems to have found primary use in optical circles,
    so the default Z should probably be 0.

    "Whilst established for microwaves and electrical circuits,
    this concept has only very recently been observed in the optical domain,
    yet is not well defined or understood."[1]

    [1]: Optical impedance of metallic nano-structures, M. Mazilu and K. Dholakia
    https://doi.org/10.1364/OE.14.007709

    [2]: http://www.gwoptics.org/learn/02_Plane_waves/01_Fabry_Perot_cavity/02_Impedance_matched.php

    /\/\-
    """

    def __init__(
        self, waveform_array: ndarray, name: str = None, impedance: float = 0.0
    ):
        """Create



        Args:
            waveform_array
        """
        self.grid = None
        self.name = name
        self.current_detector = None
        self.waveform_array = waveform_array
        self.impedance = impedance
        self.input_voltage = []  # voltage hard-imposed by the source
        self.source_voltage = []  #
        # "field" rather than "voltage" might be more meaningful
        # FIXME: these voltage time histories have a different dimensionality

    def _register_grid(self, grid: Grid, x: Number, y: Number, z: Number):
        """Register a grid for the source.

        Args:
            grid: the grid to place the source into.
            x: The x-location of the source in the grid
            y: The y-location of the source in the grid
            z: The z-location of the source in the grid

        Note:
            As its name suggests, this source is a POINT source.
            Hence it should be placed at a single coordinate tuple
            int the grid.
        """
        self.grid = grid
        self.grid.sources.append(self)
        if self.name is not None:
            if not hasattr(grid, self.name):
                setattr(grid, self.name, self)
            else:
                raise ValueError(
                    f"The grid already has an attribute with name {self.name}"
                )

        try:
            (x,), (y,), (z,) = x, y, z
        except (TypeError, ValueError):
            raise ValueError("a point source should be placed on a single grid cell.")
        self.x, self.y, self.z = grid._handle_tuple((x, y, z))

        if self.name is not None:
            detector_name += "_I"
        else:
            detector_name = None

        self.current_detector = CurrentDetector(name=detector_name)
        grid[x, y, z] = self.current_detector

    def update_E(self):

        # It is important that this step happen between the E-field update and the
        # H-field update.

        if self.grid.time_steps_passed < self.waveform_array.shape[0]:
            # check for off-by-one error here
            input_voltage = self.waveform_array[self.grid.time_steps_passed]
        else:
            input_voltage = 0.0  # one could taper the last value off smoothly instead

        if self.grid.time_steps_passed > 0:
            current = self.current_detector.I[-1][0][0][0]
        else:
            current = 0.0

        if self.impedance > 0:
            source_resistive_voltage = self.impedance * current
            output_voltage = input_voltage + source_resistive_voltage
        else:
            output_voltage = input_voltage

        # right now, this does not compensate for the cell's permittivity!

        self.grid.E[self.x, self.y, self.z, 2] += (
            output_voltage / self.grid.grid_spacing
        )

        self.input_voltage.append([[[input_voltage]]])
        self.source_voltage.append([[[output_voltage]]])

    def update_H(self):
        pass

    def __repr__(self):
        return f"{self.__class__.__name__}()"

    def __str__(self):
        s = "    " + repr(self) + "\n"
        x = f"{self.x}"
        y = f"{self.y}"
        z = f"{self.z}"
        s += f"        @ x={x}, y={y}, z={z}\n"
        return s