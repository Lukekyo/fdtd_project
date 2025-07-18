""" Python 3D FDTD Simulator """

__author__ = "Floris laporte"
__version__ = "0.3.6"

from .grid import Grid
from .sources import PointSource, LineSource, ComplexLineSource, PlaneSource, ComplexPlaneWave
from .detectors import LineDetector, BlockDetector, CurrentDetector
from .objects import Object, AbsorbingObject, AnisotropicObject
from .boundaries import PeriodicBoundary, BlochBoundary, PML
from .backend import backend
from .backend import set_backend
from .fourier import FrequencyRoutines
from .visualization import dB_map_2D, plot_detection
from .fdtd_helper import (
    um, nm, to_grid, from_grid,
    quick_diagnosis
    )
from .visualize_domain import plot_simulation_domain