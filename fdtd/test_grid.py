import unittest
from unittest.mock import MagicMock
from fdtd.grid import Grid

class TestGridSetItem(unittest.TestCase):
    def setUp(self):
        # Create a mock attribute with a _register_grid method
        self.mock_attr = MagicMock()
        self.mock_attr._register_grid = MagicMock()

        # Initialize a Grid instance
        self.grid = Grid(shape=(10, 10, 10), grid_spacing=1e-9)

    def test_setitem_single_index(self):
        # Test with a single index
        self.grid[5] = self.mock_attr
        self.mock_attr._register_grid.assert_called_once_with(
            grid=self.grid,
            x=[5],
            y=slice(None),
            z=slice(None),
        )

    def test_setitem_two_indices(self):
        # Test with two indices
        self.grid[5, 6] = self.mock_attr
        self.mock_attr._register_grid.assert_called_once_with(
            grid=self.grid,
            x=[5],
            y=[6],
            z=slice(None),
        )

    def test_setitem_three_indices(self):
        # Test with three indices
        self.grid[5, 6, 7] = self.mock_attr
        self.mock_attr._register_grid.assert_called_once_with(
            grid=self.grid,
            x=[5],
            y=[6],
            z=[7],
        )

    def test_setitem_invalid_indices(self):
        # Test with more than three indices
        with self.assertRaises(KeyError):
            self.grid[5, 6, 7, 8] = self.mock_attr

if __name__ == "__main__":
    unittest.main()