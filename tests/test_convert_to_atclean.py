#!/usr/bin/env python

import unittest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
from convert_to_atclean import (
    ConvertLightCurve,
    ConvertLoop,
    Coordinates,
    PresetColumnNames,
)


class TestConvertLightCurve(unittest.TestCase):

    def setUp(self):
        """Set up a dummy ConvertLightCurve object before each test."""
        self.colnames = MagicMock()
        self.colnames.mjd = "MJD"
        self.colnames.flux = "Flux"
        self.colnames.dflux = "DFlux"
        self.colnames.filt = "Filter"
        self.colnames.preset = "atlas"
        self.colnames.ra = "RA"
        self.colnames.dec = "Dec"

        self.lc = ConvertLightCurve("TestObject", self.colnames, control_index=0)
        self.lc.t = pd.DataFrame(
            {
                "MJD": [59000, 59001, 59002],
                "Flux": [1.2, 2.3, np.nan],
                "DFlux": [0.1, 0.2, 0.0],
                "RA": [10.5, 10.5, 10.5],
                "Dec": [-20.3, -20.3, -20.3],
                "Filter": ["g", "r", "g"],
            }
        )

    def test_move_required_cols_to_front(self):
        """Test that required columns are moved to the front."""
        self.lc.move_required_cols_to_front()
        expected_columns = ["MJD", "Flux", "DFlux", "RA", "Dec", "Filter"]
        self.assertEqual(list(self.lc.t.columns), expected_columns)

    def test_find_coords_in_t(self):
        """Test extraction of coordinates from table."""
        coords = self.lc.find_coords_in_t()
        self.assertEqual(coords.ra.angle.degree, 10.5)
        self.assertEqual(coords.dec.angle.degree, -20.3)

    def test_get_coords_with_command_line_input(self):
        """Test overriding file coordinates with command-line arguments."""
        coords = self.lc.get_coords(arg_ra=15.5, arg_dec=-25.5)
        self.assertEqual(coords.ra.angle.degree, 15.5)
        self.assertEqual(coords.dec.angle.degree, -25.5)

    def test_remove_invalid_rows(self):
        """Test removal of NaN flux values and zero dflux values."""
        self.lc.save("test_dir", overwrite=True)
        self.assertEqual(len(self.lc.t), 2)  # One row should be removed

    def test_keep_only_existing_columns(self):
        """Test keeping only existing columns when some columns are missing."""
        all_columns_to_copy = ["MJD", "Flux", "DFlux", "Filter", "NonExistentColumn"]
        self.lc.save(
            "test_dir", all_columns_to_copy=all_columns_to_copy, overwrite=True
        )
        self.assertNotIn("NonExistentColumn", self.lc.t.columns)
        self.assertNotIn("RA", self.lc.t.columns)
        self.assertNotIn("Dec", self.lc.t.columns)


class TestConvertLoop(unittest.TestCase):
    def setUp(self):
        """Set up a dummy ConvertLightCurve object before each test."""
        self.colnames = MagicMock()
        self.colnames.mjd = "MJD"
        self.colnames.flux = "Flux"
        self.colnames.dflux = "DFlux"
        self.colnames.filt = "Filter"
        self.colnames.preset = "atlas"
        self.colnames.ra = "RA"
        self.colnames.dec = "Dec"

    @patch("convert_to_atclean.SnInfoTable")
    def test_validate_args(self, mock_sninfo):
        """Test argument validation for filenames and control indices."""
        convert = ConvertLoop(self.colnames, "input_dir", "output_dir")

        # Valid case
        convert.validate_args(["file1.csv", "file2.csv"], [0, 1])

        # Mismatched list lengths
        with self.assertRaises(RuntimeError):
            convert.validate_args(["file1.csv"], [0, 1])

        # Invalid control index
        with self.assertRaises(RuntimeError):
            convert.validate_args(["file1.csv"], [-1])


if __name__ == "__main__":
    unittest.main()
