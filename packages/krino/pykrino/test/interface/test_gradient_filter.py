import unittest
import numpy as np
import pykrino

class TestGradientFilter(unittest.TestCase):

    def test_gradient_filter_import(self):
        """Test that GradientFilter class can be imported"""
        self.assertTrue(hasattr(pykrino, 'GradientFilter'))

    def test_gradient_filter_construction(self):
        """Test GradientFilter construction"""
        # Test default construction
        filter1 = pykrino.GradientFilter()
        self.assertIsNotNone(filter1)

        # Test cube boundary filter creation
        filter2 = pykrino.GradientFilter.create_cube_boundary_filter()
        self.assertIsNotNone(filter2)

        # Test custom cube boundary filter with parameters
        filter3 = pykrino.GradientFilter.create_cube_boundary_filter(-2.0, 2.0)
        self.assertIsNotNone(filter3)

    def test_filter_behavior_simulation(self):
        """Simulate filter behavior with test data"""
        # Test the logic of cube boundary filtering
        def test_cube_filter_logic(coords, min_val, max_val):
            """Test the filtering logic"""
            dir = np.array([1.0, 1.0, 1.0])

            # Apply the same logic as our cube boundary filter
            for i in range(3):
                if coords[i] <= min_val + 1e-6 or coords[i] >= max_val - 1e-6:
                    dir[i] = 0.0

            return dir

        # Test node on face (should zero out that direction)
        result = test_cube_filter_logic([1.0, 0.5, 0.5], -1.0, 1.0)
        expected = np.array([0.0, 1.0, 1.0])
        np.testing.assert_array_equal(result, expected)

        # Test node on corner (should zero out all directions)
        result = test_cube_filter_logic([1.0, 1.0, 1.0], -1.0, 1.0)
        expected = np.array([0.0, 0.0, 0.0])
        np.testing.assert_array_equal(result, expected)

        # Test interior node (should not modify direction)
        result = test_cube_filter_logic([0.5, 0.5, 0.5], -1.0, 1.0)
        expected = np.array([1.0, 1.0, 1.0])
        np.testing.assert_array_equal(result, expected)

if __name__ == "__main__":
    unittest.main()
