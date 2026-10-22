import unittest

import numpy as np
import pykrino


def generate_random_locations_in_box(bbox, num_points, seed=None):
    """
    Generate random points within a bounding box.

    Args:
        bbox: Bounding box as array with shape (2, 3) where bbox[0] is min and bbox[1] is max
        num_points: Number of random points to generate
        seed: Random seed for reproducibility

    Returns:
        Array of random points with shape (num_points, 3)
    """
    if seed is not None:
        np.random.seed(seed)

    min_vals = bbox[0]
    max_vals = bbox[1]

    # Generate random points in the range [0, 1]
    random_points = np.random.rand(num_points, 3)

    # Scale to the bounding box dimensions
    ranges = max_vals - min_vals
    random_points = random_points * ranges + min_vals

    return random_points

class TestSurfaceBase(unittest.TestCase):
    """
    Base class for testing all Krino surface wrappers.

    Derived classes should call one of the _scenario_* methods to set:

        self.surface
        self.expected_bbox
        self.known_points
        self.expected_distance_magnitudes
    """

    def setUp(self):
        self.tol = 1.0e-10
        self.num_random_points = 1000
        self.surface = None
        self.expected_bbox = None
        self.known_points = None
        self.expected_distance_magnitudes = None

    def _scenario_unit_sphere(self):
        self.surface = pykrino.Sphere(
            center=[0.0, 0.0, 0.0],
            radius=1.0,
        )

        self.expected_bbox = np.array(
            [
                [-1.0, -1.0, -1.0],
                [1.0, 1.0, 1.0],
            ],
            dtype=np.float64,
        )

        self.known_points = np.array(
            [
                [0.0, 0.0, 0.0],  # center, distance magnitude 1.0
                [1.0, 0.0, 0.0],  # on surface, distance magnitude 0.0
                [2.0, 0.0, 0.0],  # outside, distance magnitude 1.0
            ],
            dtype=np.float64,
        )

        self.expected_distance_magnitudes = np.array(
            [
                1.0,
                0.0,
                1.0,
            ],
            dtype=np.float64,
        )

    def tearDown(self):
        if hasattr(self, "temp_dir"):
            self.temp_dir.cleanup()

    def assertArrayAllClose(self, actual, expected, tol=None):
        if tol is None:
            tol = self.tol

        np.testing.assert_allclose(
            actual,
            expected,
            rtol=tol,
            atol=tol,
        )

    def check_surface_exists(self):
        self.assertIsNotNone(self.surface)

    def check_bounding_box_shape(self):
        bbox = self.surface.bounding_box()

        self.assertIsInstance(bbox, np.ndarray)
        self.assertEqual(bbox.shape, (2, 3))
        self.assertTrue(np.issubdtype(bbox.dtype, np.floating))
        self.assertTrue(np.all(np.isfinite(bbox)))

    def check_bounding_box_values(self):
        if self.expected_bbox is None:
            self.skipTest("No expected bounding box defined for this surface")

        bbox = self.surface.bounding_box()

        self.assertArrayAllClose(bbox, self.expected_bbox)

    def check_scaled_bounding_box_shape(self):
        bbox = self.surface.bounding_box(scale=1.2)

        self.assertIsInstance(bbox, np.ndarray)
        self.assertEqual(bbox.shape, (2, 3))
        self.assertTrue(np.all(np.isfinite(bbox)))

    def check_random_points_in_bounding_box(self):
        bbox = self.surface.bounding_box(scale=1.2)
        points = generate_random_locations_in_box(
            bbox,
            self.num_random_points,
            seed=123,
        )

        self.assertEqual(points.shape, (self.num_random_points, 3))
        self.assertTrue(np.all(points >= bbox[0]))
        self.assertTrue(np.all(points <= bbox[1]))

    def check_signed_distance_shape(self):
        distances = self.surface.signed_distance(self.known_points)

        self.assertIsInstance(distances, np.ndarray)
        self.assertEqual(distances.shape, (self.known_points.shape[0],))
        self.assertTrue(np.issubdtype(distances.dtype, np.floating))

    def check_signed_distance_values_are_finite(self):
        distances = self.surface.signed_distance(self.known_points)

        self.assertTrue(np.all(np.isfinite(distances)))

    def check_known_distance_magnitudes(self):
        distances = self.surface.signed_distance(self.known_points)

        self.assertArrayAllClose(
            np.abs(distances),
            self.expected_distance_magnitudes,
        )

    def check_inside_and_outside_have_opposite_signs(self):
        """
        This test assumes that the first and last known points are on opposite
        sides of the surface and neither is exactly on the surface.

        It does not assume whether inside is positive or negative.
        """

        distances = self.surface.signed_distance(self.known_points)

        first = distances[0]
        last = distances[-1]

        self.assertLess(first * last, 0.0)

    def check_bad_point_shapes_raise(self):
        bad_point_arrays = [
            np.array([0.0, 0.0, 0.0]),
            np.zeros((10, 2)),
            np.zeros((10, 4)),
            np.zeros((2, 2, 3)),
        ]

        for bad_points in bad_point_arrays:
            with self.subTest(shape=bad_points.shape):
                with self.assertRaises(Exception):
                    self.surface.signed_distance(bad_points)

    def check_float32_points_are_accepted(self):
        points = self.known_points.astype(np.float32)

        distances = self.surface.signed_distance(points)

        self.assertEqual(distances.shape, (points.shape[0],))
        self.assertTrue(np.all(np.isfinite(distances)))

    def check_repeated_calls_are_stable(self):
        distances_1 = self.surface.signed_distance(self.known_points)
        distances_2 = self.surface.signed_distance(self.known_points)

        self.assertArrayAllClose(distances_1, distances_2, tol=0.0)

    def check_batching_consistency(self):
        bbox = self.surface.bounding_box(scale=1.2)

        points = generate_random_locations_in_box(
            bbox,
            self.num_random_points,
            seed=123,
        )

        distances_all = self.surface.signed_distance(points)

        distances_a = self.surface.signed_distance(points[:400])
        distances_b = self.surface.signed_distance(points[400:])

        distances_split = np.concatenate([distances_a, distances_b])

        self.assertArrayAllClose(distances_all, distances_split)

class TestSphereSurface(TestSurfaceBase):

    def setUp(self):
        super().setUp()
        self._scenario_unit_sphere()

    def test_Sphere_construction(self):
        self.check_surface_exists()

    def test_Sphere_bounding_box_shape(self):
        self.check_bounding_box_shape()

    def test_Sphere_bounding_box_values(self):
        self.check_bounding_box_values()

    def test_Sphere_scaled_bounding_box_shape(self):
        self.check_scaled_bounding_box_shape()

    def test_Sphere_random_points_in_bounding_box(self):
        self.check_random_points_in_bounding_box()

    def test_Sphere_signed_distance_shape(self):
        self.check_signed_distance_shape()

    def test_Sphere_signed_distance_values_are_finite(self):
        self.check_signed_distance_values_are_finite()

    def test_Sphere_known_distance_magnitudes(self):
        self.check_known_distance_magnitudes()

    def test_Sphere_inside_and_outside_have_opposite_signs(self):
        self.check_inside_and_outside_have_opposite_signs()

    def test_Sphere_bad_point_shapes_raise(self):
        self.check_bad_point_shapes_raise()

    def test_Sphere_float32_points_are_accepted(self):
        self.check_float32_points_are_accepted()

    def test_Sphere_repeated_calls_are_stable(self):
        self.check_repeated_calls_are_stable()

    def test_Sphere_batching_consistency(self):
        self.check_batching_consistency()

if __name__ == "__main__":
    unittest.main()
