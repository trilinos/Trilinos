import unittest

import sys
import os
import tempfile

import pykrino
import numpy as np
from mpi4py import MPI
from pykrino.printP0 import printP0

def generate_random_locations_in_box(bbox, num_points, seed=None):
    """
    Generate random points within a bounding box.

    Args:
        bbox: Bounding box as array with shape (2, 3), where bbox[0] is min
              and bbox[1] is max.
        num_points: Number of random points to generate.
        seed: Random seed for reproducibility.

    Returns:
        Array of random points with shape (num_points, 3).
    """
    rng = np.random.default_rng(seed)

    bbox = np.asarray(bbox, dtype=np.float64)

    if bbox.shape != (2, 3):
        raise ValueError(f"bbox must have shape (2, 3), got {bbox.shape}")

    min_vals = bbox[0]
    max_vals = bbox[1]

    return rng.uniform(
        low=min_vals,
        high=max_vals,
        size=(num_points, 3),
    )


class TestComputeTrainingData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tol = 1.0e-10
        cls.envData = pykrino.EnvData(MPI.COMM_WORLD)

    def test_compute_training_data_with_sphere(self):
        """Test training-data generation using pykrino.Sphere."""

        printP0("Testing compute_training_data with PyKrino Sphere...")
            
        filename = "sphere.stl"
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        output_filename = f"training_data_{rank}.txt"
        num_query_points = 100
        pykrino.compute_training_data(filename, output_filename, num_query_points)

        self.assertTrue(os.path.exists(output_filename))

        with open(output_filename, "r", encoding="utf-8") as file:
            lines = file.readlines()

        self.assertEqual(len(lines), num_query_points + 1)

        header = lines[0].strip()
        self.assertEqual(header, "x,y,z,dist")

        data = np.loadtxt(
            output_filename,
            delimiter=",",
            skiprows=1,
        )

        self.assertIsInstance(data, np.ndarray)
        self.assertEqual(data.shape, (num_query_points, 4))
        self.assertTrue(np.all(np.isfinite(data)))

        query_points = data[:, 0:3]
        computed_distances = data[:, 3]

        self.assertEqual(query_points.shape, (num_query_points, 3))
        self.assertEqual(computed_distances.shape, (num_query_points,))
        self.assertTrue(np.all(np.isfinite(query_points)))
        self.assertTrue(np.all(np.isfinite(computed_distances)))

        sphere = pykrino.Sphere([0.0, 0.0, 0.0], 0.4)

        expected_distances = sphere.signed_distance(query_points)

        self.assertIsInstance(expected_distances, np.ndarray)
        self.assertEqual(expected_distances.shape, (num_query_points,))
        self.assertTrue(np.all(np.isfinite(expected_distances)))

        if self.envData.rank() == 0:
            max_abs_error = np.max(
                np.abs(computed_distances - expected_distances)
            )
            print(f"Maximum distance error versus pykrino.Sphere: {max_abs_error}")

        np.testing.assert_allclose(
            computed_distances,
            expected_distances,
            rtol=0,
            atol=1.e-3,
        )

    def test_sphere_known_distance_magnitudes(self):
        """Test known signed-distance magnitudes for a unit sphere."""

        distance_surface = pykrino.Sphere(
            center=[0.0, 0.0, 0.0],
            radius=1.0,
        )

        points = np.array(
            [
                [0.0, 0.0, 0.0],  # center, distance magnitude 1
                [1.0, 0.0, 0.0],  # on surface, distance magnitude 0
                [2.0, 0.0, 0.0],  # outside, distance magnitude 1
            ],
            dtype=np.float64,
        )

        distances = distance_surface.signed_distance(points)

        expected_magnitudes = np.array(
            [
                1.0,
                0.0,
                1.0,
            ],
            dtype=np.float64,
        )

        np.testing.assert_allclose(
            np.abs(distances),
            expected_magnitudes,
            rtol=self.tol,
            atol=self.tol,
        )

        # Inside and outside should have opposite signs, independent of whether
        # Krino uses positive-inside or negative-inside.
        self.assertLess(distances[0] * distances[2], 0.0)


if __name__ == "__main__":
    unittest.main()
