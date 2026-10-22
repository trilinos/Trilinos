import unittest

import pykrino


class TestEnvironment(unittest.TestCase):
    """pykrino sets up Kokkos and its communicator itself.  Clean exit is tested by test_exit_with_live_mesh and
    test_exit_after_finalize, each its own process (a python subprocess of an MPI process can hang in MPI_Init)."""

    def test_mesh_needs_no_setup(self):
        mesh = pykrino.BoundingBoxMesh("TET_4", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5, 0)
        mesh.populate_mesh()
        self.assertGreaterEqual(pykrino.parallel_rank(), 0)
        self.assertGreater(pykrino.parallel_size(), 0)

    def test_compatible_with_env_data_created_and_deleted(self):
        env_data = pykrino.EnvData(__import__("mpi4py").MPI.COMM_WORLD)
        self.assertEqual(pykrino.parallel_size(), env_data.size())
        del env_data
        mesh = pykrino.BoundingBoxMesh("TET_4", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5, 0)
        mesh.populate_mesh()


if __name__ == "__main__":
    unittest.main()
