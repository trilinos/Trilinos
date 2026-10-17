import unittest

import pykrino


class TestExitWithLiveMesh(unittest.TestCase):
    """Leaves a mesh and refinement alive at exit: pykrino.finalize (atexit) must not finalize Kokkos while they are
    alive, so their destruction during interpreter teardown does not abort (which PythonTest.cmake fails on)."""

    def test_build_mesh_and_refinement_and_keep_them(self):
        global mesh, refinement
        mesh = pykrino.BoundingBoxMesh("TET_4", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5, 0)
        refinement = pykrino.Refinement(mesh)
        mesh.populate_mesh()


if __name__ == "__main__":
    unittest.main()
