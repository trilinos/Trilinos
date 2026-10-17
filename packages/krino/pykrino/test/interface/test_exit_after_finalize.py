import unittest

import pykrino


class TestExitAfterFinalize(unittest.TestCase):
    """Releases the mesh and finalizes explicitly before exit; the atexit finalize must then be a no-op."""

    def test_release_mesh_then_finalize(self):
        mesh = pykrino.BoundingBoxMesh("TET_4", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5, 0)
        mesh.populate_mesh()
        del mesh
        pykrino.finalize()
        pykrino.finalize()  # safe to call twice


if __name__ == "__main__":
    unittest.main()
