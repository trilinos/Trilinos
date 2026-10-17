import math
import unittest
import pykrino


class TestRefinement(unittest.TestCase):
    """Unit tests for the pykrino refinement Python bindings."""

    def test_refinement_bindings_are_available(self):
        """Test that the refinement-related bindings are exported."""

        expected_bindings = [
            "Refinement",
        ]

        for binding_name in expected_bindings:
            with self.subTest(binding_name=binding_name):
                self.assertTrue(
                    hasattr(pykrino, binding_name),
                    msg=f"pykrino is missing binding: {binding_name}",
                )

    def test_refinement_methods_are_available(self):
        """Test that Refinement methods are available."""

        expected_methods = [
            "do_uniform_refinement",
            "delete_parent_elements",
            "refine_edges_too_long_in_metric",
            "refine_edges_too_long_in_size",
        ]

        for method_name in expected_methods:
            with self.subTest(method_name=method_name):
                self.assertTrue(
                    hasattr(pykrino.Refinement, method_name),
                    msg=f"Refinement is missing method: {method_name}",
                )

    def test_refinement_construction(self):
        """Test Refinement construction from a bounding box mesh, before the mesh is populated."""

        for kwargs in ({}, {"force_64bit": True}, {"assert_32bit": True}):
            with self.subTest(**kwargs):
                mesh = pykrino.BoundingBoxMesh("TET_4", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5, 0)
                refinement = pykrino.Refinement(mesh, **kwargs)  # declares its fields, so before populate_mesh
                mesh.populate_mesh()
                self.assertIsNotNone(refinement)


class TestRefineEdgesTooLong(unittest.TestCase):
    """Tests for refining to a metric or size given as a Python function of coordinates."""

    def setUp(self):
        self.mesh = pykrino.BoundingBoxMesh("TET_4", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5, 0)
        self.refinement = pykrino.Refinement(self.mesh)  # declares its fields, so before populate_mesh
        self.mesh.populate_mesh()

    @staticmethod
    def anisotropic_metric(x):
        # Size 0.1 along x and 0.5 along y and z.
        return [1.0 / 0.1**2, 0.0, 0.0,
                0.0, 1.0 / 0.5**2, 0.0,
                0.0, 0.0, 1.0 / 0.5**2]

    def test_isotropic_size_converges(self):
        self.assertEqual(0, self.refinement.refine_edges_too_long_in_size(lambda x: 0.2))
        self.assertEqual(0, self.refinement.refine_edges_too_long_in_size(lambda x: 0.2))  # nothing left

    def test_size_varying_with_position_converges(self):
        self.assertEqual(0, self.refinement.refine_edges_too_long_in_size(lambda x: 0.1 + 0.3 * x[0]))

    def test_anisotropic_metric_converges(self):
        self.assertEqual(0, self.refinement.refine_edges_too_long_in_metric(self.anisotropic_metric))

    def test_too_few_passes_returns_edges_still_too_long(self):
        self.assertGreater(self.refinement.refine_edges_too_long_in_size(lambda x: 0.02, max_num_passes=1), 0)

    def test_wrong_number_of_metric_values_raises(self):
        with self.assertRaises(RuntimeError):
            self.refinement.refine_edges_too_long_in_metric(lambda x: [1.0, 0.0, 1.0])

    def test_refine_then_smooth_in_the_same_metric(self):
        self.assertEqual(0, self.refinement.refine_edges_too_long_in_metric(self.anisotropic_metric))
        self.refinement.delete_parent_elements()

        obj = pykrino.KinematicAnisotropicMeshObjective(self.mesh, pykrino.ModifiedNeoHookeanElementObjective())
        _, coords = obj.get_node_ids_and_coordinates()
        metrics = []
        for i in range(0, len(coords), 3):
            metrics.extend(self.anisotropic_metric(coords[i:i + 3]))
        obj.set_element_metric_map_data_from_nodal_metrics(metrics)
        pykrino.SmoothMesh().smooth_mesh(self.mesh, pykrino.Optimizer(pykrino.Optimizer.Type.LBFGS), obj)


if __name__ == "__main__":
    unittest.main()
