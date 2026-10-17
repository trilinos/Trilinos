import unittest
import pykrino


class TestSmoothMesh(unittest.TestCase):
    """Unit tests for the pykrino smooth mesh Python bindings."""

    def test_smooth_mesh_bindings_are_available(self):
        """Test that the smooth-mesh-related bindings are exported."""

        expected_bindings = [
            "SmoothMesh",
            "Optimizer",
            "GradientFilter",
            "MeanRatioElementObjective",
            "MeshObjective",
            "SethHillParams",
            "SethHillElementObjective",
            "SmoothReference",
            "KinematicMeshObjective",
            "KinematicIsotropicMeshObjective",
            "KinematicConstantIsotropicMeshObjective",
            "KinematicAnisotropicMeshObjective",
            # These classes are bound with spaces in their Python names.
            "Element Objective",
            "Element Kinematic Objective",
            "Mesh Objective Base Class",
        ]

        for binding_name in expected_bindings:
            with self.subTest(binding_name=binding_name):
                self.assertTrue(
                    hasattr(pykrino, binding_name),
                    msg=f"pykrino is missing binding: {binding_name}",
                )

    def test_optimizer_enum_values_are_available(self):
        """Test that the currently-bound optimizer enum values are available."""

        self.assertTrue(hasattr(pykrino.Optimizer, "Type"))

        expected_enum_values = [
            "LBFGS",
            "SteepestDescent",
            "ROL",
        ]

        for enum_value in expected_enum_values:
            with self.subTest(enum_value=enum_value):
                self.assertTrue(hasattr(pykrino.Optimizer.Type, enum_value))

                # Because export_values() is used, these should also be available
                # directly on pykrino.Optimizer.
                self.assertTrue(hasattr(pykrino.Optimizer, enum_value))

    def test_optimizer_enum_values_are_not_callable(self):
        """Test that optimizer enum values are passed directly, not called."""

        enum_values = [
            pykrino.Optimizer.LBFGS,
            pykrino.Optimizer.SteepestDescent,
            pykrino.Optimizer.ROL,
        ]

        for enum_value in enum_values:
            with self.subTest(enum_value=enum_value):
                self.assertFalse(callable(enum_value))

    def test_optimizer_construction_with_exported_enum_values(self):
        """Test Optimizer construction using exported enum values."""

        enum_values = [
            pykrino.Optimizer.LBFGS,
            pykrino.Optimizer.SteepestDescent,
            pykrino.Optimizer.ROL,
        ]

        for enum_value in enum_values:
            with self.subTest(enum_value=enum_value):
                optimizer = pykrino.Optimizer(enum_value, 1.0e-4)
                self.assertIsNotNone(optimizer)

    def test_optimizer_construction_with_nested_enum_values(self):
        """Test Optimizer construction using Optimizer.Type enum values."""

        enum_values = [
            pykrino.Optimizer.Type.LBFGS,
            pykrino.Optimizer.Type.SteepestDescent,
            pykrino.Optimizer.Type.ROL,
        ]

        for enum_value in enum_values:
            with self.subTest(enum_value=enum_value):
                optimizer = pykrino.Optimizer(enum_value, 1.0e-4)
                self.assertIsNotNone(optimizer)

    def test_optimizer_default_tolerance(self):
        """Test Optimizer construction using the default tolerance argument."""

        optimizer = pykrino.Optimizer(pykrino.Optimizer.LBFGS)
        self.assertIsNotNone(optimizer)

    def test_element_objective_base_class_binding_name(self):
        """Test that the element objective base class is bound under its current name."""

        element_objective_base = getattr(pykrino, "Element Objective")
        self.assertIsNotNone(element_objective_base)

    def test_mean_ratio_element_objective_construction(self):
        """Test construction of the concrete mean-ratio element objective."""

        element_objective_base = getattr(pykrino, "Element Objective")

        obj = pykrino.MeanRatioElementObjective()

        self.assertIsNotNone(obj)
        self.assertIsInstance(obj, element_objective_base)

    def test_mesh_objective_base_class_binding_name(self):
        """Test that the mesh objective base class is bound under its current name."""

        mesh_objective_base = getattr(pykrino, "Mesh Objective Base Class")
        self.assertIsNotNone(mesh_objective_base)

    def test_gradient_filter_construction(self):
        """Test default GradientFilter construction."""

        gradient_filter = pykrino.GradientFilter()

        self.assertIsNotNone(gradient_filter)
        self.assertIsInstance(gradient_filter, pykrino.GradientFilter)

    def test_cube_boundary_gradient_filter_construction_with_defaults(self):
        """Test cube boundary GradientFilter construction using default bounds."""

        gradient_filter = pykrino.GradientFilter.create_cube_boundary_filter()

        self.assertIsNotNone(gradient_filter)
        self.assertIsInstance(gradient_filter, pykrino.GradientFilter)

    def test_cube_boundary_gradient_filter_construction_with_arguments(self):
        """Test cube boundary GradientFilter construction using explicit bounds."""

        # Use positional arguments unless the pybind11 binding explicitly defines
        # py::arg("min_coord") and py::arg("max_coord").
        gradient_filter = pykrino.GradientFilter.create_cube_boundary_filter(
            -2.0,
            2.0,
        )

        self.assertIsNotNone(gradient_filter)
        self.assertIsInstance(gradient_filter, pykrino.GradientFilter)

    def test_smooth_mesh_construction(self):
        """Test SmoothMesh construction."""

        smoother = pykrino.SmoothMesh()

        self.assertIsNotNone(smoother)
        self.assertTrue(hasattr(smoother, "smooth_mesh"))

    def test_mesh_objective_binding_is_available(self):
        """Test that MeshObjective class and set_gradient_filter method are bound."""

        self.assertTrue(hasattr(pykrino, "MeshObjective"))
        self.assertTrue(hasattr(pykrino.MeshObjective, "set_gradient_filter"))

    def test_seth_hill_params_import(self):
        """Test that SethHillParams class can be imported."""

        self.assertTrue(hasattr(pykrino, "SethHillParams"))

    def test_seth_hill_params_construction_with_defaults(self):
        """Test SethHillParams construction using default values."""

        params = pykrino.SethHillParams()

        self.assertIsNotNone(params)
        self.assertIsInstance(params, pykrino.SethHillParams)

    def test_seth_hill_params_construction_with_values(self):
        """Test SethHillParams construction using explicit values."""

        # Use positional arguments unless the pybind11 binding explicitly defines
        # py::arg("mu"), py::arg("kappa"), py::arg("m"), and py::arg("n").
        params = pykrino.SethHillParams(
            1.5,
            2.5,
            2,
            2,
        )

        self.assertIsNotNone(params)
        self.assertIsInstance(params, pykrino.SethHillParams)

    def test_seth_hill_objective_import(self):
        """Test that SethHillElementObjective class can be imported."""

        self.assertTrue(hasattr(pykrino, "SethHillElementObjective"))

    def test_seth_hill_objective_construction_default(self):
        """Test SethHillElementObjective construction with default parameters."""

        objective = pykrino.SethHillElementObjective()

        self.assertIsNotNone(objective)
        self.assertIsInstance(objective, pykrino.SethHillElementObjective)

    def test_seth_hill_objective_construction_with_params(self):
        """Test SethHillElementObjective construction with custom parameters."""

        params = pykrino.SethHillParams(
            1.5,
            2.5,
            2,
            2,
        )

        objective = pykrino.SethHillElementObjective(params)

        self.assertIsNotNone(objective)
        self.assertIsInstance(objective, pykrino.SethHillElementObjective)

    def test_seth_hill_objective_base_class(self):
        """Test that SethHillElementObjective is a kinematic element objective."""

        kinematic_objective_base = getattr(pykrino, "Element Kinematic Objective")
        self.assertIsNotNone(kinematic_objective_base)

        params = pykrino.SethHillParams()
        obj = pykrino.SethHillElementObjective(params)

        self.assertIsNotNone(obj)
        self.assertIsInstance(obj, kinematic_objective_base)


class TestKinematicAnisotropicMeshObjective(unittest.TestCase):
    """Tests for the anisotropic kinematic mesh objective with nodal metrics given from Python."""

    def setUp(self):
        self.mesh = pykrino.BoundingBoxMesh("TET_4", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5, 0)
        self.mesh.populate_mesh()
        self.dim = 3
        self.element_objective = pykrino.ModifiedNeoHookeanElementObjective()

    def make_objective(self):
        return pykrino.KinematicAnisotropicMeshObjective(self.mesh, self.element_objective)

    def test_node_ids_and_coordinates(self):
        obj = self.make_objective()
        node_ids, coords = obj.get_node_ids_and_coordinates()
        self.assertGreater(len(node_ids), 0)
        self.assertEqual(len(set(node_ids)), len(node_ids))
        self.assertEqual(len(coords), self.dim * len(node_ids))
        for value in coords:
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)

    def uniform_metrics(self, obj, size):
        node_ids, _ = obj.get_node_ids_and_coordinates()
        metric = [0.0] * (self.dim * self.dim)
        for i in range(self.dim):
            metric[i * self.dim + i] = 1.0 / (size * size)
        return metric * len(node_ids)

    def test_set_metrics_from_python_runs_smoothing(self):
        obj = self.make_objective()
        obj.set_element_metric_map_data_from_nodal_metrics(self.uniform_metrics(obj, 0.5))
        optimizer = pykrino.Optimizer(pykrino.Optimizer.Type.LBFGS)
        pykrino.SmoothMesh().smooth_mesh(self.mesh, optimizer, obj)

    def test_wrong_number_of_metric_values_raises(self):
        obj = self.make_objective()
        metrics = self.uniform_metrics(obj, 0.5)
        with self.assertRaises(RuntimeError):
            obj.set_element_metric_map_data_from_nodal_metrics(metrics[:-1])

    def test_non_symmetric_metric_raises(self):
        obj = self.make_objective()
        metrics = self.uniform_metrics(obj, 0.5)
        metrics[1] = 0.5  # M_01 of the first node, with M_10 still 0
        with self.assertRaises(RuntimeError):
            obj.set_element_metric_map_data_from_nodal_metrics(metrics)

    def test_missing_metric_field_raises(self):
        with self.assertRaises(RuntimeError):
            pykrino.KinematicAnisotropicMeshObjective(self.mesh, self.element_objective, "no_such_metric_field")

    def test_raw_element_metric_mapping_wrong_length_raises(self):
        obj = self.make_objective()
        with self.assertRaises(RuntimeError):
            obj.set_element_metric_map_data([1.0, 0.0])


if __name__ == "__main__":
    unittest.main()
