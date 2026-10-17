#!/usr/bin/env python3
"""
Example script demonstrating pykrino mesh smoothing functionality.

This script shows how to:
1. Load a mesh from file
2. Set up mesh smoothing parameters
3. Apply mesh smoothing using different optimizers and objectives
4. Save the smoothed mesh

Note: This is a conceptual example. Actual usage requires:
- A valid mesh file
- Proper error handling

Currently available in the API:
- Optimizers: LBFGS, SteepestDescent, ROL
- Objectives: MeanRatioElementObjective
- Gradient filters: Cube boundary filter
"""

import sys
import os

# Add pykrino to path if not installed
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    import pykrino
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure pykrino is built.")
    sys.exit(1)

def example_basic_smoothing():
    """Example of basic mesh smoothing"""
    print("=== Basic Mesh Smoothing Example ===")

    # Load mesh from file (replace with actual mesh file)
    mesh_file = "example_mesh.e"

    if not os.path.exists(mesh_file):
        print(f"Warning: Mesh file '{mesh_file}' not found. Creating dummy example.")
        # In a real scenario, you would load an actual mesh file
        # For this example, we'll just demonstrate the API
        return

    try:
        # Create MeshFromFile object
        mesh = pykrino.MeshFromFile(mesh_file)

        # Populate mesh data
        mesh.populate_mesh()

        print("Mesh loaded successfully")

        # Create mesh smoother
        smoother = pykrino.SmoothMesh()

        # Set up optimizer (LBFGS with default tolerance)
        optimizer = pykrino.Optimizer(pykrino.Optimizer.Type.LBFGS)

        # Set up objective function (Mean Ratio Norm)
        objective = pykrino.MeanRatioElementObjective()

        # Create mesh objective
        mesh_obj = pykrino.MeshObjective(mesh, objective)

        # Apply smoothing
        print("Starting mesh smoothing...")

        smoother.smooth_mesh(
            mesh,
            optimizer,
            mesh_obj,
            print_diagnostics=True
        )

        print("Mesh smoothing completed")

        # In a real application, you would now:
        # 1. Access the smoothed mesh coordinates
        # 2. Save the smoothed mesh to a new file
        # 3. Analyze mesh quality improvements

    except Exception as e:
        print(f"Error during smoothing: {e}")
        raise

def example_advanced_smoothing():
    """Example of advanced mesh smoothing with gradient filters"""
    print("\n=== Advanced Mesh Smoothing Example ===")

    # Define custom gradient filter (Python callback)
    def custom_gradient_filter(node_id, node_coords, direction):
        """
        Custom gradient filter that modifies smoothing directions.

        Args:
            node_id: Node identifier
            node_coords: Tuple of (x, y, z) coordinates
            direction: Tuple of (dx, dy, dz) current direction

        Returns:
            Modified direction as tuple (dx, dy, dz)
        """
        # Example: Zero out Z-direction for nodes on certain planes
        x, y, z = node_coords
        dx, dy, dz = direction

        # If node is on x=0 or x=1 plane, don't move in x-direction
        if abs(x) < 1e-6 or abs(x - 1.0) < 1e-6:
            dx = 0.0

        # If node is on y=0 or y=1 plane, don't move in y-direction
        if abs(y) < 1e-6 or abs(y - 1.0) < 1e-6:
            dy = 0.0

        return (dx, dy, dz)

    try:
        # Load mesh (same as basic example)
        mesh_file = "example_mesh.e"

        if not os.path.exists(mesh_file):
            print(f"Warning: Mesh file '{mesh_file}' not found. Skipping advanced example.")
            return

        mesh = pykrino.MeshFromFile(mesh_file)
        mesh.populate_mesh()

        # Create smoother
        smoother = pykrino.SmoothMesh()

        # Use different optimizer (Steepest Descent)
        optimizer = pykrino.Optimizer(pykrino.Optimizer.Type.SteepestDescent, 1e-5)

        # Use different objective (Mean Ratio Norm - only one currently available)
        objective = pykrino.MeanRatioElementObjective()

        # Create mesh objective and apply gradient filter
        mesh_obj = pykrino.MeshObjective(mesh, objective)
        mesh_obj.set_gradient_filter(pykrino.GradientFilter(custom_gradient_filter))

        # Apply smoothing
        smoother.smooth_mesh(
            mesh,
            optimizer,
            mesh_obj,
            print_diagnostics=True
        )

        print("\nAdvanced mesh smoothing completed")

    except Exception as e:
        print(f"Error during advanced smoothing: {e}")
        raise

def example_comparison():
    """Example comparing different optimizers and objectives"""
    print("\n=== Optimizer and Objective Comparison ===")

    mesh_file = "example_mesh.e"

    if not os.path.exists(mesh_file):
        print(f"Warning: Mesh file '{mesh_file}' not found. Skipping comparison.")
        return

    try:
        # Load mesh
        mesh = pykrino.MeshFromFile(mesh_file)
        mesh.populate_mesh()

        smoother = pykrino.SmoothMesh()

        # Test different combinations
        # Note: Currently only MeanRatioElementObjective is available in the bindings
        configurations = [
            ("LBFGS + MeanRatio",
             pykrino.Optimizer.Type.LBFGS,
             pykrino.MeanRatioElementObjective),
            ("SteepestDescent + MeanRatio",
             pykrino.Optimizer.Type.SteepestDescent,
             pykrino.MeanRatioElementObjective),
            ("ROL + MeanRatio",
             pykrino.Optimizer.Type.ROL,
             pykrino.MeanRatioElementObjective),
        ]

        for name, opt_type, obj_class in configurations:
            print(f"\nTesting: {name}")

            optimizer = pykrino.Optimizer(opt_type)
            objective = obj_class()

            # Create mesh objective
            mesh_obj = pykrino.MeshObjective(mesh, objective)

            # Note: In a real application, you would want to:
            # 1. Save the original mesh state
            # 2. Apply smoothing
            # 3. Compare quality metrics
            # 4. Restore original state for next test

            smoother.smooth_mesh(mesh, optimizer, mesh_obj, print_diagnostics=False)

            print(f"Completed: {name}")

    except Exception as e:
        print(f"Error during comparison: {e}")
        raise

def main():
    """Main function to run examples"""
    print("PyKrino Mesh Smoothing Examples")
    print("=" * 50)

    # Run examples
    example_basic_smoothing()
    example_advanced_smoothing()
    example_comparison()

    print("\n" + "=" * 50)
    print("All examples completed!")

if __name__ == "__main__":
    main()
