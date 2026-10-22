#!/usr/bin/env python3
"""
Example script demonstrating gradient filter functionality in pykrino.

This script shows how to:
1. Create different types of gradient filters
2. Apply gradient filters to mesh objectives
3. Use gradient filters in mesh smoothing operations

Note: This is a conceptual example. Actual usage requires:
- A valid mesh file
- Proper error handling
"""

import sys
import os
import numpy as np

# Add pykrino to path if not installed
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    import pykrino
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure pykrino is built.")
    sys.exit(1)

def example_basic_gradient_filters():
    """Example of creating basic gradient filters"""
    print("=== Basic Gradient Filter Examples ===")

    # Example 1: Create an empty gradient filter (no-op)
    empty_filter = pykrino.GradientFilter()
    print("Created empty gradient filter")

    # Example 2: Create a cube boundary filter with default bounds (-1, 1)
    cube_filter = pykrino.GradientFilter.create_cube_boundary_filter()
    print("Created cube boundary filter with default bounds [-1, 1]")

    # Example 3: Create a cube boundary filter with custom bounds
    custom_cube_filter = pykrino.GradientFilter.create_cube_boundary_filter(-2.0, 2.0)
    print("Created cube boundary filter with custom bounds [-2, 2]")

def example_custom_gradient_filter():
    """Example of creating custom gradient filters with Python callbacks"""
    print("\n=== Custom Gradient Filter Examples ===")

    # Example 1: Filter that restricts movement in Z direction
    def restrict_z_movement(node_id, node_coords, direction):
        """Zero out Z-direction movement for all nodes"""
        x, y, z = node_coords
        dx, dy, dz = direction
        dz = 0.0  # Zero out Z direction
        return (dx, dy, dz)

    z_filter = pykrino.GradientFilter(restrict_z_movement)
    print("Created filter that restricts Z-direction movement")

    # Example 2: Filter that creates a cylindrical constraint
    def cylindrical_constraint(node_id, node_coords, direction):
        """Restrict radial movement to create cylindrical constraint"""
        x, y, z = node_coords
        dx, dy, dz = direction

        # Calculate distance from Z axis
        radius = np.sqrt(x**2 + y**2)

        # If node is near the cylinder surface (radius ~ 1.0), restrict radial movement
        if abs(radius - 1.0) < 0.1:
            # Project movement onto tangential direction
            # This is a simplified example - real implementation would be more complex
            dx, dy = 0.0, 0.0

        return (dx, dy, dz)

    cylinder_filter = pykrino.GradientFilter(cylindrical_constraint)
    print("Created cylindrical constraint filter")

    # Example 3: Filter based on node location regions
    def regional_filter(node_id, node_coords, direction):
        """Apply different constraints based on spatial regions"""
        x, y, z = node_coords
        dx, dy, dz = direction

        # Region 1: Left half - restrict X movement
        if x < 0:
            dx = 0.0

        # Region 2: Right half - allow all movement
        # (no modification needed)

        # Region 3: Bottom layer - restrict Z movement
        if z < -0.5:
            dz = 0.0

        return (dx, dy, dz)

    region_filter = pykrino.GradientFilter(regional_filter)
    print("Created regional constraint filter")

def example_applying_gradient_filters():
    """Example of applying gradient filters to mesh objectives"""
    print("\n=== Applying Gradient Filters ===")

    # Note: This example shows the API usage, but requires a real mesh to work
    try:
        # Create a gradient filter
        cube_filter = pykrino.GradientFilter.create_cube_boundary_filter()
        print("Created gradient filter")

        # Method 1: Apply filter directly to mesh objective
        # mesh_obj.set_gradient_filter(cube_filter)
        print("Gradient filter can be applied directly using mesh_objective.set_gradient_filter()")

    except Exception as e:
        print(f"Note: {e}")
        print("This example requires a valid mesh setup to work completely.")

def example_smoothing_with_filters():
    """Example of using gradient filters in mesh smoothing"""
    print("\n=== Mesh Smoothing with Gradient Filters ===")

    try:
        # Load mesh (conceptual - requires actual mesh file)
        mesh_file = "example_mesh.e"

        if not os.path.exists(mesh_file):
            print(f"Note: Mesh file '{mesh_file}' not found. Showing API usage only.")
            return

        # Create mesh and populate
        # mesh = pykrino.MeshFromFile(mesh_file)
        # mesh.populate_mesh()

        # Create mesh objective
        # element_obj = pykrino.MeanRatioElementObjective()
        # mesh_obj = pykrino.MeshObjective(mesh, element_obj)

        # Create and apply gradient filter
        cube_filter = pykrino.GradientFilter.create_cube_boundary_filter()
        # mesh_obj.set_gradient_filter(cube_filter)

        print("Gradient filter applied to mesh objective")

        # Create optimizer and smoother
        # optimizer = pykrino.Optimizer(pykrino.Optimizer.Type.LBFGS)
        # smoother = pykrino.SmoothMesh()

        # Apply smoothing with gradient constraints
        # smoother.smooth_mesh(mesh, optimizer, mesh_obj, print_diagnostics=True)

        print("Mesh smoothing with gradient constraints completed")

    except Exception as e:
        print(f"Note: {e}")
        print("This example requires a valid mesh setup to work completely.")

def example_filter_composition():
    """Example of combining multiple gradient filters"""
    print("\n=== Gradient Filter Composition ===")

    # Note: Current implementation doesn't support filter composition directly,
    # but you can create a composite filter using Python callbacks

    def filter1(node_id, node_coords, direction):
        """First filter: Restrict Z movement"""
        x, y, z = node_coords
        dx, dy, dz = direction
        dz = 0.0
        return (dx, dy, dz)

    def filter2(node_id, node_coords, direction):
        """Second filter: Restrict movement on cube boundaries"""
        x, y, z = node_coords
        dx, dy, dz = direction

        # Apply cube boundary logic
        for i in range(3):
            coord = [x, y, z][i]
            if coord <= -1.0 + 1e-6 or coord >= 1.0 - 1e-6:
                [dx, dy, dz][i] = 0.0

        return (dx, dy, dz)

    def composite_filter(node_id, node_coords, direction):
        """Composite filter that applies both filters sequentially"""
        # Apply first filter
        result1 = filter1(node_id, node_coords, direction)

        # Apply second filter to the result
        result2 = filter2(node_id, node_coords, result1)

        return result2

    composite = pykrino.GradientFilter(composite_filter)
    print("Created composite gradient filter using Python callbacks")

def example_filter_visualization():
    """Example showing how gradient filters affect node movement"""
    print("\n=== Gradient Filter Effect Visualization ===")

    # Simulate the effect of a cube boundary filter
    def simulate_cube_filter(coords, direction, min_val=-1.0, max_val=1.0):
        """Simulate cube boundary filter effect"""
        result_dir = list(direction)

        for i in range(3):
            if coords[i] <= min_val + 1e-6 or coords[i] >= max_val - 1e-6:
                result_dir[i] = 0.0

        return tuple(result_dir)

    # Test cases
    test_cases = [
        ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], "Interior node"),
        ([1.0, 0.5, 0.5], [1.0, 1.0, 1.0], "Node on X+ face"),
        ([-1.0, 0.5, 0.5], [1.0, 1.0, 1.0], "Node on X- face"),
        ([0.5, 1.0, 0.5], [1.0, 1.0, 1.0], "Node on Y+ face"),
        ([0.5, -1.0, 0.5], [1.0, 1.0, 1.0], "Node on Y- face"),
        ([0.5, 0.5, 1.0], [1.0, 1.0, 1.0], "Node on Z+ face"),
        ([0.5, 0.5, -1.0], [1.0, 1.0, 1.0], "Node on Z- face"),
        ([1.0, 1.0, 1.0], [1.0, 1.0, 1.0], "Node on corner"),
    ]

    print("\nEffect of cube boundary filter on node movement:")
    print("Original Direction: (1, 1, 1)")
    print("-" * 60)

    for coords, direction, description in test_cases:
        filtered_dir = simulate_cube_filter(coords, direction)
        print(f"{description:20} @ {coords}")
        print(f"  Original: {direction}")
        print(f"  Filtered: {filtered_dir}")
        print()

def main():
    """Main function to run examples"""
    print("PyKrino Gradient Filter Examples")
    print("=" * 60)

    # Run examples
    example_basic_gradient_filters()
    example_custom_gradient_filter()
    example_applying_gradient_filters()
    example_smoothing_with_filters()
    example_filter_composition()
    example_filter_visualization()

    print("\n" + "=" * 60)
    print("All gradient filter examples completed!")
    print("\nKey Takeaways:")
    print("1. Gradient filters control node movement directions during smoothing")
    print("2. Built-in filters like cube_boundary_filter handle common cases")
    print("3. Custom filters can be created using Python callbacks")
    print("4. Filters can be applied to mesh objectives using set_gradient_filter()")
    print("5. Multiple filters can be combined using composite callbacks")

if __name__ == "__main__":
    main()
