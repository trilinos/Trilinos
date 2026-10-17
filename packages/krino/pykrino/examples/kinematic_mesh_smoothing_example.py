#!/usr/bin/env python3
"""
Example script demonstrating advanced objective functions in pykrino.

This script shows how to:
1. Use SethHillElementObjective for kinematic mesh smoothing
2. Use ModifiedNeoHookeanElementObjective for kinematic mesh smoothing
3. Configure parameters for advanced objectives

Note: These objectives are available in the bindings but may require
additional testing and validation.
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


def example_seth_hill_objective():
    """Example of using SethHillElementObjective"""
    print("=== SethHillElementObjective Example ===")

    # Method 1: Create with default parameters
    seth_hill_default = pykrino.SethHillElementObjective()
    print("Created SethHillElementObjective with default parameters")

    # Method 2: Create with custom SethHillParams
    # SethHillParams takes: mu, kappa, m, n
    params = pykrino.SethHillParams(
        mu=1.5,      # Shear modulus
        kappa=2.5,   # Bulk modulus
        m=2,         # Material parameter
        n=2          # Material parameter
    )

    seth_hill_custom = pykrino.SethHillElementObjective(params)
    print("Created SethHillElementObjective with custom parameters:")
    print(f"  mu={params.mu}, kappa={params.kappa}, m={params.m}, n={params.n}")

    # Note: SethHillElementObjective is a kinematic element objective
    # It can be used with KinematicMeshObjective
    print("\nSethHillElementObjective is a kinematic element objective")
    print("It can be used with KinematicMeshObjective for kinematic mesh smoothing")


def example_modified_neo_hookean_objective():
    """Example of using ModifiedNeoHookeanElementObjective"""
    print("\n=== ModifiedNeoHookeanElementObjective Example ===")

    # Create ModifiedNeoHookeanElementObjective
    neo_hookean = pykrino.ModifiedNeoHookeanElementObjective()
    print("Created ModifiedNeoHookeanElementObjective")

    print("\nModifiedNeoHookeanElementObjective is used for:")
    print("- Kinematic mesh smoothing")
    print("- Large deformation problems")
    print("- Hyperelastic material modeling")


def example_kinematic_mesh_objective():
    """Example of using KinematicMeshObjective with advanced objectives"""
    print("\n=== KinematicMeshObjective Example ===")

    # Note: This example shows the API usage pattern
    # Actual usage requires a valid mesh with kinematic elements

    try:
        # Create SethHill parameters
        params = pykrino.SethHillParams(
            mu=1.5,
            kappa=2.5,
            m=2,
            n=2
        )

        # Create SethHill objective
        seth_hill_obj = pykrino.SethHillElementObjective(params)
        print("Created SethHillElementObjective for kinematic mesh")

        # Create ModifiedNeoHookean objective
        neo_hookean_obj = pykrino.ModifiedNeoHookeanElementObjective()
        print("Created ModifiedNeoHookeanElementObjective for kinematic mesh")

        # Note: KinematicMeshObjective would be used like this:
        # mesh = pykrino.MeshFromFile("kinematic_mesh.e")
        # mesh.populate_mesh()
        #
        # kinematic_obj = pykrino.KinematicIsotropicMeshObjective(mesh, seth_hill_obj, pykrino.SmoothReference.EqualVolume)
        # smoother = pykrino.SmoothMesh()
        # optimizer = pykrino.Optimizer(pykrino.Optimizer.Type.LBFGS)
        # smoother.smooth_mesh(mesh, optimizer, kinematic_obj)

        print("\nKinematicMeshObjective usage pattern:")
        print("1. Load mesh with kinematic elements")
        print("2. Create kinematic objective (SethHill or ModifiedNeoHookean)")
        print("3. Create KinematicIsotropicMeshObjective or KinematicConstantIsotropicMeshObjective with mesh and objective")
        print("4. Apply smoothing with SmoothMesh")

    except Exception as e:
        print(f"Note: {e}")
        print("This example requires a valid kinematic mesh setup to work completely.")


def example_objective_comparison():
    """Example comparing different objective types"""
    print("\n=== Objective Comparison ===")

    print("Available objective types in PyKrino:")
    print()

    # MeanRatioElementObjective
    print("1. MeanRatioElementObjective:")
    print("   - Purpose: Improve mesh quality based on element mean ratio")
    print("   - Use case: General mesh smoothing")
    print("   - Status: Fully tested and available")
    mean_ratio_obj = pykrino.MeanRatioElementObjective()
    print(f"   - Created: {mean_ratio_obj}")
    print()

    # SethHillElementObjective
    print("2. SethHillElementObjective:")
    print("   - Purpose: Kinematic mesh smoothing with Seth-Hill material model")
    print("   - Use case: Large deformation problems")
    print("   - Status: Available in bindings, may need additional testing")
    params = pykrino.SethHillParams()
    seth_hill_obj = pykrino.SethHillElementObjective(params)
    print(f"   - Created: {seth_hill_obj}")
    print()

    # ModifiedNeoHookeanElementObjective
    print("3. ModifiedNeoHookeanElementObjective:")
    print("   - Purpose: Kinematic mesh smoothing with Neo-Hookean material model")
    print("   - Use case: Hyperelastic material modeling")
    print("   - Status: Available in bindings, may need additional testing")
    neo_hookean_obj = pykrino.ModifiedNeoHookeanElementObjective()
    print(f"   - Created: {neo_hookean_obj}")
    print()

    print("Note: Only MeanRatioElementObjective is currently fully tested.")
    print("SethHill and ModifiedNeoHookean objectives are available but may")
    print("require additional validation for production use.")


def example_parameter_exploration():
    """Example exploring SethHill parameter space"""
    print("\n=== SethHill Parameter Exploration ===")

    print("SethHillParams parameters:")
    print("- mu: Shear modulus (controls shear stiffness)")
    print("- kappa: Bulk modulus (controls volumetric stiffness)")
    print("- m: Material parameter (typically 2)")
    print("- n: Material parameter (typically 2)")
    print()

    # Example parameter sets
    parameter_sets = [
        {"name": "Soft material", "mu": 0.5, "kappa": 1.0, "m": 2, "n": 2},
        {"name": "Medium material", "mu": 1.5, "kappa": 2.5, "m": 2, "n": 2},
        {"name": "Stiff material", "mu": 3.0, "kappa": 5.0, "m": 2, "n": 2},
    ]

    for params_dict in parameter_sets:
        params = pykrino.SethHillParams(
            params_dict["mu"],
            params_dict["kappa"],
            params_dict["m"],
            params_dict["n"]
        )
        obj = pykrino.SethHillElementObjective(params)
        print(f"{params_dict['name']}: mu={params.mu}, kappa={params.kappa}")


def main():
    """Main function to run examples"""
    print("PyKrino Advanced Objectives Examples")
    print("=" * 60)

    # Run examples
    example_seth_hill_objective()
    example_modified_neo_hookean_objective()
    example_kinematic_mesh_objective()
    example_objective_comparison()
    example_parameter_exploration()

    print("\n" + "=" * 60)
    print("All advanced objectives examples completed!")
    print("\nKey Takeaways:")
    print("1. SethHillElementObjective: For kinematic mesh smoothing with Seth-Hill model")
    print("2. ModifiedNeoHookeanElementObjective: For hyperelastic material modeling")
    print("3. Both require SethHillParams configuration")
    print("4. Use KinematicMeshObjective for kinematic element objectives")
    print("5. MeanRatioElementObjective is recommended for general mesh smoothing")


if __name__ == "__main__":
    main()
