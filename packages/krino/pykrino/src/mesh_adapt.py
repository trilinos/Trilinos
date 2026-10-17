#!/usr/bin/env python3
"""
Mesh Adaptation Script using PyKrino

This script provides command-line mesh adaptation functionality similar to krino_mesh_adapt,
but using PyKrino's Python bindings. It supports reading input meshes, applying refinement,
and writing output meshes.

Usage:
    python mesh_adapt.py [options]

Options:
    --input_mesh FILE      Input mesh file path
    --output_mesh FILE     Output mesh file path
    --uniform_levels N     Number of uniform refinement levels (default: 1)
    --force_64bit          Force 64-bit entity IDs
    --assert_32bit         Assert 32-bit entity IDs
    --no_compose           Don't delete parent elements after refinement
    --help                 Show this help message
"""

import sys
from mpi4py import MPI
import pykrino
from pykrino.printP0 import printP0
import argparse

def parse_mesh_adapt_command_line(argv=None):
    """
    Parse mesh-adaptation command-line options into a dictionary.

    Args:
        argv (list[str], optional):
            Command-line arguments excluding the program name.
            If None, argparse uses sys.argv[1:].

    Returns:
        dict:
            Dictionary containing parsed command-line options.

    Raises:
        SystemExit:
            Raised by argparse if parsing fails or required arguments are missing.
    """

    parser = argparse.ArgumentParser(
        description="Mesh adaptation/refinement command-line parser.",
        allow_abbrev=False
    )

    # Required options
    parser.add_argument(
        "-i",
        "--input_mesh",
        required=True,
        help="Filename of input genesis mesh."
    )

    parser.add_argument(
        "-o",
        "--output_mesh",
        required=True,
        help="Filename of output genesis refined mesh to write."
    )

    # Optional integer option
    parser.add_argument(
        "-n",
        "--number_refines",
        type=int,
        default=1,
        help="Number of uniform refinement passes. Default is 1."
    )

    # Flag options
    parser.add_argument(
        "--force_64_bit",
        action="store_true",
        help="Force all new nodes and elements past 32-bit limit."
    )

    parser.add_argument(
        "--assert_32_bit",
        action="store_true",
        help="Assert that all nodes and elements are within the 32-bit limit."
    )

    parser.add_argument(
        "--no_auto_compose",
        action="store_true",
        help="Do not automatically recompose the parallel files."
    )

    args = parser.parse_args(argv)

    if args.number_refines < 1:
        parser.error("--number_refines must be >= 1")

    options = {
        "input_mesh": args.input_mesh,
        "output_mesh": args.output_mesh,
        "num_uniform_refinement_levels": args.number_refines,
        "force_64bit": args.force_64_bit,
        "assert_32bit": args.assert_32_bit,
        "auto_compose": not args.no_auto_compose,
    }
    
    return options

def adapt_mesh(options):
    print("PyKrino Mesh Adaptation Tool")
    printP0("=" * 40)
    
    print("Refinement options:")
    printP0(f"  Uniform refinement levels: {options['num_uniform_refinement_levels']}")
    printP0(f"  Force 64-bit: {options['force_64bit']}")
    printP0(f"  Assert 32-bit: {options['assert_32bit']}")
    printP0(f"  Auto compose: {options['auto_compose']}")
    printP0()
    
    mesh = pykrino.MeshFromFile(
            options['input_mesh'],
            decomp_method="RCB"
    )
        
    refinement = pykrino.Refinement(
        mesh,
        force_64bit=options['force_64bit'],
        assert_32bit=options['assert_32bit']
    )
        
    mesh.populate_mesh()

    if options['num_uniform_refinement_levels'] > 0:
        printP0(f"Applying {options['num_uniform_refinement_levels']} levels of uniform refinement...")
        refinement.do_uniform_refinement(options['num_uniform_refinement_levels'])
        printP0("Uniform refinement completed")
    
        printP0("Deleting parent elements...")
        refinement.delete_parent_elements()

    mesh.write_mesh(options['output_mesh'], compose=options['auto_compose'])

    printP0("Mesh adaptation completed successfully!")

def main():
    options = parse_mesh_adapt_command_line()
    adapt_mesh(options)

if __name__ == "__main__":
    envData = pykrino.EnvData(MPI.COMM_WORLD)
    main()
