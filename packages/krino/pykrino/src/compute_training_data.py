#!/usr/bin/env python3
"""
Python script to replicate the ComputeDistanceToSTL.computeTrainingData test
using the PyKrino interface.

This script:
1. Loads an STL surface from a file
2. Gets its bounding box and scales it
3. Generates random query points within the scaled bounding box
4. Computes signed distances for all query points
5. Writes the data to a CSV file

Usage:
    python3 compute_training_data.py [STL_FILENAME] [NUM_SAMPLES]
"""

import os
import sys
import numpy as np

import pykrino
from pykrino.printP0 import printP0


def generate_random_locations_in_box(bbox, num_points, seed=None):
    """
    Generate random points within a bounding box.

    Args:
        bbox: Bounding box as array with shape (2, 3) where bbox[0] is min and bbox[1] is max
        num_points: Number of random points to generate
        seed: Random seed for reproducibility

    Returns:
        Array of random points with shape (num_points, 3)
    """
    if seed is not None:
        np.random.seed(seed)

    min_vals = bbox[0]
    max_vals = bbox[1]

    # Generate random points in the range [0, 1]
    random_points = np.random.rand(num_points, 3)

    # Scale to the bounding box dimensions
    ranges = max_vals - min_vals
    random_points = random_points * ranges + min_vals

    return random_points


def compute_training_data(STLfilename, outputFilename, num_query_points, seed=None):
    """
    Compute training data by sampling signed distances around an STL surface.

    Args:
        STLfilename: Path to STL file
        outputFilename: Path to CSV file
        num_query_points: Number of query points to generate
        seed: Random seed for reproducibility
    """
    printP0("Loading STL surface from: {STLfilename}")
    distance_surface = pykrino.STLSurface(STLfilename)

    print("Getting bounding box and scaling by 1.2")
    bounding_box = distance_surface.bounding_box(scale=1.2)

    printP0("Generating {num_query_points} random query points")
    query_points = generate_random_locations_in_box(bounding_box, num_query_points, seed=seed)

    print("Computing signed distances")
    signed_distances = distance_surface.signed_distance(query_points)

    printP0(f"Writing data to: {outputFilename}")

    with open(outputFilename, 'w') as file:
        file.write("x,y,z,dist\n")
        for i in range(num_query_points):
            x, y, z = query_points[i]
            dist = signed_distances[i]
            file.write(f"{x},{y},{z},{dist}\n")

    printP0("Successfully wrote {num_query_points} samples to {outputFilename}")


def get_num_samples():
    """Get number of samples from environment variable or use default."""
    default_samples = 5000
    env_var = os.getenv("KRINO_SDF_NUM_SAMPLES")
    if env_var:
        try:
            return int(env_var)
        except ValueError:
            printP0("Warning: Invalid KRINO_SDF_NUM_SAMPLES value '{env_var}', using default {default_samples}")
    return default_samples


def get_stl_filename():
    """Get STL filename from environment variable or use default."""
    default_filename = "vehicle.stl"
    env_var = os.getenv("KRINO_STL_FILENAME")
    if env_var:
        return env_var
    return default_filename


def main():
    """Main function to parse arguments and run the computation."""
    # Parse command line arguments
    if len(sys.argv) > 3:
        print("Usage: python3 compute_training_data.py [STL_FILENAME] [NUM_SAMPLES]")
        sys.exit(1)

    # Get STL filename
    if len(sys.argv) >= 2:
        stl_filename = sys.argv[1]
    else:
        stl_filename = get_stl_filename()

    # Get number of samples
    if len(sys.argv) >= 3:
        try:
            num_samples = int(sys.argv[2])
        except ValueError:
            printP0("Error: Invalid number of samples '{sys.argv[2]}'")
            sys.exit(1)
    else:
        num_samples = get_num_samples()

    # Check if file exists
    if not os.path.exists(stl_filename):
        printP0("Error: STL file '{stl_filename}' does not exist")
        sys.exit(1)

    printP0("Computing {num_samples} samples from {stl_filename}")
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    output_filename = f"training_data_{rank}.txt"
    compute_training_data(stl_filename, output_filename, num_samples, seed=42)


if __name__ == "__main__":
    envData = pykrino.EnvData(MPI.COMM_WORLD)
    main()
