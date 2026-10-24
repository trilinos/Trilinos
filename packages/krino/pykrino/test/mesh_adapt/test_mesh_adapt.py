#!/usr/bin/env python3
"""
Integration test for mesh_adapt functionality using PyKrino.

This test module verifies the adapt_mesh() function from mesh_adapt.py
"""

import os
import sys
import tempfile
import shutil
import unittest

import pykrino
from mpi4py import MPI


class TestMeshAdapt(unittest.TestCase):
    """Integration tests for the mesh adaptation functionality."""

    @classmethod
    def setUpClass(cls):
        """Set up test class - initialize MPI."""
        cls.envData = pykrino.EnvData(MPI.COMM_WORLD)

    def setUp(self):
        """Set up individual test - create shared temporary directory."""
        self.create_shared_temp_dir()
        self.input_mesh = os.path.join(os.path.dirname(__file__), 'mesh.g')

        # Verify input mesh exists
        if not os.path.exists(self.input_mesh):
            raise FileNotFoundError(f"Test mesh file not found: {self.input_mesh}")

    def tearDown(self):
        self.remove_shared_temp_dir()
            
    def create_shared_temp_dir(self):
        """Create one temporary directory shared by all MPI ranks."""
        comm = MPI.COMM_WORLD
        if comm.Get_rank() == 0:
            self.test_dir = tempfile.mkdtemp()
        else:
            self.test_dir = None

        self.test_dir = comm.bcast(self.test_dir, root=0)

        # Make sure all ranks see the directory before continuing.
        comm.Barrier()

    def remove_shared_temp_dir(self):
        """Remove shared temporary directory after all ranks are done with it."""
        comm = MPI.COMM_WORLD
        comm.Barrier()

        if comm.Get_rank() == 0:
            shutil.rmtree(self.test_dir, ignore_errors=True)

        comm.Barrier()

    def test_adapt_mesh_basic(self):
        """Test basic mesh adaptation with default parameters."""
        output_mesh = os.path.join(self.test_dir, 'output.g')

        # Parse command line options
        argv = [
            '-i', self.input_mesh,
            '-o', output_mesh,
            '-n', '1'
        ]
        options = pykrino.parse_mesh_adapt_command_line(argv)

        # Verify options are parsed correctly
        self.assertEqual(options['input_mesh'], self.input_mesh)
        self.assertEqual(options['output_mesh'], output_mesh)
        self.assertEqual(options['num_uniform_refinement_levels'], 1)
        self.assertFalse(options['force_64bit'])
        self.assertFalse(options['assert_32bit'])
        self.assertTrue(options['auto_compose'])

        # Run mesh adaptation
        pykrino.adapt_mesh(options)

        # Verify output file was created
        self.assertTrue(os.path.exists(output_mesh))

        # Verify output is a valid mesh file
        self.assertGreater(os.path.getsize(output_mesh), 0)

    def test_adapt_mesh_multiple_refinement_levels(self):
        """Test mesh adaptation with multiple refinement levels."""
        output_mesh = os.path.join(self.test_dir, 'output_multi.g')

        # Parse command line options with multiple refinement levels
        argv = [
            '-i', self.input_mesh,
            '-o', output_mesh,
            '-n', '3'
        ]
        options = pykrino.parse_mesh_adapt_command_line(argv)

        # Run mesh adaptation
        pykrino.adapt_mesh(options)

        # Verify output file was created
        self.assertTrue(os.path.exists(output_mesh))
        self.assertGreater(os.path.getsize(output_mesh), 0)

    def test_adapt_mesh_64bit(self):
        """Test mesh adaptation with 64-bit entity IDs."""
        output_mesh = os.path.join(self.test_dir, 'output_64bit.g')

        # Parse command line options with 64-bit flag
        argv = [
            '-i', self.input_mesh,
            '-o', output_mesh,
            '-n', '1',
            '--force_64_bit'
        ]
        options = pykrino.parse_mesh_adapt_command_line(argv)

        # Verify 64-bit flag is set
        self.assertTrue(options['force_64bit'])

        # Run mesh adaptation
        pykrino.adapt_mesh(options)

        # Verify output file was created
        self.assertTrue(os.path.exists(output_mesh))
        self.assertGreater(os.path.getsize(output_mesh), 0)

    def test_adapt_mesh_no_compose(self):
        """Test mesh adaptation without auto-compose."""
        comm = MPI.COMM_WORLD
        num_procs = num_procs = comm.Get_size()
        rank = comm.Get_rank();
        output_mesh       = os.path.join(self.test_dir, "output_nocompose.g")
        if (num_procs == 1):
            output_mesh_split = output_mesh
        else:
            output_mesh_split = os.path.join(self.test_dir, f"output_nocompose.g.{num_procs}.{rank}")

        # Parse command line options without auto-compose
        argv = [
            '-i', self.input_mesh,
            '-o', output_mesh,
            '-n', '1',
            '--no_auto_compose'
        ]
        options = pykrino.parse_mesh_adapt_command_line(argv)

        # Verify auto-compose is disabled
        self.assertFalse(options['auto_compose'])

        # Run mesh adaptation
        pykrino.adapt_mesh(options)

        # Verify output file was created
        self.assertTrue(os.path.exists(output_mesh_split))
        self.assertGreater(os.path.getsize(output_mesh_split), 0)

    def test_adapt_mesh_32bit_assertion(self):
        """Test mesh adaptation with 32-bit assertion."""
        output_mesh = os.path.join(self.test_dir, 'output_32bit.g')

        # Parse command line options with 32-bit assertion
        argv = [
            '-i', self.input_mesh,
            '-o', output_mesh,
            '-n', '1',
            '--assert_32_bit'
        ]
        options = pykrino.parse_mesh_adapt_command_line(argv)

        # Verify 32-bit assertion is set
        self.assertTrue(options['assert_32bit'])

        # Run mesh adaptation
        pykrino.adapt_mesh(options)

        # Verify output file was created
        self.assertTrue(os.path.exists(output_mesh))
        self.assertGreater(os.path.getsize(output_mesh), 0)

    def test_parse_mesh_adapt_command_line_validation(self):
        """Test that parse_mesh_adapt_command_line validates input."""
        # Test that number_refines < 1 raises an error
        with self.assertRaises(SystemExit):
            argv = ['-i', self.input_mesh, '-o', 'output.g', '-n', '0']
            pykrino.parse_mesh_adapt_command_line(argv)


if __name__ == "__main__":
    unittest.main()
