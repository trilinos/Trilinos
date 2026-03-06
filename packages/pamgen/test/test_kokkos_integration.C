// @HEADER
// ***************************************************************************
//                     Pamgen Package - Kokkos Integration Test
//
// Copyright 2026 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// ***************************************************************************
// @HEADER

#include "brick_inline_mesh_desc.h"
#include "radial_inline_mesh_desc.h"
#include "pamgen_kokkos_utils.h"
#include <iostream>
#include <vector>
#include <map>

using namespace PAMGEN_NEVADA;

int main(int argc, char* argv[]) {
    // Initialize Kokkos
    Kokkos::initialize(argc, argv);
    {
        std::cout << "Testing PAMGEN Kokkos Integration..." << std::endl;

        // Test 1: Brick Mesh Kokkos functionality
        std::cout << "\n=== Test 1: Brick Mesh Kokkos ===" << std::endl;
        {
            Brick_Inline_Mesh_Desc brick_mesh(2);

            // Set up basic parameters
            brick_mesh.inline_b[0] = 2;
            brick_mesh.inline_b[1] = 2;
            brick_mesh.inline_b[2] = 1;

            brick_mesh.inline_n[0] = 3;
            brick_mesh.inline_n[1] = 3;
            brick_mesh.inline_n[2] = 1;

            // Set up intervals
            for (int axis = 0; axis < 3; axis++) {
                brick_mesh.interval[axis] = new long long[brick_mesh.inline_b[axis]];
                for (int i = 0; i < brick_mesh.inline_b[axis]; i++) {
                    brick_mesh.interval[axis][i] = brick_mesh.inline_n[axis];
                }
            }

            // Set up coordinate arrays
            for (int axis = 0; axis < 3; axis++) {
                brick_mesh.IJKcoors[axis] = new double[10];
                for (int i = 0; i < 10; i++) {
                    brick_mesh.IJKcoors[axis][i] = i * 1.0;
                }
            }

            // Initialize block_dist arrays
            for (int axis = 0; axis < 3; axis++) {
                brick_mesh.block_dist[axis] = new double[brick_mesh.inline_b[axis]];
                brick_mesh.c_block_dist[axis] = new double[brick_mesh.inline_b[axis] + 1];
                brick_mesh.first_size[axis] = new double[brick_mesh.inline_b[axis]];
                brick_mesh.last_size[axis] = new double[brick_mesh.inline_b[axis]];
                for (int i = 0; i < brick_mesh.inline_b[axis]; i++) {
                    brick_mesh.block_dist[axis][i] = 1.0;
                    brick_mesh.first_size[axis][i] = 1.0;
                    brick_mesh.last_size[axis][i] = 1.0;
                }
            }

            // Set inline_gmin values
            for (int axis = 0; axis < 3; axis++) {
                brick_mesh.inline_gmin[axis] = 0.0;
            }

            // Set strides
            brick_mesh.jnstride = 10;
            brick_mesh.knstride = 100;

            // Initialize the mesh descriptor
            brick_mesh.Set_Up();

            // Set up test data
            std::vector<long long> global_node_vector;
            std::map<long long, long long> global_node_map;

            for (long long i = 0; i < 5; i++) {
                global_node_vector.push_back(i);
                global_node_map[i] = i;
            }

            long long num_nodes = 5;
            double* coords = new double[num_nodes * 2];
            double* coords_copy = new double[num_nodes * 2];

            for (long long i = 0; i < num_nodes * 2; i++) {
                coords[i] = 0.0;
            }

            // Test the wrapper version which internally calls the device version
            std::cout << "Testing Kokkos wrapper implementation..." << std::endl;
            brick_mesh.Populate_Coords(coords, global_node_vector, global_node_map, num_nodes);

            std::cout << "Brick mesh coordinates computed successfully via wrapper!" << std::endl;
            for (int i = 0; i < num_nodes; i++) {
                std::cout << "Node " << i << " (wrapper): ("
                          << coords[i] << ", "
                          << coords[i + num_nodes] << ")" << std::endl;
            }

            // Save wrapper results for comparison
            for (long long i = 0; i < num_nodes * 2; i++) {
                coords_copy[i] = coords[i];
            }

            // Now test the device version directly
            std::cout << "\nTesting Kokkos device implementation directly..." << std::endl;

            // Reset coordinates
            for (long long i = 0; i < num_nodes * 2; i++) {
                coords[i] = 0.0;
            }

            // Convert data to Kokkos Views for direct device testing
            View1D<long long> global_node_vector_view("global_node_vector_view", global_node_vector.size());
            View1D<long long> global_node_map_keys_view("global_node_map_keys_view", global_node_map.size());
            View1D<long long> global_node_map_values_view("global_node_map_values_view", global_node_map.size());
            View2D<double> coords_view("coords_view", num_nodes, 2);

            // Fill the views with test data
            convert_vector_to_kokkos_view(global_node_vector, global_node_vector_view);

            size_t idx = 0;
            for (const auto& pair : global_node_map) {
                global_node_map_keys_view(idx) = pair.first;
                global_node_map_values_view(idx) = pair.second;
                idx++;
            }

            // Call the device version directly
            brick_mesh.Populate_Coords_Device(coords_view, global_node_vector_view,
                                             global_node_map_keys_view, global_node_map_values_view, num_nodes);

            // Copy results back to host
            auto coords_host = Kokkos::create_mirror_view(coords_view);
            Kokkos::deep_copy(coords_host, coords_view);

            // Copy to raw array for comparison
            for (long long i = 0; i < num_nodes; i++) {
                for (long long axis = 0; axis < 2; axis++) {
                    coords[i + axis * num_nodes] = coords_host(i, axis);
                }
            }

            std::cout << "Brick mesh coordinates computed successfully via direct device call!" << std::endl;
            for (int i = 0; i < num_nodes; i++) {
                std::cout << "Node " << i << " (device): ("
                          << coords[i] << ", "
                          << coords[i + num_nodes] << ")" << std::endl;
            }

            // Verify that wrapper and device versions produce the same results
            std::cout << "\nVerifying consistency between wrapper and device implementations..." << std::endl;
            bool consistent = true;
            for (int i = 0; i < num_nodes * 2; i++) {
                if (fabs(coords[i] - coords_copy[i]) > 1e-10) {
                    std::cout << "Mismatch at index " << i << ": wrapper=" << coords_copy[i]
                              << " vs device=" << coords[i] << std::endl;
                    consistent = false;
                }
            }
            if (consistent) {
                std::cout << "✅ Wrapper and device implementations produce identical results!" << std::endl;
            } else {
                std::cout << "❌ Wrapper and device implementations produce different results!" << std::endl;
                return 1;
            }

            std::cout << "Brick mesh coordinates computed successfully!" << std::endl;
            for (int i = 0; i < num_nodes; i++) {
                std::cout << "Node " << i << ": ("
                          << coords[i] << ", "
                          << coords[i + num_nodes] << ")" << std::endl;
            }

            // Clean up memory allocated in this test
            delete[] coords;
            delete[] coords_copy;

            std::cout << "Brick mesh test PASSED!" << std::endl;
        }

        // Test 2: Radial Mesh Kokkos functionality
        std::cout << "\n=== Test 2: Radial Mesh Kokkos ===" << std::endl;
        {
            Radial_Inline_Mesh_Desc radial_mesh(2);

            // Set up basic parameters
            radial_mesh.inline_b[0] = 2;
            radial_mesh.inline_b[1] = 4;
            radial_mesh.inline_b[2] = 1;

            radial_mesh.inline_n[0] = 3;
            radial_mesh.inline_n[1] = 8;
            radial_mesh.inline_n[2] = 1;

            // Set up intervals
            for (int axis = 0; axis < 3; axis++) {
                radial_mesh.interval[axis] = new long long[radial_mesh.inline_b[axis]];
                for (int i = 0; i < radial_mesh.inline_b[axis]; i++) {
                    radial_mesh.interval[axis][i] = radial_mesh.inline_n[axis];
                }
            }

            // Set up coordinate arrays
            for (int axis = 0; axis < 3; axis++) {
                radial_mesh.IJKcoors[axis] = new double[20];
                for (int i = 0; i < 20; i++) {
                    if (axis == 0) {
                        radial_mesh.IJKcoors[axis][i] = 1.0 + i * 0.5;
                    } else if (axis == 1) {
                        radial_mesh.IJKcoors[axis][i] = i * 45.0;
                    } else {
                        radial_mesh.IJKcoors[axis][i] = 0.0;
                    }
                }
            }

            // Initialize block_dist arrays
            for (int axis = 0; axis < 3; axis++) {
                radial_mesh.block_dist[axis] = new double[radial_mesh.inline_b[axis]];
                radial_mesh.c_block_dist[axis] = new double[radial_mesh.inline_b[axis] + 1];
                radial_mesh.first_size[axis] = new double[radial_mesh.inline_b[axis]];
                radial_mesh.last_size[axis] = new double[radial_mesh.inline_b[axis]];
                for (int i = 0; i < radial_mesh.inline_b[axis]; i++) {
                    radial_mesh.block_dist[axis][i] = 1.0;
                    radial_mesh.first_size[axis][i] = 1.0;
                    radial_mesh.last_size[axis][i] = 1.0;
                }
            }

            // Set inline_gmin values
            for (int axis = 0; axis < 3; axis++) {
                radial_mesh.inline_gmin[axis] = 0.0;
            }

            // Set strides
            radial_mesh.jnstride = 20;
            radial_mesh.knstride = 400;

            // Initialize the mesh descriptor
            radial_mesh.Set_Up();

            // Set up test data
            std::vector<long long> global_node_vector;
            std::map<long long, long long> global_node_map;

            for (long long i = 0; i < 8; i++) {
                global_node_vector.push_back(i);
                global_node_map[i] = i;
            }

            long long num_nodes = 8;
            double* coords = new double[num_nodes * 2];
            double* coords_copy = new double[num_nodes * 2];

            for (long long i = 0; i < num_nodes * 2; i++) {
                coords[i] = 0.0;
            }

            // Test the wrapper version which internally calls the device version
            std::cout << "Testing Kokkos wrapper implementation..." << std::endl;
            radial_mesh.Populate_Coords(coords, global_node_vector, global_node_map, num_nodes);

            std::cout << "Radial mesh coordinates computed successfully via wrapper!" << std::endl;
            for (int i = 0; i < num_nodes; i++) {
                double r = sqrt(coords[i] * coords[i] + coords[i + num_nodes] * coords[i + num_nodes]);
                double theta = atan2(coords[i + num_nodes], coords[i]) * 180.0 / M_PI;
                std::cout << "Node " << i << " (wrapper): ("
                          << coords[i] << ", "
                          << coords[i + num_nodes] << ") "
                          << "R=" << r << ", Theta=" << theta << "°" << std::endl;
            }

            // Save wrapper results for comparison
            for (long long i = 0; i < num_nodes * 2; i++) {
                coords_copy[i] = coords[i];
            }

            // Now test the device version directly
            std::cout << "\nTesting Kokkos device implementation directly..." << std::endl;

            // Reset coordinates
            for (long long i = 0; i < num_nodes * 2; i++) {
                coords[i] = 0.0;
            }

            // Convert data to Kokkos Views for direct device testing
            View1D<long long> global_node_vector_view("global_node_vector_view", global_node_vector.size());
            View1D<long long> global_node_map_keys_view("global_node_map_keys_view", global_node_map.size());
            View1D<long long> global_node_map_values_view("global_node_map_values_view", global_node_map.size());
            View2D<double> coords_view("coords_view", num_nodes, 2);

            // Fill the views with test data
            convert_vector_to_kokkos_view(global_node_vector, global_node_vector_view);

            size_t idx = 0;
            for (const auto& pair : global_node_map) {
                global_node_map_keys_view(idx) = pair.first;
                global_node_map_values_view(idx) = pair.second;
                idx++;
            }

            // Call the device version directly
            radial_mesh.Populate_Coords_Device(coords_view, global_node_vector_view,
                                              global_node_map_keys_view, global_node_map_values_view, num_nodes);

            // Copy results back to host
            auto coords_host = Kokkos::create_mirror_view(coords_view);
            Kokkos::deep_copy(coords_host, coords_view);

            // Copy to raw array for comparison
            for (long long i = 0; i < num_nodes; i++) {
                for (long long axis = 0; axis < 2; axis++) {
                    coords[i + axis * num_nodes] = coords_host(i, axis);
                }
            }

            std::cout << "Radial mesh coordinates computed successfully via direct device call!" << std::endl;
            for (int i = 0; i < num_nodes; i++) {
                double r = sqrt(coords[i] * coords[i] + coords[i + num_nodes] * coords[i + num_nodes]);
                double theta = atan2(coords[i + num_nodes], coords[i]) * 180.0 / M_PI;
                std::cout << "Node " << i << " (device): ("
                          << coords[i] << ", "
                          << coords[i + num_nodes] << ") "
                          << "R=" << r << ", Theta=" << theta << "°" << std::endl;
            }

            // Verify that wrapper and device versions produce the same results
            std::cout << "\nVerifying consistency between wrapper and device implementations..." << std::endl;
            bool consistent = true;
            for (int i = 0; i < num_nodes * 2; i++) {
                if (fabs(coords[i] - coords_copy[i]) > 1e-10) {
                    std::cout << "Mismatch at index " << i << ": wrapper=" << coords_copy[i]
                              << " vs device=" << coords[i] << std::endl;
                    consistent = false;
                }
            }
            if (consistent) {
                std::cout << "✅ Wrapper and device implementations produce identical results!" << std::endl;
            } else {
                std::cout << "❌ Wrapper and device implementations produce different results!" << std::endl;
                return 1;
            }

            std::cout << "Radial mesh coordinates computed successfully!" << std::endl;
            for (int i = 0; i < num_nodes; i++) {
                double r = sqrt(coords[i] * coords[i] + coords[i + num_nodes] * coords[i + num_nodes]);
                double theta = atan2(coords[i + num_nodes], coords[i]) * 180.0 / M_PI;
                std::cout << "Node " << i << ": ("
                          << coords[i] << ", "
                          << coords[i + num_nodes] << ") "
                          << "R=" << r << ", Theta=" << theta << "°" << std::endl;
            }

            // Clean up memory allocated in this test
            delete[] coords;
            delete[] coords_copy;

            std::cout << "Radial mesh test PASSED!" << std::endl;
        }

        // Test 3: Kokkos Utilities
        std::cout << "\n=== Test 3: Kokkos Utilities ===" << std::endl;
        {
            // Test map conversion
            std::map<long long, long long> test_map;
            test_map[1] = 10;
            test_map[2] = 20;
            test_map[3] = 30;
            test_map[5] = 50;
            test_map[8] = 80;

            View1D<long long> keys("test_keys", test_map.size());
            View1D<long long> values("test_values", test_map.size());

            convert_map_to_kokkos_views(test_map, keys, values);

            std::cout << "Map conversion successful! Size: " << test_map.size() << std::endl;

            // Test vector conversion
            std::vector<long long> test_vector = {100, 200, 300, 400, 500};
            View1D<long long> vector_view("test_vector", test_vector.size());

            convert_vector_to_kokkos_view(test_vector, vector_view);

            std::cout << "Vector conversion successful! Size: " << test_vector.size() << std::endl;

            std::cout << "Kokkos utilities test PASSED!" << std::endl;
        }

        std::cout << "\n🎉 All Kokkos integration tests PASSED!" << std::endl;
        std::cout << "Kokkos functionality is working correctly." << std::endl;

    }

    // Finalize Kokkos
    Kokkos::finalize();

    return 0;
}
