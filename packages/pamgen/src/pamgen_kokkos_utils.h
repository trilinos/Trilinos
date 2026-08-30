// @HEADER
// ***************************************************************************
//                     Pamgen Package - Kokkos Utilities
//
// Copyright 2026 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// ***************************************************************************
// @HEADER

#ifndef PAMGEN_KOKKOS_UTILS_H
#define PAMGEN_KOKKOS_UTILS_H

#include <Kokkos_Core.hpp>
#include <map>
#include <vector>

namespace PAMGEN_NEVADA {

// Define common Kokkos types for PAMGEN
using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;
using HostSpace = Kokkos::HostSpace;
using LayoutLeft = Kokkos::LayoutLeft;
using LayoutRight = Kokkos::LayoutRight;

// Common View types
template<typename T>
using View1D = Kokkos::View<T*, DeviceSpace>;

template<typename T>
using View2D = Kokkos::View<T**, LayoutLeft, DeviceSpace>;

template<typename T>
using HostView1D = Kokkos::View<T*, HostSpace>;

template<typename T>
using HostView2D = Kokkos::View<T**, LayoutLeft, HostSpace>;

// Device function to find entry in map represented as parallel arrays
KOKKOS_INLINE_FUNCTION
long long get_map_entry_device(const View1D<long long>& map_keys,
                              const View1D<long long>& map_values,
                              long long key,
                              long long map_size) {
    // Simple linear search for now - can be optimized later
    for (long long i = 0; i < map_size; ++i) {
        if (map_keys(i) == key) {
            return map_values(i);
        }
    }
    return -1; // Key not found
}

// Host function to convert std::map to Kokkos Views
void convert_map_to_kokkos_views(const std::map<long long, long long>& input_map,
                                View1D<long long>& keys,
                                View1D<long long>& values);

// Host function to convert std::vector to Kokkos View
void convert_vector_to_kokkos_view(const std::vector<long long>& input_vector,
                                  View1D<long long>& output_view);

// Host function to convert raw pointer to Kokkos View
template<typename T>
void convert_raw_pointer_to_kokkos_view(T* raw_ptr,
                                     HostView2D<T>& kokkos_view,
                                     long long rows,
                                     long long cols) {
    // Create host view that wraps the raw pointer
    kokkos_view = HostView2D<T>(raw_ptr, rows, cols);
}

} // namespace PAMGEN_NEVADA

#endif // PAMGEN_KOKKOS_UTILS_H