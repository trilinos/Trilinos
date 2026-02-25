// @HEADER
// ***************************************************************************
//                     Pamgen Package - Kokkos Utilities Implementation
//
// Copyright 2026 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// ***************************************************************************
// @HEADER

#include "pamgen_kokkos_utils.h"

namespace PAMGEN_NEVADA {

void convert_map_to_kokkos_views(const std::map<long long, long long>& input_map,
                                View1D<long long>& keys,
                                View1D<long long>& values) {
    // Create host mirrors first
    long long map_size = input_map.size();
    auto keys_host = Kokkos::create_mirror_view(keys);
    auto values_host = Kokkos::create_mirror_view(values);

    // Resize if needed
    if (keys_host.size() < map_size) {
        keys_host = View1D<long long>(Kokkos::view_alloc("map_keys_host", Kokkos::WithoutInitializing), map_size);
        values_host = View1D<long long>(Kokkos::view_alloc("map_values_host", Kokkos::WithoutInitializing), map_size);
    }

    // Fill the host views
    long long idx = 0;
    for (const auto& pair : input_map) {
        keys_host(idx) = pair.first;
        values_host(idx) = pair.second;
        idx++;
    }

    // Deep copy to device
    Kokkos::deep_copy(keys, keys_host);
    Kokkos::deep_copy(values, values_host);
}

void convert_vector_to_kokkos_view(const std::vector<long long>& input_vector,
                                  View1D<long long>& output_view) {
    // Create host mirror
    auto output_host = Kokkos::create_mirror_view(output_view);

    // Resize if needed
    if (output_host.size() < input_vector.size()) {
        output_host = View1D<long long>(Kokkos::view_alloc("vector_view_host", Kokkos::WithoutInitializing), input_vector.size());
    }

    // Copy data
    for (size_t i = 0; i < input_vector.size(); ++i) {
        output_host(i) = input_vector[i];
    }

    // Deep copy to device
    Kokkos::deep_copy(output_view, output_host);
}

} // namespace PAMGEN_NEVADA