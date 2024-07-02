// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_UTILITIES_HPP_
#define _COMPADRE_UTILITIES_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"

namespace Compadre {

KOKKOS_INLINE_FUNCTION
void getMidpointFromCellVertices(const member_type& teamMember, scratch_vector_type midpoint_storage, scratch_matrix_right_type cell_coordinates, const int cell_num, const int dim=3) {
Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
    auto num_nodes = cell_coordinates.extent(1)/dim;
    for (int j=0; j<dim; ++j) {
        midpoint_storage(j) = 0;
        for (size_t i=0; i<num_nodes; ++i) {
            midpoint_storage(j) += cell_coordinates(cell_num, i*dim + j) / (double)(num_nodes);
        }
    }
});
}

template <typename view_type_1, typename view_type_2>
KOKKOS_INLINE_FUNCTION
double getAreaFromVectors(const member_type& teamMember, view_type_1 v1, view_type_2 v2) {
    if (v1.extent(0)==3) {
        double area = 0;
        double val = v1[1]*v2[2] - v1[2]*v2[1];
        area += val*val;
        val = v1[2]*v2[0] - v1[0]*v2[2];
        area += val*val;
        val = v1[0]*v2[1] - v1[1]*v2[0];
        area += val*val;
        return std::sqrt(area);
    } else if (v1.extent(0)==2) {
        double area = 0;
        double val = v1[0]*v2[1] - v1[1]*v2[0];
        area += val*val;
        return std::sqrt(area);
    } else {
        compadre_kernel_assert_debug(false && "v1 in getAreaFromVectors has length != 2 or 3");
        return 0.0;
    }
}

template <typename output_memory_space, typename view_type_input_data, typename output_array_layout = typename view_type_input_data::array_layout, typename index_type=int>
Kokkos::View<int*, output_array_layout, output_memory_space> // shares layout of input by default
        filterViewByID(view_type_input_data input_data_host_or_device, index_type filtered_value) {

    // Make view on the host (does nothing if already on the host)
    auto input_data_host = Kokkos::create_mirror_view(input_data_host_or_device);
    Kokkos::deep_copy(input_data_host, input_data_host_or_device);
    Kokkos::fence();

    // Count the number of elements in the input view that match the desired value
    int num_count = 0;
    auto this_filtered_value = filtered_value;

    // call device functor here

    for (size_t i=0; i<input_data_host.extent(0); i++) {
        if (input_data_host(i) == this_filtered_value) {
             num_count++;
        }
    }
    Kokkos::fence();

    // Create a new view living on device
    Kokkos::View<int*, output_array_layout, output_memory_space> filtered_view("filterd view", num_count);
    // Gather up the indices into the new view
    int filtered_index = 0;
    for (size_t i=0; i<input_data_host.extent(0); i++) {
         if (input_data_host(i) == this_filtered_value) {
             filtered_view(filtered_index) = i;
             filtered_index++;
         }
    }
    Kokkos::fence();

    // Then copy it back out - either to host or device space based on user's request
    typedef Kokkos::View<int*, output_array_layout, output_memory_space> output_view_type;
    output_view_type filtered_view_output("output filtered view", num_count);
    Kokkos::deep_copy(filtered_view_output, filtered_view);
    Kokkos::fence();

    return filtered_view_output;
}

struct Extract {

    template <typename output_memory_space, typename view_type_input_data, typename view_type_index_data,
              enable_if_t<std::is_same<typename view_type_input_data::data_type, double**>::value
                  ||std::is_same<typename view_type_input_data::data_type, int**>::value, int> = 0>
    static Kokkos::View<typename view_type_input_data::data_type, typename view_type_input_data::array_layout, output_memory_space> // shares layout of input by default
            extractViewByIndex(view_type_input_data input_data_host_or_device, view_type_index_data index_data_host_or_device) {

        typedef typename view_type_input_data::data_type    output_data_type;
        typedef typename view_type_input_data::array_layout output_array_layout;
    
        // Make view on the host for input data (does nothing if already on the host)
        auto input_data_host = Kokkos::create_mirror_view(input_data_host_or_device);
        Kokkos::deep_copy(input_data_host, input_data_host_or_device);
        Kokkos::fence();
    
        // Make view on the host for index data (does nothing if already on the host)
        auto index_data_host = Kokkos::create_mirror_view(index_data_host_or_device);
        Kokkos::deep_copy(index_data_host, index_data_host_or_device);
        Kokkos::fence();
    
        // Create a new view to extract out the rows that belong to the filtered index
        Kokkos::View<output_data_type, output_array_layout, output_memory_space> extracted_view("extracted view", 
                index_data_host.extent(0), input_data_host.extent(1));
    
        // Loop through all the entries of index data
        for (size_t i=0; i<index_data_host.extent(0); i++) {
            for (size_t j=0; j<input_data_host.extent(1); j++) {
                extracted_view(i, j) = input_data_host(index_data_host(i), j);
            }
        }
    
        // Then copy it back out - either to host or device space based on user's request
        typedef Kokkos::View<output_data_type, output_array_layout, output_memory_space> output_view_type;
        output_view_type extracted_view_output("output extracted view", extracted_view.extent(0), extracted_view.extent(1));
        Kokkos::deep_copy(extracted_view_output, extracted_view);
        Kokkos::fence();
    
        return extracted_view_output;
    }

    template <typename output_memory_space, typename view_type_input_data, typename view_type_index_data,
              enable_if_t<std::is_same<typename view_type_input_data::data_type, double*>::value
                  ||std::is_same<typename view_type_input_data::data_type, int*>::value, int> = 0>
    static Kokkos::View<double*, typename view_type_input_data::array_layout, output_memory_space> // shares layout of input by default
            extractViewByIndex(view_type_input_data input_data_host_or_device, view_type_index_data index_data_host_or_device) {

        typedef typename view_type_input_data::data_type    output_data_type;
        typedef typename view_type_input_data::array_layout output_array_layout;
    
        // Make view on the host for input data (does nothing if already on the host)
        auto input_data_host = Kokkos::create_mirror_view(input_data_host_or_device);
        Kokkos::deep_copy(input_data_host, input_data_host_or_device);
        Kokkos::fence();
    
        // Make view on the host for index data (does nothing if already on the host)
        auto index_data_host = Kokkos::create_mirror_view(index_data_host_or_device);
        Kokkos::deep_copy(index_data_host, index_data_host_or_device);
        Kokkos::fence();
    
        // Create a new view to extract out the rows that belong to the filtered index
        Kokkos::View<output_data_type, output_array_layout, output_memory_space> extracted_view("extracted view", 
                index_data_host.extent(0));
    
        // Loop through all the entries of index data
        for (size_t i=0; i<index_data_host.extent(0); i++) {
            extracted_view(i) = input_data_host(index_data_host(i));
        }
    
        // Then copy it back out - either to host or device space based on user's request
        typedef Kokkos::View<output_data_type, output_array_layout, output_memory_space> output_view_type;
        output_view_type extracted_view_output("output extracted view", extracted_view.extent(0));
        Kokkos::deep_copy(extracted_view_output, extracted_view);
        Kokkos::fence();
    
        return extracted_view_output;
    }

};

// template <typename output_memory_space, typename view_type_input_data, typename output_array_layout = typename view_type_input_data::array_layout, typename index_type=int>
// Kokkos::View<int*, output_array_layout, output_memory_space> // shares layout of input by default
//         filterViewByID(view_type_input_data input_data_host_or_device, index_type filtered_value) {
//
//     // Make view on the device (does nothing if already on the device)
//     auto input_data_device = Kokkos::create_mirror_view(
//         device_memory_space(), input_data_host_or_device);
//     Kokkos::deep_copy(input_data_device, input_data_host_or_device);
//     Kokkos::fence();
//
//     // Count the number of elements in the input view that match the desired value
//     int num_count = 0;
//     auto this_filtered_value = filtered_value;
//
//     // call device functor here
//
//     for (int i=0; i<input_data_device.extent(0); i++) {
//         input_data_device(i)++;
//         // if (input_data_device(i) == this_filtered_value) {
//         // if (input_data_device(i) == 1) {
//         //      num_count++;
//         // }
//     }
//     Kokkos::fence();
//
//     // Create a new view living on device
//     Kokkos::View<int*, output_array_layout> filtered_view("filterd view", num_count);
//     // // Gather up the indices into the new view
//     // int filtered_index = 0;
//     // for (int i=0; i<input_data_device.extent(0); i++) {
//     //     if (input_data_device(i) == filtered_value) {
//     //         filtered_view(filtered_index) = i;
//     //         filtered_index++;
//     //     }
//     // }
//     // Kokkos::fence();
//
//     // Then copy it back out - either to host or device space based on user's request
//     typedef Kokkos::View<int*, output_array_layout, output_memory_space> output_view_type;
//     output_view_type filtered_view_output("output filtered view", num_count);
//     Kokkos::deep_copy(filtered_view_output, filtered_view);
//     Kokkos::fence();
//
//     return filtered_view_output;
// }

}  // Compadre namespace

#endif
