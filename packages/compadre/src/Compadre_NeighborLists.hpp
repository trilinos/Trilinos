// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_NEIGHBORLISTS_HPP_
#define _COMPADRE_NEIGHBORLISTS_HPP_

#include "Compadre_Typedefs.hpp"
#include <Kokkos_Core.hpp>

namespace Compadre {

//!  NeighborLists assists in accessing entries of compressed row neighborhood lists
template <typename view_type>
struct NeighborLists {

    typedef view_type internal_view_type;
    typedef Kokkos::View<global_index_type*, typename view_type::array_layout, 
            typename view_type::memory_space, typename view_type::memory_traits> internal_row_offsets_view_type;

    int _max_neighbor_list_row_storage_size;
    int _min_neighbor_list_row_storage_size;
    bool _needs_sync_to_host;
    int _number_of_targets;

    internal_row_offsets_view_type _row_offsets;
    view_type _cr_neighbor_lists;
    view_type _number_of_neighbors_list;

    typename internal_row_offsets_view_type::HostMirror _host_row_offsets;
    typename view_type::HostMirror _host_cr_neighbor_lists;
    typename view_type::HostMirror _host_number_of_neighbors_list;

/** @name Constructors
 *  Ways to initialize a NeighborLists object
 */
///@{

    //! \brief Constructor for the purpose of classes who have NeighborLists as a member object
    NeighborLists() {
        _max_neighbor_list_row_storage_size = -1;
        _min_neighbor_list_row_storage_size = -1;
        _needs_sync_to_host = true;
        _number_of_targets = 0;
    }

    /*! \brief Constructor for when compressed row `cr_neighbor_lists` is preallocated/populated, 
     *  `number_of_neighbors_list` and `neighbor_lists_row_offsets` have already been populated.
     */
    NeighborLists(view_type cr_neighbor_lists, view_type number_of_neighbors_list, 
            internal_row_offsets_view_type neighbor_lists_row_offsets, bool compute_max = true) {
        compadre_assert_release((view_type::rank==1) && 
                "cr_neighbor_lists and number_neighbors_list and neighbor_lists_row_offsets must be a 1D Kokkos view.");

        _number_of_targets = number_of_neighbors_list.extent(0);
        _number_of_neighbors_list = number_of_neighbors_list;
        _cr_neighbor_lists = cr_neighbor_lists;
        _row_offsets = neighbor_lists_row_offsets;

        _host_cr_neighbor_lists = Kokkos::create_mirror_view(_cr_neighbor_lists);
        _host_number_of_neighbors_list = Kokkos::create_mirror_view(_number_of_neighbors_list);
        _host_row_offsets = Kokkos::create_mirror_view(_row_offsets);

        Kokkos::deep_copy(_host_cr_neighbor_lists, _cr_neighbor_lists);
        Kokkos::deep_copy(_host_number_of_neighbors_list, _number_of_neighbors_list);
        Kokkos::deep_copy(_host_row_offsets, _row_offsets);
        Kokkos::fence();

        if (compute_max) {
            computeMaxNumNeighbors();
        } else {
            _max_neighbor_list_row_storage_size = -1;
        }
        _min_neighbor_list_row_storage_size = -1;

        //check neighbor_lists is large enough
        compadre_assert_release(((size_t)(this->getTotalNeighborsOverAllListsHost())<=cr_neighbor_lists.extent(0)) 
                && "neighbor_lists is not large enough to store all neighbors.");

        _needs_sync_to_host = false;
    }

    /*! \brief Constructor for when compressed row `cr_neighbor_lists` is preallocated/populated, 
     *  and `number_of_neighbors_list` is already populated, and row offsets still need to be computed.
     */
    NeighborLists(view_type cr_neighbor_lists, view_type number_of_neighbors_list) {
        compadre_assert_release((view_type::rank==1) 
                && "cr_neighbor_lists and number_neighbors_list must be a 1D Kokkos view.");

        _number_of_targets = number_of_neighbors_list.extent(0);

        _row_offsets = internal_row_offsets_view_type("row offsets", number_of_neighbors_list.extent(0));
        _number_of_neighbors_list = number_of_neighbors_list;
        _cr_neighbor_lists = cr_neighbor_lists;

        _host_cr_neighbor_lists = Kokkos::create_mirror_view(_cr_neighbor_lists);
        _host_number_of_neighbors_list = Kokkos::create_mirror_view(_number_of_neighbors_list);
        _host_row_offsets = Kokkos::create_mirror_view(_row_offsets);

        Kokkos::deep_copy(_host_cr_neighbor_lists, _cr_neighbor_lists);
        Kokkos::deep_copy(_host_number_of_neighbors_list, _number_of_neighbors_list);
        Kokkos::deep_copy(_host_row_offsets, _row_offsets);
        Kokkos::fence();

        computeRowOffsets();
        computeMaxNumNeighbors();

        //check neighbor_lists is large enough
        compadre_assert_release(((size_t)(this->getTotalNeighborsOverAllListsHost())<=cr_neighbor_lists.extent(0)) 
                && "neighbor_lists is not large enough to store all neighbors.");

        _needs_sync_to_host = false;
    }

    /*! \brief Constructor for when `number_of_neighbors_list` is already populated.
     *  Will allocate space for compressed row neighbor lists data, and will allocate then
     *  populate information for row offsets.
     */
    NeighborLists(view_type number_of_neighbors_list) {
        compadre_assert_release((view_type::rank==1) 
                && "cr_neighbor_lists and number_neighbors_list must be a 1D Kokkos view.");

        _number_of_targets = number_of_neighbors_list.extent(0);

        _row_offsets = internal_row_offsets_view_type("row offsets", number_of_neighbors_list.extent(0));
        _number_of_neighbors_list = number_of_neighbors_list;

        _host_number_of_neighbors_list = Kokkos::create_mirror_view(_number_of_neighbors_list);
        _host_row_offsets = Kokkos::create_mirror_view(_row_offsets);

        Kokkos::deep_copy(_host_number_of_neighbors_list, _number_of_neighbors_list);
        Kokkos::deep_copy(_host_row_offsets, _row_offsets);
        Kokkos::fence();

        computeRowOffsets();
        computeMaxNumNeighbors();

        _cr_neighbor_lists = view_type("compressed row neighbor lists data", this->getTotalNeighborsOverAllListsHost());
        _host_cr_neighbor_lists = Kokkos::create_mirror_view(_cr_neighbor_lists);
        Kokkos::deep_copy(_host_cr_neighbor_lists, _cr_neighbor_lists);
        Kokkos::fence();

        _needs_sync_to_host = false;
    }
///@}

/** @name Public modifiers
 */
///@{

    //! Setter function for N(i,j) indexing where N(i,j) is the index of the jth neighbor of i
    KOKKOS_INLINE_FUNCTION
    void setNeighborDevice(int target_index, int neighbor_num, int new_value) {
        _cr_neighbor_lists(_row_offsets(target_index)+neighbor_num) = new_value;
        // indicate that host view is now out of sync with device
        // but only in debug mode (notice the next line is both setting the variable and checking it was set)
        compadre_assert_debug((_needs_sync_to_host=true)==true); 
    }

    //! Calculate the maximum number of neighbors of all targets' neighborhoods (host)
    void computeMaxNumNeighbors() {
        if (_number_of_neighbors_list.extent(0)==0) {
            _max_neighbor_list_row_storage_size = 0;
        } else {
            auto number_of_neighbors_list = _number_of_neighbors_list;
            Kokkos::parallel_reduce("max number of neighbors", 
                    Kokkos::RangePolicy<typename view_type::execution_space>(0, _number_of_neighbors_list.extent(0)), 
                    KOKKOS_LAMBDA(const int i, int& t_max_num_neighbors) {
                t_max_num_neighbors = (number_of_neighbors_list(i) > t_max_num_neighbors) ? number_of_neighbors_list(i) : t_max_num_neighbors;
            }, Kokkos::Max<int>(_max_neighbor_list_row_storage_size));
            Kokkos::fence();
        }
    }

    //! Calculate the minimum number of neighbors of all targets' neighborhoods (host)
    void computeMinNumNeighbors() {
        if (_number_of_neighbors_list.extent(0)==0) {
            _min_neighbor_list_row_storage_size = 0;
        } else {
            auto number_of_neighbors_list = _number_of_neighbors_list;
            Kokkos::parallel_reduce("min number of neighbors", 
                    Kokkos::RangePolicy<typename view_type::execution_space>(0, _number_of_neighbors_list.extent(0)), 
                    KOKKOS_LAMBDA(const int i, int& t_min_num_neighbors) {
                t_min_num_neighbors = (number_of_neighbors_list(i) < t_min_num_neighbors) ? number_of_neighbors_list(i) : t_min_num_neighbors;
            }, Kokkos::Min<int>(_min_neighbor_list_row_storage_size));
            Kokkos::fence();
        }
    }

    //! Calculate the row offsets for each target's neighborhood (host)
    void computeRowOffsets() {
        auto number_of_neighbors_list = _number_of_neighbors_list;
        auto row_offsets = _row_offsets;
        Kokkos::parallel_scan("number of neighbors offsets", 
                Kokkos::RangePolicy<typename view_type::execution_space>(0, _number_of_neighbors_list.extent(0)), 
                KOKKOS_LAMBDA(const int i, global_index_type& lsum, bool final) {
            row_offsets(i) = lsum;
            lsum += number_of_neighbors_list(i);
        });
        Kokkos::deep_copy(_host_row_offsets, _row_offsets);
        Kokkos::fence();
    }

    //! Sync the host from the device (copy device data to host)
    void copyDeviceDataToHost() {
        Kokkos::deep_copy(_host_cr_neighbor_lists, _cr_neighbor_lists);
        Kokkos::fence();
        _needs_sync_to_host = false;
    }

    //! Device view into neighbor lists data (use with caution)
    view_type getNeighborLists() const {
        return _cr_neighbor_lists;
    }

    //! Device view into number of neighbors list (use with caution)
    view_type getNumberOfNeighborsList() const {
        return _number_of_neighbors_list;
    }

///@}
/** @name Public accessors
 */
///@{
    //! Get number of total targets having neighborhoods (host/device).
    KOKKOS_INLINE_FUNCTION
    int getNumberOfTargets() const {
        return _number_of_targets;
    }

    //! Get number of neighbors for a given target (host)
    int getNumberOfNeighborsHost(int target_index) const {
        compadre_assert_extreme_debug(target_index < this->getNumberOfTargets());
        return _host_number_of_neighbors_list(target_index);
    }

    //! Get number of neighbors for a given target (device)
    KOKKOS_INLINE_FUNCTION
    int getNumberOfNeighborsDevice(int target_index) const {
        compadre_kernel_assert_extreme_debug(target_index < this->getNumberOfTargets());
        return _number_of_neighbors_list(target_index);
    }

    //! Get offset into compressed row neighbor lists (host)
    global_index_type getRowOffsetHost(int target_index) const {
        compadre_assert_extreme_debug(target_index < this->getNumberOfTargets());
        return _host_row_offsets(target_index);
    }

    //! Get offset into compressed row neighbor lists (device)
    KOKKOS_INLINE_FUNCTION
    global_index_type getRowOffsetDevice(int target_index) const {
        compadre_kernel_assert_extreme_debug(target_index < this->getNumberOfTargets());
        return _row_offsets(target_index);
    }

    //! Offers N(i,j) indexing where N(i,j) is the index of the jth neighbor of i (host)
    int getNeighborHost(int target_index, int neighbor_num) const {
        compadre_assert_extreme_debug(target_index < this->getNumberOfTargets());
        compadre_assert_debug((!_needs_sync_to_host) 
                && "Stale information in host_cr_neighbor_lists. Call CopyDeviceDataToHost() to refresh.");
        compadre_assert_debug((neighbor_num<_host_number_of_neighbors_list(target_index))
                && "neighor_num exceeds number of neighbors for this target_index.");
        return _host_cr_neighbor_lists(_host_row_offsets(target_index)+neighbor_num);
    }

    //! Offers N(i,j) indexing where N(i,j) is the index of the jth neighbor of i (device)
    KOKKOS_INLINE_FUNCTION
    int getNeighborDevice(int target_index, int neighbor_num) const {
        compadre_kernel_assert_extreme_debug((neighbor_num<_number_of_neighbors_list(target_index))
                && "neighor_num exceeds number of neighbors for this target_index.");
        compadre_kernel_assert_extreme_debug(target_index < this->getNumberOfTargets());
        return _cr_neighbor_lists(_row_offsets(target_index)+neighbor_num);
    }

    //! Get the maximum number of neighbors of all targets' neighborhoods (host/device)
    KOKKOS_INLINE_FUNCTION
    int getMaxNumNeighbors() const {
        compadre_kernel_assert_debug((_max_neighbor_list_row_storage_size > -1) && "getMaxNumNeighbors() called but maximum never calculated.");
        return _max_neighbor_list_row_storage_size;
    }

    //! Get the minimum number of neighbors of all targets' neighborhoods (host/device)
    KOKKOS_INLINE_FUNCTION
    int getMinNumNeighbors() const {
        compadre_kernel_assert_debug((_min_neighbor_list_row_storage_size > -1) && "getMinNumNeighbors() called but minimum never calculated.");
        return _min_neighbor_list_row_storage_size;
    }

    //! Get the sum of the number of neighbors of all targets' neighborhoods (host)
    global_index_type getTotalNeighborsOverAllListsHost() const {
        if (this->getNumberOfTargets()==0) {
            return 0;
        } else {
            return TO_GLOBAL(this->getNumberOfNeighborsHost(this->getNumberOfTargets()-1)) + this->getRowOffsetHost(this->getNumberOfTargets()-1);
        }
    }

    //! Get the sum of the number of neighbors of all targets' neighborhoods (device)
    KOKKOS_INLINE_FUNCTION
    global_index_type getTotalNeighborsOverAllListsDevice() const {
        return TO_GLOBAL(this->getNumberOfNeighborsDevice(this->getNumberOfTargets()-1)) + this->getRowOffsetDevice(this->getNumberOfTargets()-1);
    }
///@}

}; // NeighborLists

//! CreateNeighborLists allows for the construction of an object of type NeighborLists with template deduction
template <typename view_type>
NeighborLists<view_type> CreateNeighborLists(view_type number_of_neighbors_list) {
    return NeighborLists<view_type>(number_of_neighbors_list);
}

//! CreateNeighborLists allows for the construction of an object of type NeighborLists with template deduction
template <typename view_type>
NeighborLists<view_type> CreateNeighborLists(view_type neighbor_lists, view_type number_of_neighbors_list) {
    return NeighborLists<view_type>(neighbor_lists, number_of_neighbors_list);
}

//! CreateNeighborLists allows for the construction of an object of type NeighborLists with template deduction
template <typename view_type>
NeighborLists<view_type> CreateNeighborLists(view_type neighbor_lists, view_type number_of_neighbors_list, view_type neighbor_lists_row_offsets) {
    return NeighborLists<view_type>(neighbor_lists, number_of_neighbors_list, neighbor_lists_row_offsets);
}

//! Converts 2D neighbor lists to compressed row neighbor lists
template <typename view_type_2d, typename view_type_1d = Kokkos::View<int*, typename view_type_2d::memory_space, typename view_type_2d::memory_traits> >
NeighborLists<view_type_1d> Convert2DToCompressedRowNeighborLists(view_type_2d neighbor_lists) {

    // gets total number of neighbors over all lists
    // computes calculation where the data resides (device/host)
    global_index_type total_storage_size = 0;
    Kokkos::parallel_reduce("total number of neighbors over all lists", Kokkos::RangePolicy<typename view_type_2d::execution_space>(0, neighbor_lists.extent(0)), 
            KOKKOS_LAMBDA(const int i, global_index_type& t_total_num_neighbors) {
        t_total_num_neighbors += neighbor_lists(i,0);
    }, Kokkos::Sum<global_index_type>(total_storage_size));
    Kokkos::fence();

    // view_type_1d may be on host or device, and view_type_2d may be either as well (could even be opposite)
    view_type_1d new_cr_neighbor_lists("compressed row neighbor lists", total_storage_size);
    view_type_1d new_number_of_neighbors_list("number of neighbors list", neighbor_lists.extent(0));

    // copy number of neighbors list over to view_type_1d
    // d_neighbor_lists will be accessible from view_type_1d's execution space
    auto d_neighbor_lists = create_mirror_view(typename view_type_1d::execution_space(), neighbor_lists);
    Kokkos::deep_copy(d_neighbor_lists, neighbor_lists);
    Kokkos::fence();
    Kokkos::parallel_for("copy number of neighbors to compressed row", 
            Kokkos::RangePolicy<typename view_type_1d::execution_space>(0, neighbor_lists.extent(0)), 
            KOKKOS_LAMBDA(const int i) {
        new_number_of_neighbors_list(i) = d_neighbor_lists(i,0);
    });
    Kokkos::fence();

    
    // this will calculate row offsets
    auto nla = CreateNeighborLists(new_cr_neighbor_lists, new_number_of_neighbors_list);
    auto cr_data = nla.getNeighborLists();

    // if device_execution_space can access this view, then write directly into the view
    if (Kokkos::SpaceAccessibility<device_execution_space, typename view_type_1d::memory_space>::accessible==1) {
        Kokkos::parallel_for("copy neighbor lists to compressed row", Kokkos::RangePolicy<typename view_type_1d::execution_space>(0, neighbor_lists.extent(0)), 
                KOKKOS_LAMBDA(const int i) {
            for (int j=0; j<d_neighbor_lists(i,0); ++j) {
                cr_data(nla.getRowOffsetDevice(i)+j) = d_neighbor_lists(i,j+1);
            }
        });
        Kokkos::fence();
        nla.copyDeviceDataToHost(); // has a fence at the end
    }
    // otherwise we are writing to a view that can't be seen from device (must be host space), 
    // and d_neighbor_lists was already made to be a view_type that is accessible from view_type_1d's execution_space 
    // (which we know is host) so we can do a parallel_for over the host_execution_space
    else {
        Kokkos::parallel_for("copy neighbor lists to compressed row", Kokkos::RangePolicy<host_execution_space>(0, neighbor_lists.extent(0)), 
                KOKKOS_LAMBDA(const int i) {
            for (int j=0; j<neighbor_lists(i,0); ++j) {
                cr_data(nla.getRowOffsetHost(i)+j) = d_neighbor_lists(i,j+1);
            }
        });
        Kokkos::fence();
    }

    return nla;
}

} // Compadre namespace

#endif

