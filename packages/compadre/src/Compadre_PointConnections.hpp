// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_POINTCONNECTIONS_HPP_
#define _COMPADRE_POINTCONNECTIONS_HPP_

#include "Compadre_Typedefs.hpp"
#include <Kokkos_Core.hpp>

namespace Compadre {

//!  Combines NeighborLists with the PointClouds from which it was derived
//!  Assumed that memory_space is the same as device, but it can be set to
//!  host, if desired.
template <typename view_type_1, typename view_type_2, typename nla_type, typename memory_space = device_memory_space>
struct PointConnections {

    //! source site coordinates on device
    typedef decltype(Kokkos::create_mirror_view<memory_space>(
                memory_space(), view_type_1()))
                        device_mirror_target_view_type;

    //! target site coordinates on device
    typedef decltype(Kokkos::create_mirror_view<memory_space>(
                memory_space(), view_type_2()))
                        device_mirror_source_view_type;

    device_mirror_target_view_type _target_coordinates;
    device_mirror_source_view_type _source_coordinates;
    nla_type _nla;

/** @name Constructors
 */
///@{

    //! \brief Constructor for PointConnections
    PointConnections(view_type_1 target_coordinates, 
                     view_type_2 source_coordinates,
                     nla_type nla) : _nla(nla) {

        _target_coordinates = Kokkos::create_mirror_view<memory_space>(
                memory_space(), target_coordinates);
        _source_coordinates = Kokkos::create_mirror_view<memory_space>(
                memory_space(), source_coordinates);
        Kokkos::deep_copy(_target_coordinates, target_coordinates);
        Kokkos::deep_copy(_source_coordinates, source_coordinates);

    }

    PointConnections() {}

    // copy constructor (can be used to move data from device to host or vice-versa)
    template <typename other_type_1, typename other_type_2, typename other_type_3>
    PointConnections(const PointConnections<other_type_1, other_type_2, other_type_3> &other) : 
        PointConnections(other._target_coordinates, other._source_coordinates, other._nla) {}

///@}

/** @name Public Utility
 *  
 */
///@{

    //! Returns a component of the local coordinate after transformation from global to local under the orthonormal basis V.
    KOKKOS_INLINE_FUNCTION
    static double convertGlobalToLocalCoordinate(const XYZ global_coord, const int dim, const scratch_matrix_right_type& V) {
        compadre_kernel_assert_debug(dim<3);
        // only written for 2d manifold in 3d space or 2D problem with 1D manifold
        double val = 0;
        val += global_coord.x * V(dim, 0);
        if (V.extent_int(1)>1) val += global_coord.y * V(dim, 1); 
        if (V.extent_int(1)>2) val += global_coord.z * V(dim, 2);
        return val;
    }

    //! Returns a component of the global coordinate after transformation from local to global under the orthonormal basis V^T.
    KOKKOS_INLINE_FUNCTION
    static double convertLocalToGlobalCoordinate(const XYZ local_coord, const int dim, const scratch_matrix_right_type& V) {
        double val = 0.0;
        if (dim == 0 && V.extent_int(0)==1) { // 2D problem with 1D manifold
            val = local_coord.x * V(0, dim);
        } else { // 3D problem with 2D manifold
            val = local_coord.x * V(0, dim) + local_coord.y * V(1, dim);
        }
        return val;
    }

    //! Returns Euclidean norm of a vector
    KOKKOS_INLINE_FUNCTION
    static double EuclideanVectorLength(const XYZ& delta_vector, const int dimension) {
        double inside_val = delta_vector.x*delta_vector.x;
        switch (dimension) {
        case 3:
            inside_val += delta_vector.z*delta_vector.z;
            // no break is intentional
        case 2:
            inside_val += delta_vector.y*delta_vector.y;
            // no break is intentional
        default:
            break;
        }
        return std::sqrt(inside_val);
    }


///@}

/** @name Public Modifiers
 *  Private function because information lives on the device
 */
///@{

    //! Update only target coordinates
    void setTargetCoordinates(view_type_1 target_coordinates) {
        _target_coordinates = Kokkos::create_mirror_view<memory_space>(
                memory_space(), target_coordinates);
        Kokkos::deep_copy(_target_coordinates, target_coordinates);
    }

    //! Update only source coordinates
    void setSourceCoordinates(view_type_2 source_coordinates) {
        _source_coordinates = Kokkos::create_mirror_view<memory_space>(
                memory_space(), source_coordinates);
        Kokkos::deep_copy(_source_coordinates, source_coordinates);
    }

    //! Update only target coordinates
    void setNeighborLists(nla_type nla) {
        _nla = nla;
    }


///@}

/** @name Public Accessors
 */
///@{

    //! Returns one component of the target coordinate for a particular target. Whether global or local coordinates 
    //! depends upon V being specified
    KOKKOS_INLINE_FUNCTION
    double getTargetCoordinate(const int target_index, const int dim, const scratch_matrix_right_type* V = NULL) const {
        compadre_kernel_assert_debug((_target_coordinates.extent(0) >= (size_t)target_index) && "Target index is out of range for _target_coordinates.");
        if (V==NULL) {
            return _target_coordinates(target_index, dim);
        } else {
            XYZ target_coord = XYZ(_target_coordinates(target_index, 0), 0, 0);
            if (_target_coordinates.extent_int(1)>1) target_coord[1] = _target_coordinates(target_index, 1);
            if (_target_coordinates.extent_int(1)>2) target_coord[2] = _target_coordinates(target_index, 2);
            return this->convertGlobalToLocalCoordinate(target_coord, dim, *V);
        }
    }

    //! Returns one component of the neighbor coordinate for a particular target. Whether global or local coordinates 
    //! depends upon V being specified
    KOKKOS_INLINE_FUNCTION
    double getNeighborCoordinate(const int target_index, const int neighbor_list_num, const int dim, const scratch_matrix_right_type* V = NULL) const {
        compadre_kernel_assert_debug((_source_coordinates.extent(0) >= (size_t)(this->getNeighborIndex(target_index, neighbor_list_num))) && "Source index is out of range for _source_coordinates.");
        if (V==NULL) {
            return _source_coordinates(this->getNeighborIndex(target_index, neighbor_list_num), dim);
        } else {
            XYZ neighbor_coord 
                = XYZ(_source_coordinates(this->getNeighborIndex(target_index, neighbor_list_num), 0), 0, 0);
            if (_source_coordinates.extent_int(1)>1) neighbor_coord[1] 
                = _source_coordinates(this->getNeighborIndex(target_index, neighbor_list_num), 1);
            if (_source_coordinates.extent_int(1)>2) neighbor_coord[2] 
                = _source_coordinates(this->getNeighborIndex(target_index, neighbor_list_num), 2);
            return this->convertGlobalToLocalCoordinate(neighbor_coord, dim, *V);
        }
    }

    //! Returns the relative coordinate as a vector between the target site and the neighbor site. 
    //! Whether global or local coordinates depends upon V being specified
    KOKKOS_INLINE_FUNCTION
    XYZ getRelativeCoord(const int target_index, const int neighbor_list_num, const int dimension, const scratch_matrix_right_type* V = NULL) const {
        XYZ coordinate_delta;

        coordinate_delta.x = this->getNeighborCoordinate(target_index, neighbor_list_num, 0, V) - this->getTargetCoordinate(target_index, 0, V);
        if (dimension>1) coordinate_delta.y = this->getNeighborCoordinate(target_index, neighbor_list_num, 1, V) - this->getTargetCoordinate(target_index, 1, V);
        if (dimension>2) coordinate_delta.z = this->getNeighborCoordinate(target_index, neighbor_list_num, 2, V) - this->getTargetCoordinate(target_index, 2, V);

        return coordinate_delta;
    }

    //! Mapping from [0,number of neighbors for a target] to the row that contains the source coordinates for
    //! that neighbor
    KOKKOS_INLINE_FUNCTION
    int getNeighborIndex(const int target_index, const int neighbor_list_num) const {
        return _nla.getNeighborDevice(target_index, neighbor_list_num);
    }

///@}

}; // PointConnections

} // Compadre namespace

#endif

