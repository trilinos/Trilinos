// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_SOLUTIONSET_HPP_
#define _COMPADRE_SOLUTIONSET_HPP_

#include "Compadre_Typedefs.hpp"
#include "Compadre_NeighborLists.hpp"
#include <Kokkos_Core.hpp>

namespace Compadre {

//!  All vairables and functionality related to the layout and storage of GMLS
//!  solutions (alpha values)
template <typename memory_space = device_memory_space>
struct SolutionSet {

    //! vector of user requested target operations
    Kokkos::View<TargetOperation*, memory_space> _lro; 

    //! vector containing a mapping from a target functionals enum value to the its place in the list
    //! of target functionals to be applied
    Kokkos::View<int*, memory_space> _lro_lookup; 

    //! index for where this operation begins the for _alpha coefficients
    Kokkos::View<int*, memory_space> _lro_total_offsets; 

    //! dimensions ^ rank of tensor of output for each target functional
    Kokkos::View<int*, memory_space> _lro_output_tile_size; 

    //! dimensions ^ rank of tensor of output for each sampling functional
    Kokkos::View<int*, memory_space> _lro_input_tile_size; 

    //! tensor rank of target functional (device)
    Kokkos::View<int*, memory_space> _lro_output_tensor_rank;

    //! tensor rank of sampling functional (device)
    Kokkos::View<int*, memory_space> _lro_input_tensor_rank;

    //! generated alpha coefficients (device)
    Kokkos::View<double*, layout_right, memory_space> _alphas; 

    //! additional alpha coefficients due to constraints
    int _added_alpha_size;

    //! maximum number of evaluation sites for each target (includes target site)
    int _max_evaluation_sites_per_target;

    //! used for sizing P_target_row and the _alphas view
    int _total_alpha_values;

    //! whether internal alpha values are valid (set externally on a solve)
    bool _contains_valid_alphas;

    //
    // Redundant variables (already exist in GMLS class)
    //
  
    //! Accessor to get neighbor list data, offset data, and number of neighbors per target
    NeighborLists<Kokkos::View<int*> > _neighbor_lists;

    //! generally the same as _polynomial_sampling_functional, but can differ if specified at 
    //! GMLS class instantiation
    SamplingFunctional _data_sampling_functional;

    //! dimension of the problem, set at class instantiation only
    int _dimensions;

    //! dimension of the problem, set at class instantiation only. For manifolds, generally _global_dimensions-1
    int _local_dimensions;

    //! problem type for GMLS problem, can also be set to STANDARD for normal or MANIFOLD for manifold problems
    ProblemType _problem_type;

/** @name Constructors
 */
///@{

    //! \brief Constructor for SolutionSet
    SolutionSet(SamplingFunctional data_sampling_functional,
                int dimensions, 
                int local_dimensions,
                const ProblemType problem_type) :
                    _added_alpha_size(0), 
                    _max_evaluation_sites_per_target(1),
                    _total_alpha_values(0),
                    _contains_valid_alphas(false),
                    _data_sampling_functional(data_sampling_functional),
                    _dimensions(dimensions), 
                    _local_dimensions(local_dimensions), 
                    _problem_type(problem_type) {}

    SolutionSet() : _data_sampling_functional(PointSample) {}

    //! \brief Copy constructor (can be used to move data from device to host or vice-versa)
    template <typename other_memory_space>
    SolutionSet(const SolutionSet<other_memory_space> &other) : 
            _data_sampling_functional(other._data_sampling_functional),
            _dimensions(other._dimensions),
            _local_dimensions(other._local_dimensions),
            _problem_type(other._problem_type) {

        _added_alpha_size = other._added_alpha_size;
        _max_evaluation_sites_per_target = other._max_evaluation_sites_per_target;
        _total_alpha_values = other._total_alpha_values;
        _contains_valid_alphas = false; // false until copyAlphas() is called
        _neighbor_lists = other._neighbor_lists;

        // copy from other_memory_space to memory_space (if needed)
        if (_lro.extent(0) != other._lro.extent(0)) {
            Kokkos::resize(_lro, other._lro.extent(0));
        }
        if (_lro_lookup.extent(0) != other._lro_lookup.extent(0)) {
            Kokkos::resize(_lro_lookup, other._lro_lookup.extent(0));
        }
        if (_lro_total_offsets.extent(0) != other._lro_total_offsets.extent(0)) {
            Kokkos::resize(_lro_total_offsets, other._lro_total_offsets.extent(0));
        }
        if (_lro_output_tile_size.extent(0) != other._lro_output_tile_size.extent(0)) {
            Kokkos::resize(_lro_output_tile_size, other._lro_output_tile_size.extent(0));
        }
        if (_lro_input_tile_size.extent(0) != other._lro_input_tile_size.extent(0)) {
            Kokkos::resize(_lro_input_tile_size, other._lro_input_tile_size.extent(0));
        }
        if (_lro_output_tensor_rank.extent(0) != other._lro_output_tensor_rank.extent(0)) {
            Kokkos::resize(_lro_output_tensor_rank, other._lro_output_tensor_rank.extent(0));
        }
        if (_lro_input_tensor_rank.extent(0) != other._lro_input_tensor_rank.extent(0)) {
            Kokkos::resize(_lro_input_tensor_rank, other._lro_input_tensor_rank.extent(0));
        }
        Kokkos::deep_copy(_lro, other._lro);
        Kokkos::deep_copy(_lro_lookup, other._lro_lookup);
        Kokkos::deep_copy(_lro_total_offsets, other._lro_total_offsets);
        Kokkos::deep_copy(_lro_output_tile_size, other._lro_output_tile_size);
        Kokkos::deep_copy(_lro_input_tile_size, other._lro_input_tile_size);
        Kokkos::deep_copy(_lro_output_tensor_rank, other._lro_output_tensor_rank);
        Kokkos::deep_copy(_lro_input_tensor_rank, other._lro_input_tensor_rank);

        // don't copy _alphas (expensive)
        // _alphas only copied using copyAlphas
    }

///@}

/** @name Public Accessors
 */
///@{

    // ON DEVICE

    //! Handles offset from operation input/output + extra evaluation sites
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    int getTargetOffsetIndex(const int lro_num, const int input_component, const int output_component, const int evaluation_site_local_index = 0) const {
        return ( _total_alpha_values*evaluation_site_local_index
                + _lro_total_offsets[lro_num] 
                + input_component*_lro_output_tile_size[lro_num] 
                + output_component );
    }

    //! Helper function for getting alphas for scalar reconstruction from scalar data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha0TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int evaluation_site_local_index = 0) const {
        // e.g. Dirac Delta target of a scalar field
        return getAlpha(lro, target_index, 0, 0, neighbor_index, 0, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for vector reconstruction from scalar data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha0TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int evaluation_site_local_index = 0) const {
        // e.g. gradient of a scalar field
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, 0, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for matrix reconstruction from scalar data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha0TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, 0, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for scalar reconstruction from vector data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha1TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int input_component, const int evaluation_site_local_index = 0) const {
        // e.g. divergence of a vector field
        return getAlpha(lro, target_index, 0, 0, neighbor_index, input_component, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for vector reconstruction from vector data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha1TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int input_component, const int evaluation_site_local_index = 0) const {
        // e.g. curl of a vector field
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, input_component, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for matrix reconstruction from vector data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha1TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component, const int evaluation_site_local_index = 0) const {
        // e.g. gradient of a vector field
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, input_component, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for scalar reconstruction from matrix data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha2TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, 0, 0, neighbor_index, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for vector reconstruction from matrix data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha2TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for matrix reconstruction from matrix data
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha2TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);
    }

    //! Gives index into alphas given two axes, which when incremented by the neighbor number transforms access into
    //! alphas from a rank 1 view into a rank 3 view.
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    global_index_type getAlphaIndex(const int target_index, const int alpha_column_offset) const {

        global_index_type total_neighbors_before_target = _neighbor_lists.getRowOffsetDevice(target_index);
        int total_added_alphas_before_target = target_index*_added_alpha_size;

        int alphas_per_tile_per_target = _neighbor_lists.getNumberOfNeighborsDevice(target_index) + _added_alpha_size;

        return (total_neighbors_before_target+TO_GLOBAL(total_added_alphas_before_target))
                 *TO_GLOBAL(_total_alpha_values)*TO_GLOBAL(_max_evaluation_sites_per_target)
                   + TO_GLOBAL(alpha_column_offset*alphas_per_tile_per_target);

    }

    //! Retrieves the offset for an operator based on input and output component, generic to row
    //! (but still multiplied by the number of neighbors for each row and then needs a neighbor number added 
    //! to this returned value to be meaningful)
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    int getAlphaColumnOffset(TargetOperation lro, const int output_component_axis_1, 
            const int output_component_axis_2, const int input_component_axis_1, 
            const int input_component_axis_2, const int evaluation_site_local_index = 0) const {

        const int lro_number = _lro_lookup[(int)lro];
        compadre_kernel_assert_debug((lro_number >= 0) && "getAlphaColumnOffset called for a TargetOperation that was not registered.");

        // the target functional input indexing is sized based on the output rank of the sampling
        // functional used, which can not be inferred unless a specification of target functional,
        // reconstruction space, and sampling functional are all known (as was the case at the
        // construction of this class)
        const int input_index = getSamplingOutputIndex(_data_sampling_functional, input_component_axis_1, input_component_axis_2);
        const int output_index = getTargetOutputIndex((int)lro, output_component_axis_1, output_component_axis_2, _dimensions);

        return getTargetOffsetIndex(lro_number, input_index, output_index, evaluation_site_local_index);
    }

    //! Underlying function all interface helper functions call to retrieve alpha values
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    double getAlpha(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        // lro - the operator from TargetOperations
        // target_index - the # for the target site where information is required
        // neighbor_index - the # for the neighbor of the target
        //
        // This code support up to rank 2 tensors for inputs and outputs
        //
        // scalar reconstruction from scalar data: rank 0 to rank 0
        //   provides 1 piece of information for each neighbor
        // scalar reconstruction from vector data (e.g. divergence): rank 1 to rank 0
        //   provides 'd' pieces of information for each neighbor
        // vector reconstruction from scalar data (e.g. gradient): rank 0 to rank 1
        //   provides 'd' piece of information for each neighbor
        // vector reconstruction from vector data (e.g. curl): rank 1 to rank 1
        //   provides 'd'x'd' pieces of information for each neighbor
        //
        // This function would more reasonably be called from one of the getAlphaNTensorFromNTensor
        // which is much easier to understand with respect to indexing and only requesting indices
        // that are relavent to the operator in question.
        //

        compadre_kernel_assert_debug(this->_contains_valid_alphas && 
                "getAlpha called on SolutionSet with _contains_valid_alphas=false");

        const int alpha_column_offset = this->getAlphaColumnOffset( lro, output_component_axis_1, 
                output_component_axis_2, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);

        auto alphas_index = this->getAlphaIndex(target_index, alpha_column_offset);
        return _alphas(alphas_index + neighbor_index);
    }

    //! Get the local index (internal) to GMLS for a particular TargetOperation
    //! Every TargetOperation has a global index which can be readily found in Compadre::TargetOperation
    //! but this function returns the index used inside of the GMLS class
    template<typename ms=memory_space, enable_if_t<!std::is_same<host_memory_space, ms>::value, int> = 0>
    KOKKOS_INLINE_FUNCTION
    int getTargetOperationLocalIndex(TargetOperation lro) const {
        return _lro_lookup[(int)lro];
    }

    //! Handles offset from operation input/output + extra evaluation sites
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    int getTargetOffsetIndex(const int lro_num, const int input_component, const int output_component, const int evaluation_site_local_index = 0) const {
        return ( _total_alpha_values*evaluation_site_local_index
                + _lro_total_offsets[lro_num] 
                + input_component*_lro_output_tile_size[lro_num] 
                + output_component );
    }

    //! Helper function for getting alphas for scalar reconstruction from scalar data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha0TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int evaluation_site_local_index = 0) const {
        // e.g. Dirac Delta target of a scalar field
        return getAlpha(lro, target_index, 0, 0, neighbor_index, 0, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for vector reconstruction from scalar data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha0TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int evaluation_site_local_index = 0) const {
        // e.g. gradient of a scalar field
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, 0, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for matrix reconstruction from scalar data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha0TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, 0, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for scalar reconstruction from vector data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha1TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int input_component, const int evaluation_site_local_index = 0) const {
        // e.g. divergence of a vector field
        return getAlpha(lro, target_index, 0, 0, neighbor_index, input_component, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for vector reconstruction from vector data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha1TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int input_component, const int evaluation_site_local_index = 0) const {
        // e.g. curl of a vector field
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, input_component, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for matrix reconstruction from vector data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha1TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component, const int evaluation_site_local_index = 0) const {
        // e.g. gradient of a vector field
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, input_component, 0, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for scalar reconstruction from matrix data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha2TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, 0, 0, neighbor_index, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for vector reconstruction from matrix data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha2TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);
    }

    //! Helper function for getting alphas for matrix reconstruction from matrix data
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha2TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);
    }

    //! Gives index into alphas given two axes, which when incremented by the neighbor number transforms access into
    //! alphas from a rank 1 view into a rank 3 view.
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    global_index_type getAlphaIndex(const int target_index, const int alpha_column_offset) const {

        global_index_type total_neighbors_before_target = _neighbor_lists.getRowOffsetHost(target_index);
        int total_added_alphas_before_target = target_index*_added_alpha_size;

        int alphas_per_tile_per_target = _neighbor_lists.getNumberOfNeighborsHost(target_index) + _added_alpha_size;

        return (total_neighbors_before_target+TO_GLOBAL(total_added_alphas_before_target))
                 *TO_GLOBAL(_total_alpha_values)*TO_GLOBAL(_max_evaluation_sites_per_target)
                   + TO_GLOBAL(alpha_column_offset*alphas_per_tile_per_target);

    }

    //! Retrieves the offset for an operator based on input and output component, generic to row
    //! (but still multiplied by the number of neighbors for each row and then needs a neighbor number added 
    //! to this returned value to be meaningful)
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    int getAlphaColumnOffset(TargetOperation lro, const int output_component_axis_1, 
            const int output_component_axis_2, const int input_component_axis_1, 
            const int input_component_axis_2, const int evaluation_site_local_index = 0) const {

        const int lro_number = _lro_lookup[(int)lro];
        compadre_kernel_assert_debug((lro_number >= 0) && "getAlphaColumnOffset called for a TargetOperation that was not registered.");

        // the target functional input indexing is sized based on the output rank of the sampling
        // functional used, which can not be inferred unless a specification of target functional,
        // reconstruction space, and sampling functional are all known (as was the case at the
        // construction of this class)
        const int input_index = getSamplingOutputIndex(_data_sampling_functional, input_component_axis_1, input_component_axis_2);
        const int output_index = getTargetOutputIndex((int)lro, output_component_axis_1, output_component_axis_2, _dimensions);

        return getTargetOffsetIndex(lro_number, input_index, output_index, evaluation_site_local_index);
    }

    //! Underlying function all interface helper functions call to retrieve alpha values
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    double getAlpha(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int evaluation_site_local_index = 0) const {
        // lro - the operator from TargetOperations
        // target_index - the # for the target site where information is required
        // neighbor_index - the # for the neighbor of the target
        //
        // This code support up to rank 2 tensors for inputs and outputs
        //
        // scalar reconstruction from scalar data: rank 0 to rank 0
        //   provides 1 piece of information for each neighbor
        // scalar reconstruction from vector data (e.g. divergence): rank 1 to rank 0
        //   provides 'd' pieces of information for each neighbor
        // vector reconstruction from scalar data (e.g. gradient): rank 0 to rank 1
        //   provides 'd' piece of information for each neighbor
        // vector reconstruction from vector data (e.g. curl): rank 1 to rank 1
        //   provides 'd'x'd' pieces of information for each neighbor
        //
        // This function would more reasonably be called from one of the getAlphaNTensorFromNTensor
        // which is much easier to understand with respect to indexing and only requesting indices
        // that are relavent to the operator in question.
        //
 
        compadre_assert_debug(this->_contains_valid_alphas && 
                "getAlpha called on SolutionSet with _contains_valid_alphas=false");

        const int alpha_column_offset = this->getAlphaColumnOffset( lro, output_component_axis_1, 
                output_component_axis_2, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);

        auto alphas_index = this->getAlphaIndex(target_index, alpha_column_offset);
        return _alphas(alphas_index + neighbor_index);
    }

    //! Get the local index (internal) to GMLS for a particular TargetOperation
    //! Every TargetOperation has a global index which can be readily found in Compadre::TargetOperation
    //! but this function returns the index used inside of the GMLS class
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    int getTargetOperationLocalIndex(TargetOperation lro) const {
        return _lro_lookup[(int)lro];
    }


///@}

/** @name Public Modifiers (can only call from host)
 */
///@{

    //! Copies alphas between two instances of SolutionSet
    //! Copying of alphas is intentionally omitted in copy constructor
    template <typename other_memory_space>
    void copyAlphas(SolutionSet<other_memory_space>& other) {
        if ((void*)this != (void*)&other) {
            if (_alphas.extent(0) != other._alphas.extent(0)) {
                Kokkos::resize(_alphas, other._alphas.extent(0));
            }
            Kokkos::deep_copy(_alphas, other._alphas);
            this->_contains_valid_alphas = other._contains_valid_alphas;
        }
    }


    //! Empties the vector of target functionals to apply to the reconstruction
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    void clearTargets() {
        _lro = decltype(_lro)();
        for (int i=0; i<TargetOperation::COUNT; ++i) {
            _lro_lookup[i] = -1;
        }
    }

    //! Adds a target to the vector of target functional to be applied to the reconstruction
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    void addTargets(TargetOperation lro) {
        std::vector<TargetOperation> temporary_lro_vector(1, lro);
        this->addTargets(temporary_lro_vector);
    }

    //! Adds a vector of target functionals to the vector of target functionals already to be applied to the reconstruction
    template<typename ms=memory_space, enable_if_t<std::is_same<host_memory_space, ms>::value, int> = 0>
    void addTargets(std::vector<TargetOperation> lro) {

        std::vector<TargetOperation> unique_new_lro;

        auto host_lro_lookup = create_mirror_view(_lro_lookup);
        if (_lro_lookup.extent(0) == 0) {
            _lro_lookup = decltype(_lro_lookup)("LRO Lookup", TargetOperation::COUNT);
            host_lro_lookup = create_mirror_view(_lro_lookup);
            Kokkos::deep_copy(host_lro_lookup, -1);
        } else {
            Kokkos::deep_copy(host_lro_lookup, _lro_lookup);
        }

        // loop over requested targets
        for (size_t i=0; i<lro.size(); ++i) {

            bool operation_found = false;
            // loop over existing targets registered
            for (size_t j=0; j<_lro.size(); ++j) {

                // if found
                if (_lro(j)==lro[i]) {

                    operation_found = true;

                    // the operation should now point to where the operation is stored
                    host_lro_lookup[(int)lro[i]] = j;

                    break;

                }
            }

            if (!operation_found) {
                host_lro_lookup[(int)lro[i]] = _lro.size() + unique_new_lro.size();
                unique_new_lro.push_back(lro[i]);
            }
        }

        // move unique_new_lro into _lro
        auto new_lro = decltype(_lro)("LRO", _lro.size() + unique_new_lro.size());
        for (size_t i=0; i<_lro.size(); ++i) {
            new_lro(i) = _lro(i);
        }
        for (size_t i=0; i<unique_new_lro.size(); ++i) {
            new_lro(i+_lro.size()) = unique_new_lro[i];
        }
        _lro = new_lro;

        _lro_total_offsets = decltype(_lro_total_offsets)("total offsets for alphas", _lro.size());
        _lro_output_tile_size = decltype(_lro_output_tile_size)("output tile size for each operation", _lro.size());
        _lro_input_tile_size = decltype(_lro_input_tile_size)("output tile size for each operation", _lro.size());
        _lro_output_tensor_rank = decltype(_lro_output_tensor_rank)("output tensor rank", _lro.size());
        _lro_input_tensor_rank = decltype(_lro_input_tensor_rank)("input tensor rank", _lro.size());

        auto host_lro_total_offsets = create_mirror_view(_lro_total_offsets);
        auto host_lro_output_tile_size = create_mirror_view(_lro_output_tile_size);
        auto host_lro_input_tile_size = create_mirror_view(_lro_input_tile_size);
        auto host_lro_output_tensor_rank = create_mirror_view(_lro_output_tensor_rank);
        auto host_lro_input_tensor_rank = create_mirror_view(_lro_input_tensor_rank);

        int total_offset = 0; // need total offset

        for (size_t i=0; i<_lro.size(); ++i) {
            host_lro_total_offsets(i) = total_offset;

            // allows for a tile of the product of dimension^input_tensor_rank * dimension^output_tensor_rank * the number of neighbors
            int output_tile_size = getOutputDimensionOfOperation(_lro(i), _local_dimensions);

            // the target functional input indexing is sized based on the output rank of the sampling
            // functional used
            int input_tile_size = getOutputDimensionOfSampling(_data_sampling_functional, _local_dimensions);
            host_lro_output_tile_size(i) = output_tile_size;
            host_lro_input_tile_size(i) = input_tile_size;

            total_offset += input_tile_size * output_tile_size;

            // the target functional output rank is based on the output rank of the sampling
            // functional used
            host_lro_input_tensor_rank(i) = _data_sampling_functional.output_rank;
            host_lro_output_tensor_rank(i) = getTargetOutputTensorRank((int)_lro(i));
        }

        _total_alpha_values = total_offset;

        if (_problem_type == ProblemType::MANIFOLD) {
            // if on a manifold, the total alphas values must be large enough to hold the gradient
            // of the geometry reconstruction
            _total_alpha_values = (_total_alpha_values > pown(_local_dimensions, 1)) ? 
                _total_alpha_values : pown(_local_dimensions, 1);
        }

        Kokkos::deep_copy(_lro_lookup, host_lro_lookup);
        Kokkos::deep_copy(_lro_total_offsets, host_lro_total_offsets);
        Kokkos::deep_copy(_lro_output_tile_size, host_lro_output_tile_size);
        Kokkos::deep_copy(_lro_input_tile_size, host_lro_input_tile_size);
        Kokkos::deep_copy(_lro_output_tensor_rank, host_lro_output_tensor_rank);
        Kokkos::deep_copy(_lro_input_tensor_rank, host_lro_input_tensor_rank);
    }

    //! Get a view (device) of all alphas
    decltype(_alphas) getAlphas() const { 
        compadre_assert_debug(this->_contains_valid_alphas && 
                "getAlpha called on SolutionSet with _contains_valid_alphas=false");
        return _alphas; 
    }

///@}

}; // SolutionSet

} // Compadre namespace

#endif

