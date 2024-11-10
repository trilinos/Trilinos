// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_EVALUATOR_HPP_
#define _COMPADRE_EVALUATOR_HPP_

#include "Compadre_Typedefs.hpp"
#include "Compadre_GMLS.hpp"
#include "Compadre_NeighborLists.hpp"

namespace Compadre {

//! Creates 1D subviews of data from a 2D view, generally constructed with CreateNDSliceOnDeviceView
template<typename T, typename T2, typename T3=void>
struct SubviewND { 
    
    T _data_in;
    T2 _data_original_view;
    bool _scalar_as_vector_if_needed;

    SubviewND(T data_in, T2 data_original_view, bool scalar_as_vector_if_needed) {
        _data_in = data_in;
        _data_original_view = data_original_view;
        _scalar_as_vector_if_needed = scalar_as_vector_if_needed; 
    }

    auto get1DView(const int column_num) -> decltype(Kokkos::subview(_data_in, Kokkos::ALL, column_num)) {
        if (!_scalar_as_vector_if_needed) {
            compadre_assert_debug(((size_t)column_num<_data_in.extent(1)) 
                    && "Subview asked for column > second dimension of input data.");
        }
        if ((size_t)column_num<_data_in.extent(1))
            return Kokkos::subview(_data_in, Kokkos::ALL, column_num);
        else // scalar treated as a vector (being reused for each component of the vector input that was expected)
            return Kokkos::subview(_data_in, Kokkos::ALL, 0);
    }

    auto get2DView(const int column_num, const int block_size) -> decltype(Kokkos::subview(_data_in, Kokkos::ALL, 
                Kokkos::make_pair(column_num*block_size, (column_num+1)*block_size))) {
        if (!_scalar_as_vector_if_needed) {
            compadre_assert_debug(((size_t)((column_num+1)*block_size-1)<_data_in.extent(1)) 
                    && "Subview asked for column > second dimension of input data.");
        }
        if ((size_t)((column_num+1)*block_size-1)<_data_in.extent(1)) {
            return Kokkos::subview(_data_in, Kokkos::ALL, Kokkos::make_pair(column_num*block_size, (column_num+1)*block_size));
        } else {
            compadre_assert_debug(((size_t)(block_size-1)<_data_in.extent(1)) && "Subview asked for column > second dimension of input data.");
            return Kokkos::subview(_data_in, Kokkos::ALL, Kokkos::make_pair(0,block_size));
        }
    }

    T2 copyToAndReturnOriginalView() {
        Kokkos::deep_copy(_data_original_view, _data_in);
        Kokkos::fence();
        return _data_original_view;
    }

};

//! Creates 1D subviews of data from a 1D view, generally constructed with CreateNDSliceOnDeviceView
template<typename T, typename T2>
struct SubviewND<T, T2, enable_if_t<(T::rank<2)> >
{ 

    T _data_in;
    T2 _data_original_view;
    bool _scalar_as_vector_if_needed;

    SubviewND(T data_in, T2 data_original_view, bool scalar_as_vector_if_needed) {
        _data_in = data_in;
        _data_original_view = data_original_view;
        _scalar_as_vector_if_needed = scalar_as_vector_if_needed; 
    }

    auto get1DView(const int column_num) -> decltype(Kokkos::subview(_data_in, Kokkos::ALL)) {
        // TODO: There is a valid use case for violating this assert, so in the future we may want
        // to add other logic to the evaluator function calling this so that it knows to do nothing with
        // this data.
        if (!_scalar_as_vector_if_needed) {
            compadre_assert_debug((column_num==0) && "Subview asked for column column_num!=0, but _data_in is rank 1.");
        }
        return Kokkos::subview(_data_in, Kokkos::ALL);
    }

    auto get2DView(const int column_num, const int block_size) -> decltype(Kokkos::subview(_data_in, Kokkos::ALL)) {
        compadre_assert_release((block_size==1) && "2D subview requested not compatible with one column.");
        return Kokkos::subview(_data_in, Kokkos::ALL);
    }

    T2 copyToAndReturnOriginalView() {
        Kokkos::deep_copy(_data_original_view, _data_in);
        Kokkos::fence();
        return _data_original_view;
    }

};

//! Copies data_in to the device, and then allows for access to 1D columns of data on device.
//! Handles either 2D or 1D views as input, and they can be on the host or the device.
template <typename T>
auto CreateNDSliceOnDeviceView(T sampling_input_data_host_or_device, bool scalar_as_vector_if_needed) -> SubviewND<decltype(Kokkos::create_mirror_view(
                    device_memory_space(), sampling_input_data_host_or_device)), T> {

    // makes view on the device (does nothing if already on the device)
    auto sampling_input_data_device = Kokkos::create_mirror_view(
        device_memory_space(), sampling_input_data_host_or_device);
    Kokkos::deep_copy(sampling_input_data_device, sampling_input_data_host_or_device);
    Kokkos::fence();

    return SubviewND<decltype(sampling_input_data_device),T>(sampling_input_data_device, 
            sampling_input_data_host_or_device, scalar_as_vector_if_needed);
}

//! \brief Lightweight Evaluator Helper
//! This class is a lightweight wrapper for extracting and applying all relevant data from a GMLS class
//! in order to transform data into a form that can be acted on by the GMLS operator, apply the action of
//! the GMLS operator, and then transform data again (only if on a manifold)
class Evaluator {

private:

    GMLS *_gmls;


public:

    Evaluator(GMLS *gmls) : _gmls(gmls) {
        Kokkos::fence();
    };

    ~Evaluator() {};

    //! Dot product of alphas with sampling data, FOR A SINGLE target_index,  where sampling data is in a 1D/2D Kokkos View
    //! 
    //! This function is to be used when the alpha values have already been calculated and stored for use 
    //!
    //! Only supports one output component / input component at a time. The user will need to loop over the output 
    //! components in order to fill a vector target or matrix target.
    //! 
    //! Assumptions on input data:
    //! \param sampling_input_data      [in] - 1D/2D Kokkos View (no restriction on memory space)
    //! \param column_of_input          [in] - Column of sampling_input_data to use for this input component
    //! \param lro                      [in] - Target operation from the TargetOperation enum
    //! \param target_index             [in] - Target # user wants to reconstruct target functional at, corresponds to row number of neighbor_lists
    //! \param output_component_axis_1  [in] - Row for a rank 2 tensor or rank 1 tensor, 0 for a scalar output
    //! \param output_component_axis_2  [in] - Columns for a rank 2 tensor, 0 for rank less than 2 output tensor
    //! \param input_component_axis_1   [in] - Row for a rank 2 tensor or rank 1 tensor, 0 for a scalar input
    //! \param input_component_axis_2   [in] - Columns for a rank 2 tensor, 0 for rank less than 2 input tensor
    //! \param scalar_as_vector_if_needed [in] - If a 1D view is given, where a 2D view is expected (scalar values given where a vector was expected), then the scalar will be repeated for as many components as the vector has
    template <typename view_type_data>
    double applyAlphasToDataSingleComponentSingleTargetSite(view_type_data sampling_input_data, const int column_of_input, TargetOperation lro, const int target_index, const int evaluation_site_local_index, const int output_component_axis_1, const int output_component_axis_2, const int input_component_axis_1, const int input_component_axis_2, bool scalar_as_vector_if_needed = true) const {

        double value = 0;

        const int alpha_input_output_component_index = _gmls->_h_ss.getAlphaColumnOffset(lro, output_component_axis_1, 
                output_component_axis_2, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);

        auto sampling_subview_maker = CreateNDSliceOnDeviceView(sampling_input_data, scalar_as_vector_if_needed);

        
        // gather needed information for evaluation
        auto nla = *(_gmls->getNeighborLists());
        auto alphas = _gmls->getSolutionSetDevice()->getAlphas();
        auto sampling_data_device = sampling_subview_maker.get1DView(column_of_input);
        
        auto alpha_index = _gmls->_h_ss.getAlphaIndex(target_index, alpha_input_output_component_index);
        // loop through neighbor list for this target_index
        // grabbing data from that entry of data
        Kokkos::parallel_reduce("applyAlphasToData::Device", 
                Kokkos::RangePolicy<device_execution_space>(0,nla.getNumberOfNeighborsHost(target_index)), 
                KOKKOS_LAMBDA(const int i, double& t_value) {

            t_value += sampling_data_device(nla.getNeighborDevice(target_index, i))
                *alphas(alpha_index + i);

        }, value );
        Kokkos::fence();

        return value;
    }

    //! Dot product of alphas with sampling data where sampling data is in a 1D/2D Kokkos View and output view is also 
    //! a 1D/2D Kokkos View, however THE SAMPLING DATA and OUTPUT VIEW MUST BE ON THE DEVICE!
    //! 
    //! This function is to be used when the alpha values have already been calculated and stored for use.
    //!
    //! Only supports one output component / input component at a time. The user will need to loop over the output 
    //! components in order to fill a vector target or matrix target.
    //! 
    //! Assumptions on input data:
    //! \param output_data_single_column       [out] - 1D Kokkos View (memory space must be device_memory_space())
    //! \param sampling_data_single_column      [in] - 1D Kokkos View (memory space must match output_data_single_column)
    //! \param lro                              [in] - Target operation from the TargetOperation enum
    //! \param sro                              [in] - Sampling functional from the SamplingFunctional enum
    //! \param evaluation_site_local_index      [in] - local column index of site from additional evaluation sites list or 0 for the target site
    //! \param output_component_axis_1          [in] - Row for a rank 2 tensor or rank 1 tensor, 0 for a scalar output
    //! \param output_component_axis_2          [in] - Columns for a rank 2 tensor, 0 for rank less than 2 output tensor
    //! \param input_component_axis_1           [in] - Row for a rank 2 tensor or rank 1 tensor, 0 for a scalar input
    //! \param input_component_axis_2           [in] - Columns for a rank 2 tensor, 0 for rank less than 2 input tensor
    //! \param pre_transform_local_index        [in] - For manifold problems, this is the local coordinate direction that sampling data may need to be transformed to before the application of GMLS
    //! \param pre_transform_global_index       [in] - For manifold problems, this is the global coordinate direction that sampling data can be represented in
    //! \param post_transform_local_index       [in] - For manifold problems, this is the local coordinate direction that vector output target functionals from GMLS will output into
    //! \param post_transform_global_index      [in] - For manifold problems, this is the global coordinate direction that the target functional output from GMLS will be transformed into
    //! \param transform_output_ambient         [in] - Whether or not a 1D output from GMLS is on the manifold and needs to be mapped to ambient space
    //! \param vary_on_target                   [in] - Whether the sampling functional has a tensor to act on sampling data that varies with each target site
    //! \param vary_on_neighbor                 [in] - Whether the sampling functional has a tensor to act on sampling data that varies with each neighbor site in addition to varying wit each target site
    template <typename view_type_data_out, typename view_type_data_in>
    void applyAlphasToDataSingleComponentAllTargetSitesWithPreAndPostTransform(view_type_data_out output_data_single_column, view_type_data_in sampling_data_single_column, TargetOperation lro, const SamplingFunctional sro, const int evaluation_site_local_index, const int output_component_axis_1, const int output_component_axis_2, const int input_component_axis_1, const int input_component_axis_2, const int pre_transform_local_index = -1, const int pre_transform_global_index = -1, const int post_transform_local_index = -1, const int post_transform_global_index = -1, bool vary_on_target = false, bool vary_on_neighbor = false) const {

        const int alpha_input_output_component_index = _gmls->_h_ss.getAlphaColumnOffset(lro, output_component_axis_1, 
                output_component_axis_2, input_component_axis_1, input_component_axis_2, evaluation_site_local_index);
        const int alpha_input_output_component_index2 = alpha_input_output_component_index;

        // gather needed information for evaluation
        auto nla = *(_gmls->getNeighborLists());
        auto solution_set = *(_gmls->getSolutionSetDevice());
        compadre_assert_release(solution_set._contains_valid_alphas && 
                "application of alphas called before generateAlphas() was called.");
        auto prestencil_weights = _gmls->getPrestencilWeights();

        const int num_targets = nla.getNumberOfTargets();

        // make sure input and output views have same memory space
        compadre_assert_debug((std::is_same<typename view_type_data_out::memory_space, typename view_type_data_in::memory_space>::value) && 
                "output_data_single_column view and input_data_single_column view have difference memory spaces.");

        bool weight_with_pre_T = (pre_transform_local_index>=0 && pre_transform_global_index>=0) ? true : false;
        bool target_plus_neighbor_staggered_schema = sro.use_target_site_weights;

        // loops over target indices
        Kokkos::parallel_for(team_policy(num_targets, Kokkos::AUTO),
                KOKKOS_LAMBDA(const member_type& teamMember) {

            const int target_index = teamMember.league_rank();
            teamMember.team_barrier();


            const double previous_value = output_data_single_column(target_index);

            // loops over neighbors of target_index
            auto alpha_index = solution_set.getAlphaIndex(target_index, alpha_input_output_component_index);
            double gmls_value = 0;
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember, nla.getNumberOfNeighborsDevice(target_index)), [&](const int i, double& t_value) {
                const double neighbor_varying_pre_T =  (weight_with_pre_T && vary_on_neighbor) ?
                    prestencil_weights(0, target_index, i, pre_transform_local_index, pre_transform_global_index)
                    : 1.0;

                t_value += neighbor_varying_pre_T * sampling_data_single_column(nla.getNeighborDevice(target_index, i))
                            *solution_set._alphas(alpha_index + i);

            }, gmls_value );

            // data contract for sampling functional
            double pre_T = 1.0;
            if (weight_with_pre_T) {
                if (!vary_on_neighbor && vary_on_target) {
                    pre_T = prestencil_weights(0, target_index, 0, pre_transform_local_index, 
                            pre_transform_global_index); 
                } else if (!vary_on_target) { // doesn't vary on target or neighbor
                    pre_T = prestencil_weights(0, 0, 0, pre_transform_local_index, 
                            pre_transform_global_index); 
                }
            }

            double staggered_value_from_targets = 0;
            double pre_T_staggered = 1.0;
            auto alpha_index2 = solution_set.getAlphaIndex(target_index, alpha_input_output_component_index2);
            // loops over target_index for each neighbor for staggered approaches
            if (target_plus_neighbor_staggered_schema) {
                Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember, nla.getNumberOfNeighborsDevice(target_index)), [&](const int i, double& t_value) {
                    const double neighbor_varying_pre_T_staggered =  (weight_with_pre_T && vary_on_neighbor) ?
                        prestencil_weights(1, target_index, i, pre_transform_local_index, pre_transform_global_index)
                        : 1.0;

                    t_value += neighbor_varying_pre_T_staggered * sampling_data_single_column(nla.getNeighborDevice(target_index, 0))
                                *solution_set._alphas(alpha_index2 + i);

                }, staggered_value_from_targets );

                // for staggered approaches that transform source data for the target and neighbors
                if (weight_with_pre_T) {
                    if (!vary_on_neighbor && vary_on_target) {
                        pre_T_staggered = prestencil_weights(1, target_index, 0, pre_transform_local_index, 
                                pre_transform_global_index); 
                    } else if (!vary_on_target) { // doesn't vary on target or neighbor
                        pre_T_staggered = prestencil_weights(1, 0, 0, pre_transform_local_index, 
                                pre_transform_global_index); 
                    }
                }
            }

            double added_value = pre_T*gmls_value + pre_T_staggered*staggered_value_from_targets;
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                output_data_single_column(target_index) = previous_value + added_value;
            });
        });
        Kokkos::fence();
    }

    //! Postprocessing for manifolds. Maps local chart vector solutions to ambient space.
    //! THE SAMPLING DATA and OUTPUT VIEW MUST BE ON THE DEVICE!
    //! 
    //! Only supports one output component / input component at a time. The user will need to loop over the output 
    //! components in order to transform a vector target.
    //! 
    //! Assumptions on input data:
    //! \param output_data_single_column       [out] - 1D Kokkos View (memory space must be device_memory_space())
    //! \param sampling_data_single_column      [in] - 1D Kokkos View (memory space must match output_data_single_column)
    //! \param local_dim_index                  [in] - For manifold problems, this is the local coordinate direction that sampling data may need to be transformed to before the application of GMLS
    //! \param global_dim_index                 [in] - For manifold problems, this is the global coordinate direction that sampling data can be represented in
    template <typename view_type_data_out, typename view_type_data_in>
    void applyLocalChartToAmbientSpaceTransform(view_type_data_out output_data_single_column, view_type_data_in sampling_data_single_column, const int local_dim_index, const int global_dim_index) const {

        // Does T transpose times a vector
        auto global_dimensions = _gmls->getGlobalDimensions();

        // gather needed information for evaluation
        auto nla = *(_gmls->getNeighborLists());
        const int num_targets = nla.getNumberOfTargets();

        auto tangent_directions = *(_gmls->getTangentDirections());

        // make sure input and output views have same memory space
        compadre_assert_debug((std::is_same<typename view_type_data_out::memory_space, typename view_type_data_in::memory_space>::value) && 
                "output_data_single_column view and input_data_single_column view have difference memory spaces.");

        // loops over target indices
        Kokkos::parallel_for(team_policy(num_targets, Kokkos::AUTO),
                KOKKOS_LAMBDA(const member_type& teamMember) {

            const int target_index = teamMember.league_rank();

            scratch_matrix_right_type T
                    (tangent_directions.data() + TO_GLOBAL(target_index)*TO_GLOBAL(global_dimensions)*TO_GLOBAL(global_dimensions), 
                     global_dimensions, global_dimensions);
            teamMember.team_barrier();


            const double previous_value = output_data_single_column(target_index);

            double added_value = T(local_dim_index, global_dim_index)*sampling_data_single_column(target_index);
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                output_data_single_column(target_index) = previous_value + added_value;
            });
        });
        Kokkos::fence();
    }

    //! Transformation of data under GMLS (allocates memory for output)
    //! 
    //! This function is the go-to function to be used when the alpha values have already been calculated and stored for use. The sampling functional provided instructs how a data transformation tensor is to be used on source data before it is provided to the GMLS operator. Once the sampling functional (if applicable) and the GMLS operator have been applied, this function also handles mapping the local vector back to the ambient space if working on a manifold problem and a target functional who has rank 1 output.
    //!
    //! Produces a Kokkos View as output with a Kokkos memory_space provided as a template tag by the caller. 
    //! The data type (double* or double**) must also be specified as a template type if one wish to get a 1D 
    //! Kokkos View back that can be indexed into with only one ordinal.
    //! 
    //! Assumptions on input data:
    //! \param sampling_data              [in] - 1D or 2D Kokkos View that has the layout #targets * columns of data. Memory space for data can be host or device. 
    //! \param lro                        [in] - Target operation from the TargetOperation enum
    //! \param sro_in                     [in] - Sampling functional from the SamplingFunctional enum
    //! \param scalar_as_vector_if_needed [in] - If a 1D view is given, where a 2D view is expected (scalar values given where a vector was expected), then the scalar will be repeated for as many components as the vector has
    //! \param evaluation_site_local_index [in] - 0 corresponds to evaluating at the target site itself, while a number larger than 0 indicates evaluation at a site other than the target, and specified by calling setAdditionalEvaluationSitesData on the GMLS class
    template <typename output_data_type = double**, typename output_memory_space, typename view_type_input_data, typename output_array_layout = typename view_type_input_data::array_layout>
    Kokkos::View<output_data_type, output_array_layout, output_memory_space>  // shares layout of input by default
            applyAlphasToDataAllComponentsAllTargetSites(view_type_input_data sampling_data, TargetOperation lro, const SamplingFunctional sro_in = PointSample, bool scalar_as_vector_if_needed = true, const int evaluation_site_local_index = 0) const {
        // gather needed information for evaluation
        auto nla = *(_gmls->getNeighborLists());

        // determines the number of columns needed for output after action of the target functional
        auto local_dimensions = _gmls->getLocalDimensions();
        auto global_dimensions = _gmls->getGlobalDimensions();
        int output_dimensions = getOutputDimensionOfOperation(lro, local_dimensions);
        auto problem_type = _gmls->getProblemType();

        typedef Kokkos::View<output_data_type, output_array_layout, output_memory_space> output_view_type;
        // create view on whatever memory space the user specified with their template argument when calling this function
        output_view_type target_output = createView<output_view_type>("output of target operation", 
                nla.getNumberOfTargets(), output_dimensions);

        output_view_type ambient_target_output;
        bool transform_gmls_output_to_ambient = (problem_type==MANIFOLD && getTargetOutputTensorRank(lro)==1);
        if (transform_gmls_output_to_ambient) {
            ambient_target_output = createView<output_view_type>("output of transform to ambient space", 
                nla.getNumberOfTargets(), global_dimensions);
        }

        applyAlphasToDataAllComponentsAllTargetSites(target_output, ambient_target_output, sampling_data, lro, sro_in, scalar_as_vector_if_needed, evaluation_site_local_index);

        if (transform_gmls_output_to_ambient) {
            return ambient_target_output;
        } else {
            return target_output;
        }

    }

    //! Transformation of data under GMLS (does not allocate memory for output)
    //! 
    //! If space for the output result is already allocated, this function will populate the output result view (and possibly ambient target output). The sampling functional provided instructs how a data transformation tensor is to be used on source data before it is provided to the GMLS operator. Once the sampling functional (if applicable) and the GMLS operator have been applied, this function also handles mapping the local vector back to the ambient space if working on a manifold problem and a target functional who has rank 1 output.
    //!
    //! Fills a Kokkos View of output. 
    //! 
    //! Assumptions on input data:
    //! \param target_output              [in] - 1D or 2D Kokkos View that has the resulting #targets * need output columns. Memory space for data can be host or device. 
    //! \param ambient_target_output      [in] - Same view type as target_output, but dimensions should be #targets * global_dimension if this is being filled (if not being filled, then this can be an empty view)
    //! \param sampling_data              [in] - 1D or 2D Kokkos View that has the layout #targets * columns of data. Memory space for data can be host or device. 
    //! \param lro                        [in] - Target operation from the TargetOperation enum
    //! \param sro_in                     [in] - Sampling functional from the SamplingFunctional enum
    //! \param scalar_as_vector_if_needed [in] - If a 1D view is given, where a 2D view is expected (scalar values given where a vector was expected), then the scalar will be repeated for as many components as the vector has
    //! \param evaluation_site_local_index [in] - 0 corresponds to evaluating at the target site itself, while a number larger than 0 indicates evaluation at a site other than the target, and specified by calling setAdditionalEvaluationSitesData on the GMLS class
    template <typename view_type_output_data, typename view_type_input_data, typename output_array_layout = typename view_type_input_data::array_layout>
    void applyAlphasToDataAllComponentsAllTargetSites(view_type_output_data target_output, view_type_output_data ambient_target_output, view_type_input_data sampling_data, TargetOperation lro, const SamplingFunctional sro_in = PointSample, bool scalar_as_vector_if_needed = true, const int evaluation_site_local_index = 0) const {


        // output can be device or host
        // input can be device or host
        // move everything to device and calculate there, then move back to host if necessary


        auto problem_type = _gmls->getProblemType();
        auto global_dimensions = _gmls->getGlobalDimensions();
        auto local_dimensions = _gmls->getLocalDimensions();
        int output_dimension1_of_operator = (getTargetOutputTensorRank(lro)<2) ? getOutputDimensionOfOperation(lro, local_dimensions) : std::sqrt(getOutputDimensionOfOperation(lro, local_dimensions));
        int output_dimension2_of_operator = (getTargetOutputTensorRank(lro)<2) ? 1 : std::sqrt(getOutputDimensionOfOperation(lro, local_dimensions));

        // gather needed information for evaluation
        auto nla = *(_gmls->getNeighborLists());

        // determines the number of columns needed for output after action of the target functional
        int output_dimensions = getOutputDimensionOfOperation(lro, local_dimensions);

        // special case for VectorPointSample, because if it is on a manifold it includes data transform to local charts
        auto sro = (problem_type==MANIFOLD && sro_in==VectorPointSample) ? ManifoldVectorPointSample : sro_in;
        int input_dimension_of_operator = getInputDimensionOfOperation(lro, _gmls->_data_sampling_functional, local_dimensions);

        compadre_assert_debug(target_output.extent(0)==(size_t)nla.getNumberOfTargets() 
                && "First dimension of target_output is incorrect size.\n");
        compadre_assert_debug(target_output.extent(1)==(size_t)output_dimensions 
                && "Second dimension of target_output is incorrect size.\n");

        // make sure input and output columns make sense under the target operation
        compadre_assert_debug(((output_dimensions==1 && view_type_output_data::rank==1) || view_type_output_data::rank!=1) && 
                "Output view is requested as rank 1, but the target requires a rank larger than 1. Try double** as template argument.");

        // we need to specialize a template on the rank of the output view type and the input view type
        auto sampling_subview_maker = CreateNDSliceOnDeviceView(sampling_data, scalar_as_vector_if_needed);
        auto output_subview_maker = CreateNDSliceOnDeviceView(target_output, false); // output will always be the correct dimension

        // figure out preprocessing and postprocessing
        auto prestencil_weights = _gmls->getPrestencilWeights();

        // all loop logic based on transforming data under a sampling functional
        // into something that is valid input for GMLS
        bool vary_on_target = false, vary_on_neighbor = false;
        auto sro_style = sro.transform_type;
        bool loop_global_dimensions = sro.input_rank>0 && sro_style!=Identity; 

        compadre_assert_release((_gmls->getDataSamplingFunctional()==sro || sro_style==Identity)
                && "SamplingFunctional requested for Evaluator does not match GMLS data sampling functional or is not of type 'Identity'.");

        if (sro.transform_type == Identity || sro.transform_type == SameForAll) {
            vary_on_target = false;
            vary_on_neighbor = false;
        } else if (sro.transform_type == DifferentEachTarget) {
            vary_on_target = true;
            vary_on_neighbor = false;
        } else if (sro.transform_type == DifferentEachNeighbor) {
            vary_on_target = true;
            vary_on_neighbor = true;
        }


        // only written for up to rank 1 to rank 2 (in / out)
        // loop over components of output of the target operation
        for (int axes1=0; axes1<output_dimension1_of_operator; ++axes1) {
            const int output_component_axis_1 = axes1;
            for (int axes2=0; axes2<output_dimension2_of_operator; ++axes2) {
                const int output_component_axis_2 = axes2;
                // loop over components of input of the target operation
                for (int j=0; j<input_dimension_of_operator; ++j) {
                    const int input_component_axis_1 = j;
                    const int input_component_axis_2 = 0;

                    if (loop_global_dimensions) {
                        for (int k=0; k<global_dimensions; ++k) { // loop for handling sampling functional
                            this->applyAlphasToDataSingleComponentAllTargetSitesWithPreAndPostTransform(
                                    output_subview_maker.get1DView(axes1*output_dimension2_of_operator+axes2), 
                                    sampling_subview_maker.get1DView(k), lro, sro, 
                                    evaluation_site_local_index, output_component_axis_1, output_component_axis_2, input_component_axis_1, 
                                    input_component_axis_2, j, k, -1, -1,
                                    vary_on_target, vary_on_neighbor);
                        }
                    } else if (sro_style != Identity) {
                        this->applyAlphasToDataSingleComponentAllTargetSitesWithPreAndPostTransform(
                                output_subview_maker.get1DView(axes1*output_dimension2_of_operator+axes2), 
                                sampling_subview_maker.get1DView(j), lro, sro, 
                                evaluation_site_local_index, output_component_axis_1, output_component_axis_2, input_component_axis_1, 
                                input_component_axis_2, 0, 0, -1, -1,
                                vary_on_target, vary_on_neighbor);
                    } else { // standard
                        this->applyAlphasToDataSingleComponentAllTargetSitesWithPreAndPostTransform(
                                output_subview_maker.get1DView(axes1*output_dimension2_of_operator+axes2), 
                                sampling_subview_maker.get1DView(j), lro, sro, 
                                evaluation_site_local_index, output_component_axis_1, output_component_axis_2, input_component_axis_1, 
                                input_component_axis_2);
                    }
                }
            }
        }

        bool transform_gmls_output_to_ambient = (problem_type==MANIFOLD && getTargetOutputTensorRank(lro)==1);
        if (transform_gmls_output_to_ambient) {
            Kokkos::fence();

            compadre_assert_debug(ambient_target_output.extent(0)==(size_t)nla.getNumberOfTargets() 
                    && "First dimension of target_output is incorrect size.\n");
            compadre_assert_debug(ambient_target_output.extent(1)==(size_t)global_dimensions 
                    && "Second dimension of target_output is incorrect size.\n");
            auto transformed_output_subview_maker = CreateNDSliceOnDeviceView(ambient_target_output, false); 
            // output will always be the correct dimension
            for (int i=0; i<global_dimensions; ++i) {
                for (int j=0; j<output_dimensions; ++j) {
                    this->applyLocalChartToAmbientSpaceTransform(
                            transformed_output_subview_maker.get1DView(i), output_subview_maker.get1DView(j), j, i);
                }
            }
            // copy back to whatever memory space the user requester through templating from the device
            Kokkos::deep_copy(ambient_target_output, transformed_output_subview_maker.copyToAndReturnOriginalView());
        }

        // copy back to whatever memory space the user requester through templating from the device
        Kokkos::deep_copy(target_output, output_subview_maker.copyToAndReturnOriginalView());
    }


    //! Dot product of data with full polynomial coefficient basis where sampling data is in a 1D/2D Kokkos View and output view is also 
    //! a 1D/2D Kokkos View, however THE SAMPLING DATA and OUTPUT VIEW MUST BE ON THE DEVICE!
    //! 
    //! This function is to be used when the polynomial coefficient basis has already been calculated and stored for use.
    //!
    //! Only supports one output component / input component at a time. The user will need to loop over the output 
    //! components in order to fill a vector target or matrix target.
    //! 
    //! Assumptions on input data:
    //! \param output_data_block_column       [out] - 2D Kokkos View (memory space must be device_memory_space())
    //! \param sampling_data_single_column      [in] - 1D Kokkos View (memory space must match output_data_single_column)
    //! \param sro                              [in] - Sampling functional from the SamplingFunctional enum
    //! \param target_index                     [in] - Target # user wants to reconstruct target functional at, corresponds to row number of neighbor_lists
    //! \param output_component_axis_1          [in] - Row for a rank 2 tensor or rank 1 tensor, 0 for a scalar output
    //! \param output_component_axis_2          [in] - Columns for a rank 2 tensor, 0 for rank less than 2 output tensor
    //! \param input_component_axis_1           [in] - Row for a rank 2 tensor or rank 1 tensor, 0 for a scalar input
    //! \param input_component_axis_2           [in] - Columns for a rank 2 tensor, 0 for rank less than 2 input tensor
    //! \param pre_transform_local_index        [in] - For manifold problems, this is the local coordinate direction that sampling data may need to be transformed to before the application of GMLS
    //! \param pre_transform_global_index       [in] - For manifold problems, this is the global coordinate direction that sampling data can be represented in
    //! \param post_transform_local_index       [in] - For manifold problems, this is the local coordinate direction that vector output target functionals from GMLS will output into
    //! \param post_transform_global_index      [in] - For manifold problems, this is the global coordinate direction that the target functional output from GMLS will be transformed into
    //! \param vary_on_target                   [in] - Whether the sampling functional has a tensor to act on sampling data that varies with each target site
    //! \param vary_on_neighbor                 [in] - Whether the sampling functional has a tensor to act on sampling data that varies with each neighbor site in addition to varying wit each target site
    template <typename view_type_data_out, typename view_type_data_in>
    void applyFullPolynomialCoefficientsBasisToDataSingleComponent(view_type_data_out output_data_block_column, view_type_data_in sampling_data_single_column, const SamplingFunctional sro, const int output_component_axis_1, const int output_component_axis_2, const int input_component_axis_1, const int input_component_axis_2, const int pre_transform_local_index = -1, const int pre_transform_global_index = -1, const int post_transform_local_index = -1, const int post_transform_global_index = -1, bool vary_on_target = false, bool vary_on_neighbor = false) const {

        auto nla = *(_gmls->getNeighborLists());

        auto coefficient_matrix_dims = _gmls->getPolynomialCoefficientsDomainRangeSize();
        auto coefficient_memory_layout_dims = _gmls->getPolynomialCoefficientsMemorySize();
        auto coefficient_memory_layout_dims_device = 
            Kokkos::create_mirror_view_and_copy(device_memory_space(), coefficient_memory_layout_dims);

        auto global_dimensions = _gmls->getGlobalDimensions();

        // gather needed information for evaluation
        auto coeffs         = _gmls->getFullPolynomialCoefficientsBasis();
        auto tangent_directions = *(_gmls->getTangentDirections());
        auto prestencil_weights = _gmls->getPrestencilWeights();

        const int num_targets = nla.getNumberOfTargets();

        // make sure input and output views have same memory space
        compadre_assert_debug((std::is_same<typename view_type_data_out::memory_space, typename view_type_data_in::memory_space>::value) && 
                "output_data_block_column view and input_data_single_column view have difference memory spaces.");

        bool weight_with_pre_T = (pre_transform_local_index>=0 && pre_transform_global_index>=0) ? true : false;
        bool target_plus_neighbor_staggered_schema = sro.use_target_site_weights;

        // loops over target indices
        for (int j=0; j<coefficient_matrix_dims(0); ++j) {
            Kokkos::parallel_for(team_policy(num_targets, Kokkos::AUTO),
                    KOKKOS_LAMBDA(const member_type& teamMember) {

                const int target_index = teamMember.league_rank();

                scratch_matrix_right_type T (tangent_directions.data() 
                            + TO_GLOBAL(target_index)*TO_GLOBAL(global_dimensions)*TO_GLOBAL(global_dimensions), 
                         global_dimensions, global_dimensions);
                
                scratch_matrix_right_type Coeffs;
                Coeffs = scratch_matrix_right_type(coeffs.data() 
                    + TO_GLOBAL(target_index)*TO_GLOBAL(coefficient_memory_layout_dims_device(0))
                        *TO_GLOBAL(coefficient_memory_layout_dims_device(1)),
                    coefficient_memory_layout_dims_device(0), coefficient_memory_layout_dims_device(1));

                teamMember.team_barrier();


                const double previous_value = output_data_block_column(target_index, j);

                // loops over neighbors of target_index
                double gmls_value = 0;
                Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember, nla.getNumberOfNeighborsDevice(target_index)), [&](const int i, double& t_value) {
                    const double neighbor_varying_pre_T =  (weight_with_pre_T && vary_on_neighbor) ?
                        prestencil_weights(0, target_index, i, pre_transform_local_index, pre_transform_global_index)
                        : 1.0;

                    t_value += neighbor_varying_pre_T * sampling_data_single_column(nla.getNeighborDevice(target_index, i))
                        *Coeffs(j, i+input_component_axis_1*nla.getNumberOfNeighborsDevice(target_index));

                }, gmls_value );

                // data contract for sampling functional
                double pre_T = 1.0;
                if (weight_with_pre_T) {
                    if (!vary_on_neighbor && vary_on_target) {
                        pre_T = prestencil_weights(0, target_index, 0, pre_transform_local_index, 
                                pre_transform_global_index); 
                    } else if (!vary_on_target) { // doesn't vary on target or neighbor
                        pre_T = prestencil_weights(0, 0, 0, pre_transform_local_index, 
                                pre_transform_global_index); 
                    }
                }

                double staggered_value_from_targets = 0;
                double pre_T_staggered = 1.0;
                // loops over target_index for each neighbor for staggered approaches
                if (target_plus_neighbor_staggered_schema) {
                    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember, nla.getNumberOfNeighborsDevice(target_index)), [&](const int i, double& t_value) {
                        const double neighbor_varying_pre_T_staggered =  (weight_with_pre_T && vary_on_neighbor) ?
                            prestencil_weights(1, target_index, i, pre_transform_local_index, pre_transform_global_index)
                            : 1.0;

                        t_value += neighbor_varying_pre_T_staggered * sampling_data_single_column(nla.getNeighborDevice(target_index, 0))
                            *Coeffs(j, i+input_component_axis_1*nla.getNumberOfNeighborsDevice(target_index));

                    }, staggered_value_from_targets );

                    // for staggered approaches that transform source data for the target and neighbors
                    if (weight_with_pre_T) {
                        if (!vary_on_neighbor && vary_on_target) {
                            pre_T_staggered = prestencil_weights(1, target_index, 0, pre_transform_local_index, 
                                    pre_transform_global_index); 
                        } else if (!vary_on_target) { // doesn't vary on target or neighbor
                            pre_T_staggered = prestencil_weights(1, 0, 0, pre_transform_local_index, 
                                    pre_transform_global_index); 
                        }
                    }
                }

                double added_value = (pre_T*gmls_value + pre_T_staggered*staggered_value_from_targets);
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    output_data_block_column(target_index, j) = previous_value + added_value;
                });
            });
            Kokkos::fence();
        }
    }

    //! Generation of polynomial reconstruction coefficients by applying to data in GMLS (allocates memory for output)
    //! 
    //! Polynomial reconstruction coefficients exist for each target, but there are coefficients for each neighbor (a basis for all potentional input data). This function uses a particular choice of data to contract over this basis and return the polynomial reconstructions coefficients specific to this data.
    //!
    //! Produces a Kokkos View as output with a Kokkos memory_space provided as a template tag by the caller. 
    //! The data type (double* or double**) must also be specified as a template type if one wish to get a 1D 
    //! Kokkos View back that can be indexed into with only one ordinal.
    //! 
    //! Assumptions on input data:
    //! \param sampling_data              [in] - 1D or 2D Kokkos View that has the layout #targets * columns of data. Memory space for data can be host or device. 
    //! \param sro                        [in] - Sampling functional from the SamplingFunctional enum
    //! \param scalar_as_vector_if_needed [in] - If a 1D view is given, where a 2D view is expected (scalar values given where a vector was expected), then the scalar will be repeated for as many components as the vector has
    template <typename output_data_type = double**, typename output_memory_space, typename view_type_input_data, typename output_array_layout = typename view_type_input_data::array_layout>
    Kokkos::View<output_data_type, output_array_layout, output_memory_space>  // shares layout of input by default
            applyFullPolynomialCoefficientsBasisToDataAllComponents(view_type_input_data sampling_data, bool scalar_as_vector_if_needed = true) const {

        auto local_dimensions = _gmls->getLocalDimensions();
        auto output_dimension_of_reconstruction_space = calculateBasisMultiplier(_gmls->getReconstructionSpace(), local_dimensions);
        auto coefficient_matrix_dims = _gmls->getPolynomialCoefficientsDomainRangeSize();

        // gather needed information for evaluation
        auto nla = *(_gmls->getNeighborLists());

        // determines the number of columns needed for output
        int output_dimensions = output_dimension_of_reconstruction_space;

        typedef Kokkos::View<output_data_type, output_array_layout, output_memory_space> output_view_type;
        // create view on whatever memory space the user specified with their template argument when calling this function
        output_view_type coefficient_output("output coefficients", nla.getNumberOfTargets(), 
                output_dimensions*_gmls->getPolynomialCoefficientsSize() /* number of coefficients */);

        applyFullPolynomialCoefficientsBasisToDataAllComponents(coefficient_output, sampling_data, scalar_as_vector_if_needed);

        return coefficient_output;
        
    }

    //! Generation of polynomial reconstruction coefficients by applying to data in GMLS (does not allocate memory for output)
    //! 
    //! Polynomial reconstruction coefficients exist for each target, but there are coefficients for each neighbor (a basis for all potentional input data). This function uses a particular choice of data to contract over this basis and return the polynomial reconstructions coefficients specific to this data.
    //!
    //! Produces a Kokkos View as output with a Kokkos memory_space provided as a template tag by the caller. 
    //! The data type (double* or double**) must also be specified as a template type if one wish to get a 1D 
    //! Kokkos View back that can be indexed into with only one ordinal.
    //! 
    //! Assumptions on input data:
    //! \param coefficient_output         [in] - 1D or 2D Kokkos View that has the layout #targets * #coefficients. Memory space for data can be host or device. 
    //! \param sampling_data              [in] - 1D or 2D Kokkos View that has the layout #targets * columns of data. Memory space for data can be host or device. 
    //! \param sro                        [in] - Sampling functional from the SamplingFunctional enum
    //! \param scalar_as_vector_if_needed [in] - If a 1D view is given, where a 2D view is expected (scalar values given where a vector was expected), then the scalar will be repeated for as many components as the vector has
    template <typename view_type_coefficient_output, typename view_type_input_data>
    void applyFullPolynomialCoefficientsBasisToDataAllComponents(view_type_coefficient_output coefficient_output, 
            view_type_input_data sampling_data, bool scalar_as_vector_if_needed = true) const {

        // output can be device or host
        // input can be device or host
        // move everything to device and calculate there, then move back to host if necessary

        auto global_dimensions = _gmls->getGlobalDimensions();
        auto local_dimensions = _gmls->getLocalDimensions();
        auto output_dimension_of_reconstruction_space = calculateBasisMultiplier(_gmls->getReconstructionSpace(), local_dimensions);
        auto input_dimension_of_reconstruction_space = calculateSamplingMultiplier(_gmls->getReconstructionSpace(), _gmls->getPolynomialSamplingFunctional(), local_dimensions);
        auto coefficient_matrix_dims = _gmls->getPolynomialCoefficientsDomainRangeSize();

        // gather needed information for evaluation
        auto nla = *(_gmls->getNeighborLists());

        const SamplingFunctional sro = _gmls->getDataSamplingFunctional();

        compadre_assert_debug(coefficient_output.extent(0)==(size_t)nla.getNumberOfTargets() 
                && "First dimension of coefficient_output is incorrect size.\n");
        // determines the number of columns needed for output
        compadre_assert_debug(
                coefficient_output.extent(1)==(size_t)output_dimension_of_reconstruction_space*_gmls->getPolynomialCoefficientsSize() 
                /* number of coefficients */ && "Second dimension of coefficient_output is incorrect size.\n");

        // make sure input and output columns make sense under the target operation
        compadre_assert_debug(((output_dimension_of_reconstruction_space==1 && view_type_coefficient_output::rank==1) || view_type_coefficient_output::rank!=1) && 
                "Output view is requested as rank 1, but the target requires a rank larger than 1. Try double** as template argument.");

        // we need to specialize a template on the rank of the output view type and the input view type
        auto sampling_subview_maker = CreateNDSliceOnDeviceView(sampling_data, scalar_as_vector_if_needed);
        auto output_subview_maker = CreateNDSliceOnDeviceView(coefficient_output, false); // output will always be the correct dimension

        // figure out preprocessing and postprocessing
        auto prestencil_weights = _gmls->getPrestencilWeights();

        // all loop logic based on transforming data under a sampling functional
        // into something that is valid input for GMLS
        bool vary_on_target = false, vary_on_neighbor = false;
        auto sro_style = sro.transform_type;
        bool loop_global_dimensions = sro.input_rank>0 && sro_style!=Identity; 


        if (sro.transform_type == Identity || sro.transform_type == SameForAll) {
            vary_on_target = false;
            vary_on_neighbor = false;
        } else if (sro.transform_type == DifferentEachTarget) {
            vary_on_target = true;
            vary_on_neighbor = false;
        } else if (sro.transform_type == DifferentEachNeighbor) {
            vary_on_target = true;
            vary_on_neighbor = true;
        }

        // written for up to rank 1 to rank 0 (in / out)
        // loop over components of output of the target operation
        for (int i=0; i<output_dimension_of_reconstruction_space; ++i) {
            const int output_component_axis_1 = i;
            const int output_component_axis_2 = 0;
            // loop over components of input of the target operation
            for (int j=0; j<input_dimension_of_reconstruction_space; ++j) {
                const int input_component_axis_1 = j;
                const int input_component_axis_2 = 0;

                if (loop_global_dimensions) {
                    for (int k=0; k<global_dimensions; ++k) { // loop for handling sampling functional
                        this->applyFullPolynomialCoefficientsBasisToDataSingleComponent(
                                output_subview_maker.get2DView(i,_gmls->getPolynomialCoefficientsSize()), 
                                sampling_subview_maker.get1DView(k), sro, 
                                output_component_axis_1, output_component_axis_2, input_component_axis_1, 
                                input_component_axis_2, j, k, -1, -1,
                                vary_on_target, vary_on_neighbor);
                    }
                } else if (sro_style != Identity) {
                    this->applyFullPolynomialCoefficientsBasisToDataSingleComponent(
                            output_subview_maker.get2DView(i,_gmls->getPolynomialCoefficientsSize()), 
                            sampling_subview_maker.get1DView(j), sro, 
                            output_component_axis_1, output_component_axis_2, input_component_axis_1, 
                            input_component_axis_2, 0, 0, -1, -1,
                            vary_on_target, vary_on_neighbor);
                } else { // standard
                    this->applyFullPolynomialCoefficientsBasisToDataSingleComponent(
                            output_subview_maker.get2DView(i,_gmls->getPolynomialCoefficientsSize()), 
                            sampling_subview_maker.get1DView(j), sro, 
                            output_component_axis_1, output_component_axis_2, input_component_axis_1, 
                            input_component_axis_2);
                }
            }
        }

        // copy back to whatever memory space the user requester through templating from the device
        Kokkos::deep_copy(coefficient_output, output_subview_maker.copyToAndReturnOriginalView());
    }

}; // Evaluator

} // Compadre

#endif
