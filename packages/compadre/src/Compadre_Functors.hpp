// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_FUNCTORS_HPP_
#define _COMPADRE_FUNCTORS_HPP_

#include "Compadre_Operators.hpp"
#include "Compadre_ParallelManager.hpp"
#include "Compadre_Basis.hpp"
#include "Compadre_Quadrature.hpp"
#include "Compadre_Targets.hpp"
#include "Compadre_ApplyTargetEvaluations.hpp"
#include "Compadre_Functors.hpp"
#include "Compadre_CreateConstraints.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "Compadre_GMLS.hpp"

namespace Compadre {

struct GMLSBasisData {

    Kokkos::View<double**, layout_right> _source_extra_data;
    Kokkos::View<double**, layout_right> _target_extra_data;
    Kokkos::View<double*> _epsilons; 
    Kokkos::View<double*****, layout_right> _prestencil_weights; 
    Kokkos::View<TargetOperation*> _curvature_support_operations;
    Kokkos::View<TargetOperation*> _operations;

    int _poly_order; 
    int _curvature_poly_order;
    int _NP;
    int _global_dimensions;
    int _local_dimensions;
    int _dimensions;
    int _reconstruction_space_rank;
    int _weighting_p;
    int _weighting_n;
    int _curvature_weighting_p;
    int _curvature_weighting_n;
    int _basis_multiplier;
    int _sampling_multiplier;
    int _data_sampling_multiplier;
    int _initial_index_for_batch;
    int _order_of_quadrature_points;
    int _dimension_of_quadrature_points;

    GMLS::point_connections_type _pc;
    GMLS::point_connections_type _additional_pc;
    ParallelManager _pm;
    ReconstructionSpace _reconstruction_space;
    DenseSolverType _dense_solver_type;
    ProblemType _problem_type;
    ConstraintType _constraint_type;
    SamplingFunctional _polynomial_sampling_functional;
    SamplingFunctional _data_sampling_functional;
    WeightingFunctionType _weighting_type;
    WeightingFunctionType _curvature_weighting_type;
    Quadrature _qm;
    SolutionSet<device_memory_space> _d_ss;

    // convenience variables (not from GMLS class)
    int RHS_dim_0, RHS_dim_1;
    double * RHS_data;
    int P_dim_0, P_dim_1;
    double * P_data;
    int P_target_row_dim_0, P_target_row_dim_1;
    double * P_target_row_data;
    int Coeffs_dim_0, Coeffs_dim_1;
    double * Coeffs_data;
    double * w_data;
    int max_num_neighbors;
    int max_num_rows;
    int manifold_NP;
    int max_manifold_NP;
    int this_num_cols;
    int max_poly_order;
    double * T_data;
    int ref_N_dim;
    double * ref_N_data;
    int thread_workspace_dim;
    double * manifold_curvature_coefficients_data;
    int manifold_gradient_dim;

    GMLSBasisData() : _polynomial_sampling_functional(PointSample), _data_sampling_functional(PointSample) {}
};

struct GMLSSolutionData {

    int _sampling_multiplier;
    int _initial_index_for_batch;
    SolutionSet<device_memory_space> _d_ss;

    // convenience variables (not from GMLS class)
    int this_num_cols;
    int Coeffs_dim_0, Coeffs_dim_1;
    double * Coeffs_data;
    int P_target_row_dim_0, P_target_row_dim_1;
    double * P_target_row_data;
    size_t operations_size;
    Kokkos::View<int*> number_of_neighbors_list;
    Kokkos::View<int*> additional_number_of_neighbors_list;

};

/** @name Functors
 *  Member functions that perform operations on the entire batch
 */
///@{

//! Functor to apply target evaluation to polynomial coefficients to store in _alphas
struct ApplyTargets {

    GMLSSolutionData _data;

    ApplyTargets(GMLSSolutionData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        const int local_index  = teamMember.league_rank();

        /*
         *    Data
         */

        // Coefficients for polynomial basis have overwritten _data._RHS
        scratch_matrix_right_type Coeffs = scratch_matrix_right_type(_data.Coeffs_data 
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.Coeffs_dim_0*_data.Coeffs_dim_1), 
                    _data.Coeffs_dim_0, _data.Coeffs_dim_1);
        scratch_matrix_right_type P_target_row(_data.P_target_row_data 
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_target_row_dim_0*_data.P_target_row_dim_1), 
                    _data.P_target_row_dim_0, _data.P_target_row_dim_1);

        applyTargetsToCoefficients(_data, teamMember, Coeffs, P_target_row); 
    }
};

//! Functor to evaluate targets operations on the basis
struct EvaluateStandardTargets {

    GMLSBasisData _data;

    EvaluateStandardTargets(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {
        /*
         *    Dimensions
         */

        const int local_index  = teamMember.league_rank();

        /*
         *    Data
         */

        scratch_matrix_right_type P_target_row(_data.P_target_row_data 
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_target_row_dim_0*_data.P_target_row_dim_1), 
                    _data.P_target_row_dim_0, _data.P_target_row_dim_1);

        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        /*
         *    Evaluate Standard Targets
         */

        computeTargetFunctionals(_data, teamMember, delta, thread_workspace, P_target_row);

    }
};

//! Functor to calculate prestencil weights to apply to data to transform into a format expected by a GMLS stencil
struct ComputePrestencilWeights {

    GMLSBasisData _data;

    ComputePrestencilWeights(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index  = teamMember.league_rank();
        const int dimensions   = _data._dimensions;

        /*
         *    Data
         */

        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.manifold_NP);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        scratch_matrix_right_type tangent;
        scratch_vector_type t1, t2;
        if (_data._data_sampling_functional == VaryingManifoldVectorPointSample) {
            tangent = scratch_matrix_right_type(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(1)), 
                    dimensions-1, dimensions);
            t1 = scratch_vector_type(teamMember.team_scratch(_data._pm.getTeamScratchLevel(1)), 
                    _data.max_num_neighbors);
            t2 = scratch_vector_type(teamMember.team_scratch(_data._pm.getTeamScratchLevel(1)), 
                    _data.max_num_neighbors);
            for (int j = 0; j < delta.extent_int(0); ++j) {
                delta(j) = 0;
            }
            for (int j = 0; j < thread_workspace.extent_int(0); ++j) {
                thread_workspace(j) = 0;
            }
        }


        // holds polynomial coefficients for curvature reconstruction
        scratch_matrix_right_type Q = scratch_matrix_right_type(_data.Coeffs_data 
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.Coeffs_dim_0*_data.Coeffs_dim_1),
                    _data.Coeffs_dim_0, _data.Coeffs_dim_1);

        scratch_matrix_right_type T(_data.T_data
                + TO_GLOBAL(target_index)*TO_GLOBAL(dimensions)*TO_GLOBAL(dimensions), 
                    dimensions, dimensions);

        scratch_vector_type manifold_gradient(teamMember.team_scratch(_data._pm.getTeamScratchLevel(1)), 
                _data.manifold_gradient_dim);

        /*
         *    Prestencil Weight Calculations
         */

        if (_data._data_sampling_functional == StaggeredEdgeAnalyticGradientIntegralSample) {
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    _data._prestencil_weights(0,0,0,0,0) = -1;
                    _data._prestencil_weights(1,0,0,0,0) = 1;
                });
            });
        } else if (_data._data_sampling_functional == ManifoldVectorPointSample) {
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    for (int j=0; j<dimensions; ++j) {
                        for (int k=0; k<dimensions-1; ++k) {
                            _data._prestencil_weights(0,target_index,0,k,j) =  T(k,j);
                        }
                    }
                });
            });
        } else if (_data._data_sampling_functional == StaggeredEdgeIntegralSample) {
            compadre_kernel_assert_debug(_data._problem_type==ProblemType::MANIFOLD 
                    && "StaggeredEdgeIntegralSample prestencil weight only written for manifolds.");
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,
                        _data._pc._nla.getNumberOfNeighborsDevice(target_index)), [&] (const int m) {
                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    for (int quadrature = 0; quadrature<_data._qm.getNumberOfQuadraturePoints(); ++quadrature) {
                        XYZ tangent_quadrature_coord_2d;
                        for (int j=0; j<dimensions-1; ++j) {
                            tangent_quadrature_coord_2d[j]  = _data._pc.getTargetCoordinate(target_index, j, &T);
                            tangent_quadrature_coord_2d[j] -= _data._pc.getNeighborCoordinate(target_index, m, j, &T);
                        }
                        double tangent_vector[3];
                        tangent_vector[0] = tangent_quadrature_coord_2d[0]*T(0,0) + tangent_quadrature_coord_2d[1]*T(1,0);
                        tangent_vector[1] = tangent_quadrature_coord_2d[0]*T(0,1) + tangent_quadrature_coord_2d[1]*T(1,1);
                        tangent_vector[2] = tangent_quadrature_coord_2d[0]*T(0,2) + tangent_quadrature_coord_2d[1]*T(1,2);

                        for (int j=0; j<dimensions; ++j) {
                            _data._prestencil_weights(0,target_index,m,0,j) +=  (1-_data._qm.getSite(quadrature,0))*tangent_vector[j]*_data._qm.getWeight(quadrature);
                            _data._prestencil_weights(1,target_index,m,0,j) +=  _data._qm.getSite(quadrature,0)*tangent_vector[j]*_data._qm.getWeight(quadrature);
                        }
                    }
                });
            });
        } else if (_data._data_sampling_functional == VaryingManifoldVectorPointSample) {

            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,
                        _data._pc._nla.getNumberOfNeighborsDevice(target_index)), [&] (const int m) {
                
                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    calcGradientPij(_data, teamMember, delta.data(), thread_workspace.data(), target_index, 
                        m, 0 /*alpha*/, 0 /*partial_direction*/, dimensions-1, _data._curvature_poly_order, 
                        false /*specific order only*/, &T, ReconstructionSpace::ScalarTaylorPolynomial, 
                        PointSample);
                });
                // reconstructs gradient at local neighbor index m
                double grad_xi1 = 0, grad_xi2 = 0;
                Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember,
                        _data._pc._nla.getNumberOfNeighborsDevice(target_index)), [&] (const int i, double &t_grad_xi1) {
                    double alpha_ij = 0;
                    for (int l=0; l<_data.manifold_NP; ++l) {
                        alpha_ij += delta(l)*Q(l,i);
                    }
                    XYZ rel_coord = _data._pc.getRelativeCoord(target_index, i, dimensions, &T);
                    double normal_coordinate = rel_coord[dimensions-1];

                    // apply coefficients to sample data
                    t_grad_xi1 += alpha_ij * normal_coordinate;
                }, grad_xi1);
                t1(m) = grad_xi1;

                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    calcGradientPij(_data, teamMember, delta.data(), thread_workspace.data(), target_index, 
                        m, 0 /*alpha*/, 1 /*partial_direction*/, dimensions-1, _data._curvature_poly_order, 
                        false /*specific order only*/, &T, ReconstructionSpace::ScalarTaylorPolynomial, 
                        PointSample);
                });
                Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember,
                            _data._pc._nla.getNumberOfNeighborsDevice(target_index)), [&] (const int i, double &t_grad_xi2) {
                    double alpha_ij = 0;
                    for (int l=0; l<_data.manifold_NP; ++l) {
                        alpha_ij += delta(l)*Q(l,i);
                    }
                    XYZ rel_coord = _data._pc.getRelativeCoord(target_index, i, dimensions, &T);
                    double normal_coordinate = rel_coord[dimensions-1];

                    // apply coefficients to sample data
                    if (dimensions>2) t_grad_xi2 += alpha_ij * normal_coordinate;
                }, grad_xi2);
                t2(m) = grad_xi2;
            });

            teamMember.team_barrier();

            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,
                        _data._pc._nla.getNumberOfNeighborsDevice(target_index)), [&] (const int m) {
                // constructs tangent vector at neighbor site
                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    for (int j=0; j<dimensions; ++j) {
                        tangent(0,j) = t1(m)*T(dimensions-1,j) + T(0,j);
                        tangent(1,j) = t2(m)*T(dimensions-1,j) + T(1,j);
                    }

                    // calculate norm
                    double norm = 0;
                    for (int j=0; j<dimensions; ++j) {
                        norm += tangent(0,j)*tangent(0,j);
                    }

                    // normalize first vector
                    norm = std::sqrt(norm);
                    for (int j=0; j<dimensions; ++j) {
                        tangent(0,j) /= norm;
                    }

                    // orthonormalize next vector
                    if (dimensions-1 == 2) { // 2d manifold
                        double dot_product = tangent(0,0)*tangent(1,0) 
                            + tangent(0,1)*tangent(1,1) 
                            + tangent(0,2)*tangent(1,2);
                        for (int j=0; j<dimensions; ++j) {
                            tangent(1,j) -= dot_product*tangent(0,j);
                        }
                        // normalize second vector
                        norm = 0;
                        for (int j=0; j<dimensions; ++j) {
                            norm += tangent(1,j)*tangent(1,j);
                        }
                        norm = std::sqrt(norm);
                        for (int j=0; j<dimensions; ++j) {
                            tangent(1,j) /= norm;
                        }
                    }

                    // stores matrix of tangent and normal directions as a prestencil weight
                    for (int j=0; j<dimensions; ++j) {
                        for (int k=0; k<dimensions-1; ++k) {
                            _data._prestencil_weights(0,target_index,m,k,j) =  tangent(k,j);
                        }
                    }
                });
            });
        }
        teamMember.team_barrier();
    }
};

//! Functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity
struct AssembleStandardPsqrtW {

    GMLSBasisData _data;

    AssembleStandardPsqrtW(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {
    
        /*
         *    Dimensions
         */
    
        const int target_index  = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index   = teamMember.league_rank();
        const int this_num_rows = _data._sampling_multiplier*_data._pc._nla.getNumberOfNeighborsDevice(target_index);
    
        /*
         *    Data
         */
    
        scratch_matrix_right_type PsqrtW(_data.P_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_dim_0*_data.P_dim_1), 
                    _data.P_dim_0, _data.P_dim_1);
        scratch_matrix_right_type RHS(_data.RHS_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.RHS_dim_0*_data.RHS_dim_1), 
                    _data.RHS_dim_0, _data.RHS_dim_1);
        scratch_vector_type w(_data.w_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.max_num_rows), 
                    _data.max_num_rows);
    
        // delta, used for each thread
        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);
    
        /*
         *    Assemble P*sqrt(W) and sqrt(w)*Identity
         */
    
        // creates the matrix sqrt(W)*P
        createWeightsAndP(_data, teamMember, delta, thread_workspace, PsqrtW, w, _data._dimensions, 
                _data._poly_order, true /*weight_p*/, NULL /*&V*/, _data._reconstruction_space, 
                _data._polynomial_sampling_functional);
    
        if ((_data._constraint_type == ConstraintType::NO_CONSTRAINT) && (_data._dense_solver_type != DenseSolverType::LU)) {
            // fill in RHS with Identity * sqrt(weights)
            double * rhs_data = RHS.data();
            Kokkos::parallel_for(Kokkos::TeamVectorRange(teamMember, this_num_rows), [&] (const int i) {
                rhs_data[i] = std::sqrt(w(i));
            });
        } else {
            // create global memory for matrix M = PsqrtW^T*PsqrtW
            // don't need to cast into scratch_matrix_left_type since the matrix is symmetric
            scratch_matrix_right_type M(_data.RHS_data
                    + TO_GLOBAL(local_index)*TO_GLOBAL(_data.RHS_dim_0*_data.RHS_dim_1), 
                        _data.RHS_dim_0, _data.RHS_dim_1);
            KokkosBatched::TeamVectorGemm<member_type,KokkosBatched::Trans::Transpose,KokkosBatched::Trans::NoTranspose,KokkosBatched::Algo::Gemm::Unblocked>
    	      ::invoke(teamMember,
    	    	   1.0,
    	    	   PsqrtW,
    	    	   PsqrtW,
    	    	   0.0,
    	    	   M);
            teamMember.team_barrier();
    
            // Multiply PsqrtW with sqrt(W) to get PW
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, _data.max_num_rows), [&] (const int i) {
                for (int j=0; j < _data.this_num_cols; j++) {
                     PsqrtW(i,j) = PsqrtW(i,j)*std::sqrt(w(i));
                }
            });
            teamMember.team_barrier();
    
            // conditionally fill in rows determined by constraint type
            if (_data._constraint_type == ConstraintType::NEUMANN_GRAD_SCALAR) {
                // normal vector is contained in last row of T
                scratch_matrix_right_type T(_data.T_data
                        + TO_GLOBAL(target_index)*TO_GLOBAL(_data._dimensions*_data._dimensions), 
                            _data._dimensions, _data._dimensions);
    
                // Get the number of neighbors for target index
                int num_neigh_target = _data._pc._nla.getNumberOfNeighborsDevice(target_index);
                double cutoff_p = _data._epsilons(target_index);
    
                evaluateConstraints(M, PsqrtW, _data._constraint_type, _data._reconstruction_space, 
                                    _data._NP, cutoff_p, _data._dimensions, num_neigh_target, &T);
            }
        }
        teamMember.team_barrier();
    }
};

//! Functor to create a coarse tangent approximation from a given neighborhood of points
struct ComputeCoarseTangentPlane {

    GMLSBasisData _data;

    // random number generator pool
    pool_type _random_number_pool;

    ComputeCoarseTangentPlane(GMLSBasisData data) : _data(data) {
        // seed random number generator pool
        _random_number_pool = pool_type(1);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index  = teamMember.league_rank();
        const int dimensions   = _data._dimensions;

        /*
         *    Data
         */

        scratch_matrix_right_type PsqrtW(_data.P_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_dim_0*_data.P_dim_1), 
                    _data.P_dim_0, _data.P_dim_1);
        scratch_vector_type w(_data.w_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.max_num_rows), 
                    _data.max_num_rows);
        scratch_matrix_right_type T(_data.T_data
                + TO_GLOBAL(target_index)*TO_GLOBAL(dimensions*dimensions),
                    dimensions, dimensions);

        scratch_matrix_right_type PTP(teamMember.team_scratch(_data._pm.getTeamScratchLevel(1)), 
                dimensions, dimensions);

        // delta, used for each thread
        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        /*
         *    Determine Coarse Approximation of Manifold Tangent Plane
         */

        // getting x y and z from which to derive a manifold
        createWeightsAndPForCurvature(_data, teamMember, delta, thread_workspace, PsqrtW, w, dimensions, true /* only specific order */);

        // create PsqrtW^T*PsqrtW
        KokkosBatched::TeamVectorGemm<member_type,KokkosBatched::Trans::Transpose,KokkosBatched::Trans::NoTranspose,KokkosBatched::Algo::Gemm::Unblocked>
	      ::invoke(teamMember,
	    	   1.0,
	    	   PsqrtW,
	    	   PsqrtW,
	    	   0.0,
	    	   PTP);
        teamMember.team_barrier();

        // create coarse approximation of tangent plane in first two rows of T, with normal direction in third column
        GMLS_LinearAlgebra::largestTwoEigenvectorsThreeByThreeSymmetric(teamMember, T, PTP, dimensions, 
                const_cast<pool_type&>(_random_number_pool));

        teamMember.team_barrier();
    }
};

//! Functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity for curvature
struct AssembleCurvaturePsqrtW {

    GMLSBasisData _data;

    AssembleCurvaturePsqrtW(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index  = teamMember.league_rank();
        const int this_num_neighbors = _data._pc._nla.getNumberOfNeighborsDevice(target_index);

        /*
         *    Data
         */

        scratch_matrix_right_type CurvaturePsqrtW(_data.P_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_dim_0*_data.P_dim_1), 
                    _data.P_dim_0, _data.P_dim_1);
        scratch_matrix_right_type RHS(_data.RHS_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.RHS_dim_0*_data.RHS_dim_1), 
                    _data.RHS_dim_0, _data.RHS_dim_1);
        scratch_vector_type w(_data.w_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.max_num_rows), _data.max_num_rows);
        scratch_matrix_right_type T(_data.T_data
                + TO_GLOBAL(target_index)*TO_GLOBAL(_data._dimensions*_data._dimensions), 
                    _data._dimensions, _data._dimensions);

        // delta, used for each thread
        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        //
        //  RECONSTRUCT ON THE TANGENT PLANE USING LOCAL COORDINATES
        //

        // creates the matrix sqrt(W)*P
        createWeightsAndPForCurvature(_data, teamMember, delta, thread_workspace, CurvaturePsqrtW, w, 
                _data._dimensions-1, false /* only specific order */, &T);

        // CurvaturePsqrtW is sized according to max_num_rows x this_num_cols of which in this case
        // we are only using this_num_neighbors x manifold_NP
        if (_data._dense_solver_type != DenseSolverType::LU) {
            // fill in RHS with Identity * sqrt(weights)
            double * rhs_data = RHS.data();
            Kokkos::parallel_for(Kokkos::TeamVectorRange(teamMember, this_num_neighbors), [&] (const int i) {
                rhs_data[i] = std::sqrt(w(i));
            });
        } else {
            // create global memory for matrix M = PsqrtW^T*PsqrtW
            // don't need to cast into scratch_matrix_left_type since the matrix is symmetric
            scratch_matrix_right_type M(_data.RHS_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.RHS_dim_0*_data.RHS_dim_1), 
                    _data.RHS_dim_0, _data.RHS_dim_1);
            // Assemble matrix M
            KokkosBatched::TeamVectorGemm<member_type,KokkosBatched::Trans::Transpose,
                KokkosBatched::Trans::NoTranspose,KokkosBatched::Algo::Gemm::Unblocked>
              ::invoke(teamMember,
                   1.0,
                   CurvaturePsqrtW,
                   CurvaturePsqrtW,
                   0.0,
                   M);
            teamMember.team_barrier();

            // Multiply PsqrtW with sqrt(W) to get PW
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, this_num_neighbors), [&] (const int i) {
                    for (int j=0; j < _data.manifold_NP; j++) {
                        CurvaturePsqrtW(i, j) = CurvaturePsqrtW(i, j)*std::sqrt(w(i));
                    }
            });
        }
        teamMember.team_barrier();
    }
};

//! Functor to evaluate curvature targets and construct accurate tangent direction approximation for manifolds
struct GetAccurateTangentDirections {

    GMLSBasisData _data;

    GetAccurateTangentDirections(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index  = teamMember.league_rank();
        auto dimensions = _data._dimensions;

        /*
         *    Data
         */

        scratch_matrix_right_type Q = scratch_matrix_right_type(_data.Coeffs_data 
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.Coeffs_dim_0*_data.Coeffs_dim_1),
                    _data.Coeffs_dim_0, _data.Coeffs_dim_1);
        scratch_matrix_right_type T(_data.T_data
                + TO_GLOBAL(target_index)*TO_GLOBAL(dimensions*dimensions), 
                    dimensions, dimensions);
        scratch_matrix_right_type P_target_row(_data.P_target_row_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_target_row_dim_0*_data.P_target_row_dim_1),
                    _data.P_target_row_dim_0, _data.P_target_row_dim_1);

        scratch_vector_type manifold_gradient(teamMember.team_scratch(_data._pm.getTeamScratchLevel(1)), 
                _data.manifold_gradient_dim);

        // delta, used for each thread
        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        /*
         *    Manifold
         */


        //
        //  GET TARGET COEFFICIENTS RELATED TO GRADIENT TERMS
        //
        // reconstruct grad_xi1 and grad_xi2, not used for manifold_coeffs
        computeCurvatureFunctionals(_data, teamMember, delta, thread_workspace, P_target_row, &T);
        teamMember.team_barrier();

        double grad_xi1 = 0, grad_xi2 = 0;
        for (int i=0; i<_data._pc._nla.getNumberOfNeighborsDevice(target_index); ++i) {
            for (int k=0; k<dimensions-1; ++k) {
                double alpha_ij = 0;
                int offset = _data._d_ss.getTargetOffsetIndex(0, 0, k, 0);
                Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember,
                        _data.manifold_NP), [&] (const int l, double &talpha_ij) {
                    Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                        talpha_ij += P_target_row(offset,l)*Q(l,i);
                    });
                }, alpha_ij);
                teamMember.team_barrier();
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    // stored staggered, grad_xi1, grad_xi2, grad_xi1, grad_xi2, ....
                    manifold_gradient(i*(dimensions-1) + k) = alpha_ij; 
                });
            }
            teamMember.team_barrier();

            XYZ rel_coord = _data._pc.getRelativeCoord(target_index, i, dimensions, &T);
            double normal_coordinate = rel_coord[dimensions-1];

            // apply coefficients to sample data
            grad_xi1 += manifold_gradient(i*(dimensions-1)) * normal_coordinate;
            if (dimensions>2) grad_xi2 += manifold_gradient(i*(dimensions-1)+1) * normal_coordinate;
            teamMember.team_barrier();
        }

        // Constructs high order orthonormal tangent space T and inverse of metric tensor
        Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

            double grad_xi[2] = {grad_xi1, grad_xi2};
            double T_row[3];

            // Construct T (high order approximation of orthonormal tangent vectors)
            for (int i=0; i<dimensions-1; ++i) {
                for (int j=0; j<dimensions; ++j) {
                    T_row[j] = T(i,j);
                }
                // build
                for (int j=0; j<dimensions; ++j) {
                    T(i,j) = grad_xi[i]*T(dimensions-1,j);
                    T(i,j) += T_row[j];
                }
            }

            // calculate norm
            double norm = 0;
            for (int j=0; j<dimensions; ++j) {
                norm += T(0,j)*T(0,j);
            }

            // normalize first vector
            norm = std::sqrt(norm);
            for (int j=0; j<dimensions; ++j) {
                T(0,j) /= norm;
            }

            // orthonormalize next vector
            if (dimensions-1 == 2) { // 2d manifold
                double dot_product = T(0,0)*T(1,0) + T(0,1)*T(1,1) + T(0,2)*T(1,2);
                for (int j=0; j<dimensions; ++j) {
                    T(1,j) -= dot_product*T(0,j);
                }
                // normalize second vector
                norm = 0;
                for (int j=0; j<dimensions; ++j) {
                    norm += T(1,j)*T(1,j);
                }
                norm = std::sqrt(norm);
                for (int j=0; j<dimensions; ++j) {
                    T(1,j) /= norm;
                }
            }

            // get normal vector to first two rows of T
            double norm_t_normal = 0;
            if (dimensions>2) {
                T(dimensions-1,0) = T(0,1)*T(1,2) - T(1,1)*T(0,2);
                norm_t_normal += T(dimensions-1,0)*T(dimensions-1,0);
                T(dimensions-1,1) = -(T(0,0)*T(1,2) - T(1,0)*T(0,2));
                norm_t_normal += T(dimensions-1,1)*T(dimensions-1,1);
                T(dimensions-1,2) = T(0,0)*T(1,1) - T(1,0)*T(0,1);
                norm_t_normal += T(dimensions-1,2)*T(dimensions-1,2);
            } else {
                T(dimensions-1,0) = T(1,1) - T(0,1);
                norm_t_normal += T(dimensions-1,0)*T(dimensions-1,0);
                T(dimensions-1,1) = T(0,0) - T(1,0);
                norm_t_normal += T(dimensions-1,1)*T(dimensions-1,1);
            }
            norm_t_normal = std::sqrt(norm_t_normal);
            for (int i=0; i<dimensions-1; ++i) {
                T(dimensions-1,i) /= norm_t_normal;
            }
        });
        teamMember.team_barrier();
    }
};

//! Functor to determine if tangent directions need reordered, and to reorder them if needed
//! We require that the normal is consistent with a right hand rule on the tangent vectors
struct FixTangentDirectionOrdering {

    GMLSBasisData _data;

    FixTangentDirectionOrdering(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        auto dimensions = _data._dimensions;

        /*
         *    Data
         */

        scratch_matrix_right_type T(_data.T_data + target_index*dimensions*dimensions, dimensions, dimensions);
        scratch_vector_type N(_data.ref_N_data + target_index*dimensions, dimensions);

        // take the dot product of the calculated normal in the tangent bundle with a given reference outward normal 
        // direction provided by the user. if the dot product is negative, flip the tangent vector ordering 
        // and flip the sign on the normal direction.
        Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
            compadre_kernel_assert_debug(dimensions > 1 
                && "FixTangentDirectionOrder called on manifold with a dimension of 0.");
            double dot_product = 0;
            for (int i=0; i<dimensions; ++i) {
                dot_product += T(dimensions-1,i) * N(i);

            }
            if (dot_product < 0) {
                if (dimensions==3) {
                    for (int i=0; i<dimensions; ++i) {
                        // swap tangent directions
                        double tmp = T(0,i);
                        T(0,i) = T(1,i);
                        T(1,i) = tmp;
                    }
                }
                for (int i=0; i<dimensions; ++i) {
                    // flip the sign of the normal vector
                    T(dimensions-1,i) *= -1;

                }
            }
        });
        teamMember.team_barrier();
    }
};

//! Functor to evaluate curvature targets and apply to coefficients of curvature reconstruction
struct ApplyCurvatureTargets {

    GMLSBasisData _data;

    ApplyCurvatureTargets(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index  = teamMember.league_rank();
        auto dimensions = _data._dimensions;

        /*
         *    Data
         */

        scratch_matrix_right_type Q = scratch_matrix_right_type(_data.Coeffs_data 
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.Coeffs_dim_0*_data.Coeffs_dim_1),
                    _data.Coeffs_dim_0, _data.Coeffs_dim_1);

        scratch_matrix_right_type T(_data.T_data
                + TO_GLOBAL(target_index)*TO_GLOBAL(dimensions*dimensions), 
                    dimensions, dimensions);

        scratch_vector_type manifold_coeffs(_data.manifold_curvature_coefficients_data
                + target_index*TO_GLOBAL(_data.manifold_NP), _data.manifold_NP);

        scratch_matrix_right_type P_target_row(_data.P_target_row_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_target_row_dim_0*_data.P_target_row_dim_1),
                    _data.P_target_row_dim_0, _data.P_target_row_dim_1);


        // delta, used for each thread
        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        /*
         *    Manifold
         */

        computeCurvatureFunctionals(_data, teamMember, delta, thread_workspace, P_target_row, &T);

        Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
            for (int j=0; j<_data.manifold_NP; ++j) { // set to zero
                manifold_coeffs(j) = 0;
            }
        });
        teamMember.team_barrier();
        for (int i=0; i<_data._pc._nla.getNumberOfNeighborsDevice(target_index); ++i) {
            XYZ rel_coord = _data._pc.getRelativeCoord(target_index, i, dimensions, &T);
            double normal_coordinate = rel_coord[dimensions-1];

            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    // coefficients without a target premultiplied
                    for (int j=0; j<_data.manifold_NP; ++j) {
                        manifold_coeffs(j) += Q(j,i) * normal_coordinate;
                    }
                });
            });
            teamMember.team_barrier();
        }
    }
};

//! Functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity
struct AssembleManifoldPsqrtW {

    GMLSBasisData _data;

    AssembleManifoldPsqrtW(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index  = teamMember.league_rank();
        auto dimensions = _data._dimensions;
        const int this_num_rows = _data._sampling_multiplier*_data._pc._nla.getNumberOfNeighborsDevice(target_index);

        /*
         *    Data
         */

        scratch_matrix_right_type PsqrtW(_data.P_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_dim_0)*TO_GLOBAL(_data.P_dim_1), 
                    _data.P_dim_0, _data.P_dim_1);
        scratch_matrix_right_type Q(_data.RHS_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.RHS_dim_0)*TO_GLOBAL(_data.RHS_dim_1), 
                    _data.RHS_dim_0, _data.RHS_dim_1);
        scratch_vector_type w(_data.w_data
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.max_num_rows), _data.max_num_rows);
        scratch_matrix_right_type T(_data.T_data
                + TO_GLOBAL(target_index)*TO_GLOBAL(dimensions)*TO_GLOBAL(dimensions), dimensions, dimensions);

        // delta, used for each thread
        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        /*
         *    Manifold
         */

        createWeightsAndP(_data, teamMember, delta, thread_workspace, PsqrtW, w, dimensions-1, 
                _data._poly_order, true /* weight with W*/, &T, _data._reconstruction_space, 
                _data._polynomial_sampling_functional);

        if (_data._dense_solver_type != DenseSolverType::LU) {
            // fill in RHS with Identity * sqrt(weights)
            double * Q_data = Q.data();
            Kokkos::parallel_for(Kokkos::TeamVectorRange(teamMember,this_num_rows), [&] (const int i) {
                Q_data[i] = std::sqrt(w(i));
            });
        } else {
            // create global memory for matrix M = PsqrtW^T*PsqrtW
            // don't need to cast into scratch_matrix_left_type since the matrix is symmetric
            scratch_matrix_right_type M(_data.RHS_data
                    + TO_GLOBAL(local_index)*TO_GLOBAL(_data.RHS_dim_0*_data.RHS_dim_1), 
                        _data.RHS_dim_0, _data.RHS_dim_1);

            // Assemble matrix M
            KokkosBatched::TeamVectorGemm<member_type,KokkosBatched::Trans::Transpose,KokkosBatched::Trans::NoTranspose,KokkosBatched::Algo::Gemm::Unblocked>
	          ::invoke(teamMember,
	        	   1.0,
	        	   PsqrtW,
	        	   PsqrtW,
	        	   0.0,
	        	   M);
            teamMember.team_barrier();

            // Multiply PsqrtW with sqrt(W) to get PW
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, _data.max_num_rows), [&] (const int i) {
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, _data.this_num_cols), [&] (const int j) {
                    PsqrtW(i, j) = PsqrtW(i, j)*std::sqrt(w(i));
                });
            });
        }
        teamMember.team_barrier();
    }
};

//! Functor to evaluate targets on a manifold
struct EvaluateManifoldTargets {

    GMLSBasisData _data;

    EvaluateManifoldTargets(GMLSBasisData data) : _data(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type& teamMember) const {

        /*
         *    Dimensions
         */

        const int target_index = _data._initial_index_for_batch + teamMember.league_rank();
        const int local_index  = teamMember.league_rank();
        auto dimensions = _data._dimensions;

        /*
         *    Data
         */

        scratch_matrix_right_type T(_data.T_data
                + TO_GLOBAL(target_index)*TO_GLOBAL(dimensions*dimensions), dimensions, dimensions);
        scratch_vector_type manifold_coeffs(_data.manifold_curvature_coefficients_data
                + target_index*TO_GLOBAL(_data.manifold_NP), _data.manifold_NP);
        scratch_matrix_right_type P_target_row(_data.P_target_row_data 
                + TO_GLOBAL(local_index)*TO_GLOBAL(_data.P_target_row_dim_0*_data.P_target_row_dim_1), 
                    _data.P_target_row_dim_0, _data.P_target_row_dim_1);

        // delta, used for each thread
        scratch_vector_type delta(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.this_num_cols);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_data._pm.getThreadScratchLevel(0)), 
                _data.thread_workspace_dim);

        /*
         *    Apply Standard Target Evaluations to Polynomial Coefficients
         */

        computeTargetFunctionalsOnManifold(_data, teamMember, delta, thread_workspace, P_target_row, T, manifold_coeffs);

    }
};

///@}

} // Compadre

#endif
