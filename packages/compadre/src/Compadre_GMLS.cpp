// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#include "Compadre_GMLS.hpp"
#include "Compadre_Functors.hpp"

namespace Compadre {

void GMLS::generatePolynomialCoefficients(const int number_of_batches, const bool keep_coefficients, const bool clear_cache) {

    compadre_assert_release( (keep_coefficients==false || number_of_batches==1)
                && "keep_coefficients is set to true, but number of batches exceeds 1.");

    /*
     *    Verify PointConnections are valid
     */
    this->verifyPointConnections(true);
    this->verifyAdditionalPointConnections(true);

    // ensure that solution set has neighbor list consistent with point connections
    _h_ss._neighbor_lists = _pc._nla;
    _h_ss._max_evaluation_sites_per_target = _additional_pc._nla.getMaxNumNeighbors() + 1;

    /*
     *    Generate Quadrature
     */
    this->_qm = Quadrature(_order_of_quadrature_points, _dimension_of_quadrature_points, _quadrature_type);

    /*
     *    Generate SolutionSet on device
     */
    this->_d_ss = SolutionSet<device_memory_space>(_h_ss);

    /*
     *    Operations to Device
     */

    // copy over operations
    _operations = decltype(_operations)("operations", _h_ss._lro.size());
    _host_operations = Kokkos::create_mirror_view(_operations);
    
    // make sure at least one target operation specified
    compadre_assert_release((_h_ss._lro.size() > 0) 
            && "No target operations added to GMLS class before calling generatePolynomialCoefficients().");
    
    // loop through list of linear reconstruction operations to be performed and set them on the host
    for (size_t i=0; i<_h_ss._lro.size(); ++i) _host_operations(i) = _h_ss._lro[i];

    // get copy of operations on the device
    Kokkos::deep_copy(_operations, _host_operations);

    /*
     *    Initialize Alphas and Prestencil Weights
     */

    // throw an assertion for QR solver incompatibility
    // TODO: this is a temporary location for this check, in the future the
    // constraint type could be an object that can check when given a dense_solver_type
    compadre_assert_release( (!(_dense_solver_type==DenseSolverType::QR
                && _constraint_type==ConstraintType::NEUMANN_GRAD_SCALAR))
            && "Cannot solve GMLS problems with the NEUMANN_GRAD_SCALAR constraint using QR Factorization.");

    // calculate the additional size for different constraint problems
    _d_ss._added_alpha_size = getAdditionalAlphaSizeFromConstraint(_dense_solver_type, _constraint_type);
    const int added_coeff_size = getAdditionalCoeffSizeFromConstraintAndSpace(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions);

    // initialize all alpha values to be used for taking the dot product with data to get a reconstruction 
    try {
        global_index_type total_neighbors = this->getNeighborLists()->getTotalNeighborsOverAllListsHost();
        int total_added_alphas = _pc._target_coordinates.extent(0)*_d_ss._added_alpha_size;
        _d_ss._alphas = 
            decltype(_d_ss._alphas)("alphas", (total_neighbors + TO_GLOBAL(total_added_alphas))
                    *TO_GLOBAL(_d_ss._total_alpha_values)*TO_GLOBAL(_d_ss._max_evaluation_sites_per_target));
        // this deep copy writes to all theoretically allocated memory,
        // ensuring that allocation attempted was successful
        Kokkos::deep_copy(_d_ss._alphas, 0.0);
    } catch(std::exception &e) {
       printf("Insufficient memory to store alphas: \n\n%s", e.what()); 
       throw e;
    }

    // initialize the prestencil weights that are applied to sampling data to put it into a form 
    // that the GMLS operator will be able to operate on
    auto sro = _data_sampling_functional;
    int max_num_neighbors = _pc._nla.getMaxNumNeighbors();
    try {
        _prestencil_weights = decltype(_prestencil_weights)("Prestencil weights",
                std::pow(2,sro.use_target_site_weights), 
                (sro.transform_type==DifferentEachTarget 
                        || sro.transform_type==DifferentEachNeighbor) ?
                    this->getNeighborLists()->getNumberOfTargets() : 1,
                (sro.transform_type==DifferentEachNeighbor) ?
                    max_num_neighbors : 1,
                (sro.output_rank>0) ?
                    _local_dimensions : 1,
                (sro.input_rank>0) ?
                    _global_dimensions : 1);
    } catch(std::exception &e) {
       printf("Insufficient memory to store prestencil weights: \n\n%s", e.what()); 
       throw e;
    }
    Kokkos::fence();

    /*
     *    Determine if Nonstandard Sampling Dimension or Basis Component Dimension
     */

    // calculate the dimension of the basis (a vector space on a manifold requires two components, for example)
    _basis_multiplier = calculateBasisMultiplier(_reconstruction_space, _local_dimensions);

    // calculate sampling dimension 
    _sampling_multiplier = calculateSamplingMultiplier(_reconstruction_space, _data_sampling_functional, _local_dimensions);

    // effective number of components in the basis
    _data_sampling_multiplier = getOutputDimensionOfSampling(_data_sampling_functional, _local_dimensions);

    // special case for using a higher order for sampling from a polynomial space that are gradients of a scalar polynomial
    if (_polynomial_sampling_functional == StaggeredEdgeAnalyticGradientIntegralSample) {
        // if the reconstruction is being made with a gradient of a basis, then we want that basis to be one order higher so that
        // the gradient is consistent with the convergence order expected.
        _poly_order += 1;
        _NP = this->getNP(_poly_order, _dimensions, _reconstruction_space);
    }

    /*
     *    Dimensions
     */

    // for tallying scratch space needed for device kernel calls
    int team_scratch_size_a = 0;

    // TEMPORARY, take to zero after conversion
    int team_scratch_size_b = 0;
    int thread_scratch_size_a = 0;
    int thread_scratch_size_b = 0;

    // dimensions that are relevant for each subproblem
    int max_num_rows = _sampling_multiplier*max_num_neighbors;
    int this_num_cols = _basis_multiplier*_NP;
    int manifold_NP = 0;

    // determines whether RHS will be stored implicitly
    // if true, instead of RHS storing the full sqrt(W) explicitly,
    // only the diagonal entries of sqrt(W) will be stored as a
    // 1D array beginning at entry with matrix coordinate (0,0)
    bool implicit_RHS = (_dense_solver_type != DenseSolverType::LU);

    int basis_powers_space_multiplier = (_reconstruction_space == BernsteinPolynomial) ? 2 : 1;
    if (_problem_type == ProblemType::MANIFOLD) {
        // these dimensions already calculated differ in the case of manifolds
        manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
        _NP = this->getNP(_poly_order, _dimensions-1, _reconstruction_space);
        const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
        this_num_cols = _basis_multiplier*max_manifold_NP;
        const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

        /*
         *    Calculate Scratch Space Allocations
         */

        team_scratch_size_b += scratch_matrix_right_type::shmem_size(_dimensions, _dimensions); // PTP matrix
        team_scratch_size_b += scratch_vector_type::shmem_size( (_dimensions-1)*max_num_neighbors ); // manifold_gradient

        thread_scratch_size_a += scratch_vector_type::shmem_size(this_num_cols); // delta, used for each thread
        thread_scratch_size_a += scratch_vector_type::shmem_size(
                (max_poly_order+1)*_global_dimensions*basis_powers_space_multiplier); // temporary space for powers in basis
        if (_data_sampling_functional == VaryingManifoldVectorPointSample) {
            team_scratch_size_b += scratch_vector_type::shmem_size(max_num_neighbors); // t1 work vector for prestencils
            team_scratch_size_b += scratch_vector_type::shmem_size(max_num_neighbors); // t2 work vector for prestencils
            thread_scratch_size_b += scratch_vector_type::shmem_size(_dimensions*_dimensions); // temporary tangent calculations, used for each thread
        }

        // allocate data on the device (initialized to zero)
        if (!_orthonormal_tangent_space_provided) {
            _T = Kokkos::View<double*>("tangent approximation",_pc._target_coordinates.extent(0)*_dimensions*_dimensions);
            Kokkos::deep_copy(_T, 0.0);
        } else {
            compadre_assert_release(_T.extent(0)/(_dimensions*_dimensions)==_pc._target_coordinates.extent(0) &&
                    "Provided tangent_directions has number of targets different than target_coordinates");
        }
        _manifold_curvature_coefficients = Kokkos::View<double*>("manifold curvature coefficients",
                _pc._target_coordinates.extent(0)*manifold_NP);
        Kokkos::deep_copy(_manifold_curvature_coefficients, 0.0);

    } else  { // Standard GMLS

        /*
         *    Calculate Scratch Space Allocations
         */

        thread_scratch_size_a += scratch_vector_type::shmem_size(this_num_cols); // delta, used for each thread
        thread_scratch_size_a += scratch_vector_type::shmem_size(
                (_poly_order+1)*_global_dimensions*basis_powers_space_multiplier); // temporary space for powers in basis

    }
    _pm.setTeamScratchSize(0, team_scratch_size_a);
    _pm.setTeamScratchSize(1, team_scratch_size_b);
    _pm.setThreadScratchSize(0, thread_scratch_size_a);
    _pm.setThreadScratchSize(1, thread_scratch_size_b);

    /*
     *    Calculate the size for matrix P and RHS
     */

    int RHS_dim_0, RHS_dim_1;
    getRHSDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, RHS_dim_0, RHS_dim_1);

    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Allocate Global Device Storage of Data Needed Over Multiple Calls
     */

    global_index_type max_batch_size = (_pc._target_coordinates.extent(0) + TO_GLOBAL(number_of_batches) - 1) / TO_GLOBAL(number_of_batches);
    try {
        _RHS = Kokkos::View<double*>("RHS", max_batch_size*TO_GLOBAL(RHS_dim_0)*TO_GLOBAL(RHS_dim_1));
        _P = Kokkos::View<double*>("P", max_batch_size*TO_GLOBAL(P_dim_0)*TO_GLOBAL(P_dim_1));
        _w = Kokkos::View<double*>("w", max_batch_size*TO_GLOBAL(max_num_rows));
        _Z = Kokkos::View<double*>("Z", max_batch_size*TO_GLOBAL(_d_ss._total_alpha_values*_d_ss._max_evaluation_sites_per_target*this_num_cols));

    } catch (std::exception &e) {
        printf("Failed to allocate space for RHS, P, and w. Consider increasing number_of_batches: \n\n%s", e.what());
        throw e;
    }
    Kokkos::fence();
    
    /*
     *    Calculate Optimal Threads Based On Levels of Parallelism
     */

    if (_constraint_type == ConstraintType::NEUMANN_GRAD_SCALAR) {
        compadre_assert_release( _orthonormal_tangent_space_provided
                && "Normal vectors are required for solving GMLS problems with the NEUMANN_GRAD_SCALAR constraint.");
    }

    _initial_index_for_batch = 0;
    for (int batch_num=0; batch_num<number_of_batches; ++batch_num) {

        auto this_batch_size = std::min(_pc._target_coordinates.extent(0)-_initial_index_for_batch, max_batch_size);
        Kokkos::deep_copy(_RHS, 0.0);
        Kokkos::deep_copy(_P, 0.0);
        Kokkos::deep_copy(_w, 0.0);
        Kokkos::deep_copy(_Z, 0.0);

        auto gmls_basis_data = this->extractBasisData();
        auto gmls_solution_data = this->extractSolutionData();

        
        // even kernels that should run on other # of vector lanes do not (on GPU)
        auto tp = _pm.TeamPolicyThreadsAndVectors(this_batch_size, _pm._default_threads, 1);
        //auto tp = _pm.TeamPolicyThreadsAndVectors(this_batch_size, _pm._default_threads, _pm._default_vector_lanes);
        //const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
        //const auto tp2 = Kokkos::Experimental::require(tp, work_item_property);

        if (_problem_type == ProblemType::MANIFOLD) {

            /*
             *    MANIFOLD Problems
             */

            //auto functor_name = Name(gmls_basis_data);
            //Kokkos::parallel_for("Name", tp, functor_name);
            if (!_orthonormal_tangent_space_provided) { // user did not specify orthonormal tangent directions, so we approximate them first
                // coarse tangent plane approximation construction of P^T*P
                auto functor_compute_coarse_tangent_plane = ComputeCoarseTangentPlane(gmls_basis_data);
                Kokkos::parallel_for("ComputeCoarseTangentPlane", tp, functor_compute_coarse_tangent_plane);

                // if the user provided the reference outward normal direction, then orient the computed or user provided
                // outward normal directions in the tangent bundle
                if (_reference_outward_normal_direction_provided && _use_reference_outward_normal_direction_provided_to_orient_surface) {
                    // use the reference outward normal direction provided by the user to orient
                    // the tangent bundle
                    auto functor_fix_tangent_direction_ordering = FixTangentDirectionOrdering(gmls_basis_data);
                    Kokkos::parallel_for("FixTangentDirectionOrdering", tp, functor_fix_tangent_direction_ordering);
                }

                // assembles the P*sqrt(weights) matrix and constructs sqrt(weights)*Identity for curvature
                auto functor_assemble_curvature_psqrtw = AssembleCurvaturePsqrtW(gmls_basis_data);
                Kokkos::parallel_for("AssembleCurvaturePsqrtW", tp, functor_assemble_curvature_psqrtw);

                if (_dense_solver_type == DenseSolverType::LU) {
                    // solves P^T*P against P^T*W with LU, stored in P
                    Kokkos::Profiling::pushRegion("Curvature LU Factorization");
                    // batchLU expects layout_left matrix tiles for B
                    // by giving it layout_right matrix tiles with reverse ordered ldb and ndb
                    // it effects a transpose of _P in layout_left
                    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(_pm, _RHS.data(), RHS_dim_0, RHS_dim_1, _P.data(), P_dim_1, P_dim_0, manifold_NP, manifold_NP, max_num_neighbors, this_batch_size, implicit_RHS);
                    Kokkos::Profiling::popRegion();
                } else {
                    // solves P*sqrt(weights) against sqrt(weights)*Identity with QR, stored in RHS
                    Kokkos::Profiling::pushRegion("Curvature QR+Pivoting Factorization");
                    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_dim_0, RHS_dim_1, max_num_neighbors, manifold_NP, max_num_neighbors, this_batch_size, implicit_RHS);
                    Kokkos::Profiling::popRegion();
                }

                // evaluates targets, applies target evaluation to polynomial coefficients for curvature
                auto functor_get_accurate_tangent_directions = GetAccurateTangentDirections(gmls_basis_data);
                Kokkos::parallel_for("GetAccurateTangentDirections", tp, functor_get_accurate_tangent_directions);

                // Due to converting layout, entries that are assumed zeros may become non-zeros.
                Kokkos::deep_copy(_P, 0.0);

                if (batch_num==number_of_batches-1) {
                    // copy tangent bundle from device back to host
                    _host_T = Kokkos::create_mirror_view(_T);
                    Kokkos::deep_copy(_host_T, _T);
                }
            }

            // this time assembling curvature PsqrtW matrix is using a highly accurate approximation of the tangent, previously calculated
            // assembles the P*sqrt(weights) matrix and constructs sqrt(weights)*Identity for curvature
            auto functor_assemble_curvature_psqrtw = AssembleCurvaturePsqrtW(gmls_basis_data);
            Kokkos::parallel_for("AssembleCurvaturePsqrtW", tp, functor_assemble_curvature_psqrtw);

            if (_dense_solver_type == DenseSolverType::LU) {
                // solves P^T*P against P^T*W with LU, stored in P
                Kokkos::Profiling::pushRegion("Curvature LU Factorization");
                GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(_pm, _RHS.data(), RHS_dim_0, RHS_dim_1, _P.data(), P_dim_1, P_dim_0, manifold_NP, manifold_NP, max_num_neighbors, this_batch_size, implicit_RHS);
                Kokkos::Profiling::popRegion();
            } else {
                 // solves P*sqrt(weights) against sqrt(weights)*Identity, stored in RHS
                Kokkos::Profiling::pushRegion("Curvature QR+Pivoting Factorization");
                GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_dim_0, RHS_dim_1, max_num_neighbors, manifold_NP, max_num_neighbors, this_batch_size, implicit_RHS);
                Kokkos::Profiling::popRegion();
            }

            // evaluates targets, applies target evaluation to polynomial coefficients for curvature
            auto functor_apply_curvature_targets = ApplyCurvatureTargets(gmls_basis_data);
            Kokkos::parallel_for("ApplyCurvatureTargets", tp, functor_apply_curvature_targets);
            Kokkos::fence();

            // prestencil weights calculated here. appropriate because:
            // precedes polynomial reconstruction from data (replaces contents of _RHS) 
            // follows reconstruction of geometry
            // calculate prestencil weights
            auto functor_compute_prestencil_weights = ComputePrestencilWeights(gmls_basis_data);
            Kokkos::parallel_for("ComputePrestencilWeights", tp, functor_compute_prestencil_weights);

            // Due to converting layout, entried that are assumed zeros may become non-zeros.
            Kokkos::deep_copy(_P, 0.0);

            // assembles the P*sqrt(weights) matrix and constructs sqrt(weights)*Identity
            auto functor_assemble_manifold_psqrtw = AssembleManifoldPsqrtW(gmls_basis_data);
            Kokkos::parallel_for("AssembleManifoldPsqrtW", tp, functor_assemble_manifold_psqrtw);

            // solves P*sqrt(weights) against sqrt(weights)*Identity, stored in RHS
            if (_dense_solver_type == DenseSolverType::LU) {
                Kokkos::Profiling::pushRegion("Manifold LU Factorization");
                GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(_pm, _RHS.data(), RHS_dim_0, RHS_dim_1, _P.data(), P_dim_1, P_dim_0, this_num_cols, this_num_cols, max_num_rows, this_batch_size, implicit_RHS);
                Kokkos::Profiling::popRegion();
            } else {
                Kokkos::Profiling::pushRegion("Manifold QR+Pivoting Factorization");
                GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_dim_0, RHS_dim_1, max_num_rows, this_num_cols, max_num_rows, this_batch_size, implicit_RHS);
                Kokkos::Profiling::popRegion();
            }
            Kokkos::fence();

        } else {

            /*
             *    STANDARD GMLS Problems
             */

            // assembles the P*sqrt(weights) matrix and constructs sqrt(weights)*Identity
            auto functor_assemble_standard_psqrtw = AssembleStandardPsqrtW(gmls_basis_data);
            //printf("size of assemble: %lu\n",  sizeof(functor_assemble_standard_psqrtw));
            Kokkos::parallel_for("AssembleStandardPsqrtW", tp, functor_assemble_standard_psqrtw);
            Kokkos::fence();

            // solves P*sqrt(weights) against sqrt(weights)*Identity, stored in RHS
            if (_dense_solver_type == DenseSolverType::LU) {
                    Kokkos::Profiling::pushRegion("LU Factorization");
                    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(_pm, _RHS.data(), RHS_dim_0, RHS_dim_1, _P.data(), P_dim_1, P_dim_0, this_num_cols + added_coeff_size, this_num_cols + added_coeff_size, max_num_rows + _d_ss._added_alpha_size, this_batch_size, implicit_RHS);
                    Kokkos::Profiling::popRegion();
            } else {
                Kokkos::Profiling::pushRegion("QR+Pivoting Factorization");
                if (_constraint_type != ConstraintType::NO_CONSTRAINT) {
                     GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(_pm, _RHS.data(), RHS_dim_0, RHS_dim_1, _P.data(), P_dim_1, P_dim_0, this_num_cols + added_coeff_size, this_num_cols + added_coeff_size, max_num_rows + _d_ss._added_alpha_size, this_batch_size, implicit_RHS);
                } else {
                    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_dim_0, RHS_dim_1, max_num_rows, this_num_cols, max_num_rows, this_batch_size, implicit_RHS);
                }
                Kokkos::Profiling::popRegion();
            }

            auto functor_compute_prestencil_weights = ComputePrestencilWeights(gmls_basis_data);
            Kokkos::parallel_for("ComputePrestencilWeights", tp, functor_compute_prestencil_weights);
            Kokkos::fence();
        }

        /*
         *    Calculate Optimal Threads Based On Levels of Parallelism
         */


        if (_problem_type == ProblemType::MANIFOLD) {

            /*
             *    MANIFOLD Problems
             */

            // evaluates targets, applies target evaluation to polynomial coefficients to store in _alphas
            auto functor_evaluate_manifold_targets = EvaluateManifoldTargets(gmls_basis_data);
            Kokkos::parallel_for("EvaluateManifoldTargets", tp, functor_evaluate_manifold_targets);

        } else {

            /*
             *    STANDARD GMLS Problems
             */

            // evaluates targets, applies target evaluation to polynomial coefficients to store in _alphas
            auto functor_evaluate_standard_targets = EvaluateStandardTargets(gmls_basis_data);
            Kokkos::parallel_for("EvaluateStandardTargets", tp, functor_evaluate_standard_targets);
        }

            
        // fine grain control over applying target (most expensive part after QR solve)
        ParallelManager pm;
        tp = pm.TeamPolicyThreadsAndVectors(this_batch_size, pm._default_threads, pm._default_vector_lanes);
        const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
        const auto tp2 = Kokkos::Experimental::require(tp, work_item_property);
        auto functor_apply_targets = ApplyTargets(gmls_solution_data);
        //printf("size of apply: %lu\n",  sizeof(functor_apply_targets));
        Kokkos::parallel_for("ApplyTargets", tp2, functor_apply_targets);


        _initial_index_for_batch += this_batch_size;
        if ((size_t)_initial_index_for_batch == _pc._target_coordinates.extent(0)) break;
    } // end of batch loops

    if (clear_cache) {
        // deallocate _P and _w
        _w = Kokkos::View<double*>("w",0);
        _Z = Kokkos::View<double*>("Z",0);
        if (number_of_batches > 1) { // no reason to keep coefficients if they aren't all in memory
            _RHS = Kokkos::View<double*>("RHS",0);
            _P = Kokkos::View<double*>("P",0);
            _entire_batch_computed_at_once = false;
        } else {
            if (_constraint_type != ConstraintType::NO_CONSTRAINT) {
                _RHS = Kokkos::View<double*>("RHS", 0);
                if (!keep_coefficients) _P = Kokkos::View<double*>("P", 0);
            } else {
                if (_dense_solver_type != DenseSolverType::LU) {
                    _P = Kokkos::View<double*>("P", 0);
                    if (!keep_coefficients) _RHS = Kokkos::View<double*>("RHS", 0);
                } else {
                    _RHS = Kokkos::View<double*>("RHS", 0);
                    if (!keep_coefficients) _P = Kokkos::View<double*>("P", 0);
                }
            }
            if (keep_coefficients) _store_PTWP_inv_PTW = true;
        }
    }

    /*
     *    Device to Host Copy Of Solution
     */
    // copy computed alphas back to the host
    this->_d_ss._contains_valid_alphas = true;
    this->_h_ss = SolutionSet<host_memory_space>(_d_ss);
    if (_data_sampling_functional != PointSample) {
        _host_prestencil_weights = Kokkos::create_mirror_view(_prestencil_weights);
        Kokkos::deep_copy(_host_prestencil_weights, _prestencil_weights);
    }

}

void GMLS::generateAlphas(const int number_of_batches, const bool keep_coefficients, const bool clear_cache) {

    this->generatePolynomialCoefficients(number_of_batches, keep_coefficients, clear_cache);

}

const GMLSSolutionData GMLS::extractSolutionData() const {
    auto gmls = *this;
    auto data = GMLSSolutionData();
    data._sampling_multiplier = gmls._sampling_multiplier;
    data._initial_index_for_batch = gmls._initial_index_for_batch;
    data._d_ss = gmls._d_ss;

    // store results of calculation in struct
    const int max_num_neighbors = gmls._pc._nla.getMaxNumNeighbors();
    const int max_num_rows = gmls._sampling_multiplier*max_num_neighbors;

    // applyTargetsToCoefficients currently uses data.this_num_cols for the
    // dot product range. Even for manifold problems, it is still appropriate to
    // use gmls._basis_multiplier*gmls._NP. If GMLSSolutionData was ever to be
    // used for applying solution coefficients for the curvature reconstruction,
    // the manifolf_NP would have to be used for that application (requiring an
    // extra argument to applyTargetsToCoefficients)
    data.this_num_cols = gmls._basis_multiplier*gmls._NP;

    int RHS_dim_0, RHS_dim_1;
    getRHSDims(gmls._dense_solver_type, gmls._constraint_type, gmls._reconstruction_space, 
            gmls._dimensions, max_num_rows, data.this_num_cols, RHS_dim_0, RHS_dim_1);
    int P_dim_0, P_dim_1;
    getPDims(gmls._dense_solver_type, gmls._constraint_type, gmls._reconstruction_space, 
            gmls._dimensions, max_num_rows, data.this_num_cols, P_dim_0, P_dim_1);

    if ((gmls._constraint_type == ConstraintType::NO_CONSTRAINT) && (gmls._dense_solver_type != DenseSolverType::LU)) {
        data.Coeffs_data = gmls._RHS.data();
        data.Coeffs_dim_0 = RHS_dim_0;
        data.Coeffs_dim_1 = RHS_dim_1;
    } else {
        data.Coeffs_data = gmls._P.data();
        data.Coeffs_dim_0 = P_dim_1;
        data.Coeffs_dim_1 = P_dim_0;
    }

    data.P_target_row_dim_0 = gmls._d_ss._total_alpha_values*gmls._d_ss._max_evaluation_sites_per_target;
    data.P_target_row_dim_1 = data.this_num_cols;
    data.P_target_row_data = gmls._Z.data();

    data.operations_size = gmls._operations.size();

    data.number_of_neighbors_list = gmls._pc._nla._number_of_neighbors_list;
    data.additional_number_of_neighbors_list = gmls._additional_pc._nla._number_of_neighbors_list;

    return data;
}

const GMLSBasisData GMLS::extractBasisData() const {
    auto gmls = *this;
    auto data = GMLSBasisData();
    data._source_extra_data = gmls._source_extra_data;
    data._target_extra_data = gmls._target_extra_data;
    data._pc = gmls._pc;
    data._epsilons  = gmls._epsilons ;
    data._prestencil_weights  = gmls._prestencil_weights ;
    data._additional_pc = gmls._additional_pc;
    data._poly_order  = gmls._poly_order ;
    data._curvature_poly_order = gmls._curvature_poly_order;
    data._NP = gmls._NP;
    data._global_dimensions = gmls._global_dimensions;
    data._local_dimensions = gmls._local_dimensions;
    data._dimensions = gmls._dimensions;
    data._reconstruction_space = gmls._reconstruction_space;
    data._reconstruction_space_rank = gmls._reconstruction_space_rank;
    data._dense_solver_type = gmls._dense_solver_type;
    data._problem_type = gmls._problem_type;
    data._constraint_type = gmls._constraint_type;
    data._polynomial_sampling_functional = gmls._polynomial_sampling_functional;
    data._data_sampling_functional = gmls._data_sampling_functional;
    data._curvature_support_operations = gmls._curvature_support_operations;
    data._operations = gmls._operations;
    data._weighting_type = gmls._weighting_type;
    data._curvature_weighting_type = gmls._curvature_weighting_type;
    data._weighting_p = gmls._weighting_p;
    data._weighting_n = gmls._weighting_n;
    data._curvature_weighting_p = gmls._curvature_weighting_p;
    data._curvature_weighting_n = gmls._curvature_weighting_n;
    data._basis_multiplier = gmls._basis_multiplier;
    data._sampling_multiplier = gmls._sampling_multiplier;
    data._data_sampling_multiplier = gmls._data_sampling_multiplier;
    data._initial_index_for_batch = gmls._initial_index_for_batch;
    data._pm = gmls._pm;
    data._order_of_quadrature_points = gmls._order_of_quadrature_points;
    data._dimension_of_quadrature_points = gmls._dimension_of_quadrature_points;
    data._qm = gmls._qm;
    data._d_ss = gmls._d_ss;

    data.max_num_neighbors = gmls._pc._nla.getMaxNumNeighbors();
    data.max_num_rows = gmls._sampling_multiplier*data.max_num_neighbors;
    int basis_powers_space_multiplier = (gmls._reconstruction_space == BernsteinPolynomial) ? 2 : 1;
    if (gmls._problem_type == ProblemType::MANIFOLD) {
        data.manifold_NP = GMLS::getNP(gmls._curvature_poly_order, gmls._dimensions-1, 
                ReconstructionSpace::ScalarTaylorPolynomial);
        data.max_manifold_NP = (data.manifold_NP > gmls._NP) ? data.manifold_NP : gmls._NP;
        data.this_num_cols = gmls._basis_multiplier*data.max_manifold_NP;
        data.max_poly_order = (gmls._poly_order > gmls._curvature_poly_order) ? 
                gmls._poly_order : gmls._curvature_poly_order;

        data.ref_N_data = gmls._ref_N.data();
        data.ref_N_dim = gmls._dimensions;

        data.thread_workspace_dim = (data.max_poly_order+1)*gmls._global_dimensions*basis_powers_space_multiplier;
        data.manifold_gradient_dim = (gmls._dimensions-1)*data.max_num_neighbors;

        data.manifold_curvature_coefficients_data = gmls._manifold_curvature_coefficients.data();

    } else {
        data.manifold_NP = 0;
        data.this_num_cols = gmls._basis_multiplier*gmls._NP;
        data.thread_workspace_dim = (gmls._poly_order+1)*gmls._global_dimensions*basis_powers_space_multiplier;
        data.manifold_gradient_dim = 0;
    }


    data.T_data = gmls._T.data();

    data.P_target_row_dim_0 = gmls._d_ss._total_alpha_values*gmls._d_ss._max_evaluation_sites_per_target;
    data.P_target_row_dim_1 = data.this_num_cols;
    data.P_target_row_data = gmls._Z.data();

    data.RHS_data = gmls._RHS.data();
    getRHSDims(gmls._dense_solver_type, gmls._constraint_type, gmls._reconstruction_space, 
               gmls._dimensions, data.max_num_rows, data.this_num_cols, data.RHS_dim_0, data.RHS_dim_1);

    data.P_data = gmls._P.data();
    getPDims(gmls._dense_solver_type, gmls._constraint_type, gmls._reconstruction_space, 
             gmls._dimensions, data.max_num_rows, data.this_num_cols, data.P_dim_0, data.P_dim_1);

    data.w_data = gmls._w.data();

    if ((gmls._constraint_type == ConstraintType::NO_CONSTRAINT) && (gmls._dense_solver_type != DenseSolverType::LU)) {
        data.Coeffs_data = gmls._RHS.data();
        data.Coeffs_dim_0 = data.RHS_dim_0;
        data.Coeffs_dim_1 = data.RHS_dim_1;
    } else {
        data.Coeffs_data = gmls._P.data();
        data.Coeffs_dim_0 = data.P_dim_1;
        data.Coeffs_dim_1 = data.P_dim_0;
    }

    return data;
}

} // Compadre
