#include "Compadre_GMLS.hpp"
#include "Compadre_GMLS_ApplyTargetEvaluations.hpp"
#include "Compadre_GMLS_Basis.hpp"
#include "Compadre_Quadrature.hpp"
#include "Compadre_GMLS_Targets.hpp"
#include "Compadre_Functors.hpp"
#include "Compadre_CreateConstraints.hpp"

namespace Compadre {

void GMLS::generatePolynomialCoefficients(const int number_of_batches, const bool keep_coefficients) {

    compadre_assert_release( (keep_coefficients==false || number_of_batches==1)
                && "keep_coefficients is set to true, but number of batches exceeds 1.");

    /*
     *    Generate Quadrature
     */
    this->_qm = Quadrature(_order_of_quadrature_points, _dimension_of_quadrature_points, _quadrature_type);

    /*
     *    Operations to Device
     */

    // copy over operations
    _operations = decltype(_operations)("operations", _lro.size());
    _host_operations = Kokkos::create_mirror_view(_operations);
    
    // make sure at least one target operation specified
    compadre_assert_release((_lro.size() > 0) 
            && "No target operations added to GMLS class before calling generatePolynomialCoefficients().");
    
    // loop through list of linear reconstruction operations to be performed and set them on the host
    for (size_t i=0; i<_lro.size(); ++i) _host_operations(i) = _lro[i];

    // get copy of operations on the device
    Kokkos::deep_copy(_operations, _host_operations);

    // check that if any target sites added, that neighbors_lists has equal rows
    compadre_assert_release(((size_t)_neighbor_lists.getNumberOfTargets()==_target_coordinates.extent(0)) 
            && "Neighbor lists not set in GMLS class before calling generatePolynomialCoefficients.");

    // check that if any target sites are greater than zero (could be zero), then there are more than zero source sites
    compadre_assert_release((_source_coordinates.extent(0)>0 || _target_coordinates.extent(0)==0) 
            && "Source coordinates not set in GMLS class before calling generatePolynomialCoefficients.");

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
    _added_alpha_size = getAdditionalAlphaSizeFromConstraint(_dense_solver_type, _constraint_type);
    const int added_coeff_size = getAdditionalCoeffSizeFromConstraintAndSpace(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions);

    // initialize all alpha values to be used for taking the dot product with data to get a reconstruction 
    try {
        global_index_type total_neighbors = _neighbor_lists.getTotalNeighborsOverAllListsHost();
        int total_added_alphas = _target_coordinates.extent(0)*_added_alpha_size;
        _alphas = decltype(_alphas)("alphas", (total_neighbors + TO_GLOBAL(total_added_alphas))
                    *TO_GLOBAL(_total_alpha_values)*TO_GLOBAL(_max_evaluation_sites_per_target));
    } catch(std::exception &e) {
       printf("Insufficient memory to store alphas: \n\n%s", e.what()); 
       throw e;
    }

    // initialize the prestencil weights that are applied to sampling data to put it into a form 
    // that the GMLS operator will be able to operate on
    auto sro = _data_sampling_functional;
    try {
        _prestencil_weights = decltype(_prestencil_weights)("Prestencil weights",
                std::pow(2,sro.use_target_site_weights), 
                (sro.transform_type==DifferentEachTarget 
                        || sro.transform_type==DifferentEachNeighbor) ?
                    _neighbor_lists.getNumberOfTargets() : 1,
                (sro.transform_type==DifferentEachNeighbor) ?
                    _max_num_neighbors : 1,
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
     *    Determine if Nontrivial Null Space in Solution
     */

    // check whether the sampling function acting on the basis will induce a nontrivial nullspace
    // an example would be reconstructing from gradient information, which would annihilate constants
    _nontrivial_nullspace = _polynomial_sampling_functional.nontrivial_nullspace;

    /*
     *    Determine if Nonstandard Sampling Dimension or Basis Component Dimension
     */


    // calculate the dimension of the basis (a vector space on a manifold requires two components, for example)
    _basis_multiplier = this->calculateBasisMultiplier(_reconstruction_space);

    // calculate sampling dimension 
    _sampling_multiplier = this->calculateSamplingMultiplier(_reconstruction_space, _data_sampling_functional);

    // effective number of components in the basis
    _data_sampling_multiplier = getOutputDimensionOfSampling(_data_sampling_functional);

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
    int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    int this_num_cols = _basis_multiplier*_NP;
    int manifold_NP = 0;

    if (_problem_type == ProblemType::MANIFOLD) {
        // these dimensions already calculated differ in the case of manifolds
        manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
        _NP = this->getNP(_poly_order, _dimensions-1, _reconstruction_space);
        const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
        this_num_cols = _basis_multiplier*max_manifold_NP;
        const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;
        const int max_P_row_size = ((_dimensions-1)*manifold_NP > max_manifold_NP*_total_alpha_values*_basis_multiplier) ? (_dimensions-1)*manifold_NP : max_manifold_NP*_total_alpha_values*_basis_multiplier*_max_evaluation_sites_per_target;

        /*
         *    Calculate Scratch Space Allocations
         */

        team_scratch_size_b += scratch_matrix_right_type::shmem_size(_dimensions-1, _dimensions-1); // G
        team_scratch_size_b += scratch_matrix_right_type::shmem_size(_dimensions, _dimensions); // PTP matrix
        team_scratch_size_b += scratch_vector_type::shmem_size( (_dimensions-1)*_max_num_neighbors ); // manifold_gradient

        team_scratch_size_b += scratch_vector_type::shmem_size(_max_num_neighbors*std::max(_sampling_multiplier,_basis_multiplier)); // t1 work vector for qr
        team_scratch_size_b += scratch_vector_type::shmem_size(_max_num_neighbors*std::max(_sampling_multiplier,_basis_multiplier)); // t2 work vector for qr

        team_scratch_size_b += scratch_vector_type::shmem_size(max_P_row_size); // row of P matrix, one for each operator
        thread_scratch_size_b += scratch_vector_type::shmem_size(max_manifold_NP*_basis_multiplier); // delta, used for each thread
        thread_scratch_size_b += scratch_vector_type::shmem_size((max_poly_order+1)*_global_dimensions); // temporary space for powers in basis
        if (_data_sampling_functional == VaryingManifoldVectorPointSample) {
            thread_scratch_size_b += scratch_vector_type::shmem_size(_dimensions*_dimensions); // temporary tangent calculations, used for each thread
        }


        // allocate data on the device (initialized to zero)
        _T = Kokkos::View<double*>("tangent approximation",_target_coordinates.extent(0)*_dimensions*_dimensions);
        _manifold_metric_tensor_inverse = Kokkos::View<double*>("manifold metric tensor inverse",_target_coordinates.extent(0)*(_dimensions-1)*(_dimensions-1));
        _manifold_curvature_coefficients = Kokkos::View<double*>("manifold curvature coefficients",_target_coordinates.extent(0)*manifold_NP);
        _manifold_curvature_gradient = Kokkos::View<double*>("manifold curvature gradient",_target_coordinates.extent(0)*(_dimensions-1));

    } else  { // Standard GMLS

        /*
         *    Calculate Scratch Space Allocations
         */

        team_scratch_size_a += scratch_vector_type::shmem_size(max_num_rows); // t1 work vector for qr
        team_scratch_size_a += scratch_vector_type::shmem_size(max_num_rows); // t2 work vector for qr

        // row of P matrix, one for each operator
        // +1 is for the original target site which always gets evaluated
        team_scratch_size_b += scratch_vector_type::shmem_size(this_num_cols*_total_alpha_values*_max_evaluation_sites_per_target); 

        thread_scratch_size_b += scratch_vector_type::shmem_size(this_num_cols); // delta, used for each thread
        thread_scratch_size_b += scratch_vector_type::shmem_size((_poly_order+1)*_global_dimensions); // temporary space for powers in basis
    }
#ifdef COMPADRE_USE_CUDA
    if (_basis_multiplier*_NP > 96) _pm.setThreadsPerTeam(_pm.getThreadsPerTeam() + 128);
#endif
    _pm.setTeamScratchSize(0, team_scratch_size_a);
    _pm.setTeamScratchSize(1, team_scratch_size_b);
    _pm.setThreadScratchSize(0, thread_scratch_size_a);
    _pm.setThreadScratchSize(1, thread_scratch_size_b);

    /*
     *    Calculate the size for matrix P and RHS
     */

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Allocate Global Device Storage of Data Needed Over Multiple Calls
     */

    global_index_type max_batch_size = (_target_coordinates.extent(0) + TO_GLOBAL(number_of_batches) - 1) / TO_GLOBAL(number_of_batches);
    try {
        _RHS = Kokkos::View<double*>("RHS", max_batch_size*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim));
        _P = Kokkos::View<double*>("P", max_batch_size*TO_GLOBAL(P_dim_0)*TO_GLOBAL(P_dim_1));
        _w = Kokkos::View<double*>("w", max_batch_size*TO_GLOBAL(max_num_rows));
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

        auto this_batch_size = std::min(_target_coordinates.extent(0)-_initial_index_for_batch, max_batch_size);
        Kokkos::deep_copy(_RHS, 0.0);
        Kokkos::deep_copy(_P, 0.0);
        Kokkos::deep_copy(_w, 0.0);
        
        if (_problem_type == ProblemType::MANIFOLD) {

            /*
             *    MANIFOLD Problems
             */

            if (!_orthonormal_tangent_space_provided) { // user did not specify orthonormal tangent directions, so we approximate them first
                // coarse tangent plane approximation construction of P^T*P
                _pm.CallFunctorWithTeamThreads<ComputeCoarseTangentPlane>(this_batch_size, *this);

                // if the user provided the reference outward normal direction, then orient the computed or user provided
                // outward normal directions in the tangent bundle
                if (_reference_outward_normal_direction_provided && _use_reference_outward_normal_direction_provided_to_orient_surface) {
                    // use the reference outward normal direction provided by the user to orient
                    // the tangent bundle
                    _pm.CallFunctorWithTeamThreads<FixTangentDirectionOrdering>(this_batch_size, *this); 
                }

                // assembles the P*sqrt(weights) matrix and constructs sqrt(weights)*Identity for curvature
                _pm.CallFunctorWithTeamThreads<AssembleCurvaturePsqrtW>(this_batch_size, *this);

                if (_dense_solver_type == DenseSolverType::LU) {
                    // solves P^T*P against P^T*W with LU, stored in P
                    Kokkos::Profiling::pushRegion("Curvature LU Factorization");
                    GMLS_LinearAlgebra::batchLUFactorize(_pm, _RHS.data(), RHS_square_dim, RHS_square_dim, _P.data(), P_dim_1, P_dim_0, manifold_NP, manifold_NP, _max_num_neighbors, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                    Kokkos::Profiling::popRegion();
                } else {
                    // solves P*sqrt(weights) against sqrt(weights)*Identity with QR, stored in RHS
                    Kokkos::Profiling::pushRegion("Curvature QR Factorization");
                    GMLS_LinearAlgebra::batchQRFactorize(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_square_dim, RHS_square_dim, _max_num_neighbors, manifold_NP, _max_num_neighbors, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                    Kokkos::Profiling::popRegion();
                }

                // evaluates targets, applies target evaluation to polynomial coefficients for curvature
                _pm.CallFunctorWithTeamThreads<GetAccurateTangentDirections>(this_batch_size, *this);

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
            _pm.CallFunctorWithTeamThreads<AssembleCurvaturePsqrtW>(this_batch_size, *this);

            if (_dense_solver_type == DenseSolverType::LU) {
                // solves P^T*P against P^T*W with LU, stored in P
                Kokkos::Profiling::pushRegion("Curvature LU Factorization");
                GMLS_LinearAlgebra::batchLUFactorize(_pm, _RHS.data(), RHS_square_dim, RHS_square_dim, _P.data(), P_dim_1, P_dim_0, manifold_NP, manifold_NP, _max_num_neighbors, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                Kokkos::Profiling::popRegion();
            } else {
                 // solves P*sqrt(weights) against sqrt(weights)*Identity, stored in RHS
                Kokkos::Profiling::pushRegion("Curvature QR Factorization");
                GMLS_LinearAlgebra::batchQRFactorize(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_square_dim, RHS_square_dim, _max_num_neighbors, manifold_NP, _max_num_neighbors, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                Kokkos::Profiling::popRegion();
           }

            // evaluates targets, applies target evaluation to polynomial coefficients for curvature
            _pm.CallFunctorWithTeamThreads<ApplyCurvatureTargets>(this_batch_size, *this);
            Kokkos::fence();

            // prestencil weights calculated here. appropriate because:
            // precedes polynomial reconstruction from data (replaces contents of _RHS) 
            // follows reconstruction of geometry
            // calculate prestencil weights
            _pm.CallFunctorWithTeamThreadsAndVectors<ComputePrestencilWeights>(this_batch_size, _pm.getThreadsPerTeam(_pm.getVectorLanesPerThread()), _pm.getVectorLanesPerThread(), *this);

            // Due to converting layout, entried that are assumed zeros may become non-zeros.
            Kokkos::deep_copy(_P, 0.0);

            // assembles the P*sqrt(weights) matrix and constructs sqrt(weights)*Identity
            _pm.CallFunctorWithTeamThreads<AssembleManifoldPsqrtW>(this_batch_size, *this);

            // solves P*sqrt(weights) against sqrt(weights)*Identity, stored in RHS
            // uses SVD if necessary or if explicitly asked to do so (much slower than QR)
            if (_nontrivial_nullspace || _dense_solver_type == DenseSolverType::SVD) {
                Kokkos::Profiling::pushRegion("Manifold SVD Factorization");
                GMLS_LinearAlgebra::batchSVDFactorize(_pm, true, _P.data(), P_dim_0, P_dim_1, true, _RHS.data(), RHS_square_dim, RHS_square_dim, max_num_rows, this_num_cols, max_num_rows, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                Kokkos::Profiling::popRegion();
            } else if (_dense_solver_type == DenseSolverType::LU) {
                Kokkos::Profiling::pushRegion("Manifold LU Factorization");
                GMLS_LinearAlgebra::batchLUFactorize(_pm, _RHS.data(), RHS_square_dim, RHS_square_dim, _P.data(), P_dim_1, P_dim_0, this_num_cols, this_num_cols, max_num_rows, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                Kokkos::Profiling::popRegion();
            } else {
                Kokkos::Profiling::pushRegion("Manifold QR Factorization");
                GMLS_LinearAlgebra::batchQRFactorize(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_square_dim, RHS_square_dim, max_num_rows, this_num_cols, max_num_rows, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                Kokkos::Profiling::popRegion();
            }
            Kokkos::fence();

        } else {

            /*
             *    STANDARD GMLS Problems
             */

            // assembles the P*sqrt(weights) matrix and constructs sqrt(weights)*Identity
            _pm.CallFunctorWithTeamThreads<AssembleStandardPsqrtW>(this_batch_size, *this);
            //_pm.CallFunctorWithTeamThreads<AssembleStandardPsqrtW>(this_batch_size, *this);
            Kokkos::fence();

            // solves P*sqrt(weights) against sqrt(weights)*Identity, stored in RHS
            // uses SVD if necessary or if explicitly asked to do so (much slower than QR)
            if (_constraint_type != ConstraintType::NO_CONSTRAINT) {
                if (_nontrivial_nullspace || _dense_solver_type == DenseSolverType::SVD) {
                     Kokkos::Profiling::pushRegion("SVD Factorization");
                     GMLS_LinearAlgebra::batchSVDFactorize(_pm, false, _RHS.data(), RHS_square_dim, RHS_square_dim, true, _P.data(), P_dim_1, P_dim_0, this_num_cols + added_coeff_size, this_num_cols + added_coeff_size, max_num_rows + _added_alpha_size, this_batch_size, 0, _initial_index_for_batch, NULL);
                     Kokkos::Profiling::popRegion();
                } else {
                     Kokkos::Profiling::pushRegion("LU Factorization");
                     GMLS_LinearAlgebra::batchLUFactorize(_pm, _RHS.data(), RHS_square_dim, RHS_square_dim, _P.data(), P_dim_1, P_dim_0, this_num_cols + added_coeff_size, this_num_cols + added_coeff_size, max_num_rows + _added_alpha_size, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                     Kokkos::Profiling::popRegion();
                }
            } else {
                // Solve normally for no-constraint cases
                if (_nontrivial_nullspace || _dense_solver_type == DenseSolverType::SVD) {
                    Kokkos::Profiling::pushRegion("SVD Factorization");
                    GMLS_LinearAlgebra::batchSVDFactorize(_pm, true, _P.data(), P_dim_0, P_dim_1, true, _RHS.data(), RHS_square_dim, RHS_square_dim, max_num_rows, this_num_cols, max_num_rows, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                    Kokkos::Profiling::popRegion();
                } else if (_dense_solver_type == DenseSolverType::LU) {
                    Kokkos::Profiling::pushRegion("LU Factorization");
                    GMLS_LinearAlgebra::batchLUFactorize(_pm, _RHS.data(), RHS_square_dim, RHS_square_dim, _P.data(), P_dim_1, P_dim_0, this_num_cols + added_coeff_size, this_num_cols + added_coeff_size, max_num_rows + _added_alpha_size, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                    Kokkos::Profiling::popRegion();
                } else {
                    Kokkos::Profiling::pushRegion("QR Factorization");
                    GMLS_LinearAlgebra::batchQRFactorize(_pm, _P.data(), P_dim_0, P_dim_1, _RHS.data(), RHS_square_dim, RHS_square_dim, max_num_rows, this_num_cols, max_num_rows, this_batch_size, _max_num_neighbors, _initial_index_for_batch, _host_number_of_neighbors_list.data());
                    Kokkos::Profiling::popRegion();
                }
            }

            _pm.CallFunctorWithTeamThreadsAndVectors<ComputePrestencilWeights>(this_batch_size, _pm.getThreadsPerTeam(_pm.getVectorLanesPerThread()), _pm.getVectorLanesPerThread(), *this);
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
            _pm.CallFunctorWithTeamThreads<ApplyManifoldTargets>(this_batch_size, *this);

        } else {

            /*
             *    STANDARD GMLS Problems
             */

            // evaluates targets, applies target evaluation to polynomial coefficients to store in _alphas
            _pm.CallFunctorWithTeamThreadsAndVectors<ApplyStandardTargets>(this_batch_size, _pm.getThreadsPerTeam(_pm.getVectorLanesPerThread()), _pm.getVectorLanesPerThread(), *this);

        }
        Kokkos::fence();
        _initial_index_for_batch += this_batch_size;
        if ((size_t)_initial_index_for_batch == _target_coordinates.extent(0)) break;
    } // end of batch loops

    // deallocate _P and _w
    _w = Kokkos::View<double*>("w",0);
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

    /*
     *    Device to Host Copy Of Solution
     */

    // copy computed alphas back to the host
    _host_alphas = Kokkos::create_mirror_view(_alphas);
    if (_data_sampling_functional != PointSample) {
        _host_prestencil_weights = Kokkos::create_mirror_view(_prestencil_weights);
        Kokkos::deep_copy(_host_prestencil_weights, _prestencil_weights);
    }
    Kokkos::deep_copy(_host_alphas, _alphas);
    Kokkos::fence();


}

void GMLS::generateAlphas(const int number_of_batches, const bool keep_coefficients) {

    this->generatePolynomialCoefficients(number_of_batches, keep_coefficients);

}


KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const AssembleStandardPsqrtW&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int this_num_rows = _sampling_multiplier*this->getNNeighbors(target_index);
    const int this_num_cols = _basis_multiplier*_NP;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */

    scratch_matrix_right_type PsqrtW(_P.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_0)*TO_GLOBAL(P_dim_1), P_dim_0, P_dim_1);
    scratch_matrix_right_type RHS(_RHS.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    scratch_vector_type w(_w.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(max_num_rows), max_num_rows);

    // delta, used for each thread
    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), this_num_cols);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (_poly_order+1)*_global_dimensions);

    /*
     *    Assemble P*sqrt(W) and sqrt(w)*Identity
     */

    // creates the matrix sqrt(W)*P
    this->createWeightsAndP(teamMember, delta, thread_workspace, PsqrtW, w, _dimensions, _poly_order, true /*weight_p*/, NULL /*&V*/, _reconstruction_space, _polynomial_sampling_functional);

    if ((_constraint_type == ConstraintType::NO_CONSTRAINT) && (_dense_solver_type != DenseSolverType::LU)) {
    // if (_dense_solver_type != DenseSolverType::LU) {
        // fill in RHS with Identity * sqrt(weights)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this_num_rows), [=] (const int i) {
            for(int j = 0; j < this_num_rows; ++j) {
                RHS(j,i) = (i==j) ? std::sqrt(w(i)) : 0;
            }
        });
    } else {
        // create global memory for matrix M = PsqrtW^T*PsqrtW
        // don't need to cast into scratch_matrix_left_type since the matrix is symmetric
        scratch_matrix_right_type M(_RHS.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
        GMLS_LinearAlgebra::createM(teamMember, M, PsqrtW, this_num_cols /* # of columns */, max_num_rows);

        // Multiply PsqrtW with sqrt(W) to get PW
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, max_num_rows), [=] (const int i) {
            for (int j=0; j < this_num_cols; j++) {
                 PsqrtW(i, j) = PsqrtW(i, j)*std::sqrt(w(i));
            }
        });
        teamMember.team_barrier();

        // conditionally fill in rows determined by constraint type
        if (_constraint_type == ConstraintType::NEUMANN_GRAD_SCALAR) {
            // normal vector is contained in last row of T
            scratch_matrix_right_type T(_T.data()
                + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

            // Get the number of neighbors for target index
            int num_neigh_target = this->getNNeighbors(target_index);
            double cutoff_p = _epsilons(target_index);

            evaluateConstraints(M, PsqrtW, _constraint_type, _reconstruction_space, _NP, cutoff_p, _dimensions, num_neigh_target, &T);
        }
    }
}


KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const ApplyStandardTargets&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int this_num_cols = _basis_multiplier*_NP;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */

    // Coefficients for polynomial basis have overwritten _RHS
    scratch_matrix_right_type Coeffs;
    // if (_dense_solver_type != DenseSolverType::LU) {
    if ((_constraint_type == ConstraintType::NO_CONSTRAINT) && (_dense_solver_type != DenseSolverType::LU)) {
        Coeffs = scratch_matrix_right_type(_RHS.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    } else {
        Coeffs = scratch_matrix_right_type(_P.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_1)*TO_GLOBAL(P_dim_0), P_dim_1, P_dim_0);
    }
    scratch_vector_type w(_w.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(max_num_rows), max_num_rows);

    scratch_vector_type t1(teamMember.team_scratch(_pm.getTeamScratchLevel(0)), max_num_rows);
    scratch_vector_type t2(teamMember.team_scratch(_pm.getTeamScratchLevel(0)), max_num_rows);
    scratch_matrix_right_type P_target_row(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _total_alpha_values*_max_evaluation_sites_per_target, this_num_cols);

    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), this_num_cols);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (_poly_order+1)*_global_dimensions);

    /*
     *    Apply Standard Target Evaluations to Polynomial Coefficients
     */

    // get evaluation of target functionals
    this->computeTargetFunctionals(teamMember, delta, thread_workspace, P_target_row);
    teamMember.team_barrier();

    this->applyTargetsToCoefficients(teamMember, t1, t2, Coeffs, w, P_target_row, _NP); 
    teamMember.team_barrier();
}


KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const ComputeCoarseTangentPlane&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
    const int this_num_cols = _basis_multiplier*max_manifold_NP;
    const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */

    scratch_matrix_right_type PsqrtW(_P.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_0)*TO_GLOBAL(P_dim_1), P_dim_0, P_dim_1);
    scratch_vector_type w(_w.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(max_num_rows), max_num_rows);
    scratch_matrix_right_type T(_T.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

    scratch_matrix_right_type PTP(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _dimensions, _dimensions);

    // delta, used for each thread
    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), max_manifold_NP*_basis_multiplier);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (max_poly_order+1)*_global_dimensions);

    /*
     *    Determine Coarse Approximation of Manifold Tangent Plane
     */

    // getting x y and z from which to derive a manifold
    this->createWeightsAndPForCurvature(teamMember, delta, thread_workspace, PsqrtW, w, _dimensions, true /* only specific order */);

    // create PsqrtW^T*PsqrtW
    GMLS_LinearAlgebra::createM(teamMember, PTP, PsqrtW, _dimensions /* # of columns */, this->getNNeighbors(target_index));

    // create coarse approximation of tangent plane in first two rows of T, with normal direction in third column
    GMLS_LinearAlgebra::largestTwoEigenvectorsThreeByThreeSymmetric(teamMember, T, PTP, _dimensions, 
            const_cast<pool_type&>(_random_number_pool));

    teamMember.team_barrier();

}

KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const AssembleCurvaturePsqrtW&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
    const int this_num_neighbors = this->getNNeighbors(target_index);
    const int this_num_cols = _basis_multiplier*max_manifold_NP;
    const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */

    scratch_matrix_right_type CurvaturePsqrtW(_P.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_0)*TO_GLOBAL(P_dim_1), P_dim_0, P_dim_1);
    scratch_matrix_right_type RHS(_RHS.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    scratch_vector_type w(_w.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(max_num_rows), max_num_rows);
    scratch_matrix_right_type T(_T.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

    // delta, used for each thread
    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), max_manifold_NP*_basis_multiplier);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (max_poly_order+1)*_global_dimensions);


    //
    //  RECONSTRUCT ON THE TANGENT PLANE USING LOCAL COORDINATES
    //

    // creates the matrix sqrt(W)*P
    this->createWeightsAndPForCurvature(teamMember, delta, thread_workspace, CurvaturePsqrtW, w, _dimensions-1, false /* only specific order */, &T);
    teamMember.team_barrier();

    // CurvaturePsqrtW is sized according to max_num_rows x this_num_cols of which in this case
    // we are only using this_num_neighbors x manifold_NP

    if (_dense_solver_type != DenseSolverType::LU) {
        // fill in RHS with Identity * sqrt(weights)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this_num_neighbors), [=] (const int i) {
            for(int j = 0; j < this_num_neighbors; ++j) {
                RHS(i, j) = (i==j) ? std::sqrt(w(i)) : 0;
            }
        });
    } else {
        // create global memory for matrix M = PsqrtW^T*PsqrtW
        // don't need to cast into scratch_matrix_left_type since the matrix is symmetric
        scratch_matrix_right_type M(_RHS.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
        // Assemble matrix M
        GMLS_LinearAlgebra::createM(teamMember, M, CurvaturePsqrtW, manifold_NP /* # of columns */, this_num_neighbors);

        // Multiply PsqrtW with sqrt(W) to get PW
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, this_num_neighbors), [=] (const int i) {
                for (int j=0; j < manifold_NP; j++) {
                    CurvaturePsqrtW(i, j) = CurvaturePsqrtW(i, j)*std::sqrt(w(i));
                }
        });
    }
    teamMember.team_barrier();
}

KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const GetAccurateTangentDirections&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
    const int this_num_cols = _basis_multiplier*max_manifold_NP;
    const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */
    scratch_matrix_right_type Q;
    if (_dense_solver_type != DenseSolverType::LU) {
        // Solution from QR comes from RHS
        Q = scratch_matrix_right_type(_RHS.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    } else {
        // Solution from LU comes from P
        Q = scratch_matrix_right_type(_P.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_1)*TO_GLOBAL(P_dim_0), P_dim_1, P_dim_0);
    }

    scratch_matrix_right_type T(_T.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

    scratch_vector_type manifold_gradient(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), (_dimensions-1)*_max_num_neighbors);
    scratch_matrix_right_type P_target_row(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _total_alpha_values*_max_evaluation_sites_per_target, max_manifold_NP*_basis_multiplier);

    // delta, used for each thread
    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), max_manifold_NP*_basis_multiplier);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (max_poly_order+1)*_global_dimensions);

    /*
     *    Manifold
     */


    //
    //  GET TARGET COEFFICIENTS RELATED TO GRADIENT TERMS
    //
    // reconstruct grad_xi1 and grad_xi2, not used for manifold_coeffs
    this->computeCurvatureFunctionals(teamMember, delta, thread_workspace, P_target_row, &T);
    teamMember.team_barrier();

    double grad_xi1 = 0, grad_xi2 = 0;
    for (int i=0; i<this->getNNeighbors(target_index); ++i) {
        for (int k=0; k<_dimensions-1; ++k) {
            double alpha_ij = 0;
            int offset = getTargetOffsetIndexDevice(0, 0, k, 0);
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember,
                    manifold_NP), [=] (const int l, double &talpha_ij) {
                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    talpha_ij += P_target_row(offset,l)*Q(l,i);
                });
            }, alpha_ij);
            teamMember.team_barrier();
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                manifold_gradient(i*(_dimensions-1) + k) = alpha_ij; // stored staggered, grad_xi1, grad_xi2, grad_xi1, grad_xi2, ....
            });
        }
        teamMember.team_barrier();

        XYZ rel_coord = getRelativeCoord(target_index, i, _dimensions, &T);
        double normal_coordinate = rel_coord[_dimensions-1];

        // apply coefficients to sample data
        grad_xi1 += manifold_gradient(i*(_dimensions-1)) * normal_coordinate;
        if (_dimensions>2) grad_xi2 += manifold_gradient(i*(_dimensions-1)+1) * normal_coordinate;
        teamMember.team_barrier();
    }

    // Constructs high order orthonormal tangent space T and inverse of metric tensor
    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

        double grad_xi[2] = {grad_xi1, grad_xi2};
        double T_row[3];

        // Construct T (high order approximation of orthonormal tangent vectors)
        for (int i=0; i<_dimensions-1; ++i) {
            for (int j=0; j<_dimensions; ++j) {
                T_row[j] = T(i,j);
            }
            // build
            for (int j=0; j<_dimensions; ++j) {
                T(i,j) = grad_xi[i]*T(_dimensions-1,j);
                T(i,j) += T_row[j];
            }
        }

        // calculate norm
        double norm = 0;
        for (int j=0; j<_dimensions; ++j) {
            norm += T(0,j)*T(0,j);
        }

        // normalize first vector
        norm = std::sqrt(norm);
        for (int j=0; j<_dimensions; ++j) {
            T(0,j) /= norm;
        }

        // orthonormalize next vector
        if (_dimensions-1 == 2) { // 2d manifold
            double dot_product = T(0,0)*T(1,0) + T(0,1)*T(1,1) + T(0,2)*T(1,2);
            for (int j=0; j<_dimensions; ++j) {
                T(1,j) -= dot_product*T(0,j);
            }
            // normalize second vector
            norm = 0;
            for (int j=0; j<_dimensions; ++j) {
                norm += T(1,j)*T(1,j);
            }
            norm = std::sqrt(norm);
            for (int j=0; j<_dimensions; ++j) {
                T(1,j) /= norm;
            }
        }

        // get normal vector to first two rows of T
        double norm_t_normal = 0;
        if (_dimensions>2) {
            T(_dimensions-1,0) = T(0,1)*T(1,2) - T(1,1)*T(0,2);
            norm_t_normal += T(_dimensions-1,0)*T(_dimensions-1,0);
            T(_dimensions-1,1) = -(T(0,0)*T(1,2) - T(1,0)*T(0,2));
            norm_t_normal += T(_dimensions-1,1)*T(_dimensions-1,1);
            T(_dimensions-1,2) = T(0,0)*T(1,1) - T(1,0)*T(0,1);
            norm_t_normal += T(_dimensions-1,2)*T(_dimensions-1,2);
        } else {
            T(_dimensions-1,0) = T(1,1) - T(0,1);
            norm_t_normal += T(_dimensions-1,0)*T(_dimensions-1,0);
            T(_dimensions-1,1) = T(0,0) - T(1,0);
            norm_t_normal += T(_dimensions-1,1)*T(_dimensions-1,1);
        }
        norm_t_normal = std::sqrt(norm_t_normal);
        for (int i=0; i<_dimensions-1; ++i) {
            T(_dimensions-1,i) /= norm_t_normal;
        }
    });
    teamMember.team_barrier();
}

KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const FixTangentDirectionOrdering&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();

    /*
     *    Data
     */

    scratch_matrix_right_type T(_T.data() + target_index*_dimensions*_dimensions, _dimensions, _dimensions);
    scratch_vector_type N(_ref_N.data() + target_index*_dimensions, _dimensions);

    // take the dot product of the calculated normal in the tangent bundle with a given reference outward normal 
    // direction provided by the user. if the dot product is negative, flip the tangent vector ordering 
    // and flip the sign on the normal direction.
    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
        compadre_kernel_assert_debug(_dimensions > 1 && "FixTangentDirectionOrder called on manifold with a dimension of 0.");
        double dot_product = 0;
        for (int i=0; i<_dimensions; ++i) {
            dot_product += T(_dimensions-1,i) * N(i);

        }
        if (dot_product < 0) {
            if (_dimensions==3) {
                for (int i=0; i<_dimensions; ++i) {
                    // swap tangent directions
                    double tmp = T(0,i);
                    T(0,i) = T(1,i);
                    T(1,i) = tmp;
                }
            }
            for (int i=0; i<_dimensions; ++i) {
                // flip the sign of the normal vector
                T(_dimensions-1,i) *= -1;

            }
        }
    });
    teamMember.team_barrier();
}

KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const ApplyCurvatureTargets&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
    const int this_num_cols = _basis_multiplier*max_manifold_NP;
    const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */
    scratch_matrix_right_type Q;
    if (_dense_solver_type != DenseSolverType::LU) {
        // Solution from QR comes from RHS
        Q = scratch_matrix_right_type(_RHS.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    } else {
        // Solution from LU comes from P
        Q = scratch_matrix_right_type(_P.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_1)*TO_GLOBAL(P_dim_0), P_dim_1, P_dim_0);
    }

    scratch_matrix_right_type T(_T.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

    scratch_matrix_right_type G_inv(_manifold_metric_tensor_inverse.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL((_dimensions-1))*TO_GLOBAL((_dimensions-1)), _dimensions-1, _dimensions-1);
    scratch_vector_type manifold_coeffs(_manifold_curvature_coefficients.data() + target_index*manifold_NP, manifold_NP);
    scratch_vector_type manifold_gradient_coeffs(_manifold_curvature_gradient.data() + target_index*(_dimensions-1), (_dimensions-1));

    scratch_matrix_right_type G(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _dimensions-1, _dimensions-1);
    scratch_vector_type manifold_gradient(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), (_dimensions-1)*_max_num_neighbors);
    scratch_matrix_right_type P_target_row(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _total_alpha_values*_max_evaluation_sites_per_target, max_manifold_NP*_basis_multiplier);

    // delta, used for each thread
    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), max_manifold_NP*_basis_multiplier);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (max_poly_order+1)*_global_dimensions);

    /*
     *    Manifold
     */


    //
    //  GET TARGET COEFFICIENTS RELATED TO GRADIENT TERMS
    //
    // reconstruct grad_xi1 and grad_xi2, not used for manifold_coeffs
    this->computeCurvatureFunctionals(teamMember, delta, thread_workspace, P_target_row, &T);
    teamMember.team_barrier();

    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
        for (int j=0; j<manifold_NP; ++j) { // set to zero
            manifold_coeffs(j) = 0;
        }
    });
    teamMember.team_barrier();

    double grad_xi1 = 0, grad_xi2 = 0;
    for (int i=0; i<this->getNNeighbors(target_index); ++i) {
        for (int k=0; k<_dimensions-1; ++k) {
            double alpha_ij = 0;
            int offset = getTargetOffsetIndexDevice(0, 0, k, 0);
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember,
                    manifold_NP), [=] (const int l, double &talpha_ij) {
                Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                    talpha_ij += P_target_row(offset,l)*Q(l,i);
                });
            }, alpha_ij);
            teamMember.team_barrier();
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                manifold_gradient(i*(_dimensions-1) + k) = alpha_ij; // stored staggered, grad_xi1, grad_xi2, grad_xi1, grad_xi2, ....
            });
        }
        teamMember.team_barrier();

        XYZ rel_coord = getRelativeCoord(target_index, i, _dimensions, &T);
        double normal_coordinate = rel_coord[_dimensions-1];

        // apply coefficients to sample data
        grad_xi1 += manifold_gradient(i*(_dimensions-1)) * normal_coordinate;
        if (_dimensions>2) grad_xi2 += manifold_gradient(i*(_dimensions-1)+1) * normal_coordinate;

        Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                // coefficients without a target premultiplied
                for (int j=0; j<manifold_NP; ++j) {
                    manifold_coeffs(j) += Q(j,i) * normal_coordinate;
                }
            });
        });
        teamMember.team_barrier();
    }

    // Constructs high order orthonormal tangent space T and inverse of metric tensor
    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

        manifold_gradient_coeffs(0) = grad_xi1;
        if (_dimensions>2) manifold_gradient_coeffs(1) = grad_xi2;

        // need to get 2x2 matrix of metric tensor
        G(0,0) = 1 + grad_xi1*grad_xi1;

        if (_dimensions>2) {
            G(0,1) = grad_xi1*grad_xi2;
            G(1,0) = grad_xi2*grad_xi1;
            G(1,1) = 1 + grad_xi2*grad_xi2;
        }

        double G_determinant;
        if (_dimensions==2) {
            G_determinant = G(0,0);
            compadre_kernel_assert_debug((G_determinant!=0) && "Determinant is zero.");
            G_inv(0,0) = 1/G_determinant;
        } else {
            G_determinant = G(0,0)*G(1,1) - G(0,1)*G(1,0); //std::sqrt(G_inv(0,0)*G_inv(1,1) - G_inv(0,1)*G_inv(1,0));
            compadre_kernel_assert_debug((G_determinant!=0) && "Determinant is zero.");
            {
                // inverse of 2x2
                G_inv(0,0) = G(1,1)/G_determinant;
                G_inv(1,1) = G(0,0)/G_determinant;
                G_inv(0,1) = -G(0,1)/G_determinant;
                G_inv(1,0) = -G(1,0)/G_determinant;
            }
        }

    });
    teamMember.team_barrier();
    //
    //  END OF MANIFOLD METRIC CALCULATIONS
    //
}

KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const AssembleManifoldPsqrtW&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
    const int this_num_rows = _sampling_multiplier*this->getNNeighbors(target_index);
    const int this_num_cols = _basis_multiplier*max_manifold_NP;
    const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */

    scratch_matrix_right_type PsqrtW(_P.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_0)*TO_GLOBAL(P_dim_1), P_dim_0, P_dim_1);
    scratch_matrix_right_type Q(_RHS.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    scratch_vector_type w(_w.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(max_num_rows), max_num_rows);
    scratch_matrix_right_type T(_T.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

    // delta, used for each thread
    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), max_manifold_NP*_basis_multiplier);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (max_poly_order+1)*_global_dimensions);


    /*
     *    Manifold
     */

    this->createWeightsAndP(teamMember, delta, thread_workspace, PsqrtW, w, _dimensions-1, _poly_order, true /* weight with W*/, &T, _reconstruction_space, _polynomial_sampling_functional);
    teamMember.team_barrier();

    if (_dense_solver_type != DenseSolverType::LU) {
        // fill in RHS with Identity * sqrt(weights)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this_num_rows), [=] (const int i) {
            for(int j = 0; j < this_num_rows; ++j) {
                Q(i, j) = (i==j) ? std::sqrt(w(i)) : 0;
            }
        });
    } else {
        // create global memory for matrix M = PsqrtW^T*PsqrtW
        // don't need to cast into scratch_matrix_left_type since the matrix is symmetric
        scratch_matrix_right_type M(_RHS.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
        // Assemble matrix M
        GMLS_LinearAlgebra::createM(teamMember, M, PsqrtW, this_num_cols /* # of columns */, max_num_rows);

        // Multiply PsqrtW with sqrt(W) to get PW
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, max_num_rows), [=] (const int i) {
            for (int j=0; j < this_num_cols; j++) {
                PsqrtW(i, j) = PsqrtW(i, j)*std::sqrt(w(i));
            }
        });
    }
    teamMember.team_barrier();
}

KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const ApplyManifoldTargets&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
    const int this_num_cols = _basis_multiplier*max_manifold_NP;
    const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */

    scratch_matrix_right_type Coeffs;
    if (_dense_solver_type != DenseSolverType::LU) {
        Coeffs = scratch_matrix_right_type(_RHS.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    } else {
        Coeffs = scratch_matrix_right_type(_P.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_1)*TO_GLOBAL(P_dim_0), P_dim_1, P_dim_0);
    }
    scratch_vector_type w(_w.data() 
            + TO_GLOBAL(local_index)*TO_GLOBAL(max_num_rows), max_num_rows);
    scratch_matrix_right_type T(_T.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

    scratch_matrix_right_type G_inv(_manifold_metric_tensor_inverse.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL((_dimensions-1))*TO_GLOBAL((_dimensions-1)), _dimensions-1, _dimensions-1);
    scratch_vector_type manifold_coeffs(_manifold_curvature_coefficients.data() + target_index*manifold_NP, manifold_NP);
    scratch_vector_type manifold_gradient_coeffs(_manifold_curvature_gradient.data() + target_index*(_dimensions-1), (_dimensions-1));

    scratch_vector_type t1(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _max_num_neighbors*((_sampling_multiplier>_basis_multiplier) ? _sampling_multiplier : _basis_multiplier));
    scratch_vector_type t2(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _max_num_neighbors*((_sampling_multiplier>_basis_multiplier) ? _sampling_multiplier : _basis_multiplier));
    scratch_matrix_right_type P_target_row(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _total_alpha_values*_max_evaluation_sites_per_target, max_manifold_NP*_basis_multiplier);

    // delta, used for each thread
    scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), max_manifold_NP*_basis_multiplier);
    scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (max_poly_order+1)*_global_dimensions);

    /*
     *    Apply Standard Target Evaluations to Polynomial Coefficients
     */

    this->computeTargetFunctionalsOnManifold(teamMember, delta, thread_workspace, P_target_row, T, G_inv, manifold_coeffs, manifold_gradient_coeffs);
    teamMember.team_barrier();

    this->applyTargetsToCoefficients(teamMember, t1, t2, Coeffs, w, P_target_row, _NP); 

    teamMember.team_barrier();
}


KOKKOS_INLINE_FUNCTION
void GMLS::operator()(const ComputePrestencilWeights&, const member_type& teamMember) const {

    /*
     *    Dimensions
     */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
    const int local_index  = teamMember.league_rank();

    const int max_num_rows = _sampling_multiplier*_max_num_neighbors;
    const int manifold_NP = this->getNP(_curvature_poly_order, _dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    const int max_manifold_NP = (manifold_NP > _NP) ? manifold_NP : _NP;
    const int this_num_cols = _basis_multiplier*max_manifold_NP;
    const int max_poly_order = (_poly_order > _curvature_poly_order) ? _poly_order : _curvature_poly_order;

    int RHS_square_dim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols);
    int P_dim_0, P_dim_1;
    getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, max_num_rows, this_num_cols, P_dim_0, P_dim_1);

    /*
     *    Data
     */


    scratch_vector_type t1(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _max_num_neighbors*((_sampling_multiplier>_basis_multiplier) ? _sampling_multiplier : _basis_multiplier));
    scratch_vector_type t2(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _max_num_neighbors*((_sampling_multiplier>_basis_multiplier) ? _sampling_multiplier : _basis_multiplier));

    // holds polynomial coefficients for curvature reconstruction
    scratch_matrix_right_type Q;
    if ((_constraint_type == ConstraintType::NO_CONSTRAINT) && (_dense_solver_type != DenseSolverType::LU)) {
        // Solution from QR comes from RHS
        Q = scratch_matrix_right_type(_RHS.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(RHS_square_dim)*TO_GLOBAL(RHS_square_dim), RHS_square_dim, RHS_square_dim);
    } else {
        // Solution from LU comes from P
        Q = scratch_matrix_right_type(_P.data()
            + TO_GLOBAL(local_index)*TO_GLOBAL(P_dim_1)*TO_GLOBAL(P_dim_0), P_dim_1, P_dim_0);
    }

    scratch_matrix_right_type T(_T.data() 
            + TO_GLOBAL(target_index)*TO_GLOBAL(_dimensions)*TO_GLOBAL(_dimensions), _dimensions, _dimensions);

    scratch_vector_type manifold_gradient(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), (_dimensions-1)*_max_num_neighbors);
    scratch_matrix_right_type P_target_row(teamMember.team_scratch(_pm.getTeamScratchLevel(1)), _total_alpha_values*_max_evaluation_sites_per_target, max_manifold_NP*_basis_multiplier);

    /*
     *    Prestencil Weight Calculations
     */

    if (_data_sampling_functional == StaggeredEdgeAnalyticGradientIntegralSample) {
        Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
            Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                _prestencil_weights(0,0,0,0,0) = -1;
                _prestencil_weights(1,0,0,0,0) = 1;
            });
        });
    } else if (_data_sampling_functional == ManifoldVectorPointSample) {
        Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
            Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                for (int j=0; j<_dimensions; ++j) {
                    for (int k=0; k<_dimensions-1; ++k) {
                        _prestencil_weights(0,target_index,0,k,j) =  T(k,j);
                    }
                }
            });
        });
    } else if (_data_sampling_functional == StaggeredEdgeIntegralSample) {
        compadre_kernel_assert_debug(_problem_type==ProblemType::MANIFOLD && "StaggeredEdgeIntegralSample prestencil weight only written for manifolds.");
        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this->getNNeighbors(target_index)), [=] (const int m) {
            Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                for (int quadrature = 0; quadrature<_qm.getNumberOfQuadraturePoints(); ++quadrature) {
                    XYZ tangent_quadrature_coord_2d;
                    for (int j=0; j<_dimensions-1; ++j) {
                        tangent_quadrature_coord_2d[j] = getTargetCoordinate(target_index, j, &T);
                        tangent_quadrature_coord_2d[j] -= getNeighborCoordinate(target_index, m, j, &T);
                    }
                    double tangent_vector[3];
                    tangent_vector[0] = tangent_quadrature_coord_2d[0]*T(0,0) + tangent_quadrature_coord_2d[1]*T(1,0);
                    tangent_vector[1] = tangent_quadrature_coord_2d[0]*T(0,1) + tangent_quadrature_coord_2d[1]*T(1,1);
                    tangent_vector[2] = tangent_quadrature_coord_2d[0]*T(0,2) + tangent_quadrature_coord_2d[1]*T(1,2);

                    for (int j=0; j<_dimensions; ++j) {
                        _prestencil_weights(0,target_index,m,0,j) +=  (1-_qm.getSite(quadrature,0))*tangent_vector[j]*_qm.getWeight(quadrature);
                        _prestencil_weights(1,target_index,m,0,j) +=  _qm.getSite(quadrature,0)*tangent_vector[j]*_qm.getWeight(quadrature);
                    }
                }
            });
        });
    } else if (_data_sampling_functional == VaryingManifoldVectorPointSample) {

        scratch_vector_type delta(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), manifold_NP);
        scratch_vector_type thread_workspace(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), (max_poly_order+1)*_global_dimensions);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this->getNNeighbors(target_index)), [&] (const int m) {
            
            Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                this->calcGradientPij(teamMember, delta.data(), thread_workspace.data(), target_index, m, 0 /*alpha*/, 0 /*partial_direction*/, _dimensions-1, _curvature_poly_order, false /*specific order only*/, &T, ReconstructionSpace::ScalarTaylorPolynomial, PointSample);
            });
            // reconstructs gradient at local neighbor index m
            double grad_xi1 = 0, grad_xi2 = 0;
            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember,this->getNNeighbors(target_index)), [=] (const int i, double &t_grad_xi1) {
                double alpha_ij = 0;
                for (int l=0; l<manifold_NP; ++l) {
                    alpha_ij += delta(l)*Q(l,i);
                }
                XYZ rel_coord = getRelativeCoord(target_index, i, _dimensions, &T);
                double normal_coordinate = rel_coord[_dimensions-1];

                // apply coefficients to sample data
                t_grad_xi1 += alpha_ij * normal_coordinate;
            }, grad_xi1);
            t1(m) = grad_xi1;

            Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                this->calcGradientPij(teamMember, delta.data(), thread_workspace.data(), target_index, m, 0 /*alpha*/, 1 /*partial_direction*/, _dimensions-1, _curvature_poly_order, false /*specific order only*/, &T, ReconstructionSpace::ScalarTaylorPolynomial, PointSample);
            });
            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember,this->getNNeighbors(target_index)), [=] (const int i, double &t_grad_xi2) {
                double alpha_ij = 0;
                for (int l=0; l<manifold_NP; ++l) {
                    alpha_ij += delta(l)*Q(l,i);
                }
                XYZ rel_coord = getRelativeCoord(target_index, i, _dimensions, &T);
                double normal_coordinate = rel_coord[_dimensions-1];

                // apply coefficients to sample data
                if (_dimensions>2) t_grad_xi2 += alpha_ij * normal_coordinate;
            }, grad_xi2);
            t2(m) = grad_xi2;
        });

        teamMember.team_barrier();

        scratch_matrix_right_type tangent(teamMember.thread_scratch(_pm.getThreadScratchLevel(1)), _dimensions-1, _dimensions);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this->getNNeighbors(target_index)), [=] (const int m) {
            // constructs tangent vector at neighbor site
            Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
                for (int j=0; j<_dimensions; ++j) {
                    tangent(0,j) = t1(m)*T(_dimensions-1,j) + T(0,j);
                    tangent(1,j) = t2(m)*T(_dimensions-1,j) + T(1,j);
                }

                // calculate norm
                double norm = 0;
                for (int j=0; j<_dimensions; ++j) {
                    norm += tangent(0,j)*tangent(0,j);
                }

                // normalize first vector
                norm = std::sqrt(norm);
                for (int j=0; j<_dimensions; ++j) {
                    tangent(0,j) /= norm;
                }

                // orthonormalize next vector
                if (_dimensions-1 == 2) { // 2d manifold
                    double dot_product = tangent(0,0)*tangent(1,0) + tangent(0,1)*tangent(1,1) + tangent(0,2)*tangent(1,2);
                    for (int j=0; j<_dimensions; ++j) {
                        tangent(1,j) -= dot_product*tangent(0,j);
                    }
                    // normalize second vector
                    norm = 0;
                    for (int j=0; j<_dimensions; ++j) {
                        norm += tangent(1,j)*tangent(1,j);
                    }
                    norm = std::sqrt(norm);
                    for (int j=0; j<_dimensions; ++j) {
                        tangent(1,j) /= norm;
                    }
                }

                // stores matrix of tangent and normal directions as a prestencil weight
                for (int j=0; j<_dimensions; ++j) {
                    for (int k=0; k<_dimensions-1; ++k) {
                        _prestencil_weights(0,target_index,m,k,j) =  tangent(k,j);
                    }
                }
            });
        });
    }
    teamMember.team_barrier();
}

} // Compadre
