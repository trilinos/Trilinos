#ifndef _COMPADRE_GMLS_HPP_
#define _COMPADRE_GMLS_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"

#include "Compadre_Misc.hpp"
#include "Compadre_Operators.hpp"
#include "Compadre_LinearAlgebra_Definitions.hpp"
#include "Compadre_ParallelManager.hpp"
#include "Compadre_Quadrature.hpp"
#include "Compadre_ScalarTaylorPolynomial.hpp"
#include "Compadre_DivergenceFreePolynomial.hpp"
#include "Compadre_NeighborLists.hpp"

namespace Compadre {

//!  Generalized Moving Least Squares (GMLS)
/*!
*  This class sets up a batch of GMLS problems from a given set of neighbor lists, target sites, and source sites.
*  GMLS requires a target functional, reconstruction space, and sampling functional to be specified. 
*  For a given choice of reconstruction space and sampling functional, multiple targets can be generated with very little
*  additional computation, which is why this class allows for multiple target functionals to be specified. 
*/
class GMLS {
protected:

    // random numbe generator pool
    pool_type _random_number_pool;

    // matrices that may be needed for matrix factorization on the device
    // supports batched matrix factorization dispatch

    //! contains weights for all problems
    Kokkos::View<double*> _w; 

    //! P*sqrt(w) matrix for all problems
    Kokkos::View<double*> _P;

    //! sqrt(w)*Identity matrix for all problems, later holds polynomial coefficients for all problems
    Kokkos::View<double*> _RHS;

    //! Rank 3 tensor for high order approximation of tangent vectors for all problems. First rank is
    //! for the target index, the second is for the local direction to the manifolds 0..(_dimensions-1)
    //! are tangent, _dimensions is the normal, and the third is for the spatial dimension (_dimensions)
    Kokkos::View<double*> _T;

    //! Rank 2 tensor for high order approximation of tangent vectors for all problems. First rank is
    //! for the target index, the second is for the spatial dimension (_dimensions)
    Kokkos::View<double*> _ref_N;

    //! tangent vectors information (host)
    Kokkos::View<double*>::HostMirror _host_T;

    //! reference outward normal vectors information (host)
    Kokkos::View<double*>::HostMirror _host_ref_N;

    //! metric tensor inverse for all problems
    Kokkos::View<double*> _manifold_metric_tensor_inverse;

    //! curvature polynomial coefficients for all problems
    Kokkos::View<double*> _manifold_curvature_coefficients;

    //! _dimension-1 gradient values for curvature for all problems
    Kokkos::View<double*> _manifold_curvature_gradient;

    //! Extra data available to basis functions (optional)
    Kokkos::View<double**, layout_right> _source_extra_data;

    //! Extra data available to target operations (optional)
    Kokkos::View<double**, layout_right> _target_extra_data;

    //! Accessor to get neighbor list data, offset data, and number of neighbors per target
    NeighborLists<Kokkos::View<int*> > _neighbor_lists; 
    
    //! convenient copy on host of number of neighbors
    Kokkos::View<int*, host_memory_space> _host_number_of_neighbors_list; 

    //! all coordinates for the source for which _neighbor_lists refers (device)
    Kokkos::View<double**, layout_right> _source_coordinates; 

    //! coordinates for target sites for reconstruction (device)
    Kokkos::View<double**, layout_right> _target_coordinates; 

    //! h supports determined through neighbor search (device)
    Kokkos::View<double*> _epsilons; 

    //! h supports determined through neighbor search (host)
    Kokkos::View<double*>::HostMirror _host_epsilons; 

    //! generated alpha coefficients (device)
    Kokkos::View<double*, layout_right> _alphas; 

    //! generated alpha coefficients (host)
    Kokkos::View<const double*, layout_right>::HostMirror _host_alphas;
    
    //! generated weights for nontraditional samples required to transform data into expected sampling 
    //! functional form (device). 
    Kokkos::View<double*****, layout_right> _prestencil_weights; 

    //! generated weights for nontraditional samples required to transform data into expected sampling 
    //! functional form (host)
    Kokkos::View<const double*****, layout_right>::HostMirror _host_prestencil_weights;

    //! (OPTIONAL) user provided additional coordinates for target operation evaluation (device)
    Kokkos::View<double**, layout_right> _additional_evaluation_coordinates; 

    //! (OPTIONAL) contains indices of entries in the _additional_evaluation_coordinates view (device)
    Kokkos::View<int**, layout_right> _additional_evaluation_indices; 

    //! (OPTIONAL) contains indices of entries in the _additional_evaluation_coordinates view (host)
    Kokkos::View<int**, layout_right>::HostMirror _host_additional_evaluation_indices;

    //! (OPTIONAL) contains the # of additional coordinate indices for each target
    Kokkos::View<int*> _number_of_additional_evaluation_indices; 


    //! order of basis for polynomial reconstruction
    int _poly_order; 

    //! order of basis for curvature reconstruction
    int _curvature_poly_order;

    //! dimension of basis for polynomial reconstruction
    int _NP;

    //! spatial dimension of the points, set at class instantiation only
    int _global_dimensions;

    //! dimension of the problem, set at class instantiation only. For manifolds, generally _global_dimensions-1
    int _local_dimensions;

    //! dimension of the problem, set at class instantiation only
    int _dimensions;

    //! reconstruction space for GMLS problems, set at GMLS class instantiation
    ReconstructionSpace _reconstruction_space;

    //! actual rank of reconstruction basis
    int _reconstruction_space_rank;

    //! solver type for GMLS problem - can be QR, SVD or LU
    DenseSolverType _dense_solver_type;

    //! problem type for GMLS problem, can also be set to STANDARD for normal or MANIFOLD for manifold problems
    ProblemType _problem_type;

    //! constraint type for GMLS problem
    ConstraintType _constraint_type;

    //! polynomial sampling functional used to construct P matrix, set at GMLS class instantiation
    const SamplingFunctional _polynomial_sampling_functional;

    //! generally the same as _polynomial_sampling_functional, but can differ if specified at 
    //! GMLS class instantiation
    const SamplingFunctional _data_sampling_functional;

    //! vector containing target functionals to be applied for curvature
    Kokkos::View<TargetOperation*> _curvature_support_operations;

    //! vector containing target functionals to be applied for reconstruction problem (device)
    Kokkos::View<TargetOperation*> _operations;

    //! vector containing target functionals to be applied for reconstruction problem (host)
    Kokkos::View<TargetOperation*>::HostMirror _host_operations;



    //! weighting kernel type for GMLS
    WeightingFunctionType _weighting_type;

    //! weighting kernel type for curvature problem
    WeightingFunctionType _curvature_weighting_type;

    //! power to be used for weighting kernel
    int _weighting_power;

    //! power to be used for weighting kernel for curvature
    int _curvature_weighting_power;

    //! dimension of the reconstructed function 
    //! e.g. reconstruction of vector on a 2D manifold in 3D would have _basis_multiplier of 2
    int _basis_multiplier;

    //! actual dimension of the sampling functional
    //! e.g. reconstruction of vector on a 2D manifold in 3D would have _basis_multiplier of 2
    //! e.g. in 3D, a scalar will be 1, a vector will be 3, and a vector of reused scalars will be 1
    int _sampling_multiplier;

    //! effective dimension of the data sampling functional
    //! e.g. in 3D, a scalar will be 1, a vector will be 3, and a vector of reused scalars will be 3
    int _data_sampling_multiplier;

    //! whether or not operator to be inverted for GMLS problem has a nontrivial nullspace (requiring SVD)
    bool _nontrivial_nullspace;

    //! whether or not the orthonormal tangent directions were provided by the user. If they are not,
    //! then for the case of calculations on manifolds, a GMLS approximation of the tangent space will
    //! be made and stored for use.
    bool _orthonormal_tangent_space_provided; 

    //! whether or not the reference outward normal directions were provided by the user. 
    bool _reference_outward_normal_direction_provided;

    //! whether or not to use reference outward normal directions to orient the surface in a manifold problem. 
    bool _use_reference_outward_normal_direction_provided_to_orient_surface;

    //! whether entire calculation was computed at once
    //! the alternative is that it was broken up over many smaller groups, in which case
    //! this is false, and so the _RHS matrix can not be stored or requested
    bool _entire_batch_computed_at_once;

    //! whether polynomial coefficients were requested to be stored (in a state not yet applied to data)
    bool _store_PTWP_inv_PTW;

    //! initial index for current batch
    int _initial_index_for_batch;

    //! maximum number of neighbors over all target sites
    int _max_num_neighbors;

    //! maximum number of evaluation sites for each target (includes target site)
    int _max_evaluation_sites_per_target;

    //! vector of user requested target operations
    std::vector<TargetOperation> _lro; 

    //! vector containing a mapping from a target functionals enum value to the its place in the list
    //! of target functionals to be applied
    std::vector<int> _lro_lookup; 

    //! index for where this operation begins the for _alpha coefficients (device)
    Kokkos::View<int*> _lro_total_offsets; 

    //! index for where this operation begins the for _alpha coefficients (host)
    Kokkos::View<int*>::HostMirror _host_lro_total_offsets; 

    //! dimensions ^ rank of tensor of output for each target functional (device)
    Kokkos::View<int*> _lro_output_tile_size; 

    //! dimensions ^ rank of tensor of output for each target functional (host)
    Kokkos::View<int*>::HostMirror _host_lro_output_tile_size; 

    //! dimensions ^ rank of tensor of output for each sampling functional (device)
    Kokkos::View<int*> _lro_input_tile_size; 

    //! dimensions ^ rank of tensor of output for each sampling functional (host)
    Kokkos::View<int*>::HostMirror _host_lro_input_tile_size; 

    //! tensor rank of target functional (device)
    Kokkos::View<int*> _lro_output_tensor_rank;

    //! tensor rank of target functional (host)
    Kokkos::View<int*>::HostMirror _host_lro_output_tensor_rank;

    //! tensor rank of sampling functional (device)
    Kokkos::View<int*> _lro_input_tensor_rank;

    //! tensor rank of sampling functional (host)
    Kokkos::View<int*>::HostMirror _host_lro_input_tensor_rank;

    //! used for sizing P_target_row and the _alphas view
    int _total_alpha_values;

    //! additional alpha coefficients due to constraints
    int _added_alpha_size;

    //! determines scratch level spaces and is used to call kernels
    ParallelManager _pm;

    //! order of exact polynomial integration for quadrature rule
    int _order_of_quadrature_points;

    //! dimension of quadrature rule
    int _dimension_of_quadrature_points;

    //! quadrature rule type 
    std::string _quadrature_type;

    //! manages and calculates quadrature
    Quadrature _qm;


/** @name Private Modifiers
 *  Private function because information lives on the device
 */
///@{

    /*! \brief Evaluates the polynomial basis under a particular sampling function. Generally used to fill a row of P.
        \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param target_index         [in] - target number
        \param neighbor_index       [in] - index of neighbor for this target with respect to local numbering [0,...,number of neighbors for target]
        \param alpha                [in] - double to determine convex combination of target and neighbor site at which to evaluate polynomials. (1-alpha)*neighbor + alpha*target
        \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
        \param poly_order           [in] - polynomial basis degree
        \param specific_order_only  [in] - boolean for only evaluating one degree of polynomial when true
        \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
        \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
        \param sampling_strategy    [in] - sampling functional specification
        \param additional_evaluation_local_index [in] - local index for evaluation sites 
    */
    KOKKOS_INLINE_FUNCTION
    void calcPij(const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int dimension, const int poly_order, bool specific_order_only = false, const scratch_matrix_right_type* V = NULL, const ReconstructionSpace reconstruction_space = ReconstructionSpace::ScalarTaylorPolynomial, const SamplingFunctional sampling_strategy = PointSample, const int additional_evaluation_local_index = 0) const;

    /*! \brief Evaluates the gradient of a polynomial basis under the Dirac Delta (pointwise) sampling function.
        \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param target_index         [in] - target number
        \param neighbor_index       [in] - index of neighbor for this target with respect to local numbering [0,...,number of neighbors for target]
        \param alpha                [in] - double to determine convex combination of target and neighbor site at which to evaluate polynomials. (1-alpha)*neighbor + alpha*target
        \param partial_direction    [in] - direction that partial is taken with respect to, e.g. 0 is x direction, 1 is y direction
        \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
        \param poly_order           [in] - polynomial basis degree
        \param specific_order_only  [in] - boolean for only evaluating one degree of polynomial when true
        \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
        \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
        \param sampling_strategy    [in] - sampling functional specification
        \param additional_evaluation_local_index [in] - local index for evaluation sites 
    */
    KOKKOS_INLINE_FUNCTION
    void calcGradientPij(const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int partial_direction, const int dimension, const int poly_order, bool specific_order_only, const scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional sampling_strategy, const int additional_evaluation_local_index = 0) const;

    /*! \brief Evaluates the Hessian of a polynomial basis under the Dirac Delta (pointwise) sampling function.
        \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param target_index         [in] - target number
        \param neighbor_index       [in] - index of neighbor for this target with respect to local numbering [0,...,number of neighbors for target]
        \param alpha                [in] - double to determine convex combination of target and neighbor site at which to evaluate polynomials. (1-alpha)*neighbor + alpha*target
        \param partial_direction_1  [in] - first direction that partial is taken with respect to, e.g. 0 is x direction, 1 is y direction
        \param partial_direction_2  [in] - second direction that partial is taken with respect to, e.g. 0 is x direction, 1 is y direction
        \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
        \param poly_order           [in] - polynomial basis degree
        \param specific_order_only  [in] - boolean for only evaluating one degree of polynomial when true
        \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
        \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
        \param sampling_strategy    [in] - sampling functional specification
        \param additional_evaluation_local_index [in] - local index for evaluation sites 
    */
    KOKKOS_INLINE_FUNCTION
    void calcHessianPij(const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int partial_direction_1, const int partial_direction_2, const int dimension, const int poly_order, bool specific_order_only, const scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional sampling_strategy, const int additional_evaluation_local_index = 0) const;

    /*! \brief Fills the _P matrix with either P or P*sqrt(w)
        \param teamMember           [in] - Kokkos::TeamPolicy member type (created by parallel_for)
        \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param P                   [out] - 2D Kokkos View which will contain evaluation of sampling functional on polynomial basis for each neighbor the target has (stored column major)
        \param w                   [out] - 1D Kokkos View which will contain weighting kernel values for the target with each neighbor if weight_p = true
        \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
        \param polynomial_order     [in] - polynomial basis degree
        \param weight_p             [in] - boolean whether to fill w with kernel weights
        \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
        \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
        \param sampling_strategy    [in] - sampling functional specification
    */
    KOKKOS_INLINE_FUNCTION
    void createWeightsAndP(const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P, scratch_vector_type w, const int dimension, int polynomial_order, bool weight_p = false, scratch_matrix_right_type* V = NULL, const ReconstructionSpace reconstruction_space = ReconstructionSpace::ScalarTaylorPolynomial, const SamplingFunctional sampling_strategy = PointSample) const;

    /*! \brief Fills the _P matrix with P*sqrt(w) for use in solving for curvature

         Uses _curvature_poly_order as the polynomial order of the basis

        \param teamMember           [in] - Kokkos::TeamPolicy member type (created by parallel_for)
        \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param P                   [out] - 2D Kokkos View which will contain evaluation of sampling functional on polynomial basis for each neighbor the target has (stored column major)
        \param w                   [out] - 1D Kokkos View which will contain weighting kernel values for the target with each neighbor if weight_p = true
        \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
        \param only_specific_order  [in] - boolean for only evaluating one degree of polynomial when true
        \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
    */
    KOKKOS_INLINE_FUNCTION
    void createWeightsAndPForCurvature(const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P, scratch_vector_type w, const int dimension, bool only_specific_order, scratch_matrix_right_type* V = NULL) const;

    /*! \brief Evaluates a polynomial basis with a target functional applied to each member of the basis
        \param teamMember                   [in] - Kokkos::TeamPolicy member type (created by parallel_for)
        \param delta                    [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace         [in/out] - scratch space that is allocated so that each team has its own copy. Must be at least as large is the _poly_order*_global_dimensions.
        \param P_target_row                [out] - 1D Kokkos View where the evaluation of the polynomial basis is stored
    */
    KOKKOS_INLINE_FUNCTION
    void computeTargetFunctionals(const member_type& teamMember, scratch_vector_type t1, scratch_vector_type t2, scratch_matrix_right_type P_target_row) const;

    /*! \brief Evaluates a polynomial basis for the curvature with a gradient target functional applied

        _operations is used by this function which is set through a modifier function

        \param teamMember                   [in] - Kokkos::TeamPolicy member type (created by parallel_for)
        \param delta                    [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace         [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _curvature_poly_order*the spatial dimension of the polynomial basis.
        \param P_target_row                [out] - 1D Kokkos View where the evaluation of the polynomial basis is stored
        \param V                            [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
    */
    KOKKOS_INLINE_FUNCTION
    void computeCurvatureFunctionals(const member_type& teamMember, scratch_vector_type t1, scratch_vector_type t2, scratch_matrix_right_type P_target_row, const scratch_matrix_right_type* V, const local_index_type local_neighbor_index = -1) const;

    /*! \brief Evaluates a polynomial basis with a target functional applied, using information from the manifold curvature

         _operations is used by this function which is set through a modifier function

        \param teamMember                   [in] - Kokkos::TeamPolicy member type (created by parallel_for)
        \param delta                    [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
        \param thread_workspace         [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _curvature_poly_order*the spatial dimension of the polynomial basis.
        \param P_target_row                [out] - 1D Kokkos View where the evaluation of the polynomial basis is stored
        \param V                            [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
        \param G_inv                        [in] - (_dimensions-1)*(_dimensions-1) Kokkos View containing inverse of metric tensor
        \param curvature_coefficients       [in] - polynomial coefficients for curvature
        \param curvature_gradients          [in] - approximation of gradient of curvature, Kokkos View of size (_dimensions-1)
    */
    KOKKOS_INLINE_FUNCTION
    void computeTargetFunctionalsOnManifold(const member_type& teamMember, scratch_vector_type t1, scratch_vector_type t2, scratch_matrix_right_type P_target_row, scratch_matrix_right_type V, scratch_matrix_right_type G_inv, scratch_vector_type curvature_coefficients, scratch_vector_type curvature_gradients) const;

    //! Helper function for applying the evaluations from a target functional to the polynomial coefficients
    KOKKOS_INLINE_FUNCTION
    void applyTargetsToCoefficients(const member_type& teamMember, scratch_vector_type t1, scratch_vector_type t2, scratch_matrix_right_type Q, scratch_vector_type w, scratch_matrix_right_type P_target_row, const int target_NP) const;

///@}



/** @name Private Accessors
 *  Private function because information lives on the device
 */
///@{

    //! Returns number of neighbors for a particular target
    KOKKOS_INLINE_FUNCTION
    int getNNeighbors(const int target_index) const {
        return _neighbor_lists.getNumberOfNeighborsDevice(target_index);
    }

    //! Mapping from [0,number of neighbors for a target] to the row that contains the source coordinates for
    //! that neighbor
    KOKKOS_INLINE_FUNCTION
    int getNeighborIndex(const int target_index, const int neighbor_list_num) const {
        return _neighbor_lists.getNeighborDevice(target_index, neighbor_list_num);
    }

    //! Returns the maximum neighbor lists size over all target sites
    KOKKOS_INLINE_FUNCTION
    int getMaxNNeighbors() const {
        return _neighbor_lists.getMaxNumNeighbors();
    }

    //! Returns the maximum number of evaluation sites over all target sites (target sites are included in total)
    KOKKOS_INLINE_FUNCTION
    int getMaxEvaluationSitesPerTarget() const {
        return _max_evaluation_sites_per_target;
    }

    //! (OPTIONAL)
    //! Returns number of additional evaluation sites for a particular target
    KOKKOS_INLINE_FUNCTION
    int getNEvaluationSitesPerTarget(const int target_index) const {
        return _number_of_additional_evaluation_indices(target_index)+1;
    }

    //! (OPTIONAL)
    //! Mapping from [0,number of additional evaluation sites for a target] to the row that contains the coordinates for
    //! that evaluation
    KOKKOS_INLINE_FUNCTION
    int getAdditionalEvaluationIndex(const int target_index, const int additional_list_num) const {
        compadre_kernel_assert_debug((additional_list_num >= 1) 
            && "additional_list_num must be greater than or equal to 1, unlike neighbor lists which begin indexing at 0.");
        return _additional_evaluation_indices(target_index, additional_list_num);
    }

    //! Returns one component of the target coordinate for a particular target. Whether global or local coordinates 
    //! depends upon V being specified
    KOKKOS_INLINE_FUNCTION
    double getTargetCoordinate(const int target_index, const int dim, const scratch_matrix_right_type* V = NULL) const {
        compadre_kernel_assert_debug((_target_coordinates.extent(0) >= (size_t)target_index) && "Target index is out of range for _target_coordinates.");
        if (V==NULL) {
            return _target_coordinates(target_index, dim);
        } else {
            XYZ target_coord = XYZ( _target_coordinates(target_index, 0), 
                                    _target_coordinates(target_index, 1), 
                                    _target_coordinates(target_index, 2));
            return this->convertGlobalToLocalCoordinate(target_coord, dim, V);
        }
    }

    //! (OPTIONAL)
    //! Returns one component of the additional evaluation coordinates. Whether global or local coordinates 
    //! depends upon V being specified
    KOKKOS_INLINE_FUNCTION
    double getTargetAuxiliaryCoordinate(const int target_index, const int additional_list_num, const int dim, const scratch_matrix_right_type* V = NULL) const {
        auto additional_evaluation_index = getAdditionalEvaluationIndex(target_index, additional_list_num);
        compadre_kernel_assert_debug((_additional_evaluation_coordinates.extent(0) >= (size_t)additional_evaluation_index) && "Additional evaluation index is out of range for _additional_evaluation_coordinates.");
        if (V==NULL) {
            return _additional_evaluation_coordinates(additional_evaluation_index, dim);
        } else {
            XYZ additional_target_coord = XYZ( _additional_evaluation_coordinates(additional_evaluation_index, 0),
                                               _additional_evaluation_coordinates(additional_evaluation_index, 1),
                                               _additional_evaluation_coordinates(additional_evaluation_index, 2));
            return this->convertGlobalToLocalCoordinate(additional_target_coord, dim, V);
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
            XYZ neighbor_coord = XYZ(_source_coordinates(this->getNeighborIndex(target_index, neighbor_list_num), 0), _source_coordinates(this->getNeighborIndex(target_index, neighbor_list_num), 1), _source_coordinates(this->getNeighborIndex(target_index, neighbor_list_num), 2));
            return this->convertGlobalToLocalCoordinate(neighbor_coord, dim, V);
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

    //! Returns a component of the local coordinate after transformation from global to local under the orthonormal basis V.
    KOKKOS_INLINE_FUNCTION
    double convertGlobalToLocalCoordinate(const XYZ global_coord, const int dim, const scratch_matrix_right_type* V) const {
        // only written for 2d manifold in 3d space
        double val = 0;
        val += global_coord.x * (*V)(dim, 0);
        val += global_coord.y * (*V)(dim, 1); // can't be called from dimension 1 problem
        if (_dimensions>2) val += global_coord.z * (*V)(dim, 2);
        return val;
    }

    //! Returns a component of the global coordinate after transformation from local to global under the orthonormal basis V^T.
    KOKKOS_INLINE_FUNCTION
    double convertLocalToGlobalCoordinate(const XYZ local_coord, const int dim, const scratch_matrix_right_type* V) const {
        // only written for 2d manifold in 3d space
        double val;
        if (dim == 0 && _dimensions==2) { // 2D problem with 1D manifold
            val = local_coord.x * (*V)(0, dim);
        } else if (dim == 0) { // 3D problem with 2D manifold
            val = local_coord.x * ((*V)(0, dim) + (*V)(1, dim));
        } else if (dim == 1) { // 3D problem with 2D manifold
            val = local_coord.y * ((*V)(0, dim) + (*V)(1, dim));
        }
        return val;
    }

    //! Handles offset from operation input/output + extra evaluation sites
    int getTargetOffsetIndexHost(const int lro_num, const int input_component, const int output_component, const int additional_evaluation_local_index = 0) const {
        return ( _total_alpha_values*additional_evaluation_local_index 
                + _host_lro_total_offsets[lro_num] 
                + input_component*_host_lro_output_tile_size[lro_num] 
                + output_component );
    }

    //! Handles offset from operation input/output + extra evaluation sites
    KOKKOS_INLINE_FUNCTION
    int getTargetOffsetIndexDevice(const int lro_num, const int input_component, const int output_component, const int additional_evaluation_local_index = 0) const {
        return ( _total_alpha_values*additional_evaluation_local_index 
                + _lro_total_offsets[lro_num] 
                + input_component*_lro_output_tile_size[lro_num] 
                + output_component );
    }

///@}

/** @name Private Utility
 *  
 */
///@{

    //! Parses a string to determine solver type
    static DenseSolverType parseSolverType(const std::string& dense_solver_type) {
        std::string solver_type_to_lower = dense_solver_type;
        transform(solver_type_to_lower.begin(), solver_type_to_lower.end(), solver_type_to_lower.begin(), ::tolower);
        if (solver_type_to_lower == "svd") {
            return DenseSolverType::SVD;
        } else if (solver_type_to_lower == "lu") {
            return DenseSolverType::LU;
        } else {
            return DenseSolverType::QR;
        }
    }

    //! Parses a string to determine problem type
    static ProblemType parseProblemType(const std::string& problem_type) {
        std::string problem_type_to_lower = problem_type;
        transform(problem_type_to_lower.begin(), problem_type_to_lower.end(), problem_type_to_lower.begin(), ::tolower);
        if (problem_type_to_lower == "standard") {
            return ProblemType::STANDARD;
        } else if (problem_type_to_lower == "manifold") {
            return ProblemType::MANIFOLD;
        } else {
            return ProblemType::STANDARD;
        }
    }

    //! Parses a string to determine constraint type
    static ConstraintType parseConstraintType(const std::string& constraint_type) {
        std::string constraint_type_to_lower = constraint_type;
        transform(constraint_type_to_lower.begin(), constraint_type_to_lower.end(), constraint_type_to_lower.begin(), ::tolower);
        if (constraint_type_to_lower == "none") {
            return ConstraintType::NO_CONSTRAINT;
        } else if (constraint_type_to_lower == "neumann_grad_scalar") {
            return ConstraintType::NEUMANN_GRAD_SCALAR;
        } else {
            return ConstraintType::NO_CONSTRAINT;
        }
    }

///@}

public:

/** @name Instantiation / Destruction
 *  
 */
///@{

    //! Minimal constructor providing no data (neighbor lists, source sites, target sites) 
    GMLS(ReconstructionSpace reconstruction_space,
        const SamplingFunctional polynomial_sampling_strategy,
        const SamplingFunctional data_sampling_strategy,
        const int poly_order,
        const int dimensions = 3,
        const std::string dense_solver_type = std::string("QR"),
        const std::string problem_type = std::string("STANDARD"),
        const std::string constraint_type = std::string("NO_CONSTRAINT"),
        const int manifold_curvature_poly_order = 2) : 
            _poly_order(poly_order),
            _curvature_poly_order(manifold_curvature_poly_order),
            _dimensions(dimensions),
            _reconstruction_space(reconstruction_space),
            _dense_solver_type(parseSolverType(dense_solver_type)),
            _problem_type(parseProblemType(problem_type)),
            _constraint_type(parseConstraintType(constraint_type)),
            _polynomial_sampling_functional(((_problem_type == ProblemType::MANIFOLD) 
                        && (polynomial_sampling_strategy == VectorPointSample)) ? ManifoldVectorPointSample : polynomial_sampling_strategy),
            _data_sampling_functional(((_problem_type == ProblemType::MANIFOLD) 
                        && (data_sampling_strategy == VectorPointSample)) ? ManifoldVectorPointSample : data_sampling_strategy)
            {

        // seed random number generator pool
        _random_number_pool = pool_type(1);

        _NP = this->getNP(_poly_order, dimensions, _reconstruction_space);
        Kokkos::fence();

        // register curvature operations for manifold problems
        if (_problem_type == ProblemType::MANIFOLD) {
            _curvature_support_operations = Kokkos::View<TargetOperation*>
                ("operations needed for manifold gradient reconstruction", 1);
            auto curvature_support_operations_mirror = 
                Kokkos::create_mirror_view(_curvature_support_operations);
            curvature_support_operations_mirror(0) = 
                TargetOperation::GradientOfScalarPointEvaluation;
            Kokkos::deep_copy(_curvature_support_operations, curvature_support_operations_mirror);
        }

        _lro_lookup = std::vector<int>(TargetOperation::COUNT,-1); // hard coded against number of operations defined
        _lro = std::vector<TargetOperation>();

        // various initializations
        _total_alpha_values = 0;

        _weighting_type = WeightingFunctionType::Power;
        _curvature_weighting_type = WeightingFunctionType::Power;
        _weighting_power = 2;
        _curvature_weighting_power = 2;

        _reconstruction_space_rank = ActualReconstructionSpaceRank[_reconstruction_space];

        _basis_multiplier = 1;
        _sampling_multiplier = 1;

        _nontrivial_nullspace = false;
        _orthonormal_tangent_space_provided = false; 
        _reference_outward_normal_direction_provided = false;
        _use_reference_outward_normal_direction_provided_to_orient_surface = false;
        _entire_batch_computed_at_once = true;
        _store_PTWP_inv_PTW = false;

        _initial_index_for_batch = 0;

        _max_num_neighbors = 0;
        _max_evaluation_sites_per_target = 1;

        _global_dimensions = dimensions;
        if (_problem_type == ProblemType::MANIFOLD) {
            _local_dimensions = dimensions-1;
        } else {
            _local_dimensions = dimensions;
        }

        _order_of_quadrature_points = 0;
        _dimension_of_quadrature_points = 0;
    }

    //! Constructor for the case when the data sampling functional does not match the polynomial
    //! sampling functional. Only case anticipated is staggered Laplacian.
    GMLS(const int poly_order,
         const int dimensions = 3,
         const std::string dense_solver_type = std::string("QR"),
         const std::string problem_type = std::string("STANDARD"),
         const std::string constraint_type = std::string("NO_CONSTRAINT"),
         const int manifold_curvature_poly_order = 2)
      : GMLS(ReconstructionSpace::VectorOfScalarClonesTaylorPolynomial, VectorPointSample, VectorPointSample, poly_order, dimensions, dense_solver_type, problem_type, constraint_type, manifold_curvature_poly_order) {}

    //! Constructor for the case when nonstandard sampling functionals or reconstruction spaces
    //! are to be used. Reconstruction space and sampling strategy can only be set at instantiation.
    GMLS(ReconstructionSpace reconstruction_space,
         SamplingFunctional dual_sampling_strategy,
         const int poly_order,
         const int dimensions = 3,
         const std::string dense_solver_type = std::string("QR"),
         const std::string problem_type = std::string("STANDARD"),
         const std::string constraint_type = std::string("NO_CONSTRAINT"),
         const int manifold_curvature_poly_order = 2)
      : GMLS(reconstruction_space, dual_sampling_strategy, dual_sampling_strategy, poly_order, dimensions, dense_solver_type, problem_type, constraint_type, manifold_curvature_poly_order) {}

    //! Destructor
    ~GMLS(){
    };

///@}

/** @name Functors
 *  Member functions that perform operations on the entire batch
 */
///@{

    //! Tag for functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity
    struct AssembleStandardPsqrtW{};

    //! Tag for functor to evaluate targets, apply target evaluation to polynomial coefficients to
    //! store in _alphas
    struct ApplyStandardTargets{};

    //! Tag for functor to create a coarse tangent approximation from a given neighborhood of points
    struct ComputeCoarseTangentPlane{};

    //! Tag for functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity for curvature
    struct AssembleCurvaturePsqrtW{};

    //! Tag for functor to evaluate curvature targets and construct accurate tangent direction approximation for manifolds
    struct GetAccurateTangentDirections{};

    //! Tag for functor to determine if tangent directions need reordered, and to reorder them if needed
    struct FixTangentDirectionOrdering{};

    //! Tag for functor to evaluate curvature targets and apply to coefficients of curvature reconstruction
    struct ApplyCurvatureTargets{};

    //! Tag for functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity
    struct AssembleManifoldPsqrtW{};

    //! Tag for functor to evaluate targets, apply target evaluation to polynomial coefficients to store in _alphas
    struct ApplyManifoldTargets{};

    //! Tag for functor to calculate prestencil weights to apply to data to transform into a format expected by a GMLS stencil
    struct ComputePrestencilWeights{};



    //! Functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity
    KOKKOS_INLINE_FUNCTION
    void operator() (const AssembleStandardPsqrtW&, const member_type& teamMember) const;

    //! Functor to evaluate targets, apply target evaluation to polynomial coefficients to store in _alphas
    KOKKOS_INLINE_FUNCTION
    void operator() (const ApplyStandardTargets&, const member_type& teamMember) const;

    //! Functor to create a coarse tangent approximation from a given neighborhood of points
    KOKKOS_INLINE_FUNCTION
    void operator() (const ComputeCoarseTangentPlane&, const member_type& teamMember) const;

    //! Functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity for curvature
    KOKKOS_INLINE_FUNCTION
    void operator() (const AssembleCurvaturePsqrtW&, const member_type& teamMember) const;

    //! Functor to evaluate curvature targets and construct accurate tangent direction approximation for manifolds
    KOKKOS_INLINE_FUNCTION
    void operator() (const GetAccurateTangentDirections&, const member_type& teamMember) const;

    //! Functor to determine if tangent directions need reordered, and to reorder them if needed
    //! We require that the normal is consistent with a right hand rule on the tangent vectors
    KOKKOS_INLINE_FUNCTION
    void operator() (const FixTangentDirectionOrdering&, const member_type& teamMember) const;

    //! Functor to evaluate curvature targets and apply to coefficients of curvature reconstruction
    KOKKOS_INLINE_FUNCTION
    void operator() (const ApplyCurvatureTargets&, const member_type& teamMember) const;

    //! Functor to assemble the P*sqrt(weights) matrix and construct sqrt(weights)*Identity
    KOKKOS_INLINE_FUNCTION
    void operator() (const AssembleManifoldPsqrtW&, const member_type& teamMember) const;

    //! Functor to evaluate targets, apply target evaluation to polynomial coefficients to store in _alphas
    KOKKOS_INLINE_FUNCTION
    void operator() (const ApplyManifoldTargets&, const member_type& teamMember) const;

    //! Functor to calculate prestencil weights to apply to data to transform into a format expected by a GMLS stencil
    KOKKOS_INLINE_FUNCTION
    void operator() (const ComputePrestencilWeights&, const member_type& teamMember) const;

///@}

/** @name Public Utility
 *  
 */
///@{

    //! Returns size of the basis for a given polynomial order and dimension
    //! General to dimension 1..3 and polynomial order m
    //! The divfree options will return the divergence-free basis if true
    KOKKOS_INLINE_FUNCTION
    static int getNP(const int m, const int dimension = 3, const ReconstructionSpace r_space = ReconstructionSpace::ScalarTaylorPolynomial) {
        if (r_space != ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial) {
            return ScalarTaylorPolynomialBasis::getSize(m, dimension);
        } else {
            return DivergenceFreePolynomialBasis::getSize(m, dimension);
        }
    }

    //! Returns number of neighbors needed for unisolvency for a given basis order and dimension
    KOKKOS_INLINE_FUNCTION
    static int getNN(const int m, const int dimension = 3, const ReconstructionSpace r_space = ReconstructionSpace::ScalarTaylorPolynomial) {
        // may need div-free argument in the future
        const int np = getNP(m, dimension, r_space);
        int nn = np;
        switch (dimension) {
            case 3:
                nn = np * (1.7 + m*0.1);
                break;
            case 2:
                nn = np * (1.4 + m*0.03);
                break;
            case 1:
                nn = np * 1.1;
        }
        return nn;
    }

    /*! \brief Evaluates the weighting kernel
        \param r                [in] - Euclidean distance of relative vector. Euclidean distance of (target - neighbor) in some basis.
        \param h                [in] - window size. Kernel is guaranteed to take on a value of zero if it exceeds h.
        \param weighting_type   [in] - weighting type to be evaluated as the kernel. e,g. power, Gaussian, etc..
        \param power            [in] - power parameter to be given to the kernel.
    */
    KOKKOS_INLINE_FUNCTION
    static double Wab(const double r, const double h, const WeightingFunctionType& weighting_type, const int power) {
        if (weighting_type == WeightingFunctionType::Power) {
            return std::pow(1.0-std::abs(r/h), power) * double(1.0-std::abs(r/h)>0.0);
        } else if (weighting_type == WeightingFunctionType::Gaussian) {
            // 2.5066282746310002416124 = sqrt(2*pi)
            double h_over_3 = h/3.0;
            return 1./( h_over_3 * 2.5066282746310002416124 ) * std::exp(-.5*r*r/(h_over_3*h_over_3));
        } else if (weighting_type == WeightingFunctionType::CubicSpline) {
            double x = r/h;
            return ((1-x)+x*(1-x)*(1-2*x)) * double(r<=h);
        } else { // no weighting power specified
            compadre_kernel_assert_release(false && "Invalid WeightingFunctionType selected.");
            return 0; 
        }
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


    //! Helper function for finding alpha coefficients
    static int getTargetOutputIndex(const int operation_num, const int output_component_axis_1, const int output_component_axis_2, const int dimensions) {
        const int axis_1_size = (TargetOutputTensorRank[operation_num] > 1) ? dimensions : 1;
        return axis_1_size*output_component_axis_1 + output_component_axis_2; // 0 for scalar, 0 for vector;
    }

    //! Helper function for finding alpha coefficients
    static int getSamplingOutputIndex(const SamplingFunctional sf, const int output_component_axis_1, const int output_component_axis_2) {
        const int axis_1_size = (sf.output_rank > 1) ? sf.output_rank : 1;
        return axis_1_size*output_component_axis_1 + output_component_axis_2; // 0 for scalar, 0 for vector;
    }

    //! Output rank for sampling operation
    static int getOutputRankOfSampling(SamplingFunctional sro) {
        return sro.output_rank;
    }

    //! Input rank for sampling operation
    static int getInputRankOfSampling(SamplingFunctional sro) {
        return sro.input_rank;
    }
    

///@}

/** @name Accessors
 *  Retrieve member variables through public member functions
 */
///@{


    //! Returns (size of the basis used in instance's polynomial reconstruction) x (data input dimension)
    host_managed_local_index_type getPolynomialCoefficientsDomainRangeSize() const { 
        host_managed_local_index_type sizes("sizes", 2);
        sizes(0) = _basis_multiplier*_NP;
        sizes(1) = _sampling_multiplier*getMaxNNeighbors();
        return sizes;
    }

    //! Returns size of the basis used in instance's polynomial reconstruction
    int getPolynomialCoefficientsSize() const {
        auto sizes = this->getPolynomialCoefficientsDomainRangeSize();
        return sizes(0);
    }

    //! Returns 2D array size in memory on which coefficients are stored
    host_managed_local_index_type getPolynomialCoefficientsMemorySize() const {
        auto M_by_N = this->getPolynomialCoefficientsDomainRangeSize();
        compadre_assert_release(_entire_batch_computed_at_once 
                && "Entire batch not computed at once, so getFullPolynomialCoefficientsBasis() can not be called.");
        compadre_assert_release(_store_PTWP_inv_PTW
                && "generateAlphas() called with keep_coefficients set to false.");
        host_managed_local_index_type sizes("sizes", 2);
        if ((_constraint_type == ConstraintType::NO_CONSTRAINT) && (_dense_solver_type != DenseSolverType::LU)) {
            int rhsdim = getRHSSquareDim(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, M_by_N[1], M_by_N[0]);
            sizes(0) = rhsdim;
            sizes(1) = rhsdim;
        } else {
            getPDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, M_by_N[1], M_by_N[0], sizes(1), sizes(0));
        }
        return sizes;
    }

    //! Dimension of the GMLS problem, set only at class instantiation
    int getDimensions() const { return _dimensions; }

    //! Dimension of the GMLS problem's point data (spatial description of points in ambient space), set only at class instantiation
    int getGlobalDimensions() const { return _global_dimensions; }

    //! Local dimension of the GMLS problem (less than global dimension if on a manifold), set only at class instantiation
    int getLocalDimensions() const { return _local_dimensions; }

    //! Get dense solver type
    DenseSolverType getDenseSolverType() { return _dense_solver_type; }

    //! Get problem type
    ProblemType getProblemType() { return _problem_type; }

    //! Get constraint type
    ConstraintType getConstraintType() { return _constraint_type; }

    //! Type for weighting kernel for GMLS problem
    WeightingFunctionType getWeightingType() const { return _weighting_type; }

    //! Type for weighting kernel for curvature 
    WeightingFunctionType getManifoldWeightingType() const { return _curvature_weighting_type; }

    //! Power for weighting kernel for GMLS problem
    int getWeightingPower() const { return _weighting_power; }

    //! Power for weighting kernel for curvature
    int getManifoldWeightingPower() const { return _curvature_weighting_power; }

    //! Number of quadrature points
    int getNumberOfQuadraturePoints() const { return _qm.getNumberOfQuadraturePoints(); }

    //! Order of quadrature points
    int getOrderOfQuadraturePoints() const { return _order_of_quadrature_points; }

    //! Dimensions of quadrature points
    int getDimensionOfQuadraturePoints() const { return _dimension_of_quadrature_points; }

    //! Type of quadrature points
    std::string getQuadratureType() const { return _quadrature_type; }

    //! Get neighbor list accessor
    decltype(_neighbor_lists)* getNeighborLists() { return &_neighbor_lists; }

    //! Get a view (device) of all tangent direction bundles.
    decltype(_T) getTangentDirections() const { return _T; }

    //! Get a view (device) of all reference outward normal directions.
    decltype(_ref_N) getReferenceNormalDirections() const { return _ref_N; }

    //! Get component of tangent or normal directions for manifold problems
    double getTangentBundle(const int target_index, const int direction, const int component) const {
        // Component index 0.._dimensions-2 will return tangent direction
        // Component index _dimensions-1 will return the normal direction
        scratch_matrix_right_type::HostMirror 
                T(_host_T.data() + target_index*_dimensions*_dimensions, _dimensions, _dimensions);
        return T(direction, component);
    }

    //! Get component of tangent or normal directions for manifold problems
    double getReferenceNormalDirection(const int target_index, const int component) const {
        compadre_assert_debug(_reference_outward_normal_direction_provided && 
                "getRefenceNormalDirection called, but reference outwrad normal directions were never provided.");
        scratch_vector_type::HostMirror 
                ref_N(_host_ref_N.data() + target_index*_dimensions, _dimensions);
        return ref_N(component);
    }

    //! Get the local index (internal) to GMLS for a particular TargetOperation
    //! Every TargetOperation has a global index which can be readily found in Compadre::TargetOperation
    //! but this function returns the index used inside of the GMLS class
    int getTargetOperationLocalIndex(TargetOperation lro) const {
        return _lro_lookup[(int)lro];
    }

    //! Get a view (device) of all rank 2 preprocessing tensors
    //! This is a rank 5 tensor that is able to provide data transformation
    //! into a form that GMLS is able to operate on. The ranks are as follows:
    //!
    //! 1 - Either size 2 if it operates on the target site and neighbor site (staggered schemes)
    //!     or 1 if it operates only on the neighbor sites (almost every scheme)
    //!
    //! 2 - If the data transform varies with each target site (but could be the same for each neighbor of that target site), then this is the number of target sites
    //!
    //! 3 - If the data transform varies with each neighbor of each target site, then this is the number of neighbors for each respective target (max number of neighbors for all target sites is its uniform size)
    //!
    //! 4 - Data transform resulting in rank 1 data for the GMLS operator will have size _local_dimensions, otherwise 1
    //!
    //! 5 - Data transform taking in rank 1 data will have size _global_dimensions, otherwise 1
    decltype(_prestencil_weights) getPrestencilWeights() const { 
        return _prestencil_weights;
    }

    //! Retrieves the offset for an operator based on input and output component, generic to row
    //! (but still multiplied by the number of neighbors for each row and then needs a neighbor number added 
    //! to this returned value to be meaningful)
    int getAlphaColumnOffset(TargetOperation lro, const int output_component_axis_1, 
            const int output_component_axis_2, const int input_component_axis_1, 
            const int input_component_axis_2, const int additional_evaluation_local_index = 0) const {

        const int lro_number = _lro_lookup[(int)lro];
        compadre_assert_debug((lro_number >= 0) && "getAlphaColumnOffset called for a TargetOperation that was not registered.");

        // the target functional input indexing is sized based on the output rank of the sampling
        // functional used, which can not be inferred unless a specification of target functional,
        // reconstruction space, and sampling functional are all known (as was the case at the
        // construction of this class)
        const int input_index = getSamplingOutputIndex(_data_sampling_functional, input_component_axis_1, input_component_axis_2);
        const int output_index = getTargetOutputIndex((int)lro, output_component_axis_1, output_component_axis_2, _dimensions);

        return getTargetOffsetIndexHost(lro_number, input_index, output_index, additional_evaluation_local_index);
    }

    //! Get a view (device) of all alphas
    decltype(_alphas) getAlphas() const { return _alphas; }

    //! Get a view (device) of all polynomial coefficients basis
    decltype(_RHS) getFullPolynomialCoefficientsBasis() const { 
        compadre_assert_release(_entire_batch_computed_at_once 
                && "Entire batch not computed at once, so getFullPolynomialCoefficientsBasis() can not be called.");
        compadre_assert_release(_store_PTWP_inv_PTW
                && "generateAlphas() called with keep_coefficients set to false.");
        if ((_constraint_type == ConstraintType::NO_CONSTRAINT) && (_dense_solver_type != DenseSolverType::LU)) {
            return _RHS; 
        } else {
            return _P; 
        }
    }

    //! Get the polynomial sampling functional specified at instantiation
    SamplingFunctional getPolynomialSamplingFunctional() const { return _polynomial_sampling_functional; }
 
    //! Get the data sampling functional specified at instantiation (often the same as the polynomial sampling functional)
    SamplingFunctional getDataSamplingFunctional() const { return _data_sampling_functional; }

    //! Get the reconstruction space specified at instantiation
    ReconstructionSpace getReconstructionSpace() const { return _reconstruction_space; }

    //! Helper function for getting alphas for scalar reconstruction from scalar data
    double getAlpha0TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int additional_evaluation_site = 0) const {
        // e.g. Dirac Delta target of a scalar field
        return getAlpha(lro, target_index, 0, 0, neighbor_index, 0, 0, additional_evaluation_site);
    }

    //! Helper function for getting alphas for vector reconstruction from scalar data
    double getAlpha0TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int additional_evaluation_site = 0) const {
        // e.g. gradient of a scalar field
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, 0, 0, additional_evaluation_site);
    }

    //! Helper function for getting alphas for matrix reconstruction from scalar data
    double getAlpha0TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int additional_evaluation_site = 0) const {
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, 0, 0, additional_evaluation_site);
    }

    //! Helper function for getting alphas for scalar reconstruction from vector data
    double getAlpha1TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int input_component, const int additional_evaluation_site = 0) const {
        // e.g. divergence of a vector field
        return getAlpha(lro, target_index, 0, 0, neighbor_index, input_component, 0, additional_evaluation_site);
    }

    //! Helper function for getting alphas for vector reconstruction from vector data
    double getAlpha1TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int input_component, const int additional_evaluation_site = 0) const {
        // e.g. curl of a vector field
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, input_component, 0, additional_evaluation_site);
    }

    //! Helper function for getting alphas for matrix reconstruction from vector data
    double getAlpha1TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component, const int additional_evaluation_site = 0) const {
        // e.g. gradient of a vector field
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, input_component, 0, additional_evaluation_site);
    }

    //! Helper function for getting alphas for scalar reconstruction from matrix data
    double getAlpha2TensorTo0Tensor(TargetOperation lro, const int target_index, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int additional_evaluation_site = 0) const {
        return getAlpha(lro, target_index, 0, 0, neighbor_index, input_component_axis_1, input_component_axis_2, additional_evaluation_site);
    }

    //! Helper function for getting alphas for vector reconstruction from matrix data
    double getAlpha2TensorTo1Tensor(TargetOperation lro, const int target_index, const int output_component, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int additional_evaluation_site = 0) const {
        return getAlpha(lro, target_index, output_component, 0, neighbor_index, input_component_axis_1, input_component_axis_2, additional_evaluation_site);
    }

    //! Helper function for getting alphas for matrix reconstruction from matrix data
    double getAlpha2TensorTo2Tensor(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int additional_evaluation_site = 0) const {
        return getAlpha(lro, target_index, output_component_axis_1, output_component_axis_2, neighbor_index, input_component_axis_1, input_component_axis_2, additional_evaluation_site);
    }

    //! Gives index into alphas given two axes, which when incremented by the neighbor number transforms access into
    //! alphas from a rank 1 view into a rank 3 view.
    KOKKOS_INLINE_FUNCTION
    global_index_type getAlphaIndexDevice(const int target_index, const int alpha_column_offset) const {

        global_index_type total_neighbors_before_target = _neighbor_lists.getRowOffsetDevice(target_index);
        int total_added_alphas_before_target = target_index*_added_alpha_size;

        int alphas_per_tile_per_target = _neighbor_lists.getNumberOfNeighborsDevice(target_index) + _added_alpha_size;

        return (total_neighbors_before_target+TO_GLOBAL(total_added_alphas_before_target))
                 *TO_GLOBAL(_total_alpha_values)*TO_GLOBAL(_max_evaluation_sites_per_target)
                   + TO_GLOBAL(alpha_column_offset*alphas_per_tile_per_target);

    }

    //! Gives index into alphas given two axes, which when incremented by the neighbor number transforms access into
    //! alphas from a rank 1 view into a rank 3 view.
    global_index_type getAlphaIndexHost(const int target_index, const int alpha_column_offset) const {

        global_index_type total_neighbors_before_target = _neighbor_lists.getRowOffsetHost(target_index);
        int total_added_alphas_before_target = target_index*_added_alpha_size;

        int alphas_per_tile_per_target = _neighbor_lists.getNumberOfNeighborsHost(target_index) + _added_alpha_size;

        return (total_neighbors_before_target+TO_GLOBAL(total_added_alphas_before_target))
                 *TO_GLOBAL(_total_alpha_values)*TO_GLOBAL(_max_evaluation_sites_per_target)
                   + TO_GLOBAL(alpha_column_offset*alphas_per_tile_per_target);

    }

    //! Underlying function all interface helper functions call to retrieve alpha values
    double getAlpha(TargetOperation lro, const int target_index, const int output_component_axis_1, const int output_component_axis_2, const int neighbor_index, const int input_component_axis_1, const int input_component_axis_2, const int additional_evaluation_site = 0) const {
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

        const int alpha_column_offset = this->getAlphaColumnOffset( lro, output_component_axis_1, 
                output_component_axis_2, input_component_axis_1, input_component_axis_2, additional_evaluation_site);

        auto alphas_index = this->getAlphaIndexHost(target_index, alpha_column_offset);
        return _host_alphas(alphas_index + neighbor_index);
    }

    //! Returns a stencil to transform data from its existing state into the input expected 
    //! for some sampling functionals.
    double getPreStencilWeight(SamplingFunctional sro, const int target_index, const int neighbor_index, bool for_target, const int output_component = 0, const int input_component = 0) const {
        // for certain sampling strategies, linear combinations of the neighbor and target value are needed
        // for the traditional PointSample, this value is 1 for the neighbor and 0 for the target
        if (sro == PointSample ) {
            if (for_target) return 0; else return 1;
        }

        // these check conditions on the sampling operator and change indexing on target and neighbors
        // in order to reuse information, such as if the same data transformation is used, regardless
        // of target site or neighbor site
        const int target_index_in_weights = 
            (sro.transform_type==DifferentEachTarget 
                    || sro.transform_type==DifferentEachNeighbor) ?
                target_index : 0;
        const int neighbor_index_in_weights = 
            (sro.transform_type==DifferentEachNeighbor) ?
                neighbor_index : 0;

        return _host_prestencil_weights((int)for_target, target_index_in_weights, neighbor_index_in_weights, 
                    output_component, input_component);
    }

    //! Dimensions ^ output rank for target operation
    int getOutputDimensionOfOperation(TargetOperation lro, bool ambient = false) const {
        int return_val;
        if (ambient) return_val = std::pow(_global_dimensions, TargetOutputTensorRank[(int)lro]);
        else return_val = std::pow(_local_dimensions, TargetOutputTensorRank[(int)lro]);
        return return_val;
    }

    //! Dimensions ^ input rank for target operation (always in local chart if on a manifold, never ambient space)
    int getInputDimensionOfOperation(TargetOperation lro) const {
        return this->_host_lro_input_tile_size[_lro_lookup[(int)lro]];
        // this is the same return values as the OutputDimensionOfSampling for the GMLS class's SamplingFunctional
    }

    //! Dimensions ^ output rank for sampling operation 
    //! (always in local chart if on a manifold, never ambient space)
    int getOutputDimensionOfSampling(SamplingFunctional sro) const {
        return std::pow(_local_dimensions, sro.output_rank);
    }

    //! Dimensions ^ output rank for sampling operation 
    //! (always in ambient space, never local chart on a manifold)
    int getInputDimensionOfSampling(SamplingFunctional sro) const {
        return std::pow(_global_dimensions, sro.input_rank);
    }

    //! Calculate basis_multiplier
    int calculateBasisMultiplier(const ReconstructionSpace rs) const {
        // calculate the dimension of the basis 
        // (a vector space on a manifold requires two components, for example)
        return std::pow(_local_dimensions, ActualReconstructionSpaceRank[(int)rs]);
    }

    //! Calculate sampling_multiplier
    int calculateSamplingMultiplier(const ReconstructionSpace rs, const SamplingFunctional sro) const {
        // this would normally be SamplingOutputTensorRank[_data_sampling_functional], but we also want to handle the
        // case of reconstructions where a scalar basis is reused as a vector, and this handles everything
        // this handles scalars, vectors, and scalars that are reused as vectors
        int bm = this->calculateBasisMultiplier(rs);
        int sm = this->getOutputDimensionOfSampling(sro);
        if (rs == ReconstructionSpace::VectorOfScalarClonesTaylorPolynomial) {
            // storage and computational efficiency by reusing solution to scalar problem for 
            // a vector problem (in 3d, 27x cheaper computation, 9x cheaper storage)
            sm = std::min(bm,sm);
        }
        return sm;
    }


///@}


/** @name Modifiers
 *  Changed member variables through public member functions
 */
///@{

    void resetCoefficientData() {
        if (_RHS.extent(0) > 0)
            _RHS = Kokkos::View<double*>("RHS",0);
    }

    //! Sets basic problem data (neighbor lists, source coordinates, and target coordinates)
    template<typename view_type_1, typename view_type_2, typename view_type_3, typename view_type_4>
    void setProblemData(
            view_type_1 neighbor_lists,
            view_type_2 source_coordinates,
            view_type_3 target_coordinates,
            view_type_4 epsilons) {
        this->setNeighborLists<view_type_1>(neighbor_lists);
        this->setSourceSites<view_type_2>(source_coordinates);
        this->setTargetSites<view_type_3>(target_coordinates);
        this->setWindowSizes<view_type_4>(epsilons);
    }

    //! Sets basic problem data (neighbor lists data, number of neighbors list, source coordinates, and target coordinates)
    template<typename view_type_1, typename view_type_2, typename view_type_3, typename view_type_4>
    void setProblemData(
            view_type_1 cr_neighbor_lists,
            view_type_1 number_of_neighbors_list,
            view_type_2 source_coordinates,
            view_type_3 target_coordinates,
            view_type_4 epsilons) {
        this->setNeighborLists<view_type_1>(cr_neighbor_lists, number_of_neighbors_list);
        this->setSourceSites<view_type_2>(source_coordinates);
        this->setTargetSites<view_type_3>(target_coordinates);
        this->setWindowSizes<view_type_4>(epsilons);
    }

    //! (OPTIONAL) Sets additional evaluation sites for each target site
    template<typename view_type_1, typename view_type_2>
    void setAdditionalEvaluationSitesData(
            view_type_1 additional_evaluation_indices,
            view_type_2 additional_evaluation_coordinates) {
        this->setAuxiliaryEvaluationIndicesLists<view_type_1>(additional_evaluation_indices);
        this->setAuxiliaryEvaluationCoordinates<view_type_2>(additional_evaluation_coordinates);
    }

    //! Sets neighbor list information from compressed row neighborhood lists data (if same view_type).
    template <typename view_type>
    typename std::enable_if<view_type::rank==1&&std::is_same<decltype(_neighbor_lists)::internal_view_type,view_type>::value==1, void>::type 
            setNeighborLists(view_type neighbor_lists, view_type number_of_neighbors_list) {

        _neighbor_lists = NeighborLists<view_type>(neighbor_lists, number_of_neighbors_list);
        _max_num_neighbors = _neighbor_lists.getMaxNumNeighbors();
        _host_number_of_neighbors_list = decltype(_host_number_of_neighbors_list)("host number of neighbors list", _neighbor_lists.getNumberOfTargets());
        Kokkos::parallel_for("copy neighbor list sizes", Kokkos::RangePolicy<host_execution_space>(0, _host_number_of_neighbors_list.extent(0)), KOKKOS_LAMBDA(const int i) {
            _host_number_of_neighbors_list(i) = _neighbor_lists.getNumberOfNeighborsHost(i);
        });
        Kokkos::fence();
        this->resetCoefficientData();

    }

    //! Sets neighbor list information from compressed row neighborhood lists data (if different view_type).
    template <typename view_type>
    typename std::enable_if<view_type::rank==1&&std::is_same<decltype(_neighbor_lists)::internal_view_type,view_type>::value==0, void>::type 
            setNeighborLists(view_type neighbor_lists, view_type number_of_neighbors_list) {

        typedef decltype(_neighbor_lists)::internal_view_type gmls_view_type;
        gmls_view_type d_neighbor_lists("compressed row neighbor lists data", neighbor_lists.extent(0));
        gmls_view_type d_number_of_neighbors_list("number of neighbors list", number_of_neighbors_list.extent(0));
        Kokkos::deep_copy(d_neighbor_lists, neighbor_lists);
        Kokkos::deep_copy(d_number_of_neighbors_list, number_of_neighbors_list);
        Kokkos::fence();
        _neighbor_lists = NeighborLists<gmls_view_type>(d_neighbor_lists, d_number_of_neighbors_list);
        _max_num_neighbors = _neighbor_lists.getMaxNumNeighbors();
        _host_number_of_neighbors_list = decltype(_host_number_of_neighbors_list)("host number of neighbors list", _neighbor_lists.getNumberOfTargets());
        Kokkos::parallel_for("copy neighbor list sizes", Kokkos::RangePolicy<host_execution_space>(0, _host_number_of_neighbors_list.extent(0)), KOKKOS_LAMBDA(const int i) {
            _host_number_of_neighbors_list(i) = _neighbor_lists.getNumberOfNeighborsHost(i);
        });
        Kokkos::fence();
        this->resetCoefficientData();
            
    }

    //! Sets neighbor list information. Should be # targets x maximum number of neighbors for any target + 1.
    //! first entry in ever row should be the number of neighbors for the corresponding target.
    template <typename view_type>
    typename std::enable_if<view_type::rank==2, void>::type setNeighborLists(view_type neighbor_lists) {
    
        _neighbor_lists = Convert2DToCompressedRowNeighborLists<decltype(neighbor_lists), Kokkos::View<int*> >(neighbor_lists);
        _max_num_neighbors = _neighbor_lists.getMaxNumNeighbors();
        _host_number_of_neighbors_list = decltype(_host_number_of_neighbors_list)("host number of neighbors list", _neighbor_lists.getNumberOfTargets());
        Kokkos::parallel_for("copy neighbor list sizes", Kokkos::RangePolicy<host_execution_space>(0, _host_number_of_neighbors_list.extent(0)), KOKKOS_LAMBDA(const int i) {
            _host_number_of_neighbors_list(i) = _neighbor_lists.getNumberOfNeighborsHost(i);
        });
        Kokkos::fence();
        this->resetCoefficientData();

    }

    //! Sets source coordinate information. Rows of this 2D-array should correspond to neighbor IDs contained in the entries
    //! of the neighbor lists 2D array.
    template<typename view_type>
    void setSourceSites(view_type source_coordinates) {

        // allocate memory on device
        _source_coordinates = decltype(_source_coordinates)("device neighbor coordinates",
                source_coordinates.extent(0), source_coordinates.extent(1));

        typedef typename view_type::memory_space input_array_memory_space;
        if (std::is_same<input_array_memory_space, device_memory_space>::value) {
            // check if on the device, then copy directly
            // if it is, then it doesn't match the internal layout we use
            // then copy to the host mirror
            // switches potential layout mismatches
            Kokkos::deep_copy(_source_coordinates, source_coordinates);
        } else {
            // if is on the host, copy to the host mirror
            // then copy to the device
            // switches potential layout mismatches
            auto host_source_coordinates = Kokkos::create_mirror_view(_source_coordinates);
            Kokkos::deep_copy(host_source_coordinates, source_coordinates);
            // switches memory spaces
            Kokkos::deep_copy(_source_coordinates, host_source_coordinates);
        }
        this->resetCoefficientData();
    }

    //! Sets source coordinate information. Rows of this 2D-array should correspond to neighbor IDs contained in the entries
    //! of the neighbor lists 2D array.
    template<typename view_type>
    void setSourceSites(decltype(_source_coordinates) source_coordinates) {
        // allocate memory on device
        _source_coordinates = source_coordinates;
        this->resetCoefficientData();
    }

    //! Sets target coordinate information. Rows of this 2D-array should correspond to rows of the neighbor lists.
    template<typename view_type>
    void setTargetSites(view_type target_coordinates) {
        // allocate memory on device
        _target_coordinates = decltype(_target_coordinates)("device target coordinates",
                target_coordinates.extent(0), target_coordinates.extent(1));

        typedef typename view_type::memory_space input_array_memory_space;
        if (std::is_same<input_array_memory_space, device_memory_space>::value) {
            // check if on the device, then copy directly
            // if it is, then it doesn't match the internal layout we use
            // then copy to the host mirror
            // switches potential layout mismatches
            Kokkos::deep_copy(_target_coordinates, target_coordinates);
        } else {
            // if is on the host, copy to the host mirror
            // then copy to the device
            // switches potential layout mismatches
            auto host_target_coordinates = Kokkos::create_mirror_view(_target_coordinates);
            Kokkos::deep_copy(host_target_coordinates, target_coordinates);
            // switches memory spaces
            Kokkos::deep_copy(_target_coordinates, host_target_coordinates);
        }
        _number_of_additional_evaluation_indices 
            = decltype(_number_of_additional_evaluation_indices)("number of additional evaluation indices", target_coordinates.extent(0));
        Kokkos::deep_copy(_number_of_additional_evaluation_indices, 0);
        this->resetCoefficientData();
    }

    //! Sets target coordinate information. Rows of this 2D-array should correspond to rows of the neighbor lists.
    template<typename view_type>
    void setTargetSites(decltype(_target_coordinates) target_coordinates) {
        // allocate memory on device
        _target_coordinates = target_coordinates;
        _number_of_additional_evaluation_indices 
            = decltype(_number_of_additional_evaluation_indices)("number of additional evaluation indices", target_coordinates.extent(0));
        Kokkos::deep_copy(_number_of_additional_evaluation_indices, 0);
        this->resetCoefficientData();
    }

    //! Sets window sizes, also called the support of the kernel
    template<typename view_type>
    void setWindowSizes(view_type epsilons) {

        // allocate memory on device
        _epsilons = decltype(_epsilons)("device epsilons", epsilons.extent(0));

        _host_epsilons = Kokkos::create_mirror_view(_epsilons);
        Kokkos::deep_copy(_host_epsilons, epsilons);
        // copy data from host to device
        Kokkos::deep_copy(_epsilons, _host_epsilons);
        this->resetCoefficientData();
    }

    //! Sets window sizes, also called the support of the kernel (device)
    template<typename view_type>
    void setWindowSizes(decltype(_epsilons) epsilons) {
        // allocate memory on device
        _epsilons = epsilons;
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets orthonormal tangent directions for reconstruction on a manifold. The first rank of this 2D array 
    //! corresponds to the target indices, i.e., rows of the neighbor lists 2D array. The second rank is the 
    //! ordinal of the tangent direction (spatial dimensions-1 are tangent, last one is normal), and the third 
    //! rank is indices into the spatial dimension.
    template<typename view_type>
    void setTangentBundle(view_type tangent_directions) {
        // accept input from user as a rank 3 tensor
        // but convert data to a rank 2 tensor with the last rank of dimension = _dimensions x _dimensions
        // this allows for nonstrided views on the device later

        // allocate memory on device
        _T = decltype(_T)("device tangent directions", _target_coordinates.extent(0)*_dimensions*_dimensions);

        compadre_assert_release( (std::is_same<decltype(_T)::memory_space, typename view_type::memory_space>::value) &&
                "Memory space does not match between _T and tangent_directions");

        auto this_dimensions = _dimensions;
        auto this_T = _T;
        // rearrange data on device from data given on host
        Kokkos::parallel_for("copy tangent vectors", Kokkos::RangePolicy<device_execution_space>(0, _target_coordinates.extent(0)), KOKKOS_LAMBDA(const int i) {
            scratch_matrix_right_type T(this_T.data() + i*this_dimensions*this_dimensions, this_dimensions, this_dimensions);
            for (int j=0; j<this_dimensions; ++j) {
                for (int k=0; k<this_dimensions; ++k) {
                    T(j,k) = tangent_directions(i, j, k);
                }
            }
        });
        _orthonormal_tangent_space_provided = true;

        // copy data from device back to host in rearranged format
        _host_T = Kokkos::create_mirror_view(_T);
        Kokkos::deep_copy(_host_T, _T);
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets outward normal direction. For manifolds this may be used for orienting surface. It is also accessible
    //! for sampling operators that require a normal direction.
    template<typename view_type>
    void setReferenceOutwardNormalDirection(view_type outward_normal_directions, bool use_to_orient_surface = true) {
        // accept input from user as a rank 2 tensor
        
        // allocate memory on device
        _ref_N = decltype(_ref_N)("device normal directions", _target_coordinates.extent(0)*_dimensions);
        // to assist LAMBDA capture
        auto this_ref_N = this->_ref_N;
        auto this_dimensions = this->_dimensions;

        // rearrange data on device from data given on host
        Kokkos::parallel_for("copy normal vectors", Kokkos::RangePolicy<device_execution_space>(0, _target_coordinates.extent(0)), KOKKOS_LAMBDA(const int i) {
            for (int j=0; j<this_dimensions; ++j) {
                this_ref_N(i*this_dimensions + j) = outward_normal_directions(i, j);
            }
        });
        Kokkos::fence();
        _reference_outward_normal_direction_provided = true;
        _use_reference_outward_normal_direction_provided_to_orient_surface = use_to_orient_surface;

        // copy data from device back to host in rearranged format
        _host_ref_N = Kokkos::create_mirror_view(_ref_N);
        Kokkos::deep_copy(_host_ref_N, _ref_N);
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets extra data to be used by sampling functionals in certain instances.
    template<typename view_type>
    void setSourceExtraData(view_type extra_data) {

        // allocate memory on device
        _source_extra_data = decltype(_source_extra_data)("device source extra data", extra_data.extent(0), extra_data.extent(1));

        auto host_extra_data = Kokkos::create_mirror_view(_source_extra_data);
        Kokkos::deep_copy(host_extra_data, extra_data);
        // copy data from host to device
        Kokkos::deep_copy(_source_extra_data, host_extra_data);
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets extra data to be used by sampling functionals in certain instances. (device)
    template<typename view_type>
    void setSourceExtraData(decltype(_source_extra_data) extra_data) {
        // allocate memory on device
        _source_extra_data = extra_data;
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets extra data to be used by target operations in certain instances.
    template<typename view_type>
    void setTargetExtraData(view_type extra_data) {

        // allocate memory on device
        _target_extra_data = decltype(_target_extra_data)("device target extra data", extra_data.extent(0), extra_data.extent(1));

        auto host_extra_data = Kokkos::create_mirror_view(_target_extra_data);
        Kokkos::deep_copy(host_extra_data, extra_data);
        // copy data from host to device
        Kokkos::deep_copy(_target_extra_data, host_extra_data);
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets extra data to be used by target operations in certain instances. (device)
    template<typename view_type>
    void setTargetExtraData(decltype(_target_extra_data) extra_data) {
        // allocate memory on device
        _target_extra_data = extra_data;
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets additional points for evaluation of target operation on polynomial reconstruction.
    //! If this is never called, then the target sites are the only locations where the target
    //! operations will be evaluated and applied to polynomial reconstructions.
    template <typename view_type>
    void setAuxiliaryEvaluationCoordinates(view_type evaluation_coordinates) {
        // allocate memory on device
        _additional_evaluation_coordinates = decltype(_additional_evaluation_coordinates)("device additional evaluation coordinates",
            evaluation_coordinates.extent(0), evaluation_coordinates.extent(1));

        typedef typename view_type::memory_space input_array_memory_space;
        if (std::is_same<input_array_memory_space, device_memory_space>::value) {
            // check if on the device, then copy directly
            // if it is, then it doesn't match the internal layout we use
            // then copy to the host mirror
            // switches potential layout mismatches
            Kokkos::deep_copy(_additional_evaluation_coordinates, evaluation_coordinates);
        } else {
            // if is on the host, copy to the host mirror
            // then copy to the device
            // switches potential layout mismatches
            auto host_additional_evaluation_coordinates = Kokkos::create_mirror_view(_additional_evaluation_coordinates);
            Kokkos::deep_copy(host_additional_evaluation_coordinates, evaluation_coordinates);
            // switches memory spaces
            Kokkos::deep_copy(_additional_evaluation_coordinates, host_additional_evaluation_coordinates);
        }
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets additional points for evaluation of target operation on polynomial reconstruction.
    //! If this is never called, then the target sites are the only locations where the target
    //! operations will be evaluated and applied to polynomial reconstructions. (device)
    template <typename view_type>
    void setAuxiliaryEvaluationCoordinates(decltype(_additional_evaluation_coordinates) evaluation_coordinates) {
        _additional_evaluation_coordinates = evaluation_coordinates;
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets the additional target evaluation coordinate indices list information. Should be # targets x maximum number of indices
    //! evaluation indices for any target + 1. first entry in every row should be the number of indices for the corresponding target.
    template <typename view_type>
    void setAuxiliaryEvaluationIndicesLists(view_type indices_lists) {
        // allocate memory on device
        _additional_evaluation_indices = decltype(_additional_evaluation_indices)("device additional evaluation indices",
            indices_lists.extent(0), indices_lists.extent(1));

        _host_additional_evaluation_indices = Kokkos::create_mirror_view(_additional_evaluation_indices);

        typedef typename view_type::memory_space input_array_memory_space;
        if (std::is_same<input_array_memory_space, device_memory_space>::value) {
            // check if on the device, then copy directly
            // if it is, then it doesn't match the internal layout we use
            // then copy to the host mirror
            // switches potential layout mismatches
            Kokkos::deep_copy(_additional_evaluation_indices, indices_lists);
            Kokkos::deep_copy(_host_additional_evaluation_indices, _additional_evaluation_indices);
        } else {
            // if is on the host, copy to the host mirror
            // then copy to the device
            // switches potential layout mismatches
            Kokkos::deep_copy(_host_additional_evaluation_indices, indices_lists);
            // copy data from host to device
            Kokkos::deep_copy(_additional_evaluation_indices, _host_additional_evaluation_indices);
        }

        _max_evaluation_sites_per_target = 1;
        auto number_of_additional_evaluation_indices = _number_of_additional_evaluation_indices;
        auto additional_evaluation_indices = _additional_evaluation_indices;
        Kokkos::parallel_reduce("additional evaluation indices", 
                Kokkos::RangePolicy<device_execution_space>(0, _additional_evaluation_indices.extent(0)), 
                KOKKOS_LAMBDA(const int i, int& t_max_evaluation_sites_per_target) {
            number_of_additional_evaluation_indices(i) = additional_evaluation_indices(i,0);
            t_max_evaluation_sites_per_target = (t_max_evaluation_sites_per_target > number_of_additional_evaluation_indices(i)+1) 
                                                ? t_max_evaluation_sites_per_target : number_of_additional_evaluation_indices(i)+1;
        }, Kokkos::Max<int>(_max_evaluation_sites_per_target));
        Kokkos::fence();
        this->resetCoefficientData();
    }

    //! (OPTIONAL)
    //! Sets the additional target evaluation coordinate indices list information. Should be # targets x maximum number of indices
    //! evaluation indices for any target + 1. first entry in every row should be the number of indices for the corresponding target.
    template <typename view_type>
    void setAuxiliaryEvaluationIndicesLists(decltype(_additional_evaluation_indices) indices_lists) {
        // allocate memory on device
        _additional_evaluation_indices = indices_lists;

        _host_additional_evaluation_indices = Kokkos::create_mirror_view(_additional_evaluation_indices);
        // copy data from host to device
        Kokkos::deep_copy(_host_additional_evaluation_indices, _additional_evaluation_indices);

        _max_evaluation_sites_per_target = 1;
        auto number_of_additional_evaluation_indices = _number_of_additional_evaluation_indices;
        auto additional_evaluation_indices = _additional_evaluation_indices;
        Kokkos::parallel_reduce("additional evaluation indices", 
                Kokkos::RangePolicy<device_execution_space>(0, _additional_evaluation_indices.extent(0)), 
                KOKKOS_LAMBDA(const int i, int& t_max_evaluation_sites_per_target) {
            number_of_additional_evaluation_indices(i) = additional_evaluation_indices(i,0);
            t_max_evaluation_sites_per_target = (t_max_evaluation_sites_per_target > number_of_additional_evaluation_indices(i)+1) 
                                                ? t_max_evaluation_sites_per_target : number_of_additional_evaluation_indices(i)+1;
        }, Kokkos::Max<int>(_max_evaluation_sites_per_target));
        Kokkos::fence();
        this->resetCoefficientData();
    }


    //! Type for weighting kernel for GMLS problem
    void setWeightingType( const std::string &wt) {
        std::string wt_to_lower = wt;
        transform(wt_to_lower.begin(), wt_to_lower.end(), wt_to_lower.begin(), ::tolower);
        if (wt_to_lower == "power") {
            _weighting_type = WeightingFunctionType::Power;
        } else if (wt_to_lower == "gaussian") {
            _weighting_type = WeightingFunctionType::Gaussian;
        } else if (wt_to_lower == "cubicspline") {
            _weighting_type = WeightingFunctionType::CubicSpline;
        } else {
            // Power is default
            _weighting_type = WeightingFunctionType::Power;
        }
        this->resetCoefficientData();
    }

    //! Type for weighting kernel for GMLS problem
    void setWeightingType( const WeightingFunctionType wt) {
        _weighting_type = wt;
        this->resetCoefficientData();
    }

    //! Type for weighting kernel for curvature 
    void setCurvatureWeightingType( const std::string &wt) {
        std::string wt_to_lower = wt;
        transform(wt_to_lower.begin(), wt_to_lower.end(), wt_to_lower.begin(), ::tolower);
        if (wt_to_lower == "power") {
            _curvature_weighting_type = WeightingFunctionType::Power;
        } else if (wt_to_lower == "gaussian") {
            _curvature_weighting_type = WeightingFunctionType::Gaussian;
        } else if (wt_to_lower == "cubicspline") {
            _curvature_weighting_type = WeightingFunctionType::CubicSpline;
        } else {
            // Power is default
            _curvature_weighting_type = WeightingFunctionType::Power;
        }
        this->resetCoefficientData();
    }

    //! Type for weighting kernel for curvature
    void setCurvatureWeightingType( const WeightingFunctionType wt) {
        _curvature_weighting_type = wt;
        this->resetCoefficientData();
    }

    //! Sets basis order to be used when reoncstructing any function
    void setPolynomialOrder(const int poly_order) {
        _poly_order = poly_order;
        _NP = this->getNP(_poly_order, _dimensions, _reconstruction_space);
        this->resetCoefficientData();
    }

    //! Sets basis order to be used when reoncstructing curvature
    void setCurvaturePolynomialOrder(const int manifold_poly_order) {
        _curvature_poly_order = manifold_poly_order;
        this->resetCoefficientData();
    }

    //! Power for weighting kernel for GMLS problem
    void setWeightingPower(int wp) { 
        _weighting_power = wp;
        this->resetCoefficientData();
    }

    //! Power for weighting kernel for curvature
    void setCurvatureWeightingPower(int wp) { 
        _curvature_weighting_power = wp;
        this->resetCoefficientData();
    }

    //! Number quadrature points to use
    void setOrderOfQuadraturePoints(int order) { 
        _order_of_quadrature_points = order;
        this->resetCoefficientData();
    }

    //! Dimensions of quadrature points to use
    void setDimensionOfQuadraturePoints(int dim) { 
        _dimension_of_quadrature_points = dim;
        this->resetCoefficientData();
    }

    //! Type of quadrature points
    void setQuadratureType(std::string quadrature_type) { 
        _quadrature_type = quadrature_type;
        this->resetCoefficientData();
    }

    //! Adds a target to the vector of target functional to be applied to the reconstruction
    void addTargets(TargetOperation lro) {
        std::vector<TargetOperation> temporary_lro_vector(1, lro);
        this->addTargets(temporary_lro_vector);
        this->resetCoefficientData();
    }

    //! Adds a vector of target functionals to the vector of target functionals already to be applied to the reconstruction
    void addTargets(std::vector<TargetOperation> lro) {
        // if called multiple times with different dimensions, only the last
        // dimension called with is used for all

        // loop over requested targets
        for (size_t i=0; i<lro.size(); ++i) {

            bool operation_found = false;
            // loop over existing targets registered
            for (size_t j=0; j<_lro.size(); ++j) {

                // if found
                if (_lro[j]==lro[i]) {

                    operation_found = true;

                    // the operation should now point to where the operation is stored
                    _lro_lookup[(int)lro[i]] = j;

                    break;

                }
            }

            if (!operation_found) {
                _lro_lookup[(int)lro[i]] = _lro.size();
                _lro.push_back(lro[i]);
            }
        }

        _lro_total_offsets = Kokkos::View<int*>("total offsets for alphas", _lro.size());
        _lro_output_tile_size = Kokkos::View<int*>("output tile size for each operation", _lro.size());
        _lro_input_tile_size = Kokkos::View<int*>("output tile size for each operation", _lro.size());
        _lro_output_tensor_rank = Kokkos::View<int*>("output tensor rank", _lro.size());
        _lro_input_tensor_rank = Kokkos::View<int*>("input tensor rank", _lro.size());

        _host_lro_total_offsets = create_mirror_view(_lro_total_offsets);
        _host_lro_output_tile_size = create_mirror_view(_lro_output_tile_size);
        _host_lro_input_tile_size = create_mirror_view(_lro_input_tile_size);
        _host_lro_output_tensor_rank = create_mirror_view(_lro_output_tensor_rank);
        _host_lro_input_tensor_rank = create_mirror_view(_lro_input_tensor_rank);

        int total_offset = 0; // need total offset
        int output_offset = 0;
        int input_offset = 0;
        for (size_t i=0; i<_lro.size(); ++i) {
            _host_lro_total_offsets(i) = total_offset;

            // allows for a tile of the product of dimension^input_tensor_rank * dimension^output_tensor_rank * the number of neighbors
            int output_tile_size = getOutputDimensionOfOperation(_lro[i], false /* uses _local_dimension */);

            // the target functional input indexing is sized based on the output rank of the sampling
            // functional used
            int input_tile_size = getOutputDimensionOfSampling(_data_sampling_functional);
            _host_lro_output_tile_size(i) = output_tile_size;
            _host_lro_input_tile_size(i) = input_tile_size;

            total_offset += input_tile_size * output_tile_size;
            output_offset += output_tile_size;
            input_offset += input_tile_size;

            // the target functional output rank is based on the output rank of the sampling
            // functional used
            _host_lro_input_tensor_rank(i) = _data_sampling_functional.output_rank;
            _host_lro_output_tensor_rank(i) = TargetOutputTensorRank[(int)_lro[i]];
        }

        _total_alpha_values = total_offset;

        if (_problem_type == ProblemType::MANIFOLD) {
            // if on a manifold, the total alphas values must be large enough to hold the gradient
            // of the geometry reconstruction
            _total_alpha_values = (_total_alpha_values > std::pow(_local_dimensions, 1)) ? 
                _total_alpha_values : std::pow(_local_dimensions, 1);
        }

        Kokkos::deep_copy(_lro_total_offsets, _host_lro_total_offsets);
        Kokkos::deep_copy(_lro_output_tile_size, _host_lro_output_tile_size);
        Kokkos::deep_copy(_lro_input_tile_size, _host_lro_input_tile_size);
        Kokkos::deep_copy(_lro_output_tensor_rank, _host_lro_output_tensor_rank);
        Kokkos::deep_copy(_lro_input_tensor_rank, _host_lro_input_tensor_rank);
        this->resetCoefficientData();
    }

    //! Empties the vector of target functionals to apply to the reconstruction
    void clearTargets() {
        _lro.clear();
        for (int i=0; i<TargetOperation::COUNT; ++i) {
            _lro_lookup[i] = -1;
        }
        this->resetCoefficientData();
    }

    /*! \brief Generates polynomial coefficients by setting up and solving least squares problems
    //! Sets up the batch of GMLS problems to be solved for. Provides alpha values
    //! that can later be contracted against data or degrees of freedom to form a
    //! global linear system.
    //! \param number_of_batches    [in] - how many batches to break up the total workload into (for storage)
    //! \param keep_coefficients    [in] - whether to store (P^T W P)^-1 * P^T * W
    */
    void generatePolynomialCoefficients(const int number_of_batches = 1, const bool keep_coefficients = false);

    /*! \brief Meant to calculate target operations and apply the evaluations to the previously 
    //! constructed polynomial coefficients. But now that is inside of generatePolynomialCoefficients because
    //! it must be to handle number_of_batches>1. Effectively, this just calls generatePolynomialCoefficients.
    //! \param number_of_batches    [in] - how many batches to break up the total workload into (for storage)
    //! \param keep_coefficients    [in] - whether to store (P^T W P)^-1 * P^T * W
    */
    void generateAlphas(const int number_of_batches = 1, const bool keep_coefficients = false);

///@}




}; // GMLS Class
} // Compadre

#endif


