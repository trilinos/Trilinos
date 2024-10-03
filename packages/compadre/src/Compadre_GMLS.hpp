// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_GMLS_HPP_
#define _COMPADRE_GMLS_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"

#include "Compadre_Misc.hpp"
#include "Compadre_Operators.hpp"
#include "Compadre_LinearAlgebra_Definitions.hpp"
#include "Compadre_ParallelManager.hpp"
#include "Compadre_Quadrature.hpp"
#include "Compadre_SolutionSet.hpp"
#include "Compadre_ScalarTaylorPolynomial.hpp"
#include "Compadre_DivergenceFreePolynomial.hpp"
#include "Compadre_BernsteinPolynomial.hpp"
#include "Compadre_NeighborLists.hpp"
#include "Compadre_PointConnections.hpp"

namespace Compadre {

class  Evaluator;
struct GMLSBasisData;
struct GMLSSolutionData;

//!  Generalized Moving Least Squares (GMLS)
/*!
*  This class sets up a batch of GMLS problems from a given set of neighbor lists, target sites, and source sites.
*  GMLS requires a target functional, reconstruction space, and sampling functional to be specified. 
*  For a given choice of reconstruction space and sampling functional, multiple targets can be generated with very little
*  additional computation, which is why this class allows for multiple target functionals to be specified. 
*/
class GMLS {

friend class Evaluator;

public:

    typedef PointConnections<Kokkos::View<double**, layout_right>, 
            Kokkos::View<double**, layout_right>, 
            NeighborLists<Kokkos::View<int*> > > 
                point_connections_type;

    typedef NeighborLists<Kokkos::View<int*> > neighbor_lists_type;

    typedef Kokkos::View<double**, layout_right> coordinates_type; 

private:

    // matrices that may be needed for matrix factorization on the device
    // supports batched matrix factorization dispatch

    //! contains weights for all problems
    Kokkos::View<double*> _w; 

    //! P*sqrt(w) matrix for all problems
    Kokkos::View<double*> _P;

    //! sqrt(w)*Identity matrix for all problems, later holds polynomial coefficients for all problems
    Kokkos::View<double*> _RHS;

    //! stores evaluations of targets applied to basis
    Kokkos::View<double*> _Z;

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

    //! connections between points and neighbors
    point_connections_type _pc;

    //! h supports determined through neighbor search (device)
    Kokkos::View<double*> _epsilons; 
    
    //! generated weights for nontraditional samples required to transform data into expected sampling 
    //! functional form (device). 
    Kokkos::View<double*****, layout_right> _prestencil_weights; 

    //! generated weights for nontraditional samples required to transform data into expected sampling 
    //! functional form (host)
    Kokkos::View<const double*****, layout_right>::HostMirror _host_prestencil_weights;

    //! (OPTIONAL) connections between additional points and neighbors
    point_connections_type _additional_pc;

    //! Solution Set (contains all alpha values from solution and alpha layout methods)
    // _h_ss is private so that getSolutionSetHost() must be called
    // which ensures that the copy of the solution to device is necessary
    SolutionSet<host_memory_space> _h_ss;
    SolutionSet<device_memory_space> _d_ss;

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

    //! problem type for GMLS problem, can also be set to STANDARD for normal or MANIFOLD for manifold problems <br>
    //! <b>NOTE: can only be set at object instantiation</b>
    ProblemType _problem_type;

    //! constraint type for GMLS problem
    ConstraintType _constraint_type;

    //! polynomial sampling functional used to construct P matrix, set at GMLS class instantiation <br>
    //! <b>NOTE: can only be set at object instantiation</b>
    SamplingFunctional _polynomial_sampling_functional;

    //! generally the same as _polynomial_sampling_functional, but can differ if specified at 
    //  can only be set at object instantiation
    //! GMLS class instantiation
    SamplingFunctional _data_sampling_functional;

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

    //! first parameter to be used for weighting kernel
    int _weighting_p;

    //! second parameter to be used for weighting kernel
    int _weighting_n;

    //! first parameter to be used for weighting kernel for curvature
    int _curvature_weighting_p;

    //! second parameter to be used for weighting kernel for curvature
    int _curvature_weighting_n;

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

private:

/** @name Private Modifiers
 *  Private function because information lives on the device
 */
///@{

    //! (OPTIONAL)
    //! Sets additional points for evaluation of target operation on polynomial reconstruction.
    //! If this is never called, then the target sites are the only locations where the target
    //! operations will be evaluated and applied to polynomial reconstructions.
    template <typename view_type>
    void setAuxiliaryEvaluationCoordinates(view_type evaluation_coordinates) {
        // allocate memory on device
        auto additional_evaluation_coordinates = coordinates_type("device additional evaluation coordinates",
            evaluation_coordinates.extent(0), evaluation_coordinates.extent(1));

        typedef typename view_type::memory_space input_array_memory_space;
        if (std::is_same<input_array_memory_space, device_memory_space>::value) {
            // check if on the device, then copy directly
            // if it is, then it doesn't match the internal layout we use
            // then copy to the host mirror
            // switches potential layout mismatches
            Kokkos::deep_copy(additional_evaluation_coordinates, evaluation_coordinates);
        } else {
            // if is on the host, copy to the host mirror
            // then copy to the device
            // switches potential layout mismatches
            auto host_additional_evaluation_coordinates = Kokkos::create_mirror_view(additional_evaluation_coordinates);
            Kokkos::deep_copy(host_additional_evaluation_coordinates, evaluation_coordinates);
            // switches memory spaces
            Kokkos::deep_copy(additional_evaluation_coordinates, host_additional_evaluation_coordinates);
        }
        this->resetCoefficientData();
        _additional_pc.setSourceCoordinates(additional_evaluation_coordinates);
    }

    //! (OPTIONAL)
    //! Sets additional points for evaluation of target operation on polynomial reconstruction.
    //! If this is never called, then the target sites are the only locations where the target
    //! operations will be evaluated and applied to polynomial reconstructions. (device)
    template <typename view_type>
    void setAuxiliaryEvaluationCoordinates(coordinates_type evaluation_coordinates) {
        this->resetCoefficientData();
        _additional_pc.setSourceCoordinates(evaluation_coordinates);
    }

    //! (OPTIONAL)
    //! Sets the additional target evaluation indices list information from compressed row format (if same view_type)
    template <typename view_type>
    typename std::enable_if<view_type::rank==1&&std::is_same<neighbor_lists_type::internal_view_type,view_type>::value==1, void>::type 
            setAuxiliaryEvaluationIndicesLists(view_type additional_evaluation_indices, view_type number_of_neighbors_list) {

        auto additional_nla = NeighborLists<view_type>(additional_evaluation_indices, number_of_neighbors_list);
        this->resetCoefficientData();
        _additional_pc.setNeighborLists(additional_nla);

    }

    //! (OPTIONAL)
    //! Sets the additional target evaluation indices list information from compressed row format (if different view_type)
    template <typename view_type>
    typename std::enable_if<view_type::rank==1&&std::is_same<neighbor_lists_type::internal_view_type,view_type>::value==0, void>::type 
            setAuxiliaryEvaluationIndicesLists(view_type additional_evaluation_indices, view_type number_of_neighbors_list) {

        typedef neighbor_lists_type::internal_view_type gmls_view_type;
        gmls_view_type d_additional_evaluation_indices("compressed row additional evaluation indices lists data", additional_evaluation_indices.extent(0));
        gmls_view_type d_number_of_neighbors_list("number of additional evaluation indices", number_of_neighbors_list.extent(0));
        Kokkos::deep_copy(d_additional_evaluation_indices, additional_evaluation_indices);
        Kokkos::deep_copy(d_number_of_neighbors_list, number_of_neighbors_list);
        Kokkos::fence();
        auto additional_nla = NeighborLists<gmls_view_type>(d_additional_evaluation_indices, d_number_of_neighbors_list);
        this->resetCoefficientData();
        _additional_pc.setNeighborLists(additional_nla);

    }

    //! (OPTIONAL)
    //! Sets the additional target evaluation indices list information. Should be # targets x maximum number of indices
    //! evaluation indices for any target + 1. first entry in every row should be the number of indices for the corresponding target.
    template <typename view_type>
    typename std::enable_if<view_type::rank==2, void>::type setAuxiliaryEvaluationIndicesLists(view_type additional_evaluation_indices) {
    
        auto additional_nla = Convert2DToCompressedRowNeighborLists<decltype(additional_evaluation_indices), Kokkos::View<int*> >(additional_evaluation_indices);
        this->resetCoefficientData();
        _additional_pc.setNeighborLists(additional_nla);

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
        if (solver_type_to_lower == "lu") {
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

    //! Maximal constructor, not intended for users
    GMLS(ReconstructionSpace reconstruction_space,
        const SamplingFunctional polynomial_sampling_strategy,
        const SamplingFunctional data_sampling_strategy,
        const int poly_order,
        const int dimensions,
        const DenseSolverType dense_solver_type,
        const ProblemType problem_type,
        const ConstraintType constraint_type,
        const int manifold_curvature_poly_order) : 
            _poly_order(poly_order),
            _curvature_poly_order(manifold_curvature_poly_order),
            _dimensions(dimensions),
            _reconstruction_space(reconstruction_space),
            _dense_solver_type(dense_solver_type),
            _problem_type(problem_type),
            _constraint_type(constraint_type),
            _polynomial_sampling_functional(((_problem_type == ProblemType::MANIFOLD) 
                        && (polynomial_sampling_strategy == VectorPointSample)) ? ManifoldVectorPointSample : polynomial_sampling_strategy),
            _data_sampling_functional(((_problem_type == ProblemType::MANIFOLD) 
                        && (data_sampling_strategy == VectorPointSample)) ? ManifoldVectorPointSample : data_sampling_strategy)
            {

        compadre_assert_release(poly_order<11 && "Unsupported polynomial order (>=11).");
        compadre_assert_release(manifold_curvature_poly_order<11 && "Unsupported curvature polynomial order (>=11).");
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

        // various initializations

        _weighting_type = WeightingFunctionType::Power;
        _curvature_weighting_type = WeightingFunctionType::Power;
        _weighting_p = 2;
        _weighting_n = 1;
        _curvature_weighting_p = 2;
        _curvature_weighting_n = 1;

        _reconstruction_space_rank = getActualReconstructionSpaceRank(_reconstruction_space);

        _basis_multiplier = 1;
        _sampling_multiplier = 1;

        _orthonormal_tangent_space_provided = false; 
        _reference_outward_normal_direction_provided = false;
        _use_reference_outward_normal_direction_provided_to_orient_surface = false;
        _entire_batch_computed_at_once = true;
        _store_PTWP_inv_PTW = false;

        _initial_index_for_batch = 0;

        _global_dimensions = dimensions;
        if (_problem_type == ProblemType::MANIFOLD) {
            _local_dimensions = dimensions-1;
        } else {
            _local_dimensions = dimensions;
        }

        _order_of_quadrature_points = 0;
        _dimension_of_quadrature_points = 0;

        _h_ss = SolutionSet<host_memory_space>(
                _data_sampling_functional,
                _dimensions, 
                _local_dimensions,
                _problem_type);
        _h_ss._max_evaluation_sites_per_target = 1;
    }

    //! Maximal constructor, but with string arguments
    GMLS(ReconstructionSpace reconstruction_space,
        const SamplingFunctional polynomial_sampling_strategy,
        const SamplingFunctional data_sampling_strategy,
        const int poly_order,
        const int dimensions = 3,
        const std::string dense_solver_type = std::string("QR"),
        const std::string problem_type = std::string("STANDARD"),
        const std::string constraint_type = std::string("NO_CONSTRAINT"),
        const int manifold_curvature_poly_order = 2)
      : GMLS(reconstruction_space, polynomial_sampling_strategy, data_sampling_strategy, poly_order, dimensions, parseSolverType(dense_solver_type), parseProblemType(problem_type), parseConstraintType(constraint_type), manifold_curvature_poly_order) {}

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
        if (r_space == ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial) {
            return DivergenceFreePolynomialBasis::getSize(m, dimension);
        } else if (r_space == ReconstructionSpace::BernsteinPolynomial) {
            return BernsteinPolynomialBasis::getSize(m, dimension);
        } else {
            return ScalarTaylorPolynomialBasis::getSize(m, dimension);
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
        \param p                [in] - parameter to be given to the kernel (see Wab definition for context).
        \param n                [in] - parameter to be given to the kernel (see Wab definition for context).
    */
    KOKKOS_INLINE_FUNCTION
    static double Wab(const double r, const double h, const WeightingFunctionType& weighting_type, const int p, const int n) {
        if (weighting_type == WeightingFunctionType::Power) {
            // compactly supported on [0,h]
            // (1 - |r/h|^n)^p
            // p=0,n=1 -> Uniform, boxcar
            // p=1,n=1 -> triangular
            // p=1,n=2 -> Epanechnikov, parabolic
            // p=2,n=2 -> Quartic, biweight
            // p=3,n=2 -> Triweight
            // p=3,n=3 -> Tricube
            double abs_r_over_h_to_n = std::abs(r/h);
            if (n>1) abs_r_over_h_to_n = std::pow(abs_r_over_h_to_n, n);
            return std::pow(1.0-abs_r_over_h_to_n, p) * double(1.0-abs_r_over_h_to_n>0.0);
        } else if (weighting_type == WeightingFunctionType::CubicSpline) {
            // compactly supported on [0,h]
            // invariant to p and n
            double x = std::abs(r/h);
            return ((1-x)+x*(1-x)*(1-2*x)) * double(x<=1.0);
        } else if (weighting_type == WeightingFunctionType::CardinalCubicBSpline) {
            // compactly supported on [0,h]
            // invariant to p and n
            // Calculate the value using a cardinal cubic b-spline kernel (often just called cubic b spline)
            double x = std::abs(r/h);
            if (x < 0.5) return 1.0 + 6.0 * x * x * (-1.0 + x);
            if (x < 1.0) return 2.0 * (1.0 + x * (-3.0 + 3.0 * x - 1.0 * x * x));
            return 0.0;
        } else if (weighting_type == WeightingFunctionType::Cosine) {
            // compactly supported on [0,h]
            double pi = 3.14159265358979323846;
            double abs_r_over_h_to_n = std::abs(r/h);
            return std::cos(0.5*pi*r/h) * double(1.0-abs_r_over_h_to_n>0.0);
        } else if (weighting_type == WeightingFunctionType::Gaussian) {
            // NOT compactly supported on [0,h], but approximately 0 at h with >> p
            // invariant to n, p is number of standard deviations at distance h
            // 2.5066282746310002416124 = sqrt(2*pi)
            double h_over_p = h/p;
            double abs_r_over_h_to_n = std::abs(r/h);
            return double(1.0-abs_r_over_h_to_n>0.0)/( h_over_p * 2.5066282746310002416124 ) * std::exp(-.5*r*r/(h_over_p*h_over_p));
        } else if (weighting_type == WeightingFunctionType::Sigmoid) {
            // NOT compactly supported on [0,h], but approximately 0 at h with >> p
            // n=0 is sigmoid, n==2 is logistic, with larger p making Wab decay more quickly
            double abs_r_over_h_to_n = std::abs(r/h);
            return  double(1.0-abs_r_over_h_to_n>0.0) / (std::exp(p*r) + std::exp(-p*r) + n);
        } else { // unsupported type
            compadre_kernel_assert_release(false && "Invalid WeightingFunctionType selected.");
            return 0; 
        }
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
        sizes(1) = _sampling_multiplier*this->getNeighborLists()->getMaxNumNeighbors();
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
            getRHSDims(_dense_solver_type, _constraint_type, _reconstruction_space, _dimensions, M_by_N[1], M_by_N[0], sizes(0), sizes(1));
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
    DenseSolverType getDenseSolverType() const { return _dense_solver_type; }

    //! Get problem type
    ProblemType getProblemType() const { return _problem_type; }

    //! Get constraint type
    ConstraintType getConstraintType() const { return _constraint_type; }

    //! Get basis order used for reconstruction
    int getPolynomialOrder() const { return _poly_order; }

    //! Get basis order used for curvature reconstruction
    int getCurvaturePolynomialOrder() const { return _curvature_poly_order; }

    //! Type for weighting kernel for GMLS problem
    WeightingFunctionType getWeightingType() const { return _weighting_type; }

    //! Type for weighting kernel for curvature 
    WeightingFunctionType getManifoldWeightingType() const { return _curvature_weighting_type; }

    //! Get parameter for weighting kernel for GMLS problem
    int getWeightingParameter(const int index = 0) const { 
        if (index==1) {
            return _weighting_n;
        } else {
            return _weighting_p; 
        }
    }

    //! Get parameter for weighting kernel for curvature
    int getManifoldWeightingParameter(const int index = 0) const {
        if (index==1) {
            return _curvature_weighting_n;
        } else {
            return _curvature_weighting_p; 
        }
    }

    //! Number of quadrature points
    int getNumberOfQuadraturePoints() const { return _qm.getNumberOfQuadraturePoints(); }

    //! Order of quadrature points
    int getOrderOfQuadraturePoints() const { return _order_of_quadrature_points; }

    //! Dimensions of quadrature points
    int getDimensionOfQuadraturePoints() const { return _dimension_of_quadrature_points; }

    //! Type of quadrature points
    std::string getQuadratureType() const { return _quadrature_type; }

    //! Get neighbor list accessor
    neighbor_lists_type* getNeighborLists() const { 
        return const_cast<neighbor_lists_type*>(&_pc._nla); 
    }

    //! Get a view (device) of all point connection info
    decltype(_pc)* getPointConnections() { return &_pc; }

    //! (OPTIONAL) Get additional evaluation sites neighbor list-like accessor
    neighbor_lists_type* getAdditionalEvaluationIndices() const { 
        return const_cast<neighbor_lists_type*>(&_additional_pc._nla);
    }

    //! (OPTIONAL) Get a view (device) of all additional evaluation point connection info
    decltype(_additional_pc)* getAdditionalPointConnections() { return &_additional_pc; }

    //! Get a view (device) of all window sizes
    decltype(_epsilons)* getWindowSizes() { return &_epsilons; }

    //! Get a view (device) of all tangent direction bundles.
    decltype(_T)* getTangentDirections() { return &_T; }

    //! Get a view (device) of all reference outward normal directions.
    decltype(_ref_N)* getReferenceNormalDirections() { return &_ref_N; }

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

    //! Get solution set on host
    decltype(_h_ss)* getSolutionSetHost(bool alpha_validity_check=true) { 
        if (!_h_ss._contains_valid_alphas && _d_ss._contains_valid_alphas) {
            // solution solved for on device, but now solution
            // requested on the host
            _h_ss.copyAlphas(_d_ss);
        }
        compadre_assert_release((!alpha_validity_check || _h_ss._contains_valid_alphas) &&
                "getSolutionSetHost() called with invalid alpha values.");
        return &_h_ss; 
    }

    //! Get solution set on device
    decltype(_d_ss)* getSolutionSetDevice(bool alpha_validity_check=true) { 
        compadre_assert_release((!alpha_validity_check || _d_ss._contains_valid_alphas) &&
                "getSolutionSetDevice() called with invalid alpha values.");
        return &_d_ss; 
    }

    //! Check if GMLS solution set contains valid alpha values (has generateAlphas been called)
    bool containsValidAlphas() const { return this->_d_ss._contains_valid_alphas; }

    //! Get GMLS solution data
    const GMLSSolutionData extractSolutionData() const;

    //! Get GMLS basis data
    const GMLSBasisData extractBasisData() const;

///@}


/** @name Modifiers
 *  Changed member variables through public member functions
 */
///@{

    void resetCoefficientData() {
        if (_RHS.extent(0) > 0)
            _RHS = Kokkos::View<double*>("RHS",0);
        _h_ss._contains_valid_alphas = false;
        _d_ss._contains_valid_alphas = false;
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

    //! (OPTIONAL) Sets additional evaluation sites for each target site
    template<typename view_type_1, typename view_type_2>
    void setAdditionalEvaluationSitesData(
            view_type_1 cr_additional_evaluation_indices,
            view_type_1 number_of_additional_evaluation_indices,
            view_type_2 additional_evaluation_coordinates) {
        this->setAuxiliaryEvaluationIndicesLists<view_type_1>(cr_additional_evaluation_indices, 
                number_of_additional_evaluation_indices);
        this->setAuxiliaryEvaluationCoordinates<view_type_2>(additional_evaluation_coordinates);
    }

    //! Sets neighbor list information from compressed row neighborhood lists data (if same view_type).
    template <typename view_type>
    typename std::enable_if<view_type::rank==1&&std::is_same<neighbor_lists_type::internal_view_type,view_type>::value==1, void>::type 
            setNeighborLists(view_type neighbor_lists, view_type number_of_neighbors_list) {

        auto nla = NeighborLists<view_type>(neighbor_lists, number_of_neighbors_list);
        this->resetCoefficientData();
        _pc.setNeighborLists(nla);
    }

    //! Sets neighbor list information from compressed row neighborhood lists data (if different view_type).
    template <typename view_type>
    typename std::enable_if<view_type::rank==1&&std::is_same<neighbor_lists_type::internal_view_type,view_type>::value==0, void>::type 
            setNeighborLists(view_type neighbor_lists, view_type number_of_neighbors_list) {

        typedef neighbor_lists_type::internal_view_type gmls_view_type;
        gmls_view_type d_neighbor_lists("compressed row neighbor lists data", neighbor_lists.extent(0));
        gmls_view_type d_number_of_neighbors_list("number of neighbors list", number_of_neighbors_list.extent(0));
        Kokkos::deep_copy(d_neighbor_lists, neighbor_lists);
        Kokkos::deep_copy(d_number_of_neighbors_list, number_of_neighbors_list);
        Kokkos::fence();
        auto nla = NeighborLists<gmls_view_type>(d_neighbor_lists, d_number_of_neighbors_list);
        this->resetCoefficientData();
        _pc.setNeighborLists(nla);
    }

    //! Sets neighbor list information. Should be # targets x maximum number of neighbors for any target + 1.
    //! first entry in ever row should be the number of neighbors for the corresponding target.
    template <typename view_type>
    typename std::enable_if<view_type::rank==2, void>::type setNeighborLists(view_type neighbor_lists) {
    
        auto nla = Convert2DToCompressedRowNeighborLists<decltype(neighbor_lists), Kokkos::View<int*> >(neighbor_lists);
        this->resetCoefficientData();
        _pc.setNeighborLists(nla);
    }

    //! Sets source coordinate information. Rows of this 2D-array should correspond to neighbor IDs contained in the entries
    //! of the neighbor lists 2D array.
    template<typename view_type>
    void setSourceSites(view_type source_coordinates) {

        // allocate memory on device
        auto sc = coordinates_type("device neighbor coordinates",
                source_coordinates.extent(0), source_coordinates.extent(1));

        typedef typename view_type::memory_space input_array_memory_space;
        if (std::is_same<input_array_memory_space, device_memory_space>::value) {
            // check if on the device, then copy directly
            // if it is, then it doesn't match the internal layout we use
            // then copy to the host mirror
            // switches potential layout mismatches
            Kokkos::deep_copy(sc, source_coordinates);
        } else {
            // if is on the host, copy to the host mirror
            // then copy to the device
            // switches potential layout mismatches
            auto host_source_coordinates = Kokkos::create_mirror_view(sc);
            Kokkos::deep_copy(host_source_coordinates, source_coordinates);
            // switches memory spaces
            Kokkos::deep_copy(sc, host_source_coordinates);
        }
        this->resetCoefficientData();
        _pc.setSourceCoordinates(sc);
    }

    //! Sets source coordinate information. Rows of this 2D-array should correspond to neighbor IDs contained in the entries
    //! of the neighbor lists 2D array.
    template<typename view_type>
    void setSourceSites(coordinates_type source_coordinates) {
        // allocate memory on device
        this->resetCoefficientData();
        _pc.setSourceCoordinates(source_coordinates);
    }

    //! Sets target coordinate information. Rows of this 2D-array should correspond to rows of the neighbor lists.
    template<typename view_type>
    void setTargetSites(view_type target_coordinates) {
        // allocate memory on device
        auto tc = coordinates_type("device target coordinates",
                target_coordinates.extent(0), target_coordinates.extent(1));

        typedef typename view_type::memory_space input_array_memory_space;
        if (std::is_same<input_array_memory_space, device_memory_space>::value) {
            // check if on the device, then copy directly
            // if it is, then it doesn't match the internal layout we use
            // then copy to the host mirror
            // switches potential layout mismatches
            Kokkos::deep_copy(tc, target_coordinates);
        } else {
            // if is on the host, copy to the host mirror
            // then copy to the device
            // switches potential layout mismatches
            auto host_target_coordinates = Kokkos::create_mirror_view(tc);
            Kokkos::deep_copy(host_target_coordinates, target_coordinates);
            // switches memory spaces
            Kokkos::deep_copy(tc, host_target_coordinates);
        }
        if (this->getAdditionalEvaluationIndices()->getNumberOfTargets() != target_coordinates.extent_int(0)) {
            this->setAuxiliaryEvaluationIndicesLists(
                    Kokkos::View<int*>(),
                    Kokkos::View<int*>("number of additional evaluation indices", 
                                       target_coordinates.extent(0))
            );
        }
        this->resetCoefficientData();
        _pc.setTargetCoordinates(tc);
        if (_additional_pc._target_coordinates.extent(0) != _pc._target_coordinates.extent(0)) {
            _additional_pc.setTargetCoordinates(tc);
        }
    }

    //! Sets target coordinate information. Rows of this 2D-array should correspond to rows of the neighbor lists.
    template<typename view_type>
    void setTargetSites(coordinates_type target_coordinates) {
        if (this->getAdditionalEvaluationIndices()->getNumberOfTargets() != target_coordinates.extent(0)) {
            this->setAuxiliaryEvaluationIndicesLists(
                    Kokkos::View<int*>(),
                    Kokkos::View<int*>("number of additional evaluation indices", 
                                       target_coordinates.extent(0))
            );
        }
        this->resetCoefficientData();
        _pc.setTargetCoordinates(target_coordinates);
        if (_additional_pc._target_coordinates.extent(0) != _pc._target_coordinates.extent(0)) {
            _additional_pc.setTargetCoordinates(target_coordinates);
        }
    }

    //! Sets window sizes, also called the support of the kernel
    template<typename view_type>
    void setWindowSizes(view_type epsilons) {

        // allocate memory on device
        _epsilons = decltype(_epsilons)("device epsilons", epsilons.extent(0));

        // copy data from host to device
        Kokkos::deep_copy(_epsilons, epsilons);
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
        compadre_assert_release( tangent_directions.size() % (_dimensions*_dimensions) == 0 &&
                "tangent_directions must have cardinality #number of targets * #dimensions * #dimensions");
        auto num_targets = tangent_directions.size() / (_dimensions*_dimensions);
        _T = decltype(_T)("device tangent directions", num_targets*_dimensions*_dimensions);

        compadre_assert_release( (std::is_same<decltype(_T)::memory_space, typename view_type::memory_space>::value) &&
                "Memory space does not match between _T and tangent_directions");

        auto this_dimensions = _dimensions;
        auto this_T = _T;
        // rearrange data on device from data given on host
        Kokkos::parallel_for("copy tangent vectors", Kokkos::RangePolicy<device_execution_space>(0, num_targets), KOKKOS_LAMBDA(const int i) {
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
        _ref_N = decltype(_ref_N)("device normal directions", _pc._target_coordinates.extent(0)*_dimensions);
        // to assist LAMBDA capture
        auto this_ref_N = this->_ref_N;
        auto this_dimensions = this->_dimensions;

        // rearrange data on device from data given on host
        Kokkos::parallel_for("copy normal vectors", Kokkos::RangePolicy<device_execution_space>(0, _pc._target_coordinates.extent(0)), KOKKOS_LAMBDA(const int i) {
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

    //! Set dense solver type type
    void setDenseSolverType(const DenseSolverType dst) {
        _dense_solver_type = dst;
        this->resetCoefficientData();
    }

    //! Set constraint type
    void setConstraintType(const ConstraintType ct) {
        _constraint_type = ct;
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
        } else if (wt_to_lower == "cardinalcubicbspline") {
            _weighting_type = WeightingFunctionType::CardinalCubicBSpline;
        } else if (wt_to_lower == "cosine") {
            _weighting_type = WeightingFunctionType::Cosine;
        } else if (wt_to_lower == "sigmoid") {
            _weighting_type = WeightingFunctionType::Sigmoid;
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
        } else if (wt_to_lower == "cardinalcubicbspline") {
            _curvature_weighting_type = WeightingFunctionType::CardinalCubicBSpline;
        } else if (wt_to_lower == "cosine") {
            _curvature_weighting_type = WeightingFunctionType::Cosine;
        } else if (wt_to_lower == "sigmoid") {
            _curvature_weighting_type = WeightingFunctionType::Sigmoid;
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

    //! Sets basis order to be used when reconstructing any function
    void setPolynomialOrder(const int poly_order) {
        compadre_assert_release(poly_order<11 && "Unsupported polynomial order (>=11).");
        _poly_order = poly_order;
        _NP = this->getNP(_poly_order, _dimensions, _reconstruction_space);
        this->resetCoefficientData();
    }

    //! Sets basis order to be used when reconstruction curvature
    void setCurvaturePolynomialOrder(const int curvature_poly_order) {
        compadre_assert_release(curvature_poly_order<11 && "Unsupported curvature polynomial order (>=11).");
        _curvature_poly_order = curvature_poly_order;
        this->resetCoefficientData();
    }

    //! Parameter for weighting kernel for GMLS problem
    //! index = 0 sets p paramater for weighting kernel
    //! index = 1 sets n paramater for weighting kernel
    void setWeightingParameter(int wp, int index = 0) { 
        if (index==1) {
            _weighting_n = wp;
        } else {
            _weighting_p = wp;
        }
        this->resetCoefficientData();
    }

    //! Parameter for weighting kernel for curvature
    //! index = 0 sets p paramater for weighting kernel
    //! index = 1 sets n paramater for weighting kernel
    void setCurvatureWeightingParameter(int wp, int index = 0) { 
        if (index==1) {
            _curvature_weighting_n = wp;
        } else {
            _curvature_weighting_p = wp;
        }
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
        _h_ss.addTargets(lro);
        this->resetCoefficientData();
    }

    //! Adds a vector of target functionals to the vector of target functionals already to be applied to the reconstruction
    void addTargets(std::vector<TargetOperation> lro) {
        _h_ss.addTargets(lro);
        this->resetCoefficientData();
    }

    //! Empties the vector of target functionals to apply to the reconstruction
    void clearTargets() {
        _h_ss.clearTargets();
        this->resetCoefficientData();
    }

    //! Verify whether _pc is valid
    bool verifyPointConnections(bool assert_valid = false) {
        bool result = (_pc._target_coordinates.extent_int(0)==_pc._nla.getNumberOfTargets());
        compadre_assert_release((!assert_valid || result) &&
                "Target coordinates and neighbor lists have different size.");
        
        result &= (_pc._source_coordinates.extent(1)==_pc._target_coordinates.extent(1));
        compadre_assert_release((!assert_valid || result) &&
                "Source coordinates and target coordinates have different dimensions.");

        result &= (_pc._source_coordinates.extent(0)>0||_pc._target_coordinates.extent(0)==0);
        compadre_assert_release((!assert_valid || result) &&
                "Source coordinates not set in GMLS class before calling generatePolynomialCoefficients.");
        return result;
    }

    //! Verify whether _additional_pc is valid
    bool verifyAdditionalPointConnections(bool assert_valid = false) {
        bool result = (_additional_pc._target_coordinates.extent_int(0)==_additional_pc._nla.getNumberOfTargets());
        compadre_assert_release((!assert_valid || result) &&
                "Target coordinates and additional evaluation indices have different size.");
       
        result &= (_pc._target_coordinates.extent(0)==_additional_pc._target_coordinates.extent(0));
        compadre_assert_release((!assert_valid || result) &&
                "Target coordinates and additional evaluation indices have different size.");

        if (_additional_pc._source_coordinates.extent(0)>0) {
            result &= (_additional_pc._target_coordinates.extent(1)==_additional_pc._source_coordinates.extent(1));
            compadre_assert_release((!assert_valid || result) &&
                    "Target coordinates and additional evaluation coordinates have different dimensions.");
        }
        return result;
    }

    /*! \brief Generates polynomial coefficients by setting up and solving least squares problems
        Sets up the batch of GMLS problems to be solved for. Provides alpha values
        that can later be contracted against data or degrees of freedom to form a
        global linear system.

        \param number_of_batches    [in] - how many batches to break up the total workload into (for storage)
        \param keep_coefficients    [in] - whether to store (P^T W P)^-1 * P^T * W
        \param clear_cache          [in] - clear whatever temporary memory was used for calculations that  
                                           \a keep_coefficients hasn't indicated needs to be saved
    
        \a clear_cache should be \a false to be used for debugging as it will leave all data structures used to 
        compute the solution intact. If \a clear_cache is set \a true and \a keep_coefficients is \a true, it will 
        allow the polynomial coefficients to persist after calculation.
    */
    void generatePolynomialCoefficients(const int number_of_batches = 1, const bool keep_coefficients = false, 
            const bool clear_cache = true);

    /*! \brief Meant to calculate target operations and apply the evaluations to the previously 
        constructed polynomial coefficients. But now that is inside of generatePolynomialCoefficients because
        it must be to handle number_of_batches>1. Effectively, this just calls generatePolynomialCoefficients.

        \param number_of_batches    [in] - how many batches to break up the total workload into (for storage)
        \param keep_coefficients    [in] - whether to store (P^T W P)^-1 * P^T * W
        \param clear_cache          [in] - clear whatever temporary memory was used for calculations that  
                                           \a keep_coefficients hasn't indicated needs to be saved
    
        \a clear_cache should be \a false to be used for debugging as it will leave all data structures used to 
        compute the solution intact. If \a clear_cache is set \a true and \a keep_coefficients is \a true, it will 
        allow the polynomial coefficients to persist after calculation.
    */
    void generateAlphas(const int number_of_batches = 1, const bool keep_coefficients = false, 
            const bool clear_cache = true);

///@}


}; // GMLS Class
} // Compadre

#endif


