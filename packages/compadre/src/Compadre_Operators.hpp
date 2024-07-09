// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_OPERATORS_HPP_
#define _COMPADRE_OPERATORS_HPP_

#include "Compadre_Typedefs.hpp"

#define make_sampling_functional(input, output, targets, nontrivial, transform) SamplingFunctional(input, output, targets, nontrivial, transform, __COUNTER__)

namespace Compadre {

    //! Available target functionals
    enum TargetOperation : int {
        //! Point evaluation of a scalar
        ScalarPointEvaluation,
        //! Point evaluation of a vector (reconstructs entire vector at once, requiring a 
        //! ReconstructionSpace having a sufficient number of components in the basis)
        VectorPointEvaluation, 
        //! Point evaluation of the laplacian of a scalar (could be on a manifold or not)
        LaplacianOfScalarPointEvaluation,
        //! Point evaluation of the laplacian of each component of a vector
        VectorLaplacianPointEvaluation,
        //! Point evaluation of the gradient of a scalar
        GradientOfScalarPointEvaluation,
        //! Point evaluation of the gradient of a vector (results in a matrix, NOT CURRENTLY IMPLEMENTED)
        GradientOfVectorPointEvaluation,
        //! Point evaluation of the divergence of a vector (results in a scalar)
        DivergenceOfVectorPointEvaluation,
        //! Point evaluation of the curl of a vector (results in a vector)
        CurlOfVectorPointEvaluation,
        //! Point evaluation of the curl of a curl of a vector (results in a vector)
        CurlCurlOfVectorPointEvaluation,
        //! Point evaluation of the partial with respect to x of a scalar
        PartialXOfScalarPointEvaluation,
        //! Point evaluation of the partial with respect to y of a scalar
        PartialYOfScalarPointEvaluation,
        //! Point evaluation of the partial with respect to z of a scalar
        PartialZOfScalarPointEvaluation,
        //! Point evaluation of the chained staggered Laplacian acting on VectorTaylorPolynomial 
        //! basis + StaggeredEdgeIntegralSample sampling functional
        ChainedStaggeredLaplacianOfScalarPointEvaluation,
        //! Point evaluation of Gaussian curvature
        GaussianCurvaturePointEvaluation,
        //! Average values of a cell using quadrature
        //! Supported on 2D faces in 3D problems (manifold) and 2D cells in 2D problems
        CellAverageEvaluation,
        //! Integral values over cell using quadrature
        //! Supported on 2D faces in 3D problems (manifold) and 2D cells in 2D problems
        CellIntegralEvaluation,
        //! Average value of vector dotted with normal direction
        //! Supported on 1D edges in 3D problems (2D-manifold) and 1D edges on 2D cells
        FaceNormalAverageEvaluation,
        //! Integral value of vector dotted with normal direction
        //! Supported on 1D edges in 3D problems (2D-manifold) and 1D edges on 2D cells
        FaceNormalIntegralEvaluation,
        //! Average value of vector dotted with tangent directions
        //! Supported on 1D edges in 3D problems (2D-manifold) and 1D edges on 2D cells
        EdgeTangentAverageEvaluation,
        //! Integral value of vector dotted with tangent directions
        //! Supported on 1D edges in 3D problems (2D-manifold) and 1D edges on 2D cells
        EdgeTangentIntegralEvaluation,
        //! Should be the total count of all available target functionals
        COUNT=20,
    };

    //! Rank of target functional output for each TargetOperation 
    //! Rank of target functional input for each TargetOperation is based on the output
    //! rank of the SamplingFunctional used on the polynomial basis
    KOKKOS_INLINE_FUNCTION
    int getTargetOutputTensorRank(const int& index) {
        constexpr int TargetOutputTensorRank[] {
            0, ///< PointEvaluation
            1, ///< VectorPointEvaluation
            0, ///< LaplacianOfScalarPointEvaluation
            1, ///< VectorLaplacianPointEvaluation
            1, ///< GradientOfScalarPointEvaluation
            2, ///< GradientOfVectorPointEvaluation
            0, ///< DivergenceOfVectorPointEvaluation
            1, ///< CurlOfVectorPointEvaluation
            1, ///< CurlCurlOfVectorPointEvaluation
            0, ///< PartialXOfScalarPointEvaluation
            0, ///< PartialYOfScalarPointEvaluation
            0, ///< PartialZOfScalarPointEvaluation
            0, ///< ChainedStaggeredLaplacianOfScalarPointEvaluation
            0, ///< GaussianCurvaturePointEvaluation
            0, ///< CellAverageEvaluation
            0, ///< CellIntegralEvaluation
            0, ///< FaceNormalAverageEvaluation
            0, ///< FaceNormalIntegralEvaluation
            0, ///< EdgeTangentAverageEvaluation
            0, ///< EdgeTangentIntegralEvaluation
        };
        return TargetOutputTensorRank[index];
    }

    //! Space in which to reconstruct polynomial
    enum ReconstructionSpace : int {
        //! Scalar polynomial basis centered at the target site and scaled by sum of basis powers 
        //! e.g. \f$(x-x_t)^2*(y-y_t)*(z-z_t)^3/factorial(2+1+3)\f$ would be a member of 3rd order in 3D, where 
        //! \f$(x_t,y_t,z_t)\f$ is the coordinate of the target site in 3D coordinates.
        ScalarTaylorPolynomial,
        //! Vector polynomial basis having # of components _dimensions, or (_dimensions-1) in the case of manifolds)
        VectorTaylorPolynomial,
        //! Scalar basis reused as many times as there are components in the vector
        //! resulting in a much cheaper polynomial reconstruction
        VectorOfScalarClonesTaylorPolynomial,
        //! Divergence-free vector polynomial basis
        DivergenceFreeVectorTaylorPolynomial,
        //! Bernstein polynomial basis
        BernsteinPolynomial,
    };

    //! Number of actual components in the ReconstructionSpace
    KOKKOS_INLINE_FUNCTION
    int getActualReconstructionSpaceRank(const int& index) {
        constexpr int ActualReconstructionSpaceRank[] = {
            0, ///< ScalarTaylorPolynomial
            1, ///< VectorTaylorPolynomial
            0, ///< VectorOfScalarClonesTaylorPolynomial
            0, ///< DivergenceFreeVectorTaylorPolynomial
            0, ///< BernsteinPolynomial
        };
        return ActualReconstructionSpaceRank[index];
    }

    //! Describes the SamplingFunction relationship to targets, neighbors
    enum SamplingTransformType : int {
        Identity,           ///< No action performed on data before GMLS target operation
        SameForAll,         ///< Each neighbor for each target all apply the same transform
        DifferentEachTarget,   ///< Each target applies a different data transform, but the same to each neighbor
        DifferentEachNeighbor, ///< Each target applies a different transform for each neighbor
    };

    struct SamplingFunctional {
        //! for uniqueness
        size_t id;
        //! Rank of sampling functional input for each SamplingFunctional
        int input_rank;
        //! Rank of sampling functional output for each SamplingFunctional
        int output_rank;
        //! Whether or not the SamplingTensor acts on the target site as well as the neighbors.
        //! This makes sense only in staggered schemes, when each target site is also a source site
        bool use_target_site_weights;
        //! Whether the SamplingFunctional + ReconstructionSpace results in a nontrivial nullspace
        bool nontrivial_nullspace;
        //! Describes the SamplingFunction relationship to targets, neighbors
        int transform_type;

        //SamplingFunctional(std::string name, const int input_rank_, const int output_rank_,
        KOKKOS_INLINE_FUNCTION
        constexpr SamplingFunctional(const int input_rank_, const int output_rank_,
                const bool use_target_site_weights_, const bool nontrivial_nullspace_,
                const int transform_type_, const int id_) : 
                id(id_), input_rank(input_rank_), output_rank(output_rank_),
                use_target_site_weights(use_target_site_weights_), nontrivial_nullspace(nontrivial_nullspace_),
                transform_type(transform_type_) {}

        KOKKOS_INLINE_FUNCTION
        constexpr bool operator == (const SamplingFunctional sf) const {
            return id == sf.id;
        }

        KOKKOS_INLINE_FUNCTION
        constexpr bool operator != (const SamplingFunctional sf) const {
            return id != sf.id;
        }

    };

    //! Available sampling functionals
    constexpr SamplingFunctional 

        //! Point evaluations of the scalar source function
        PointSample = make_sampling_functional(0,0,false,false,(int)Identity),

        //! Point evaluations of the entire vector source function
        VectorPointSample = make_sampling_functional(1,1,false,false,(int)Identity),

        //! Point evaluations of the entire vector source function 
        //! (but on a manifold, so it includes a transform into local coordinates)
        ManifoldVectorPointSample = make_sampling_functional(1,1,false,false,(int)DifferentEachTarget),

        //! Analytical integral of a gradient source vector is just a difference of the scalar source at neighbor and target
        StaggeredEdgeAnalyticGradientIntegralSample = make_sampling_functional(0,0,true,true,(int)SameForAll),

        //! Samples consist of the result of integrals of a vector dotted with the tangent along edges between neighbor and target
        StaggeredEdgeIntegralSample = make_sampling_functional(1,0,true,true,(int)DifferentEachNeighbor),

        //! For integrating polynomial dotted with normal over an edge
        VaryingManifoldVectorPointSample = make_sampling_functional(1,1,false,false,(int)DifferentEachNeighbor),

        //! For integrating polynomial dotted with normal over an edge
        FaceNormalIntegralSample = make_sampling_functional(1,0,false,false,(int)Identity),

        //! For polynomial dotted with normal on edge
        FaceNormalAverageSample = make_sampling_functional(1,0,false,false,(int)Identity),

        //! For integrating polynomial dotted with tangent over an edge
        EdgeTangentIntegralSample = make_sampling_functional(1,0,false,false,(int)Identity),

        //! For polynomial dotted with tangent
        EdgeTangentAverageSample = make_sampling_functional(1,0,false,false,(int)Identity),

        //! For polynomial integrated on cells
        CellAverageSample = make_sampling_functional(0,0,false,false,(int)DifferentEachNeighbor),

        //! For polynomial integrated on cells
        CellIntegralSample = make_sampling_functional(0,0,false,false,(int)DifferentEachNeighbor);

    //! Dense solver type
    enum DenseSolverType : int {
        //! QR+Pivoting factorization performed on P*sqrt(w) matrix
        QR, 
        //! LU factorization performed on P^T*W*P matrix
        LU, 
    };

    //! Problem type, that optionally can handle manifolds
    enum ProblemType : int {
        //! Standard GMLS problem type
        STANDARD, 
        //! Solve GMLS problem on a manifold (will use QR or SVD to solve the resultant GMLS 
        //! problem dependent on SamplingNontrivialNullspace
        MANIFOLD, 
    };

    //! Constraint type
    enum ConstraintType : int {
        //! No constraint
        NO_CONSTRAINT,
        //! Neumann Gradient Scalar Type
        NEUMANN_GRAD_SCALAR,
    };

    //! Available weighting kernel function types
    enum WeightingFunctionType : int {
        Power,
        Gaussian,
        CubicSpline,
        CardinalCubicBSpline,
        Cosine,
        Sigmoid
    };

    //! Coordinate type for input and output format of vector data on manifold problems.
    //! Anything without a manifold is always Ambient.
    enum CoordinatesType : int {
        Ambient, ///< a 2D manifold in 3D in ambient coordinates would have 3 components for a vector
        Local,   ///< a 2D manifold in 3D in local coordinates would have 2 components for a vector
    };

    /*
       \page pageGMLSOperators GMLS Target Operations, Spaces, and Sampling Functionals
       \section sectionGMLSTarget Target Operations

       \section sectionGMLSSpace Reconstruction Spaces

       \section sectionGMLSSampling Sampling Functionals 
    */


} // Compadre

#endif
