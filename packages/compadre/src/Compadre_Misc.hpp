// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_MISC_HPP_
#define _COMPADRE_MISC_HPP_

#include "Compadre_Operators.hpp"

namespace Compadre {

struct XYZ {

    double x;
    double y;
    double z;

    KOKKOS_INLINE_FUNCTION
    XYZ(double _x = 0.0, double _y = 0.0, double _z = 0.0) : x(_x), y(_y), z(_z) {}

    KOKKOS_INLINE_FUNCTION
    XYZ(const scalar_type* arry) : x(arry[0]), y(arry[1]), z(arry[2]) {};
    
    XYZ(const std::vector<scalar_type>& vec) : x(vec[0]), y(vec[1]), z(vec[2]) {};
    
    KOKKOS_INLINE_FUNCTION
    XYZ operator += (const XYZ& other)
        { x += other.x;
          y += other.y;
          z += other.z;
          return *this;    }
    KOKKOS_INLINE_FUNCTION
    XYZ operator += (XYZ& other)
        { x += other.x;
          y += other.y;
          z += other.z;
          return *this;    }
    KOKKOS_INLINE_FUNCTION
    XYZ operator -= (const XYZ& other)
        { x -= other.x;
          y -= other.y;
          z -= other.z;
          return *this;    }
    KOKKOS_INLINE_FUNCTION
    XYZ operator -= (XYZ& other)
        { x -= other.x;
          y -= other.y;
          z -= other.z;
          return *this;    }
    KOKKOS_INLINE_FUNCTION
    XYZ operator *= (const double& scaling)
        { x *= scaling;
          y *= scaling;
          z *= scaling;
          return *this;    }
    KOKKOS_INLINE_FUNCTION
    scalar_type& operator [](const int i) {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                return z;
        }
    }
    KOKKOS_INLINE_FUNCTION
    scalar_type operator [](const int i) const {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                return z;
        }
    }
    KOKKOS_INLINE_FUNCTION
    XYZ operator *(double scalar) {
        XYZ result;
        result.x = scalar*x;
        result.y = scalar*y;
        result.z = scalar*z;
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    size_t extent(int comp = 0) {
        return 3;
    }

    KOKKOS_INLINE_FUNCTION
    int extent_int(int comp = 0) {
        return 3;
    }

}; // XYZ


KOKKOS_INLINE_FUNCTION
XYZ operator + ( const XYZ& vecA, const XYZ& vecB ) {
    return XYZ( vecA.x + vecB.x, vecA.y + vecB.y, vecA.z + vecB.z); }

    KOKKOS_INLINE_FUNCTION
XYZ operator - ( const XYZ& vecA, const XYZ& vecB ) {
    return XYZ( vecA.x - vecB.x, vecA.y - vecB.y, vecA.z - vecB.z); }
    
KOKKOS_INLINE_FUNCTION
XYZ operator * ( const XYZ& vecA, const XYZ& vecB ) {
    return XYZ( vecA.x * vecB.x, vecA.y * vecB.y, vecA.z * vecB.z); }

KOKKOS_INLINE_FUNCTION
XYZ operator + ( const XYZ& vecA, const scalar_type& constant ) {
    return XYZ( vecA.x + constant, vecA.y + constant, vecA.z + constant); }

KOKKOS_INLINE_FUNCTION
XYZ operator + ( const scalar_type& constant, const XYZ& vecA ) {
    return XYZ( vecA.x + constant, vecA.y + constant, vecA.z + constant); }

KOKKOS_INLINE_FUNCTION
XYZ operator - ( const XYZ& vecA, const scalar_type& constant ) {
    return XYZ( vecA.x - constant, vecA.y - constant, vecA.z - constant); }

KOKKOS_INLINE_FUNCTION
XYZ operator - ( const scalar_type& constant,  const XYZ& vecA ) {
    return XYZ( constant - vecA.x, constant - vecA.y, constant - vecA.z); }

KOKKOS_INLINE_FUNCTION
XYZ operator * ( const XYZ& vecA, const scalar_type& constant ) {
    return XYZ( vecA.x * constant, vecA.y * constant, vecA.z * constant); }

KOKKOS_INLINE_FUNCTION
XYZ operator * (const scalar_type& constant, const XYZ& vecA) {
    return XYZ( vecA.x * constant, vecA.y * constant, vecA.z * constant); }

KOKKOS_INLINE_FUNCTION
XYZ operator / ( const XYZ& vecA, const scalar_type& constant ) {
    return XYZ( vecA.x / constant, vecA.y / constant, vecA.z / constant); }

inline std::ostream& operator << ( std::ostream& os, const XYZ& vec ) {
    os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" ; return os; }

//! n^p (n^0 returns 1, regardless of n)
KOKKOS_INLINE_FUNCTION
int pown(int n, unsigned p) {
    // O(p) implementation
    int y = 1;
    for (unsigned i=0; i<p; i++) {
        y *= n;
    }
    return y;
    // O(log(p)) implementation
    //int result = 1;
    //while (p) {
    //    if (p & 0x1) {
    //        result *= n;
    //    }
    //    n *= n;
    //    p >>= 1;
    //}
    //return result;
}

KOKKOS_INLINE_FUNCTION
int getAdditionalAlphaSizeFromConstraint(DenseSolverType dense_solver_type, ConstraintType constraint_type) {
    // Return the additional constraint size
    if (constraint_type == NEUMANN_GRAD_SCALAR) {
        return 1;
    } else {
        return 0;
    }
}

KOKKOS_INLINE_FUNCTION
int getAdditionalCoeffSizeFromConstraintAndSpace(DenseSolverType dense_solver_type, ConstraintType constraint_type, ReconstructionSpace reconstruction_space, const int dimension) {
    // Return the additional constraint size
    int added_alpha_size = getAdditionalAlphaSizeFromConstraint(dense_solver_type, constraint_type);
    if (reconstruction_space == VectorTaylorPolynomial) {
        return dimension*added_alpha_size;
    } else {
        return added_alpha_size;
    }
}

KOKKOS_INLINE_FUNCTION
void getRHSDims(DenseSolverType dense_solver_type, ConstraintType constraint_type, ReconstructionSpace reconstruction_space, const int dimension, const int M, const int N, int &RHS_row, int &RHS_col) {
    // Return the appropriate size for _RHS. Since in LU, the system solves P^T*P against P^T*W.
    // We store P^T*P in the RHS space, which means RHS can be much smaller compared to the
    // case for QR/SVD where the system solves PsqrtW against sqrtW*Identity

    int added_coeff_size = getAdditionalCoeffSizeFromConstraintAndSpace(dense_solver_type, constraint_type, reconstruction_space, dimension);

    if (constraint_type == NEUMANN_GRAD_SCALAR) {
        RHS_row = RHS_col = N + added_coeff_size;
    } else {
        if (dense_solver_type != LU) {
            RHS_row = N;
            RHS_col = M;
        } else {
            RHS_row = RHS_col = N + added_coeff_size;
        }
    }
}

KOKKOS_INLINE_FUNCTION
void getPDims(DenseSolverType dense_solver_type, ConstraintType constraint_type, ReconstructionSpace reconstruction_space, const int dimension, const int M, const int N, int &out_row, int &out_col) {
    // Return the appropriate size for _P.
    // In the case of solving with LU and additional constraint is used, _P needs
    // to be resized to include additional row(s) based on the type of constraint.

    int added_coeff_size = getAdditionalCoeffSizeFromConstraintAndSpace(dense_solver_type, constraint_type, reconstruction_space, dimension);
    int added_alpha_size = getAdditionalAlphaSizeFromConstraint(dense_solver_type, constraint_type);

    if (constraint_type == NEUMANN_GRAD_SCALAR) {
        out_row = M + added_alpha_size;
        out_col = N + added_coeff_size;
    } else {
        if (dense_solver_type == LU) {
            out_row = M + added_alpha_size;
            out_col = N + added_coeff_size;
        } else {
            out_row = M;
            out_col = N;
        }
    }
}

//! Helper function for finding alpha coefficients
KOKKOS_INLINE_FUNCTION
int getTargetOutputIndex(const int operation_num, const int output_component_axis_1, const int output_component_axis_2, const int dimensions) {
    const int axis_1_size = (getTargetOutputTensorRank(operation_num) > 1) ? dimensions : 1;
    return axis_1_size*output_component_axis_1 + output_component_axis_2; // 0 for scalar, 0 for vector;
}

//! Helper function for finding alpha coefficients
KOKKOS_INLINE_FUNCTION
int getSamplingOutputIndex(const SamplingFunctional sf, const int output_component_axis_1, const int output_component_axis_2) {
    const int axis_1_size = (sf.output_rank > 1) ? sf.output_rank : 1;
    return axis_1_size*output_component_axis_1 + output_component_axis_2; // 0 for scalar, 0 for vector;
}

//! Input rank for sampling operation
KOKKOS_INLINE_FUNCTION
int getInputRankOfSampling(SamplingFunctional sro) {
    return sro.input_rank;
}

//! Dimensions ^ output rank for sampling operation 
//! (always in local chart if on a manifold, never ambient space)
KOKKOS_INLINE_FUNCTION
int getOutputDimensionOfSampling(SamplingFunctional sro, const int local_dimensions) {
    return pown(local_dimensions, sro.output_rank);
}

//! Dimensions ^ output rank for sampling operation 
//! (always in ambient space, never local chart on a manifold)
KOKKOS_INLINE_FUNCTION
int getInputDimensionOfSampling(SamplingFunctional sro, const int global_dimensions) {
    return pown(global_dimensions, sro.input_rank);
}

//! Calculate basis_multiplier
KOKKOS_INLINE_FUNCTION
int calculateBasisMultiplier(const ReconstructionSpace rs, const int local_dimensions) {
    // calculate the dimension of the basis 
    // (a vector space on a manifold requires two components, for example)
    return pown(local_dimensions, getActualReconstructionSpaceRank((int)rs));
}

//! Calculate sampling_multiplier
KOKKOS_INLINE_FUNCTION
int calculateSamplingMultiplier(const ReconstructionSpace rs, const SamplingFunctional sro, const int local_dimensions) {
    // this would normally be SamplingOutputTensorRank[_data_sampling_functional], but we also want to handle the
    // case of reconstructions where a scalar basis is reused as a vector, and this handles everything
    // this handles scalars, vectors, and scalars that are reused as vectors
    int bm = calculateBasisMultiplier(rs, local_dimensions);
    int sm = getOutputDimensionOfSampling(sro, local_dimensions);
    if (rs == ReconstructionSpace::VectorOfScalarClonesTaylorPolynomial) {
        // storage and computational efficiency by reusing solution to scalar problem for 
        // a vector problem (in 3d, 27x cheaper computation, 9x cheaper storage)
        sm =(bm < sm) ? bm : sm;
    }
    return sm;
}

//! Dimensions ^ output rank for target operation
KOKKOS_INLINE_FUNCTION
int getOutputDimensionOfOperation(TargetOperation lro, const int local_dimensions) {
    return pown(local_dimensions, getTargetOutputTensorRank(lro));
}

//! Dimensions ^ input rank for target operation (always in local chart if on a manifold, never ambient space)
KOKKOS_INLINE_FUNCTION
int getInputDimensionOfOperation(TargetOperation lro, SamplingFunctional sro, const int local_dimensions) {
    // this is the same return values as the OutputDimensionOfSampling for the GMLS class's SamplingFunctional
    return getOutputDimensionOfSampling(sro, local_dimensions);
}

} // Compadre

#endif
