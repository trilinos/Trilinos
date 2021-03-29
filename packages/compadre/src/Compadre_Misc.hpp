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

} // Compadre

#endif
