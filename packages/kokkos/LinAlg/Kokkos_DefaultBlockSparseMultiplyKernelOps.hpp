//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DEFAULTBLOCKSPARSEMULTIPLYKERNELOPS_HPP
#define KOKKOS_DEFAULTBLOCKSPARSEMULTIPLYKERNELOPS_HPP

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

#ifdef __CUDACC__
#include <Teuchos_ScalarTraitsCUDA.hpp>
#else
#include <Teuchos_ScalarTraits.hpp>
#endif

namespace Kokkos {

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class Ordinal,class DomainScalar,class RangeScalar>
void densematvec(Ordinal Nrows, Ordinal Ncols,
                 Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  unsigned offset = 0;
  for(Ordinal c=0; c<Ncols; ++c) {
    for(Ordinal r=0; r<Nrows; ++r) {
      y[r] += alpha*A[offset++]*x[c];
    }
  }
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
void dense_matvec_1x1( Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*A[0]*x[0];
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
void dense_matvec_2x2(Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[2]*x[1]);
  y[1] += alpha*(A[1]*x[0] + A[3]*x[1]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
void dense_matvec_3x3(Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[3]*x[1] + A[6]*x[2]);
  y[1] += alpha*(A[1]*x[0] + A[4]*x[1] + A[7]*x[2]);
  y[2] += alpha*(A[2]*x[0] + A[5]*x[1] + A[8]*x[2]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
void dense_matvec_4x4(Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[4]*x[1] + A[8]*x[2] + A[12]*x[3]);
  y[1] += alpha*(A[1]*x[0] + A[5]*x[1] + A[9]*x[2] + A[13]*x[3]);
  y[2] += alpha*(A[2]*x[0] + A[6]*x[1] + A[10]*x[2] + A[14]*x[3]);
  y[3] += alpha*(A[3]*x[0] + A[7]*x[1] + A[11]*x[2] + A[15]*x[3]);
}

template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
struct DefaultBlockSparseMultiplyOp1 {
  // mat data
  const Ordinal *rptr, *cptr, *bptr;
  const Ordinal *bindx, *indx;
  const Scalar  *vals;
  // matvec params
  RangeScalar        alpha, beta;
  size_t numBlockRows;
  // mv data
  const DomainScalar  *x;
  RangeScalar         *y;
  size_t xstride, ystride;
  size_t numVecs;

  inline KERNEL_PREFIX void execute(size_t i) {
    const size_t row = i;
    for(size_t v=0; v<numVecs; ++v) {
      const Ordinal Nrows = rptr[row+1]-rptr[row];
      const DomainScalar* xvec = x+v*xstride;
      RangeScalar* yvec = y+v*ystride;
      RangeScalar* yy = &yvec[rptr[row]];
  
      for(Ordinal i=0; i<Nrows; ++i) yy[i] = beta*yy[i];
  
      for (Ordinal b=bptr[row]; b<bptr[row+1]; ++b) {
        const Ordinal col = bindx[b];
        const Ordinal Ncols = cptr[col+1]-cptr[col];
  
        const Scalar* A = &vals[indx[b]];
        const Scalar* xx = &xvec[cptr[col]];
  
        if (Nrows == Ncols) {
          switch(Nrows) {
          case 1: dense_matvec_1x1(alpha, A, xx, yy); break;
          case 2: dense_matvec_2x2(alpha, A, xx, yy); break;
          case 3: dense_matvec_3x3(alpha, A, xx, yy); break;
          case 4: dense_matvec_4x4(alpha, A, xx, yy); break;
          default: densematvec(Nrows, Ncols, alpha, A, xx, yy);
          }
        }
        else {
          densematvec(Nrows,Ncols,alpha,A,xx,yy);
        }
      }
    }
  }
};

} // namespace Kokkos

#endif
