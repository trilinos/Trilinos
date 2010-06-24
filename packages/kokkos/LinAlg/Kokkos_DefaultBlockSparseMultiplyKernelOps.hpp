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


// FINISH TODO : This code is missing the calls to ScalarTraits::conjugate() 
//               necessary for the transposed kernels to work for complex-valued scalars.
//               Note, this will require also the addition of conjugate functionality to the 1x1 kernel

namespace Kokkos {

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class Ordinal,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
densematvec(Ordinal Nrows, Ordinal Ncols,
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

/** Form dense transpose-matrix-vector product y = A*x */
template<class Scalar,class Ordinal,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
densematvec_trans(Ordinal Nrows, Ordinal Ncols,
                  Scalar alpha, const Scalar* A,
                  const DomainScalar* x, RangeScalar* y)
{
  unsigned offset = 0;
  for(Ordinal c=0; c<Ncols; ++c) {
    for(Ordinal r=0; r<Nrows; ++r) {
      y[c] += alpha*A[offset++]*x[r];
    }
  }
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_1x1( Scalar alpha, const Scalar* A,
                  const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*A[0]*x[0];
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_2x2(Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[2]*x[1]);
  y[1] += alpha*(A[1]*x[0] + A[3]*x[1]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_3x3(Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[3]*x[1] + A[6]*x[2]);
  y[1] += alpha*(A[1]*x[0] + A[4]*x[1] + A[7]*x[2]);
  y[2] += alpha*(A[2]*x[0] + A[5]*x[1] + A[8]*x[2]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_4x4(Scalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[4]*x[1] + A[8]*x[2] + A[12]*x[3]);
  y[1] += alpha*(A[1]*x[0] + A[5]*x[1] + A[9]*x[2] + A[13]*x[3]);
  y[2] += alpha*(A[2]*x[0] + A[6]*x[1] + A[10]*x[2] + A[14]*x[3]);
  y[3] += alpha*(A[3]*x[0] + A[7]*x[1] + A[11]*x[2] + A[15]*x[3]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_trans_2x2(Scalar alpha, const Scalar* A,
                       const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[1]*x[1]);
  y[1] += alpha*(A[2]*x[0] + A[3]*x[1]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_trans_3x3(Scalar alpha, const Scalar* A,
                       const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[1]*x[1] + A[2]*x[2]);
  y[1] += alpha*(A[3]*x[0] + A[4]*x[1] + A[5]*x[2]);
  y[2] += alpha*(A[6]*x[0] + A[7]*x[1] + A[8]*x[2]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_trans_4x4(Scalar alpha, const Scalar* A,
                       const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[1]*x[1] + A[2]*x[2] + A[3]*x[3]);
  y[1] += alpha*(A[4]*x[0] + A[5]*x[1] + A[6]*x[2] + A[7]*x[3]);
  y[2] += alpha*(A[8]*x[0] + A[9]*x[1] + A[10]*x[2] + A[11]*x[3]);
  y[3] += alpha*(A[12]*x[0] + A[13]*x[1] + A[14]*x[2] + A[15]*x[3]);
}

template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
struct DefaultBlockSparseMultiplyOp1 {
  // mat data
  const Ordinal *rptr, *cptr, *bptr;
  const Ordinal *bindx, *indx;
  const Scalar  *vals;
  // matvec params
  Scalar        alpha, beta;
  size_t numBlockRows;
  // mv data
  const DomainScalar  *x;
  RangeScalar         *y;
  size_t xstride, ystride;

  inline KERNEL_PREFIX void execute(size_t myIter) {
    const size_t myBlockRow = myIter % numBlockRows;
    const size_t myRHS = (myIter - myBlockRow) / numBlockRows;
    const Ordinal numRowsInBlock = rptr[myBlockRow+1]-rptr[myBlockRow];
    const DomainScalar* xvec = x+myRHS*xstride;
    RangeScalar*        yvec = y+myRHS*ystride;
    RangeScalar* yy = &yvec[rptr[myBlockRow]];
    // init my block of y to zero
    for(Ordinal i=0; i<numRowsInBlock; ++i) yy[i] = beta*yy[i];
    // loop over the block in my row and do the multiplication
    for (Ordinal b=bptr[myBlockRow]; b<bptr[myBlockRow+1]; ++b) {
      // get pointers into A and x
      const Ordinal col = bindx[b];
      const Ordinal numColsInBlock = cptr[col+1]-cptr[col];
      const Scalar* A = &vals[indx[b]];
      const DomainScalar* xx = &xvec[cptr[col]];
      // do the GEMM
      if (numRowsInBlock == numColsInBlock) {
        switch(numRowsInBlock) {
          case 1: dense_matvec_1x1(alpha, A, xx, yy); break;
          case 2: dense_matvec_2x2(alpha, A, xx, yy); break;
          case 3: dense_matvec_3x3(alpha, A, xx, yy); break;
          case 4: dense_matvec_4x4(alpha, A, xx, yy); break;
          default: densematvec(numRowsInBlock, numColsInBlock, alpha, A, xx, yy);
        }
      }
      else {
        densematvec(numRowsInBlock,numColsInBlock,alpha,A,xx,yy);
      }
    }
  }
};

template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
struct DefaultBlockSparseMultiplyOp1Transpose {
  // mat data
  const Ordinal *rptr, *cptr, *bptr;
  const Ordinal *bindx, *indx;
  const Scalar  *vals;
  // matvec params
  Scalar        alpha;
  size_t numBlockRows;
  // mv data
  const DomainScalar  *x;
  RangeScalar         *y;
  size_t xstride, ystride;

  inline KERNEL_PREFIX void execute(const size_t myRHS) {
    // get pointers into X and Y for my assigned RHS
    const DomainScalar* xvec = x+myRHS*xstride;
    RangeScalar*        yvec = y+myRHS*ystride;
    for (size_t bR=0; bR<numBlockRows; ++bR) {
      // accumulate bR-th block bR (transposed) times xvec[bR] block entry into yvec
      const DomainScalar* xx = &xvec[rptr[bR]];
      const Ordinal Nrows = rptr[bR+1]-rptr[bR];
  
      for (Ordinal b=bptr[bR]; b<bptr[bR+1]; ++b) {
        const Ordinal col = bindx[b];
        const Ordinal Ncols = cptr[col+1]-cptr[col];
  
        const Scalar* A = &vals[indx[b]];
        RangeScalar* yy = &yvec[cptr[col]];
  
        if (Nrows == Ncols) {
          switch(Nrows) {
          case 1: dense_matvec_1x1(alpha, A, xx, yy); break;
          case 2: dense_matvec_trans_2x2(alpha, A, xx, yy); break;
          case 3: dense_matvec_trans_3x3(alpha, A, xx, yy); break;
          case 4: dense_matvec_trans_4x4(alpha, A, xx, yy); break;
          default: densematvec_trans(Nrows, Ncols, alpha, A, xx, yy);
          }
        }
        else {
          densematvec_trans(Nrows,Ncols,alpha,A,xx,yy);
        }
      }
    }
  }
};

} // namespace Kokkos

#endif
