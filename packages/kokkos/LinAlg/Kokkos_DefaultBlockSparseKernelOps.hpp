//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DEFAULTBLOCKSPARSEKERNELOPS_HPP
#define KOKKOS_DEFAULTBLOCKSPARSEKERNELOPS_HPP

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
            RangeScalar alpha, const Scalar* A,
            const DomainScalar* x, RangeScalar* y)
{
  unsigned offset = 0;
  for(Ordinal c=0; c<Ncols; ++c) {
    for(Ordinal r=0; r<Nrows; ++r) {
      y[r] += alpha*A[offset++]*x[c];
    }
  }
}

/** Form dense matrix-vector product y = A*x but skip the diagonal of A.
(Used by the solve operation.) */
template<class Scalar,class Ordinal,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
densematvec_skipdiag(Ordinal N,
                     RangeScalar alpha, const Scalar* A,
                     const DomainScalar* x, RangeScalar* y)
{
  unsigned offset = 0;
  for(Ordinal c=0; c<N; ++c) {
    for(Ordinal r=0; r<N; ++r) {
      if (r != c) {
        y[r] += alpha*A[offset]*x[c];
      }
      ++offset;
    }
  }
}

/** Solve square dense upper-triangular matrix.
    On entry y is the rhs, on exit y is the solution.
*/
template<class Scalar,class Ordinal,class RangeScalar>
inline KERNEL_PREFIX void 
solve_block_upper(Ordinal N, const Scalar* A, RangeScalar* y)
{
  int n = N;
  for(int r=n-1; r>=0; --r) {
    for(int c=n-1; c>r; --c) {
      y[r] -= A[c*n+r]*y[c];
    }
    y[r] /= A[r*n+r];
  }
}

/** Transpose-Solve square dense upper-triangular matrix.
    On entry y is the rhs, on exit y is the solution.
*/
template<class Scalar,class Ordinal,class RangeScalar>
inline KERNEL_PREFIX void 
solve_block_upper_transpose(Ordinal N, const Scalar* A, RangeScalar* y)
{
  int n = N;
  for(int c=0; c<N; ++c) {
    for(int r=0; r<c; ++r) {
      y[c] -= A[c*n+r]*y[r];
    }
    y[c] /= A[c*n+c];
  }
}

/** Solve square dense lower-triangular matrix.
    On entry y is the rhs, on exit y is the solution.
*/
template<class Scalar,class Ordinal,class RangeScalar>
inline KERNEL_PREFIX void 
solve_block_lower(Ordinal N, const Scalar* A, RangeScalar* y)
{
  int n = N;
  for(int r=0; r<n; ++r) {
    for(int c=0; c<r; ++c) {
      y[r] -= A[c*n+r]*y[c];
    }
    y[r] /= A[r*n+r];
  }
}

/** Transpose-Solve square dense lower-triangular matrix.
    On entry y is the rhs, on exit y is the solution.
*/
template<class Scalar,class Ordinal,class RangeScalar>
inline KERNEL_PREFIX void 
solve_block_lower_transpose(Ordinal N, const Scalar* A, RangeScalar* y)
{
  int n = N;
  for(int c=n-1; c>=0; --c) {
    for(int r=n-1; r>c; --r) {
      y[c] -= A[c*n+r]*y[r];
    }
    y[c] /= A[c*n+c];
  }
}

/** Form dense transpose-matrix-vector product y = A*x */
template<class Scalar,class Ordinal,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void
densematvec_trans(Ordinal Nrows, Ordinal Ncols,
                  RangeScalar alpha, const Scalar* A,
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
dense_matvec_1x1( RangeScalar alpha, const Scalar* A,
                  const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*A[0]*x[0];
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_2x2(RangeScalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[2]*x[1]);
  y[1] += alpha*(A[1]*x[0] + A[3]*x[1]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_3x3(RangeScalar alpha, const Scalar* A,
                 const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[3]*x[1] + A[6]*x[2]);
  y[1] += alpha*(A[1]*x[0] + A[4]*x[1] + A[7]*x[2]);
  y[2] += alpha*(A[2]*x[0] + A[5]*x[1] + A[8]*x[2]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_4x4(RangeScalar alpha, const Scalar* A,
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
dense_matvec_trans_2x2(RangeScalar alpha, const Scalar* A,
                       const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[1]*x[1]);
  y[1] += alpha*(A[2]*x[0] + A[3]*x[1]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_trans_3x3(RangeScalar alpha, const Scalar* A,
                       const DomainScalar* x, RangeScalar* y)
{
  y[0] += alpha*(A[0]*x[0] + A[1]*x[1] + A[2]*x[2]);
  y[1] += alpha*(A[3]*x[0] + A[4]*x[1] + A[5]*x[2]);
  y[2] += alpha*(A[6]*x[0] + A[7]*x[1] + A[8]*x[2]);
}

/** Form dense matrix-vector product y = A*x */
template<class Scalar,class DomainScalar,class RangeScalar>
inline KERNEL_PREFIX void 
dense_matvec_trans_4x4(RangeScalar alpha, const Scalar* A,
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
  const Ordinal *rptr, *cptr;
  const size_t *bptr;
  const Ordinal *bindx, *indx;
  const Scalar  *vals;
  // matvec params
  RangeScalar        alpha, beta;
  size_t numBlockRows;
  // mv data
  const DomainScalar  *x;
  RangeScalar         *y;
  size_t xstride, ystride;

  inline KERNEL_PREFIX void execute(size_t myIter) {
    RangeScalar zero = Teuchos::ScalarTraits<RangeScalar>::zero();
    const size_t myBlockRow = myIter % numBlockRows;
    const size_t myRHS = (myIter - myBlockRow) / numBlockRows;
    const Ordinal numRowsInBlock = rptr[myBlockRow+1]-rptr[myBlockRow];
    const DomainScalar* xvec = x+myRHS*xstride;
    RangeScalar*        yvec = y+myRHS*ystride;
    RangeScalar* yy = &yvec[rptr[myBlockRow]];

    // init my block of y to beta*y
    if (beta == zero) for(Ordinal i=0; i<numRowsInBlock; ++i) yy[i] = zero;
    else for(Ordinal i=0; i<numRowsInBlock; ++i) yy[i] = beta*yy[i];

    // loop over the block in my row and do the multiplication
    for (size_t b=bptr[myBlockRow]; b<bptr[myBlockRow+1]; ++b) {
      // get pointers into A and x
      const Ordinal col = bindx[b];
      const Ordinal numColsInBlock = cptr[col+1]-cptr[col];
      const Scalar* A = &vals[indx[b]];
      const DomainScalar* xx = &xvec[cptr[col]];
      // do the GEMV
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
  const Ordinal *rptr, *cptr;
  const size_t *bptr;
  const Ordinal *bindx, *indx;
  const Scalar  *vals;
  // matvec params
  RangeScalar        alpha;
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
  
      for (size_t b=bptr[bR]; b<bptr[bR+1]; ++b) {
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

template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
struct DefaultBlockSparseSolveOp1 {
  // mat data
  const Ordinal *rptr, *cptr;
  const size_t *bptr;
  const Ordinal *bindx, *indx;
  const Scalar  *vals;
  // solve params
  bool unitDiag, upper;
  size_t numBlockRows;
  // multivector data
  const RangeScalar  *y;
  DomainScalar       *x;
  size_t xstride, ystride;

  //find x such that A*x = y

  inline KERNEL_PREFIX void execute(size_t i) {
    //solve for i-th vector of multivector
    const RangeScalar* yvec = y+i*ystride;
    DomainScalar*      xvec = x+i*xstride;
    RangeScalar one = Teuchos::ScalarTraits<RangeScalar>::one();

    //copy y into x:
    Ordinal numPointRows = rptr[numBlockRows];
    for(Ordinal ptrow=0; ptrow<numPointRows; ++ptrow) xvec[ptrow] = yvec[ptrow];

    for(size_t r = 0; r < numBlockRows; ++r) {
      //if upper-triangular then loop from row N to 0 (bottom to top), but
      //if lower-triangular then loop from 0 to N (top to bottom):
      Ordinal row = upper ? numBlockRows-r-1 : r;

      const Ordinal numRowsInBlock = rptr[row+1]-rptr[row];
      DomainScalar* xx = &xvec[rptr[row]];

      // loop over the blocks in this row and do the multiplication
      for (size_t b=bptr[row]; b<bptr[row+1]; ++b) {
        //if upper-triangular then loop right-to-left (backwards) across
        //the row; but if lower-triangular then loop left-to-right:
        size_t blk = upper ? bptr[row+1]-(b-bptr[row])-1 : b;

        // get pointers into A and y
        const Ordinal col = bindx[blk];
        const Ordinal numColsInBlock = cptr[col+1]-cptr[col];
        const Scalar* A = &vals[indx[blk]];
        DomainScalar* yy = &xvec[cptr[col]];
        if (col==row) {
          if (unitDiag) {
            densematvec_skipdiag(numRowsInBlock, -one, A, yy, xx);
          }
          else {
            if (upper) solve_block_upper(numRowsInBlock, A, xx);
            else solve_block_lower(numRowsInBlock, A, xx);
          }
        }
        else {
          // do the GEMV
          if (numRowsInBlock == numColsInBlock) {
            switch(numRowsInBlock) {
              case 1: dense_matvec_1x1(-one, A, yy, xx); break;
              case 2: dense_matvec_2x2(-one, A, yy, xx); break;
              case 3: dense_matvec_3x3(-one, A, yy, xx); break;
              case 4: dense_matvec_4x4(-one, A, yy, xx); break;
              default: densematvec(numRowsInBlock, numColsInBlock, -one, A, yy, xx);
            }
          }
          else {
            densematvec(numRowsInBlock,numColsInBlock,-one,A,yy,xx);
          }
        }
      }
    }
  }
};

template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
struct DefaultBlockSparseTransposeSolveOp1 {
  // mat data
  const Ordinal *rptr, *cptr;
  const size_t *bptr;
  const Ordinal *bindx, *indx;
  const Scalar  *vals;
  // solve params
  bool unitDiag, upper;
  size_t numBlockRows;
  // multivector data
  const RangeScalar  *y;
  DomainScalar       *x;
  size_t xstride, ystride;

  //find x such that A^T*x = y

  inline KERNEL_PREFIX void execute(size_t i) {
    //solve for i-th vector of multivector
    const RangeScalar* yvec = y+i*ystride;
    DomainScalar*      xvec = x+i*xstride;
    RangeScalar one = Teuchos::ScalarTraits<RangeScalar>::one();

    //copy y into x:
    Ordinal numPointRows = rptr[numBlockRows];
    for(Ordinal ptrow=0; ptrow<numPointRows; ++ptrow) xvec[ptrow] = yvec[ptrow];

    for(size_t r = 0; r < numBlockRows; ++r) {
      //for the transpose solve, if upper-triangular then
      //loop from row 0 to row N, but if lower-triangular then
      //loop from row N to row 0:
      Ordinal row = upper ? r : numBlockRows-r-1;

      const Ordinal numRowsInBlock = rptr[row+1]-rptr[row];
      DomainScalar* xx = &xvec[rptr[row]];

      // loop over the blocks in this row and do the multiplication
      for (size_t b=bptr[row]; b<bptr[row+1]; ++b) {
        //for the tranpose solve, if upper-triangular then loop
        //left-to-right across the row, but if lower-triangular then
        //loop right-to-left across the row:
        size_t blk = upper ? b : bptr[row+1]-(b-bptr[row])-1;

        // get pointers into A and y
        const Ordinal col = bindx[blk];
        const Ordinal numColsInBlock = cptr[col+1]-cptr[col];
        const Scalar* A = &vals[indx[blk]];
        DomainScalar* yy = &xvec[cptr[col]];
        if (col==row) {
          if (unitDiag) {
            densematvec_skipdiag(numRowsInBlock, -one, A, yy, xx);
          }
          else {
            if (upper) solve_block_upper_transpose(numRowsInBlock, A, xx);
            else solve_block_lower_transpose(numRowsInBlock, A, xx);
          }
        }
        else {
          // do the GEMV
          if (numRowsInBlock == numColsInBlock) {
            switch(numRowsInBlock) {
              case 1: dense_matvec_1x1(-one, A, xx, yy); break;
              case 2: dense_matvec_2x2(-one, A, xx, yy); break;
              case 3: dense_matvec_3x3(-one, A, xx, yy); break;
              case 4: dense_matvec_4x4(-one, A, xx, yy); break;
              default: densematvec(numRowsInBlock, numColsInBlock, -one, A, xx, yy);
            }
          }
          else {
            densematvec(numRowsInBlock,numColsInBlock,-one,A,xx,yy);
          }
        }
      }
    }
  }
};

} // namespace Kokkos

#endif
