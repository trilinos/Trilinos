//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef KOKKOS_DEFAULTSPARSEMULTIPLY_KERNELOPS_HPP
#define KOKKOS_DEFAULTSPARSEMULTIPLY_KERNELOPS_HPP

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

#ifdef __CUDACC__
#include <Teuchos_ScalarTraitsCUDA.hpp>
#else
#include <Teuchos_ScalarTraits.hpp>
#endif

// TODO: begs/ends instead of offsets/rowsizes seems clever, as it allows the same code to be
// use for packed and non-packed storage, by setting ends=begs+1. 
// however, in the packed scenario, it requires an additional (redundant) N reads
// this may become significant as NNZ -> N
// this should be revisited; we may have to add DefaultSparseMultiplyPacked and DefaultSparseTransposeMultiplyPacked kernels

namespace Kokkos {

  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseMultiplyOp1 {
    // mat data
    const size_t  *begs;
    const size_t  *ends;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      const size_t rhs = (i - row) / numRows;
      RangeScalar tmp = Teuchos::ScalarTraits<RangeScalar>::zero();
      const DomainScalar *xj = x + rhs * xstride;
      RangeScalar        *yj = y + rhs * ystride;
      for (size_t c=begs[row]; c != ends[row]; ++c) {
        tmp += (RangeScalar)vals[c] * (RangeScalar)xj[inds[c]];
      }
      if (NO_BETA_AND_OVERWRITE) {
        yj[row] = (RangeScalar)alpha * tmp;
      }
      else {
        RangeScalar tmp2 = beta * yj[row];
        yj[row] = (RangeScalar)(alpha * tmp + tmp2);
      }
    }
  };

  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseTransposeMultiplyOp1 {
    // mat data
    const size_t  *begs;
    const size_t  *ends;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // multiply entire matrix for rhs i
      const size_t rhs = i;
      const DomainScalar *xj = x + rhs * xstride;
      const RangeScalar RANGE_ZERO = Teuchos::ScalarTraits<RangeScalar>::zero();
      RangeScalar        *yj = y + rhs * ystride;
      for (size_t row=0; row < numCols; ++row) {
        if (NO_BETA_AND_OVERWRITE) {
          yj[row] = RANGE_ZERO;
        }
        else {
          yj[row] = (RangeScalar)(yj[row] * beta);
        }
      }
      for (size_t row=0; row < numRows; ++row) {
        for (size_t c=begs[row]; c != ends[row]; ++c) {
          yj[inds[c]] += (RangeScalar)(alpha * Teuchos::ScalarTraits<RangeScalar>::conjugate(vals[c]) * (RangeScalar)xj[row]);
        }
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseMultiplyOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      const size_t rhs = (i - row) / numRows;
      RangeScalar tmp = Teuchos::ScalarTraits<RangeScalar>::zero();
      const DomainScalar *xj = x + rhs * xstride;
      RangeScalar        *yj = y + rhs * ystride;
      const Scalar  *curval = vals_beg[row];
      const Ordinal *curind = inds_beg[row];
      for (size_t j=0; j != numEntries[row]; ++j) {
        tmp += (RangeScalar)curval[j] * (RangeScalar)xj[curind[j]];
      }
      if (NO_BETA_AND_OVERWRITE) {
        yj[row] = (RangeScalar)alpha * tmp;
      }
      else {
        RangeScalar tmp2 = beta * yj[row];
        yj[row] = (RangeScalar)(alpha * tmp + tmp2);
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseTransposeMultiplyOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // multiply entire matrix for rhs i
      const size_t rhs = i;
      const RangeScalar RANGE_ZERO = Teuchos::ScalarTraits<RangeScalar>::zero();
      const DomainScalar *xj = x + rhs * xstride;
      RangeScalar        *yj = y + rhs * ystride;
      for (size_t row=0; row < numCols; ++row) {
        if (NO_BETA_AND_OVERWRITE) {
          yj[row] = RANGE_ZERO;
        }
        else {
          yj[row] = (RangeScalar)(yj[row] * beta);
        }
      }
      for (size_t row=0; row < numRows; ++row) {
        const Scalar  *rowval = vals_beg[row];
        const Ordinal *rowind = inds_beg[row];
        for (size_t j=0; j != numEntries[row]; ++j) {
          yj[rowind[j]] += (RangeScalar)(alpha * Teuchos::ScalarTraits<RangeScalar>::conjugate(rowval[j]) * (RangeScalar)xj[row]);
        }
      }
    }
  };

} // namespace Kokkos

#endif
