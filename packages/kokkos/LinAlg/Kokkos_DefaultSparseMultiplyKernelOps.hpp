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

namespace Kokkos {

  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseMultiplyOp1 {
    // mat data
    const size_t  *offsets;
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
      for (size_t c=offsets[row]; c != offsets[row+1]; ++c) {
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
    const size_t  *offsets;
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
        for (size_t c=offsets[row]; c != offsets[row+1]; ++c) {
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
