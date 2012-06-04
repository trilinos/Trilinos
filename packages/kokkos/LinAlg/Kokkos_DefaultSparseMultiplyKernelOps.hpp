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
  struct DefaultSparseMultiplyOp {
    // mat data
    const size_t  *ptrs;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t numRHS, xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t row) {
      const Scalar  *v = vals + ptrs[row];
      const Ordinal *i = inds + ptrs[row],
                   *ie = inds + ptrs[row+1];
      if (NO_BETA_AND_OVERWRITE) {
        for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] = Teuchos::ScalarTraits<RangeScalar>::zero();
      }
      else {
        for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] *= beta;
      }
      // save the extra multiplication if possible
      if (alpha == Teuchos::ScalarTraits<RangeScalar>::one()) {
        while (i != ie) 
        {
          const  Scalar val = *v++;
          const Ordinal ind = *i++;
          for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] += (RangeScalar)val * (RangeScalar)x[j*xstride+ind];
        }
      }
      else { // alpha != one
        while (i != ie) 
        {
          const  Scalar val = *v++;
          const Ordinal ind = *i++;
          for (size_t j=0; j<numRHS; ++j) y[j*ystride+row] += alpha * (RangeScalar)val * (RangeScalar)x[j*xstride+ind];
        }
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar, int NO_BETA_AND_OVERWRITE>
  struct DefaultSparseTransposeMultiplyOp {
    // mat data
    const size_t  *ptrs;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    RangeScalar        alpha, beta;
    size_t numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t numRHS, xstride, ystride;

    inline void execute() {
      if (NO_BETA_AND_OVERWRITE) {
        for (size_t j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (size_t row=0; row<numCols; ++row) {
            yp[row] = Teuchos::ScalarTraits<RangeScalar>::zero();
          }
        }
      }
      else {
        for (size_t j=0; j<numRHS; ++j) {
          RangeScalar *yp = y+j*ystride;
          for (size_t row=0; row<numCols; ++row) {
            yp[row] *= beta;
          }
        }
      }
      // save the extra multiplication if possible
      if (alpha == Teuchos::ScalarTraits<RangeScalar>::one()) {
        for (size_t colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v = vals + ptrs[colAt];
          const Ordinal *i = inds + ptrs[colAt],
                       *ie = inds + ptrs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind
          while (i != ie) {
            const  Scalar val = Teuchos::ScalarTraits<RangeScalar>::conjugate(*v++);
            const Ordinal ind = *i++;
            for (size_t j=0; j<numRHS; ++j) y[j*ystride+ind] += (RangeScalar)val * (RangeScalar)x[j*xstride+colAt];
          }
        }
      }
      else { // alpha != one
        for (size_t colAt=0; colAt < numRows; ++colAt) {
          const Scalar  *v = vals + ptrs[colAt];
          const Ordinal *i = inds + ptrs[colAt],
                       *ie = inds + ptrs[colAt+1];
          // sparse outer product: AT[:,colAt] * X[ind
          while (i != ie) {
            const  Scalar val = Teuchos::ScalarTraits<RangeScalar>::conjugate(*v++);
            const Ordinal ind = *i++;
            for (size_t j=0; j<numRHS; ++j) y[j*ystride+ind] += alpha * (RangeScalar)val * (RangeScalar)x[j*xstride+colAt];
          }
        }
      }
    }
  };

} // namespace Kokkos

#endif
