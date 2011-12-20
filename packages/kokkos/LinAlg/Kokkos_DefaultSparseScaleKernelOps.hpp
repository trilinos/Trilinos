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

#ifndef KOKKOS_DEFAULTSPARSESCALE_KERNELOPS_HPP
#define KOKKOS_DEFAULTSPARSESCALE_KERNELOPS_HPP

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

#ifdef __CUDACC__
#include <Teuchos_ScalarTraitsCUDA.hpp>
#else
#include <Teuchos_ScalarTraits.hpp>
#endif

namespace Kokkos {

  template <class Scalar, class Ordinal, class VectorScalar>
  struct DefaultSparseScaleOp1 {
    // mat data
    const size_t  *begs;
    const size_t  *ends;
    const Ordinal *inds;
    Scalar  *vals;
    // matvec params
    size_t numRows;
    // mv data
    const VectorScalar         *x;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      for (size_t c=begs[row]; c != ends[row]; ++c) {
	vals[c]*=(Scalar)x[row];
      }
    }
  };


  template <class Scalar, class Ordinal, class VectorScalar>
  struct DefaultSparseScaleOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    // matvec params
    size_t numRows;
    // mv data
    const VectorScalar  *x;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      Scalar  *curval = vals_beg[row];
      for (size_t j=0; j != numEntries[row]; ++j) {
	curval[j]*=(Scalar)x[row];
      }
    }
  };

} // namespace Kokkos

#endif
