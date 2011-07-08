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
      const size_t rhs = (i - row) / numRows;
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
      const size_t rhs = (i - row) / numRows;
      Scalar  *curval = vals_beg[row];
      const Ordinal *curind = inds_beg[row];
      for (size_t j=0; j != numEntries[row]; ++j) {
	curval[j]*=(Scalar)x[row];
      }
    }
  };

} // namespace Kokkos

#endif
