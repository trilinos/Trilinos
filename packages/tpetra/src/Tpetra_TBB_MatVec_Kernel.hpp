
#ifndef TPETRA_TBB_MATVEC_KERNEL_HPP
#define TPETRA_TBB_MATVEC_KERNEL_HPP

// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp"

#ifdef HAVE_TPETRA_TBB

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace Tpetra {

template<typename OrdinalType, typename ScalarType>
struct crsmv_op {
  crsmv_op(OrdinalType numrows,
           const OrdinalType* rowptr,
           const OrdinalType* colindices,
           const ScalarType* coefs,
           const ScalarType* xvec,
           ScalarType* yvec)
   : nrows(numrows), rptr(rowptr), colinds(colindices), vals(coefs),
     y(yvec), x(xvec) {}

  OrdinalType nrows;
  const OrdinalType* rptr;
  const OrdinalType* colinds;
  const ScalarType* vals;
  ScalarType* y;
  const ScalarType* x;

  void operator()(const tbb::blocked_range<int>& range) const {

    for(int i=range.begin(); i!=range.end(); ++i) {
      OrdinalType jbeg = rptr[i];
      OrdinalType jend = rptr[i+1];
      ScalarType sum = 0.0;

      for(OrdinalType j=jbeg; j<jend; ++j) {
        sum += vals[j]*x[colinds[j]];
      }
      y[i] = sum;
    }

  }
};

template<typename OrdinalType, typename ScalarType>
void threaded_crsmv(OrdinalType numrows,
                    const OrdinalType* rowptr,
                    const OrdinalType* colindices,
                    const ScalarType* coefs,
                    const ScalarType* xvec,
                    ScalarType* yvec)
{
  crsmv_op<OrdinalType,ScalarType> op(numrows,
                                      rowptr, colindices, coefs,
                                      xvec, yvec);

  OrdinalType grain_size = numrows>8 ? numrows/8 : 1;

  tbb::parallel_for(tbb::blocked_range<OrdinalType>(0, numrows, grain_size), op);
}

}//namespace Tpetra

#endif

#endif

