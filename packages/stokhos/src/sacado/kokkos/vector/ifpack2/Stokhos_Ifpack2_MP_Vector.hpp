// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_IFPACK2_MP_VECTOR_HPP
#define STOKHOS_IFPACK2_MP_VECTOR_HPP

// This header file should be included whenever compiling any Ifpack2
// code with Stokhos scalar types

// MP includes and specializations
#include "Stokhos_Tpetra_MP_Vector.hpp"
#include "Ifpack2_Krylov_MP_Vector.hpp"

// Specialization of LocalReciprocalThreshold functor in
// Ifpack2_Details_Chebyshev_def.hpp
namespace Ifpack2 {
namespace Details {

template<class XV, class SizeType>
struct LocalReciprocalThreshold;

template<class S, class ... P, class SizeType>
struct LocalReciprocalThreshold<
  Kokkos::View< Sacado::MP::Vector<S>*,P...>, SizeType > {
  typedef Kokkos::View<Sacado::MP::Vector<S>*,P...> XV;

  static void
  compute (const XV& X,
           const typename XV::non_const_value_type& minVal)
  {
    if (!Sacado::is_constant(minVal)) {
      Kokkos::Impl::raise_error(
        "LocalReciprocalThreshold not implemented for non-constant minVal");
    }

    typedef typename Kokkos::FlatArrayType<XV>::type Flat_XV;
    Flat_XV flat_X = X;
    LocalReciprocalThreshold< Flat_XV, SizeType >::compute( flat_X,
                                                            minVal.coeff(0) );
  }
};

}
}

#endif // STOKHOS_IFPACK2_MP_VECTOR_HPP
