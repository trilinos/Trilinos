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

#ifndef AMESOS2_BASKER_MP_VECTOR_HPP
#define AMESOS2_BASKER_MP_VECTOR_HPP

#include "Amesos2_config.h"
#ifdef HAVE_AMESOS2_BASKER

#include "Amesos2_Basker.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Amesos2_Solver_MP_Vector.hpp"

// Specialization of BASKER_ScalarTraits for MP::Vector
template <class T> struct BASKER_ScalarTraits;
template <class S>
struct BASKER_ScalarTraits< Sacado::MP::Vector<S> > {
  typedef Sacado::MP::Vector<S> val_type;
  typedef Kokkos::Details::ArithTraits<val_type> KAT;
  typedef typename KAT::mag_type magnitudeType;
  static inline val_type reciprocal(val_type c){ return 1.0/c; }
  static inline val_type divide(val_type a, val_type b){ return a/b; }
  static inline magnitudeType approxABS(val_type a) { return KAT::abs(a); }
  static inline magnitudeType abs(val_type a) { return KAT::abs(a); }
  static inline bool gt (val_type a, val_type b){ return (a>b); }
};

namespace Amesos2 {

  // Enable MP::Vector as a valid Scalar type for Basker
  template <class ST>
  struct TypeMap< Basker,Sacado::MP::Vector<ST> > {
    static Sacado::MP::Vector<ST> dtype;
    typedef Sacado::MP::Vector<ST> type;
    typedef typename Kokkos::Details::ArithTraits< Sacado::MP::Vector<ST> >::mag_type magnitude_type;
  };

  // Specialize our specialization for create_solver_with_supported_type
  // to pass the scalar type directly to Basker
  template < class ST, class LO, class GO, class NO >
  struct create_mp_vector_solver_impl < Basker, ST, LO, GO, NO > {
    typedef Sacado::MP::Vector<ST> SC;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> Matrix;
    typedef Tpetra::MultiVector<SC,LO,GO,NO> Vector;
    static Teuchos::RCP<Solver<Matrix,Vector> >
    apply(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B )
    {
      ctassert<
        Meta::is_same<
          typename MatrixTraits<Matrix>::scalar_t,
          typename MultiVecAdapter<Vector>::scalar_t
        >::value
      > same_scalar_assertion;
      (void)same_scalar_assertion; // This stops the compiler from warning about unused declared variables

      // If our assertion did not fail, then create and return a new solver
      typedef Tpetra::CrsMatrix<Sacado::MP::Vector<ST>,LO,GO,NO> Matrix;
      typedef Tpetra::MultiVector<Sacado::MP::Vector<ST>,LO,GO,NO> Vector;
      return Teuchos::rcp( new Basker<Matrix,Vector>(A, X, B) );
    }
  };
}
#endif // HAVE_AMESOS2_BASKER

#endif // AMESOS2_BASKER_MP_VECTOR_HPP
