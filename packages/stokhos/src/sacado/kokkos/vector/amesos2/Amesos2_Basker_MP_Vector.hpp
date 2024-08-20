// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  typedef Kokkos::ArithTraits<val_type> KAT;
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
    typedef Sacado::MP::Vector<ST> dtype;
    typedef Sacado::MP::Vector<ST> type;
    typedef typename Kokkos::ArithTraits< Sacado::MP::Vector<ST> >::mag_type magnitude_type;
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
        std::is_same_v<
          typename MatrixTraits<Matrix>::scalar_t,
          typename MultiVecAdapter<Vector>::scalar_t
        >
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
