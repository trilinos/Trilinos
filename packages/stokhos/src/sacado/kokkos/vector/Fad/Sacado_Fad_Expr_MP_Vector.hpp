// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_FAD_EXPR_MP_VECTOR_HPP
#define SACADO_FAD_EXPR_MP_VECTOR_HPP

#include "Sacado_Fad_Expression.hpp"

namespace Stokhos {
  template <typename Ord, typename Val, int Num, typename Dev>
  class StaticFixedStorage;
}

namespace Sacado {

  namespace MP {
    template <typename S> class Vector;
  }

  namespace Fad {

    //! Constant expression template
    /*!
     * This template class represents a constant expression.
     */
    template <typename Ord, typename Val, int VecNum, typename Dev>
    class ConstExpr< Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> > > {

    public:

      typedef Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> > ConstT;
      typedef typename ConstT::value_type val_type;

      //! Typename of argument values
      typedef ConstT value_type;

      //! Typename of scalar's (which may be different from ConstT)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Typename of base-expressions
      typedef ConstT base_expr_type;

      //! Constructor
      KOKKOS_INLINE_FUNCTION
      ConstExpr(const ConstT& constant) : constant_(constant) {}

      //! Return value of operation
      KOKKOS_INLINE_FUNCTION
      const ConstT& val() const { return constant_; }

      //! Return value of operation
      KOKKOS_INLINE_FUNCTION
      const val_type& val(int j) const { return constant_.fastAccessCoeff(j); }

    protected:

      //! The constant
      const ConstT& constant_;

    }; // class ConstExpr

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXPR_MP_VECTOR_HPP
