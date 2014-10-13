// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
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

#ifndef SACADO_CACHEFAD_EXPRESSION_HPP
#define SACADO_CACHEFAD_EXPRESSION_HPP

#include "Sacado_Traits.hpp"

namespace Sacado {

  namespace CacheFad {

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all Fad expression
     * template classes.
     */
    template <typename ExprT>
    class Expr {};

    //! Constant expression template
    /*!
     * This template class represents a constant expression.
     */
    template <typename ConstT>
    class ConstExpr {

    public:

      //! Typename of argument values
      typedef ConstT value_type;

      //! Typename of scalar values
      typedef ConstT scalar_type;

      //! Constructor
      KOKKOS_INLINE_FUNCTION
      ConstExpr(const ConstT& constant) : constant_(constant) {}

      //! Return size of the derivative array of the operation
      KOKKOS_INLINE_FUNCTION
      int size() const { return 0; }

      //! Return if operation has fast access
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const { return 1; }

      //! Cache values
      KOKKOS_INLINE_FUNCTION
      void cache() const {}

      //! Return value of operation
      KOKKOS_INLINE_FUNCTION
      value_type val() const { return constant_; }

      //! Return derivative component \c i of operation
      KOKKOS_INLINE_FUNCTION
      value_type dx(int i) const { return value_type(0); }

      //! Return derivative component \c i of operation
      KOKKOS_INLINE_FUNCTION
      value_type fastAccessDx(int i) const { return value_type(0); }

    protected:

      //! The constant
      const ConstT& constant_;

    }; // class ConstExpr

  } // namespace CacheFad

} // namespace Sacado

#endif // SACADO_CACHEFAD_EXPRESSION_HPP
