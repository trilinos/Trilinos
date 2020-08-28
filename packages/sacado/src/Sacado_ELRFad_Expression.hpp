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

#ifndef SACADO_ELRFAD_EXPRESSION_HPP
#define SACADO_ELRFAD_EXPRESSION_HPP

#include "Sacado_Traits.hpp"

namespace Sacado {

  namespace ELRFad {

    //! Meta-function for determining concrete base expression
    /*!
     * This determines the concrete base expression type of each leaf in
     * an expression tree.  The Promote meta-function is then used to promote
     * all of the leaves to a single expression type that the whole expression
     * can be assigned/promoted to.  This allows Promote to operate on
     * expressions as well as AD types.
     */
    template <typename> struct BaseExpr {};

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all Fad expression
     * template classes.
     */
    template <typename ExprT>
    class Expr {};

    //! Meta-function for determining nesting with an expression
    /*!
     * This determines the level of nesting within nested Fad types.
     * The default implementation works for any type that isn't a Fad type
     * or an expression of Fad types.
     */
    template <typename T>
    struct ExprLevel {
      static const unsigned value = 0;
    };

    template <typename T>
    struct ExprLevel< Expr<T> > {
      static const unsigned value =
        ExprLevel< typename Expr<T>::value_type >::value + 1;
    };

    //! Determine whether a given type is an expression
    template <typename T>
    struct IsFadExpr {
      static const bool value = false;
    };

    template <typename T>
    struct IsFadExpr< Expr<T> > {
      static const bool value = true;
    };

    //! Constant expression template
    /*!
     * This template class represents a constant expression.
     */
    template <typename ConstT>
    class ConstExpr {

    public:

      //! Typename of argument values
      typedef ConstT value_type;

      //! Typename of scalar's (which may be different from ConstT)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Typename of base-expressions
      typedef ConstT base_expr_type;

      //! Number of arguments
      static const int num_args = 0;

      //! Constructor
      SACADO_INLINE_FUNCTION
      ConstExpr(const ConstT& constant) : constant_(constant) {}

      //! Return size of the derivative array of the operation
      SACADO_INLINE_FUNCTION
      int size() const { return 0; }

      //! Return whether value should be updated
      SACADO_INLINE_FUNCTION
      bool updateValue() const { return true; }

      //! Return value of operation
      SACADO_INLINE_FUNCTION
      value_type val() const { return constant_; }

      //! Return partials w.r.t. arguments
      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {}

      //! Rturn tangent component \c i of arguments
      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {}

      //! Return tangent component \c i of argument \c Arg
      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const { return 0.0; }

      //! Return whether argument is active
      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const { return false; }

      //! Return whether expression is linear
      SACADO_INLINE_FUNCTION
      bool isLinear() const { return true; }

    protected:

      //! The constant
      const ConstT& constant_;

    }; // class ConstExpr

    template <typename T> struct ExprConstRef {
      typedef const T& type;
    };
    template <typename T> struct ExprConstRef< ConstExpr<T> > {
      typedef const ConstExpr<T> type;
    };

  } // namespace ELRFad

  template <typename T>
  struct IsExpr< ELRFad::Expr<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct BaseExprType< ELRFad::Expr<T> > {
    typedef typename ELRFad::Expr<T>::base_expr_type type;
  };

} // namespace Sacado

#include "Sacado_SFINAE_Macros.hpp"

#endif // SACADO_ELRFAD_EXPRESSION_HPP
