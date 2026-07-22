// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
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

#ifndef SACADO_FAD_EXPRESSION_HPP
#define SACADO_FAD_EXPRESSION_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_Fad_ExpressionFwd.hpp"

namespace Sacado {

  namespace Fad {

    //! Meta-function for determining concrete base expression
    /*!
     * This determines the concrete base expression type of each leaf in
     * an expression tree.  The Promote meta-function is then used to promote
     * all of the leaves to a single expression type that the whole expression
     * can be assigned/promoted to.  This allows Promote to operate on
     * expressions as well as AD types.
     */
    template <typename> struct BaseExpr {};

    struct ExprSpecDefault {};
    template <typename ExprT> struct ExprSpec {
      typedef ExprSpecDefault type;
    };

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all Fad expression
     * template classes.
     */
    template <typename ExprT, typename Spec>
    class Expr {
    public:
      typedef ExprT value_type;
    };

    template <typename ExprT, typename Spec>
    struct ExprSpec< Expr<ExprT,Spec> > {
      typedef Spec type;
    };

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

      //! Constructor
      SACADO_INLINE_FUNCTION
      ConstExpr(const ConstT& constant) : constant_(constant) {}

      //! Return value of operation
      SACADO_INLINE_FUNCTION
      const ConstT& val() const { return constant_; }

      //! Return value of operation
      SACADO_INLINE_FUNCTION
      const ConstT& val(int j) const { return constant_; }

    protected:

      //! The constant
      const ConstT& constant_;

    }; // class ConstExpr

  } // namespace Fad

  template <typename T>
  struct IsExpr< Fad::Expr<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct BaseExprType< Fad::Expr<T> > {
    typedef typename Fad::Expr<T>::base_expr_type type;
  };

  template <typename T>
  struct ValueType< Fad::ConstExpr<T> > {
    typedef typename Fad::ConstExpr<T>::value_type type;
  };

} // namespace Sacado

#include "Sacado_SFINAE_Macros.hpp"

#endif // SACADO_FAD_EXPRESSION_HPP
