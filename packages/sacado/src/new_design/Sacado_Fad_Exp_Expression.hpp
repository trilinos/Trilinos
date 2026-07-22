// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_EXPRESSION_HPP
#define SACADO_FAD_EXP_EXPRESSION_HPP

#include "Sacado_Traits.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    //! Wrapper for a generic expression template
    /*!
     * This class is used to limit the overload set for building up
     * expressions.  Each expression object should derive from this
     * using CRTP:
     *
     * \code
     * class T : public Expr<T> { ... };
     * \endcode
     *
     * In this case the default implementation here should be correct for
     * any expression class.  If not, an expression class is free to change
     * the implementation through partial specialization.
     */
    template <typename T>
    class Expr {
    public:

      //! Typename of derived object, returned by derived()
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>
       */
      typedef T derived_type;

      //! Return derived object
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>.  This will only compile if this infact the case.
       */
      SACADO_INLINE_FUNCTION
      const derived_type& derived() const {
        return static_cast<const derived_type&>(*this);
      }

      //! Return derived object
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>.  This will only compile if this infact the case.
       */
      SACADO_INLINE_FUNCTION
      const volatile derived_type& derived() const volatile {
        return static_cast<const volatile derived_type&>(*this);
      }
    };

    //! Meta-function for determining nesting with an expression
    /*!
     * This determines the level of nesting within nested Fad types.
     * The default implementation works for any type that isn't a Fad type
     * or an expression of Fad types.
     */
    template <typename T>
    struct ExprLevel {
      static constexpr unsigned value = 0;
    };

    template <typename T>
    struct ExprLevel< Expr<T> > {
      static constexpr unsigned value =
        ExprLevel< typename Expr<T>::value_type >::value + 1;
    };

    //! Determine whether a given type is an expression
    template <typename T>
    struct IsFadExpr {
      static constexpr bool value = false;
    };

    template <typename T>
    struct IsFadExpr< Expr<T> > {
      static constexpr bool value = true;
    };

    // Tag for delegating expression template specializations
    class ExprSpecDefault {};

  } // namespace Exp
  } // namespace Fad

  template <typename T>
  struct IsExpr< Fad::Exp::Expr<T> > {
    static constexpr bool value = true;
  };

  template <typename T>
  struct BaseExprType< Fad::Exp::Expr<T> > {
    typedef typename BaseExprType<T>::type type;
  };

} // namespace Sacado

#include "Sacado_SFINAE_Macros.hpp"

#endif // SACADO_FAD_EXP_EXPRESSION_HPP
