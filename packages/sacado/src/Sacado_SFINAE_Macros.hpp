// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_SFINAE_MACROS_H
#define SACADO_SFINAE_MACROS_H

#include <type_traits>

#include "Sacado_mpl_enable_if.hpp"
#include "Sacado_mpl_disable_if.hpp"
#include "Sacado_mpl_type_wrap.hpp"

/* Define some macros useful for disabling template function overloads */
#define SACADO_ENABLE_IF_SAME(TYPE1, TYPE2, RETURN_TYPE)              \
  typename mpl::enable_if_c<std::is_convertible< TYPE1 , TYPE2 >::value && ExprLevel<TYPE1>::value == ExprLevel<TYPE2>::value, RETURN_TYPE >::type
#define SACADO_ENABLE_EXPR_FUNC(RETURN_TYPE) \
  SACADO_ENABLE_IF_SAME(typename Expr<S>::value_type, value_type, RETURN_TYPE)
#define SACADO_ENABLE_EXPR_CTOR_DEF SACADO_ENABLE_EXPR_FUNC(void*)
#define SACADO_ENABLE_EXPR_CTOR_DECL SACADO_ENABLE_EXPR_CTOR_DEF = 0
#define SACADO_FAD_ENABLE_EXPR_FUNC \
  SACADO_ENABLE_IF_SAME(typename Expr<S>::value_type, typename FAD::value_type, FAD&)

#define SACADO_EXP_ENABLE_EXPR_FUNC(RETURN_TYPE) \
  SACADO_ENABLE_IF_SAME(typename Expr<S>::derived_type::value_type, value_type, RETURN_TYPE)
#define SACADO_EXP_ENABLE_EXPR_CTOR_DEF \
  typename mpl::enable_if_c<                                            \
    std::is_convertible< typename Expr<S>::derived_type::value_type ,   \
                         value_type >::value &&                         \
    (ExprLevel< typename Expr<S>::derived_type::value_type >::value ==  \
     ExprLevel< value_type >::value) &&                                 \
    !is_view, void* >::type
#define SACADO_EXP_ENABLE_EXPR_CTOR_DECL SACADO_EXP_ENABLE_EXPR_CTOR_DEF = 0
#define SACADO_FAD_EXP_ENABLE_EXPR_FUNC \
  SACADO_ENABLE_IF_SAME(typename Expr<S>::derived_type::value_type, typename FAD::value_type, FAD&)

#define SACADO_ENABLE_IF_CONVERTIBLE(TYPE1, TYPE2, RETURN_TYPE)              \
  typename Sacado::mpl::enable_if<std::is_convertible< TYPE1 , TYPE2 >, RETURN_TYPE >::type
#define SACADO_ENABLE_VALUE_FUNC(RETURN_TYPE) \
  SACADO_ENABLE_IF_CONVERTIBLE(S, value_type, RETURN_TYPE)
#define SACADO_ENABLE_VALUE_CTOR_DEF SACADO_ENABLE_VALUE_FUNC(void*)
#define SACADO_ENABLE_VALUE_CTOR_DECL SACADO_ENABLE_VALUE_CTOR_DEF = 0

#define SACADO_EXP_ENABLE_VALUE_CTOR_DEF \
  typename mpl::enable_if_c<                                            \
    std::is_convertible< S , value_type >::value &&             \
    !is_view, void* >::type
#define SACADO_EXP_ENABLE_VALUE_CTOR_DECL SACADO_EXP_ENABLE_VALUE_CTOR_DEF = 0

#define SACADO_FAD_OP_ENABLE_EXPR_EXPR(OP)                              \
  typename mpl::enable_if_c< IsFadExpr<T1>::value && IsFadExpr<T2>::value && \
                             ExprLevel<T1>::value == ExprLevel<T2>::value, \
                             Expr< OP< T1, T2 > >                       \
                           >::type
#define SACADO_FAD_EXP_OP_ENABLE_EXPR_EXPR(OP)                              \
  typename mpl::enable_if_c< IsFadExpr<T1>::value && IsFadExpr<T2>::value && \
                             ExprLevel<T1>::value == ExprLevel<T2>::value, \
                             OP< typename Expr<T1>::derived_type, typename Expr<T2>::derived_type, false, false, typename T1::expr_spec_type > \
                           >::type
#define SACADO_FAD_OP_ENABLE_SCALAR_EXPR(OP)                            \
  typename mpl::disable_if<std::is_same< typename Expr<T>::value_type, typename Expr<T>::scalar_type>, Expr< OP< ConstExpr<typename Expr<T>::scalar_type>, Expr<T> > > >::type
#define SACADO_FAD_OP_ENABLE_EXPR_SCALAR(OP)                            \
  typename mpl::disable_if<std::is_same< typename Expr<T>::value_type, typename Expr<T>::scalar_type>, Expr< OP< Expr<T>, ConstExpr<typename Expr<T>::scalar_type> > > >::type
#define SACADO_FAD_EXP_OP_ENABLE_SCALAR_EXPR(OP)                            \
  typename mpl::disable_if<std::is_same< typename T::value_type, typename T::scalar_type>, OP< typename T::scalar_type, typename Expr<T>::derived_type, true, false, typename T::expr_spec_type > >::type
#define SACADO_FAD_EXP_OP_ENABLE_EXPR_SCALAR(OP)                            \
  typename mpl::disable_if<std::is_same< typename T::value_type, typename T::scalar_type>, OP< typename Expr<T>::derived_type, typename T::scalar_type, false, true, typename T::expr_spec_type > >::type

#endif // SACADO_SFINAE_MACROS_H
