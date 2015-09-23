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
// @HEADER

#ifndef SACADO_SFINAE_MACROS_H
#define SACADO_SFINAE_MACROS_H

#include "Sacado_mpl_enable_if.hpp"
#include "Sacado_mpl_disable_if.hpp"
#include "Sacado_mpl_is_same.hpp"
#include "Sacado_mpl_is_convertible.hpp"
#include "Sacado_mpl_type_wrap.hpp"

/* Define some macros useful for disabling template function overloads */
#define SACADO_ENABLE_IF_SAME(TYPE1, TYPE2, RETURN_TYPE)              \
  typename mpl::enable_if_c<mpl::is_convertible< TYPE1 , TYPE2 >::value && ExprLevel<TYPE1>::value == ExprLevel<TYPE2>::value, RETURN_TYPE >::type
#define SACADO_ENABLE_EXPR_FUNC(RETURN_TYPE) \
  SACADO_ENABLE_IF_SAME(typename Expr<S>::value_type, value_type, RETURN_TYPE)
#define SACADO_ENABLE_EXPR_CTOR_DEF SACADO_ENABLE_EXPR_FUNC(void*)
#define SACADO_ENABLE_EXPR_CTOR_DECL SACADO_ENABLE_EXPR_CTOR_DEF = 0
#define SACADO_FAD_ENABLE_EXPR_FUNC \
  SACADO_ENABLE_IF_SAME(typename Expr<S>::value_type, typename FAD::value_type, FAD&)

#define SACADO_ENABLE_IF_CONVERTIBLE(TYPE1, TYPE2, RETURN_TYPE)              \
  typename Sacado::mpl::enable_if<Sacado::mpl::is_convertible< TYPE1 , TYPE2 >, RETURN_TYPE >::type
#define SACADO_ENABLE_VALUE_FUNC(RETURN_TYPE) \
  SACADO_ENABLE_IF_CONVERTIBLE(S, value_type, RETURN_TYPE)
#define SACADO_ENABLE_VALUE_CTOR_DEF SACADO_ENABLE_VALUE_FUNC(void*)
#define SACADO_ENABLE_VALUE_CTOR_DECL SACADO_ENABLE_VALUE_CTOR_DEF = 0

#define SACADO_FAD_OP_ENABLE_EXPR_EXPR(OP)                              \
  typename mpl::enable_if_c< IsFadExpr<T1>::value && IsFadExpr<T2>::value && \
                             ExprLevel<T1>::value == ExprLevel<T2>::value, \
                             Expr< OP< T1, T2 > >                       \
                           >::type
#define SACADO_FAD_OP_ENABLE_SCALAR_EXPR(OP)                            \
  typename mpl::disable_if<mpl::is_same< typename Expr<T>::value_type, typename Expr<T>::scalar_type>, Expr< OP< ConstExpr<typename Expr<T>::scalar_type>, Expr<T> > > >::type
#define SACADO_FAD_OP_ENABLE_EXPR_SCALAR(OP)                            \
  typename mpl::disable_if<mpl::is_same< typename Expr<T>::value_type, typename Expr<T>::scalar_type>, Expr< OP< Expr<T>, ConstExpr<typename Expr<T>::scalar_type> > > >::type

#endif // SACADO_SFINAE_MACROS_H
