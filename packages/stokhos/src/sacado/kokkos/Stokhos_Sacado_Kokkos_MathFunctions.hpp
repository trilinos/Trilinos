// $Id$ 
// $Source$ 
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SACADO_KOKKOS_MATHFUNCTIONS_HPP
#define STOKHOS_SACADO_KOKKOS_MATHFUNCTIONS_HPP

#include "KokkosArray_Macros.hpp"

#define UNARYFUNC_MACRO(OP,FADOP)					\
namespace Sacado {							\
									\
  namespace MP {							\
    template <typename T, typename N> class FADOP;			\
    template <typename T, typename N> class Expr;			\
									\
    template <typename T, typename N>					\
    FADOP< T,N >							\
    OP (const Expr<T,N>&);						\
  }									\
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::MP::OP;							\
}

UNARYFUNC_MACRO(exp, ExpOp)
UNARYFUNC_MACRO(log, LogOp)
UNARYFUNC_MACRO(log10, Log10Op)
UNARYFUNC_MACRO(sqrt, SqrtOp)
UNARYFUNC_MACRO(cos, CosOp)
UNARYFUNC_MACRO(sin, SinOp)
UNARYFUNC_MACRO(tan, TanOp)
UNARYFUNC_MACRO(acos, ACosOp)
UNARYFUNC_MACRO(asin, ASinOp)
UNARYFUNC_MACRO(atan, ATanOp)
UNARYFUNC_MACRO(cosh, CoshOp)
UNARYFUNC_MACRO(sinh, SinhOp)
UNARYFUNC_MACRO(tanh, TanhOp)
UNARYFUNC_MACRO(acosh, ACoshOp)
UNARYFUNC_MACRO(asinh, ASinhOp)
UNARYFUNC_MACRO(atanh, ATanhOp)
UNARYFUNC_MACRO(abs, AbsOp)
UNARYFUNC_MACRO(fabs, FAbsOp)

#undef UNARYFUNC_MACRO

#define BINARYFUNC_MACRO(OP,FADOP)					\
namespace Sacado {							\
									\
  namespace MP {							\
    template <typename T1, typename T2, typename N> class FADOP;	\
    template <typename T, typename N> class Expr;			\
									\
    template <typename T1, typename T2, typename N>			\
    FADOP< T1, T2, N >							\
    OP (const Expr<T1,N>&,						\
	const Expr<T2,N>&);						\
									\
    template <typename T, typename N>					\
    FADOP< typename T::value_type, T, N >				\
    OP (const typename T::value_type&,					\
	const Expr<T,N>&);						\
    									\
    template <typename T, typename N>					\
    FADOP< T, typename T::value_type, N >				\
    OP (const Expr<T,N>&,						\
	const typename T::value_type&);					\
  }									\
  									\
}									\
                                                                        \
namespace std {                                                         \
  using Sacado::MP::OP;							\
}

BINARYFUNC_MACRO(atan2, Atan2Op)
BINARYFUNC_MACRO(pow, PowerOp)
BINARYFUNC_MACRO(max, MaxOp)
BINARYFUNC_MACRO(min, MinOp)

#undef BINARYFUNC_MACRO

#endif // STOKHOS_STATIC_ARRAY_TRAITS_HPP
