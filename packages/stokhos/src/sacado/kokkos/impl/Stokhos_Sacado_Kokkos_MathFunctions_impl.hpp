// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Stokhos_Sacado_Kokkos_MathFucntions_impl.hpp> without macros defined"

#else

#define UNARYFUNC_MACRO(OP,FADOP)					\
namespace Sacado {							\
									\
  namespace MP {							\
    template <typename T, typename N> class FADOP;			\
    template <typename T, typename N> class Expr;			\
									\
    template <typename T>						\
    FADOP< T,KOKKOSARRAY_MACRO_DEVICE >					\
    OP (const Expr<T,KOKKOSARRAY_MACRO_DEVICE>&);			\
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
    template <typename T1, typename T2>					\
    FADOP< T1, T2, KOKKOSARRAY_MACRO_DEVICE >				\
    OP (const Expr<T1,KOKKOSARRAY_MACRO_DEVICE>&,			\
	const Expr<T2,KOKKOSARRAY_MACRO_DEVICE>&);			\
    									\
    template <typename T>						\
    FADOP< typename T::value_type, T, KOKKOSARRAY_MACRO_DEVICE >	\
    OP (const typename T::value_type&,					\
	const Expr<T,KOKKOSARRAY_MACRO_DEVICE>&);			\
    									\
    template <typename T>						\
    FADOP< T, typename T::value_type, KOKKOSARRAY_MACRO_DEVICE >	\
    OP (const Expr<T,KOKKOSARRAY_MACRO_DEVICE>&,			\
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

#endif
