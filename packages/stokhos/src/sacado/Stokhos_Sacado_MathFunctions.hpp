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

#ifndef STOKHOS_SACADO_MATH_FUNCTIONS_HPP
#define STOKHOS_SACADO_MATH_FUNCTIONS_HPP

#define UNARYFUNC_MACRO(OP,FADOP)					\
namespace Sacado {							\
									\
  namespace PCE {							\
    template <typename T, typename S> class OrthogPoly;			\
    template <typename T, typename S>					\
    OrthogPoly<T,S> OP (const OrthogPoly<T,S>&);			\
  }									\
									\
  namespace ETPCE {							\
    template <typename T> class FADOP;					\
    template <typename T> class Expr;					\
    template <typename T>						\
    Expr< FADOP< Expr<T> > > OP (const Expr<T>&);			\
  }									\
									\
  namespace ETV {							\
    template <typename T> class FADOP;					\
    template <typename T> class Expr;					\
    template <typename T>						\
    Expr< FADOP< Expr<T> > > OP (const Expr<T>&);			\
  }									\
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::PCE::OP;						\
  using Sacado::ETPCE::OP;						\
  using Sacado::ETV::OP;						\
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
  namespace PCE {							\
    template <typename T, typename S> class OrthogPoly;			\
    template <typename T, typename S>					\
    OrthogPoly<T,S> OP (const OrthogPoly<T,S>&,				\
			const OrthogPoly<T,S>&);			\
    template <typename T, typename S>					\
    OrthogPoly<T,S> OP (const T&,					\
			const OrthogPoly<T,S>&);			\
    template <typename T, typename S>					\
    OrthogPoly<T,S> OP (const OrthogPoly<T,S>&,				\
			const T&);					\
  }									\
									\
  namespace ETPCE {							\
    template <typename T1, typename T2> class FADOP;			\
    template <typename T> class Expr;					\
    template <typename T> class ConstExpr;				\
    template <typename T1, typename T2>					\
    Expr< FADOP< Expr<T1>, Expr<T2> > >					\
    OP (const Expr<T1>&, const Expr<T2>&);				\
									\
    template <typename T>						\
    Expr< FADOP< Expr<T>, Expr<T> > >					\
    OP (const Expr<T>&, const Expr<T>&);				\
									\
    template <typename T>						\
    Expr< FADOP< typename Expr<T>::value_type, Expr<T> > >		\
    OP (const typename Expr<T>::value_type&, const Expr<T>&);		\
									\
    template <typename T>						\
    Expr< FADOP< Expr<T>, typename Expr<T>::value_type > >		\
    OP (const Expr<T>&, const typename Expr<T>::value_type&);		\
  }									\
									\
  namespace ETV {							\
    template <typename T1, typename T2> class FADOP;			\
    template <typename T> class Expr;					\
    template <typename T1, typename T2>					\
    Expr< FADOP< Expr<T1>, Expr<T2> > >					\
    OP (const Expr<T1>&, const Expr<T2>&);				\
									\
    template <typename T>						\
    Expr< FADOP< Expr<T>, Expr<T> > >					\
    OP (const Expr<T>&, const Expr<T>&);				\
									\
    template <typename T>						\
    Expr< FADOP< typename Expr<T>::value_type, Expr<T> > >		\
    OP (const typename Expr<T>::value_type&, const Expr<T>&);		\
									\
    template <typename T>						\
    Expr< FADOP< Expr<T>, typename Expr<T>::value_type > >		\
    OP (const Expr<T>&, const typename Expr<T>::value_type&);		\
  }									\
									\
}									\
                                                                        \
namespace std {                                                         \
  using Sacado::PCE::OP;						\
  using Sacado::ETPCE::OP;						\
  using Sacado::ETV::OP;						\
}

BINARYFUNC_MACRO(atan2, Atan2Op)
BINARYFUNC_MACRO(pow, PowerOp)
BINARYFUNC_MACRO(max, MaxOp)
BINARYFUNC_MACRO(min, MinOp)

#undef BINARYFUNC_MACRO

#endif // STOKHOS_SACADO_MATH_FUNCTIONS_HPP 
