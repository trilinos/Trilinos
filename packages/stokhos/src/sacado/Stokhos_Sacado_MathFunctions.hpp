// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
UNARYFUNC_MACRO(cbrt, CbrtOp)
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
