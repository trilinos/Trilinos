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

#ifndef SACADO_ETV_VECTOR_OPS_HPP
#define SACADO_ETV_VECTOR_OPS_HPP

#include "Sacado_cmath.hpp"
#include <ostream>	// for std::ostream

#define ETV_UNARYOP_MACRO(OPNAME,OP,OPER)				\
namespace Sacado {							\
  namespace ETV {							\
									\
    template <typename ExprT>						\
    class OP {};							\
									\
    template <typename ExprT>						\
    class Expr< OP<ExprT> > {						\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
      typedef typename ExprT::storage_type storage_type;		\
      typedef typename ExprT::base_expr_type base_expr_type;            \
									\
      Expr(const ExprT& expr_) : expr(expr_)  {}			\
									\
      std::string name() const {					\
	return std::string(#OPER) + expr.name();			\
      }									\
									\
      int size() const {						\
	return expr.size();						\
      }									\
									\
      bool hasFastAccess(int sz) const {				\
	return expr.hasFastAccess(sz);					\
      }									\
									\
      value_type val() const {						\
	return OPER(expr.val());					\
      }									\
      									\
      value_type coeff(int i) const {					\
	return OPER(expr.coeff(i));					\
      }									\
									\
      value_type fastAccessCoeff(int i) const {				\
	return OPER(expr.fastAccessCoeff(i));				\
      }									\
									\
    protected:								\
									\
      const ExprT& expr;						\
									\
    };									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T> > >					\
    OPNAME (const Expr<T>& expr)					\
    {									\
      typedef OP< Expr<T> > expr_t;					\
      									\
      return Expr<expr_t>(expr);					\
    }									\
  }									\
}

ETV_UNARYOP_MACRO(operator+, UnaryPlusOp, +)
ETV_UNARYOP_MACRO(operator-, UnaryMinusOp, -)
ETV_UNARYOP_MACRO(exp, ExpOp, std::exp)
ETV_UNARYOP_MACRO(log, LogOp, std::log)
ETV_UNARYOP_MACRO(log10, Log10Op, std::log10)
ETV_UNARYOP_MACRO(sqrt, SqrtOp, std::sqrt)
ETV_UNARYOP_MACRO(cos, CosOp, std::cos)
ETV_UNARYOP_MACRO(sin, SinOp, std::sin)
ETV_UNARYOP_MACRO(tan, TanOp, std::tan)
ETV_UNARYOP_MACRO(acos, ACosOp, std::acos)
ETV_UNARYOP_MACRO(asin, ASinOp, std::asin)
ETV_UNARYOP_MACRO(atan, ATanOp, std::atan)
ETV_UNARYOP_MACRO(cosh, CoshOp, std::cosh)
ETV_UNARYOP_MACRO(sinh, SinhOp, std::sinh)
ETV_UNARYOP_MACRO(tanh, TanhOp, std::tanh)
ETV_UNARYOP_MACRO(acosh, ACoshOp, std::acosh)
ETV_UNARYOP_MACRO(asinh, ASinhOp, std::asinh)
ETV_UNARYOP_MACRO(atanh, ATanhOp, std::atanh)
ETV_UNARYOP_MACRO(abs, AbsOp, std::abs)
ETV_UNARYOP_MACRO(fabs, FAbsOp, std::fabs)

#undef ETV_UNARYOP_MACRO

#define ETV_BINARYOP_MACRO(OPNAME,OP,OPER)				\
namespace Sacado {							\
  namespace ETV {							\
									\
    template <typename ExprT1, typename ExprT2>				\
    class OP {};							\
									\
    template <typename T1, typename T2>					\
    class Expr< OP< Expr<T1>, Expr<T2> > > {				\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT1::value_type value_type_1;			\
      typedef typename ExprT2::value_type value_type_2;			\
      typedef typename Sacado::Promote<value_type_1,			\
				       value_type_2>::type value_type;  \
									\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
      									\
      std::string name() const {					\
	return expr1.name() + std::string(#OPER) + expr2.name();	\
      }									\
									\
      int size() const {						\
	int sz1 = expr1.size(), sz2 = expr2.size();			\
	return sz1 > sz2 ? sz1 : sz2;					\
      }									\
									\
      bool hasFastAccess(int sz) const {				\
	return expr1.hasFastAccess(sz) && expr2.hasFastAccess(sz);	\
      }									\
									\
      value_type val() const {						\
	return expr1.val() OPER expr2.val();				\
      }									\
      									\
      value_type coeff(int i) const {					\
	return (expr1.coeff(i) OPER expr2.coeff(i));			\
      }									\
									\
      value_type fastAccessCoeff(int i) const {				\
	return (expr1.fastAccessCoeff(i) OPER expr2.fastAccessCoeff(i)); \
      }									\
      									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ExprT2& expr2;						\
									\
    };									\
									\
    template <typename T1>						\
    class Expr< OP< Expr<T1>, typename Expr<T1>::value_type> > {	\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef typename ExprT1::value_type value_type;			\
      typedef typename ExprT1::value_type ConstT;			\
									\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
      Expr(const ExprT1& expr1_, const ConstT& c_) :			\
	expr1(expr1_), c(c_) {}						\
									\
      std::string name() const {					\
	return expr1.name() + std::string(#OPER) + std::string("c");	\
      }									\
									\
      int size() const { return expr1.size(); }				\
									\
      bool hasFastAccess(int sz) const {				\
	return expr1.hasFastAccess(sz);					\
      }									\
									\
      value_type val() const {						\
	return (expr1.val() OPER c);					\
      }									\
      									\
      value_type coeff(int i) const {					\
	return (expr1.coeff(i) OPER c);					\
      }									\
									\
      value_type fastAccessCoeff(int i) const {				\
	return (expr1.fastAccessCoeff(i) OPER c);			\
      }									\
									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ConstT& c;							\
    };									\
									\
    template <typename T2>						\
    class Expr< OP< typename Expr<T2>::value_type, Expr<T2> > > {	\
									\
    public:								\
									\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT2::value_type value_type;			\
      typedef typename ExprT2::value_type ConstT;			\
									\
      typedef typename ExprT2::storage_type storage_type;		\
      typedef typename ExprT2::base_expr_type base_expr_type;           \
									\
      Expr(const ConstT& c_, const ExprT2& expr2_) :			\
	c(c_), expr2(expr2_) {}						\
									\
      std::string name() const {					\
	return std::string("c") + std::string(#OPER) + expr2.name();	\
      }									\
									\
      int size() const { return expr2.size(); }				\
									\
      bool hasFastAccess(int sz) const {				\
	return expr2.hasFastAccess(sz);					\
      }									\
									\
      value_type val() const {						\
	return (c OPER expr2.val());					\
      }									\
      									\
      value_type coeff(int i) const {					\
	return (c OPER expr2.coeff(i));					\
      }									\
									\
      value_type fastAccessCoeff(int i) const {				\
	return (c OPER expr2.fastAccessCoeff(i));			\
      }									\
      									\
    protected:								\
									\
      const ConstT& c;							\
      const ExprT2& expr2;						\
    };									\
									\
    template <typename T1, typename T2>					\
    inline Expr< OP< Expr<T1>, Expr<T2> > >				\
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)		\
    {									\
      typedef OP< Expr<T1>, Expr<T2> > expr_t;				\
    									\
      return Expr<expr_t>(expr1, expr2);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>, Expr<T> > >				\
    OPNAME (const Expr<T>& expr1, const Expr<T>& expr2)			\
    {									\
      typedef OP< Expr<T>, Expr<T> > expr_t;				\
    									\
      return Expr<expr_t>(expr1, expr2);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< typename Expr<T>::value_type, Expr<T> > >		\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< ConstT, Expr<T> > expr_t;				\
									\
      return Expr<expr_t>(c, expr);					\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>, typename Expr<T>::value_type > >		\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< Expr<T>, ConstT > expr_t;				\
									\
      return Expr<expr_t>(expr, c);					\
    }									\
  }									\
}

ETV_BINARYOP_MACRO(operator+, AdditionOp, +)
ETV_BINARYOP_MACRO(operator-, SubtractionOp, -)
ETV_BINARYOP_MACRO(operator*, MultiplicationOp, *)
ETV_BINARYOP_MACRO(operator/, DivisionOp, /)

#undef ETV_BINARYOP_MACRO

#define ETV_BINARYOP_MACRO(OPNAME,OP,OPER)				\
namespace Sacado {							\
  namespace ETV {							\
									\
    template <typename ExprT1, typename ExprT2>				\
    class OP {};							\
									\
    template <typename T1, typename T2>					\
    class Expr< OP< Expr<T1>, Expr<T2> > > {				\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT1::value_type value_type_1;			\
      typedef typename ExprT2::value_type value_type_2;			\
      typedef typename Sacado::Promote<value_type_1,			\
				       value_type_2>::type value_type;  \
									\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
      									\
      std::string name() const {					\
	return expr1.name() + std::string(#OPER) + expr2.name();	\
      }									\
									\
      int size() const {						\
	int sz1 = expr1.size(), sz2 = expr2.size();			\
	return sz1 > sz2 ? sz1 : sz2;					\
      }									\
									\
      bool hasFastAccess(int sz) const {				\
	return expr1.hasFastAccess(sz) && expr2.hasFastAccess(sz);	\
      }									\
									\
      value_type val() const {						\
	return OPER(expr1.val(), expr2.val());				\
      }									\
      									\
      value_type coeff(int i) const {					\
	return OPER(expr1.coeff(i), expr2.coeff(i));			\
      }									\
									\
      value_type fastAccessCoeff(int i) const {				\
	return OPER(expr1.fastAccessCoeff(i), expr2.fastAccessCoeff(i)); \
      }									\
      									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ExprT2& expr2;						\
									\
    };									\
									\
    template <typename T1>						\
    class Expr< OP< Expr<T1>, typename Expr<T1>::value_type> > {	\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef typename ExprT1::value_type value_type;			\
      typedef typename ExprT1::value_type ConstT;			\
									\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
      Expr(const ExprT1& expr1_, const ConstT& c_) :			\
	expr1(expr1_), c(c_) {}						\
									\
      std::string name() const {					\
	return expr1.name() + std::string(#OPER) + std::string("c");	\
      }									\
									\
      int size() const { return expr1.size(); }				\
									\
      bool hasFastAccess(int sz) const {				\
	return expr1.hasFastAccess(sz);					\
      }									\
									\
      value_type val() const {						\
	return OPER(expr1.val(), c);					\
      }									\
      									\
      value_type coeff(int i) const {					\
	return OPER(expr1.coeff(i), c);					\
      }									\
									\
      value_type fastAccessCoeff(int i) const {				\
	return OPER(expr1.fastAccessCoeff(i), c);			\
      }									\
									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ConstT& c;							\
    };									\
									\
    template <typename T2>						\
    class Expr< OP< typename Expr<T2>::value_type, Expr<T2> > > {	\
									\
    public:								\
									\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT2::value_type value_type;			\
      typedef typename ExprT2::value_type ConstT;			\
									\
      typedef typename ExprT2::storage_type storage_type;		\
      typedef typename ExprT2::base_expr_type base_expr_type;           \
									\
      Expr(const ConstT& c_, const ExprT2& expr2_) :			\
	c(c_), expr2(expr2_) {}						\
									\
      std::string name() const {					\
	return std::string("c") + std::string(#OPER) + expr2.name();	\
      }									\
									\
      int size() const { return expr2.size(); }				\
									\
      bool hasFastAccess(int sz) const {				\
	return expr2.hasFastAccess(sz);					\
      }									\
									\
      value_type val() const {						\
	return OPER(c, expr2.val());					\
      }									\
      									\
      value_type coeff(int i) const {					\
	return OPER(c, expr2.coeff(i));					\
      }									\
									\
      value_type fastAccessCoeff(int i) const {				\
	return OPER(c, expr2.fastAccessCoeff(i));			\
      }									\
      									\
    protected:								\
									\
      const ConstT& c;							\
      const ExprT2& expr2;						\
    };									\
									\
    template <typename T1, typename T2>					\
    inline Expr< OP< Expr<T1>, Expr<T2> > >				\
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)		\
    {									\
      typedef OP< Expr<T1>, Expr<T2> > expr_t;				\
    									\
      return Expr<expr_t>(expr1, expr2);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>, Expr<T> > >				\
    OPNAME (const Expr<T>& expr1, const Expr<T>& expr2)			\
    {									\
      typedef OP< Expr<T>, Expr<T> > expr_t;				\
    									\
      return Expr<expr_t>(expr1, expr2);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< typename Expr<T>::value_type, Expr<T> > >		\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< ConstT, Expr<T> > expr_t;				\
									\
      return Expr<expr_t>(c, expr);					\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>, typename Expr<T>::value_type > >		\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< Expr<T>, ConstT > expr_t;				\
									\
      return Expr<expr_t>(expr, c);					\
    }									\
  }									\
}

ETV_BINARYOP_MACRO(atan2, Atan2Op, std::atan2)
ETV_BINARYOP_MACRO(pow, PowerOp, std::pow)
ETV_BINARYOP_MACRO(max, MaxOp, std::max)
ETV_BINARYOP_MACRO(min, MinOp, std::min)

#undef ETV_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define ETV_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace ETV {							\
    template <typename ExprT1, typename ExprT2>				\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return expr1.val() OP expr2.val();				\
    }									\
									\
    template <typename ExprT2>						\
    inline bool								\
    operator OP (const typename Expr<ExprT2>::value_type& a,		\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return a OP expr2.val();						\
    }									\
									\
    template <typename ExprT1>						\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const typename Expr<ExprT1>::value_type& b)		\
    {									\
      return expr1.val() OP b;						\
    }									\
  }									\
}

ETV_RELOP_MACRO(==)
ETV_RELOP_MACRO(!=)
ETV_RELOP_MACRO(<)
ETV_RELOP_MACRO(>)
ETV_RELOP_MACRO(<=)
ETV_RELOP_MACRO(>=)
ETV_RELOP_MACRO(<<=)
ETV_RELOP_MACRO(>>=)
ETV_RELOP_MACRO(&)
ETV_RELOP_MACRO(|)

#undef ETV_RELOP_MACRO

namespace Sacado {

  namespace ETV {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace ETV

} // namespace Sacado


//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace ETV {

    template <typename ExprT>
    bool toBool(const Expr<ExprT>& x) {
      bool is_zero = true;
      for (int i=0; i<x.size(); i++)
	is_zero = is_zero && (x.coeff(i) == 0.0);
      return !is_zero;
    }

  } // namespace ETV

} // namespace Sacado

#define PCE_BOOL_MACRO(OP)						\
namespace Sacado {							\
  namespace ETV {							\
    template <typename ExprT1, typename ExprT2>				\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return toBool(expr1) OP toBool(expr2);				\
    }									\
									\
    template <typename ExprT2>						\
    inline bool								\
    operator OP (const typename Expr<ExprT2>::value_type& a,		\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return a OP toBool(expr2);					\
    }									\
									\
    template <typename ExprT1>						\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const typename Expr<ExprT1>::value_type& b)		\
    {									\
      return toBool(expr1) OP b;					\
    }									\
  }									\
}

PCE_BOOL_MACRO(&&)
PCE_BOOL_MACRO(||)

#undef PCE_BOOL_MACRO


//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace ETV {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      typedef typename Expr<ExprT>::value_type value_type;
      typedef typename Expr<ExprT>::storage_type storage_type;
      Vector<value_type, storage_type> a(x);
      os << a;
      return os;
    }

  } // namespace ETV

} // namespace Sacado

#endif // SACADO_ETV_VECTOR_OPS_HPP
