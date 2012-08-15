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

#ifndef SACADO_CACHEFAD_OPS_HPP
#define SACADO_CACHEFAD_OPS_HPP

#include "Sacado_CacheFad_Expression.hpp"
#include "Sacado_cmath.hpp"
#include "Sacado_dummy_arg.hpp"
#include <ostream>	// for std::ostream

namespace Sacado {
  namespace CacheFad {

    //
    // UnaryPlusOp
    //

    template <typename ExprT>
    class UnaryPlusOp {};

    template <typename ExprT>
    class Expr< UnaryPlusOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;

      typedef typename ExprT::scalar_type scalar_type;

      Expr(const ExprT& expr_) : expr(expr_)  {}

      int size() const { return expr.size(); }

      bool updateValue() const { return expr.updateValue(); }

      void cache() const {
	expr.cache();
      }

      value_type val() const {
	return expr.val();
      }

      bool isLinear() const {
        return expr.isLinear();
      }

      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr.dx(i);
      }
 
      const value_type fastAccessDx(int i) const {
        return expr.fastAccessDx(i);
      }

    protected:

      const ExprT& expr;
    };

    template <typename T>
    inline Expr< UnaryPlusOp< Expr<T> > >
    operator+ (const Expr<T>& expr)
    {
      typedef UnaryPlusOp< Expr<T> > expr_t;
      
      return Expr<expr_t>(expr);
    }

    //
    // UnaryMinusOp
    //
    template <typename ExprT>
    class UnaryMinusOp {};

    template <typename ExprT>
    class Expr< UnaryMinusOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;

      typedef typename ExprT::scalar_type scalar_type;

      Expr(const ExprT& expr_) : expr(expr_)  {}

      int size() const { return expr.size(); }

      bool updateValue() const { return expr.updateValue(); }

      void cache() const {
	expr.cache();
      }

      value_type val() const {
	return -expr.val();
      }

      bool isLinear() const {
        return expr.isLinear();
      }

      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      const value_type dx(int i) const {
        return -expr.dx(i);
      }
 
      const value_type fastAccessDx(int i) const {
        return -expr.fastAccessDx(i);
      }

    protected:

      const ExprT& expr;
    };

    template <typename T>
    inline Expr< UnaryMinusOp< Expr<T> > >
    operator- (const Expr<T>& expr)
    {
      typedef UnaryMinusOp< Expr<T> > expr_t;
      
      return Expr<expr_t>(expr);
    }

    //
    // AbsOp
    //

    template <typename ExprT>
    class AbsOp {};

    template <typename ExprT>
    class Expr< AbsOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;

      typedef typename ExprT::scalar_type scalar_type;

      Expr(const ExprT& expr_) : expr(expr_)  {}

      int size() const { return expr.size(); }

      bool updateValue() const { return expr.updateValue(); }

      void cache() const {
	expr.cache();
	v = expr.val();
	v_pos = (v >= 0);
      }

      value_type val() const {
	return std::abs(v);
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      const value_type dx(int i) const {
        return v_pos ? expr.dx(i) : value_type(-expr.dx(i));
      }
 
      const value_type fastAccessDx(int i) const {
        return v_pos ? expr.fastAccessDx(i) : value_type(-expr.fastAccessDx(i));
      }

    protected:

      const ExprT& expr;
      mutable value_type v;
      mutable bool v_pos;
    };

    template <typename T>
    inline Expr< AbsOp< Expr<T> > >
    abs (const Expr<T>& expr)
    {
      typedef AbsOp< Expr<T> > expr_t;
      
      return Expr<expr_t>(expr);
    }

    //
    // FAbsOp
    //

    template <typename ExprT>
    class FAbsOp {};

    template <typename ExprT>
    class Expr< FAbsOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;

      typedef typename ExprT::scalar_type scalar_type;

      Expr(const ExprT& expr_) : expr(expr_)  {}

      int size() const { return expr.size(); }

      bool updateValue() const { return expr.updateValue(); }

      void cache() const {
	expr.cache();
	v = expr.val();
	v_pos = (v >= 0);
      }

      value_type val() const {
	return std::fabs(v);
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      const value_type dx(int i) const {
        return v_pos ? expr.dx(i) : value_type(-expr.dx(i));
      }
 
      const value_type fastAccessDx(int i) const {
        return v_pos ? expr.fastAccessDx(i) : value_type(-expr.fastAccessDx(i));
      }

    protected:

      const ExprT& expr;
      mutable value_type v;
      mutable bool v_pos;
    };

    template <typename T>
    inline Expr< FAbsOp< Expr<T> > >
    fabs (const Expr<T>& expr)
    {
      typedef FAbsOp< Expr<T> > expr_t;
      
      return Expr<expr_t>(expr);
    }

  }
}

#define FAD_UNARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE)			\
namespace Sacado {							\
  namespace CacheFad {							\
									\
    template <typename ExprT>						\
    class OP {};							\
									\
    template <typename ExprT>						\
    class Expr< OP<ExprT> > {						\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
      typedef typename ExprT::scalar_type scalar_type;			\
									\
      Expr(const ExprT& expr_) : expr(expr_)  {}			\
									\
      int size() const { return expr.size(); }				\
									\
      bool hasFastAccess() const { return expr.hasFastAccess(); }	\
									\
      bool isPassive() const { return expr.isPassive();}		\
									\
      bool updateValue() const { return expr.updateValue(); }		\
									\
      void cache() const {						\
	expr.cache();							\
        v = expr.val();							\
	PARTIAL;							\
      }									\
									\
      value_type val() const {						\
	return VALUE;							\
      }									\
									\
      value_type dx(int i) const {					\
	return expr.dx(i)*a;						\
      }									\
									\
      value_type fastAccessDx(int i) const {				\
	return expr.fastAccessDx(i)*a;					\
      }									\
									\
    protected:								\
									\
      const ExprT& expr;						\
      mutable value_type v;						\
      mutable value_type a;						\
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

FAD_UNARYOP_MACRO(exp,
		  ExpOp, 
		  a = std::exp(v),
		  a)
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  a=value_type(1)/v,
		  std::log(v))
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  a = value_type(1)/(std::log(value_type(10))*v),
		  std::log10(v))
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp, 
		  a = value_type(1)/(value_type(2)*std::sqrt(v)),
		  std::sqrt(v))
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  a = -std::sin(v),
		  std::cos(v))
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  a = std::cos(v),
		  std::sin(v))
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  a = value_type(1)+std::tan(v)*std::tan(v),
		  std::tan(v))
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  a = value_type(-1)/std::sqrt(value_type(1)-v*v),
		  std::acos(v))
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  a = value_type(1)/std::sqrt(value_type(1)-v*v),
		  std::asin(v))
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  a = value_type(1)/(value_type(1)+v*v),
		  std::atan(v))
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  a = std::sinh(v),
		  std::cosh(v))
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  a = std::cosh(v),
		  std::sinh(v))
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  a = value_type(1)/(std::cosh(v)*std::cosh(v)),
		  std::tanh(v))
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  a = value_type(1)/std::sqrt((v-value_type(1))*(v+value_type(1))),
		  std::acosh(v))
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  a = value_type(1)/std::sqrt(value_type(1)+v*v),
		  std::asinh(v))
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  a = value_type(1)/(value_type(1)-v*v),
		  std::atanh(v))

#undef FAD_UNARYOP_MACRO

//
// Binary operators
//
namespace Sacado {
  namespace CacheFad {

    //
    // AdditionOp
    //

    template <typename ExprT1, typename ExprT2>
    class AdditionOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< AdditionOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
      }

      value_type val() const {
	return expr1.val()+expr2.val();
      }

      bool isLinear() const {
        return expr1.isLinear() && expr2.isLinear();
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i) + expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i) + expr2.fastAccessDx(i);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< AdditionOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
      }

      value_type val() const {
	return expr1.val() + expr2.val();
      }

      bool isLinear() const {
        return expr1.isLinear();
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;

    };

    template <typename T1, typename ExprT2>
    class Expr< AdditionOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
      }

      value_type val() const {
	return expr1.val() + expr2.val();
      }

      bool isLinear() const {
        return expr2.isLinear();
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;

    };

    //
    // SubtractionOp
    //

    template <typename ExprT1, typename ExprT2>
    class SubtractionOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< SubtractionOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
      }

      value_type val() const {
	return expr1.val()-expr2.val();
      }

      bool isLinear() const {
        return expr1.isLinear() && expr2.isLinear();
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i) - expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i) - expr2.fastAccessDx(i);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< SubtractionOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
      }

      value_type val() const {
	return expr1.val() - expr2.val();
      }

      bool isLinear() const {
        return expr1.isLinear();
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;

    };

    template <typename T1, typename ExprT2>
    class Expr< SubtractionOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
      }

      value_type val() const {
	return expr1.val() - expr2.val();
      }

      bool isLinear() const {
        return expr2.isLinear();
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return -expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return -expr2.fastAccessDx(i);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;

    };

    //
    // MultiplicationOp
    //

    template <typename ExprT1, typename ExprT2>
    class MultiplicationOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< MultiplicationOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
        v1 = expr1.val();
	v2 = expr2.val();
      }

      value_type val() const {
	return v1*v2;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
	  return v1*expr2.dx(i) + expr1.dx(i)*v2;
	else if (expr1.size() > 0)
	  return expr1.dx(i)*v2;
	else
	  return v1*expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return v1*expr2.fastAccessDx(i) + expr1.fastAccessDx(i)*v2;
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type_1 v1;
      mutable value_type_2 v2;

    };

    template <typename ExprT1, typename T2>
    class Expr< MultiplicationOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
      }

      value_type val() const {
	return expr1.val()*expr2.val();
      }

      bool isLinear() const {
        return expr1.isLinear();
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i)*expr2.val();
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*expr2.val();
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;

    };

    template <typename T1, typename ExprT2>
    class Expr< MultiplicationOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
      }

      value_type val() const {
	return expr1.val()*expr2.val();
      }

      bool isLinear() const {
        return expr2.isLinear();
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.val()*expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return expr1.val()*expr2.fastAccessDx(i);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;

    };

    //
    // DivisionOp
    //

    template <typename ExprT1, typename ExprT2>
    class DivisionOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< DivisionOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
        const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	a = value_type(1)/v2;
	v = v1*a;
	b = -v/v2;
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
	if (expr1.size() > 0 && expr2.size() > 0)
	  return expr1.dx(i)*a + expr2.dx(i)*b;
	else if (expr1.size() > 0)
	  return expr1.dx(i)*a;
	else
	  return expr1.val()*b;
      }

      const value_type fastAccessDx(int i) const {
        return  expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b;
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type a;
      mutable value_type b;

    };

    template <typename ExprT1, typename T2>
    class Expr< DivisionOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
	const value_type_1 v1 = expr1.val();
	a = value_type(1)/expr2.val();
	v = v1*a;
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return expr1.isLinear();
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i)*a;
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a;
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable value_type a;

    };

    template <typename T1, typename ExprT2>
    class Expr< DivisionOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
	const value_type_2 v2 = expr2.val();
	v = expr1.val()/v2;
	b = -v/v2;
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr2.dx(i)*b;
      }

      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i)*b;
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type b;

    };

    //
    // Atan2Op
    //

    template <typename ExprT1, typename ExprT2>
    class Atan2Op {};

    template <typename ExprT1, typename ExprT2>
    class Expr< Atan2Op<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
        const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	a = value_type(1)/(v1*v1 + v2*v2);
	b = -v1*a;
	a = v2*a;
	v = std::atan2(v1,v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
	  return expr1.dx(i)*a + expr2.dx(i)*b;
	else if (expr1.size() > 0)
	  return expr1.dx(i)*a;
	else
	  return expr1.val()*b;
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b;
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type a;
      mutable value_type b;

    };

    template <typename ExprT1, typename T2>
    class Expr< Atan2Op<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	a = v2/(v1*v1 + v2*v2);
	v = std::atan2(v1,v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i)*a;
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a;
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable value_type a;

    };

    template <typename T1, typename ExprT2>
    class Expr< Atan2Op< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	b = -v1/(v1*v1 + v2*v2);
	v = std::atan2(v1,v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr2.dx(i)*b;
      }

      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i)*b;
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type b;

    };

    //
    // PowerOp
    //

    template <typename ExprT1, typename ExprT2>
    class PowerOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< PowerOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
        const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::pow(v1,v2);
	if (v1 == value_type(0)) {
	  a = value_type(0);
	  b = value_type(0);
	}
	else {
	  a = v*v2/v1;
	  b = v*std::log(v1);
	}
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
	  return expr1.dx(i)*a + expr2.dx(i)*b;
	else if (expr1.size() > 0)
	  return expr1.dx(i)*a;
	else
	  return expr1.val()*b;
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b;
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type a;
      mutable value_type b;

    };

    template <typename ExprT1, typename T2>
    class Expr< PowerOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::pow(v1,v2);
	if (v1 == value_type(0)) {
	  a = value_type(0);
	}
	else {
	  a = v*v2/v1;
	}
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr1.dx(i)*a;
      }

      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a;
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable value_type a;

    };

    template <typename T1, typename ExprT2>
    class Expr< PowerOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::pow(v1,v2);
	if (v1 == value_type(0)) {
	  b = value_type(0);
	}
	else {
	  b = v*std::log(v1);
	}
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return expr2.dx(i)*b;
      }

      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i)*b;
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type b;

    };

    //
    // MaxOp
    //

    template <typename ExprT1, typename ExprT2>
    class MaxOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< MaxOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
        const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::max(v1,v2);
	max_v1 = (v1 >= v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return max_v1 ? expr1.dx(i) : expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return max_v1 ? expr1.fastAccessDx(i) : expr2.fastAccessDx(i);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool max_v1;

    };

    template <typename ExprT1, typename T2>
    class Expr< MaxOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::max(v1,v2);
	max_v1 = (v1 >= v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return max_v1 ? expr1.dx(i) : value_type(0);
      }

      const value_type fastAccessDx(int i) const {
        return max_v1 ? expr1.fastAccessDx(i) : value_type(0);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable bool max_v1;

    };

    template <typename T1, typename ExprT2>
    class Expr< MaxOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::max(v1,v2);
	max_v1 = (v1 >= v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return max_v1 ? value_type(0) : expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return max_v1 ? value_type(0) : expr2.fastAccessDx(i);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool max_v1;

    };

    //
    // MinOp
    //

    template <typename ExprT1, typename ExprT2>
    class MinOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< MinOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	int sz1 = expr1.size(), sz2 = expr2.size();
	return sz1 > sz2 ? sz1 : sz2;
      }

      bool updateValue() const {
	return expr1.updateValue() && expr2.updateValue();
      }

      void cache() const {
	expr1.cache();
	expr2.cache();
        const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::min(v1,v2);
	min_v1 = (v1 <= v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return min_v1 ? expr1.dx(i) : expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return min_v1 ? expr1.fastAccessDx(i) : expr2.fastAccessDx(i);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool min_v1;

    };

    template <typename ExprT1, typename T2>
    class Expr< MinOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr1.size();
      }

      bool updateValue() const {
	return expr1.updateValue();
      }

      void cache() const {
	expr1.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::min(v1,v2);
	min_v1 = (v1 <= v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      const value_type dx(int i) const {
        return min_v1 ? expr1.dx(i) : value_type(0);
      }

      const value_type fastAccessDx(int i) const {
        return min_v1 ? expr1.fastAccessDx(i) : value_type(0);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable bool min_v1;

    };

    template <typename T1, typename ExprT2>
    class Expr< MinOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
				   scalar_type_2>::type scalar_type;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      int size() const {
	return expr2.size();
      }

      bool updateValue() const {
	return expr2.updateValue();
      }

      void cache() const {
	expr2.cache();
	const value_type_1 v1 = expr1.val();
	const value_type_2 v2 = expr2.val();
	v = std::min(v1,v2);
	min_v1 = (v1 <= v2);
      }

      value_type val() const {
	return v;
      }

      bool isLinear() const {
        return false;
      }

      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      const value_type dx(int i) const {
        return min_v1 ? value_type(0) : expr2.dx(i);
      }

      const value_type fastAccessDx(int i) const {
        return min_v1 ? value_type(0) : expr2.fastAccessDx(i);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool min_v1;

    };

  }

}

#define FAD_BINARYOP_MACRO(OPNAME,OP)					\
namespace Sacado {							\
  namespace CacheFad {							\
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
    inline Expr< OP< ConstExpr<typename Expr<T>::value_type>,		\
		     Expr<T> > >					\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef OP< ConstT, Expr<T> > expr_t;				\
									\
      return Expr<expr_t>(ConstT(c), expr);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>,						\
		     ConstExpr<typename Expr<T>::value_type> > >	\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef OP< Expr<T>, ConstT > expr_t;				\
									\
      return Expr<expr_t>(expr, ConstT(c));				\
    }									\
  }									\
}


FAD_BINARYOP_MACRO(operator+, AdditionOp)
FAD_BINARYOP_MACRO(operator-, SubtractionOp)
FAD_BINARYOP_MACRO(operator*, MultiplicationOp)
FAD_BINARYOP_MACRO(operator/, DivisionOp)
FAD_BINARYOP_MACRO(atan2, Atan2Op)
FAD_BINARYOP_MACRO(pow, PowerOp)
FAD_BINARYOP_MACRO(max, MaxOp)
FAD_BINARYOP_MACRO(min, MinOp)

#undef FAD_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace CacheFad {							\
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

FAD_RELOP_MACRO(==)
FAD_RELOP_MACRO(!=)
FAD_RELOP_MACRO(<)
FAD_RELOP_MACRO(>)
FAD_RELOP_MACRO(<=)
FAD_RELOP_MACRO(>=)
FAD_RELOP_MACRO(<<=)
FAD_RELOP_MACRO(>>=)
FAD_RELOP_MACRO(&)
FAD_RELOP_MACRO(|)

#undef FAD_RELOP_MACRO

namespace Sacado {

  namespace CacheFad {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace CacheFad

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace CacheFad {

    template <typename ExprT>
    bool toBool(const Expr<ExprT>& x) {
      bool is_zero = (x.val() == 0.0);
      for (int i=0; i<x.size(); i++)
	is_zero = is_zero && (x.dx(i) == 0.0);
      return !is_zero;
    }

  } // namespace Fad

} // namespace Sacado

#define FAD_BOOL_MACRO(OP)						\
namespace Sacado {							\
  namespace CacheFad {							\
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

FAD_BOOL_MACRO(&&)
FAD_BOOL_MACRO(||)

#undef FAD_BOOL_MACRO

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace CacheFad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os << x.val() << " [";
      
      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

  } // namespace CacheFad

} // namespace Sacado


#endif // SACADO_CACHEFAD_OPS_HPP
