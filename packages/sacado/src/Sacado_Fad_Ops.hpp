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

#ifndef SACADO_FAD_OPS_HPP
#define SACADO_FAD_OPS_HPP

#include "Sacado_Fad_Expression.hpp"
#include "Sacado_cmath.hpp"
#include <ostream>	// for std::ostream

#define FAD_UNARYOP_MACRO(OPNAME,OP,VALUE,DX,FASTACCESSDX)		\
namespace Sacado {							\
  namespace Fad {							\
									\
    template <typename ExprT>						\
    class OP {};							\
									\
    template <typename ExprT>						\
    class Expr< OP<ExprT> > {						\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
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
      value_type val() const {						\
	return VALUE;							\
      }									\
									\
      value_type dx(int i) const {					\
	return DX;							\
      }									\
									\
      value_type fastAccessDx(int i) const {				\
	return FASTACCESSDX;						\
      }									\
									\
    protected:								\
									\
      const ExprT& expr;						\
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

FAD_UNARYOP_MACRO(operator+,
		  UnaryPlusOp, 
		  expr.val(),
		  expr.dx(i),
		  expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(operator-,
		  UnaryMinusOp, 
		  -expr.val(),
		  -expr.dx(i),
		  -expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(exp,
		  ExpOp, 
		  std::exp(expr.val()),
		  std::exp(expr.val())*expr.dx(i),
		  std::exp(expr.val())*expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  std::log(expr.val()),
		  expr.dx(i)/expr.val(),
		  expr.fastAccessDx(i)/expr.val())
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  std::log10(expr.val()),
		  expr.dx(i)/( std::log(value_type(10))*expr.val()),
		  expr.fastAccessDx(i) / ( std::log(value_type(10))*expr.val()))
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp, 
		  std::sqrt(expr.val()),
		  expr.dx(i)/(value_type(2)* std::sqrt(expr.val())),
		  expr.fastAccessDx(i)/(value_type(2)* std::sqrt(expr.val())))
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  std::cos(expr.val()),
		  -expr.dx(i)* std::sin(expr.val()),
		  -expr.fastAccessDx(i)* std::sin(expr.val()))
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  std::sin(expr.val()),
		  expr.dx(i)* std::cos(expr.val()),
		  expr.fastAccessDx(i)* std::cos(expr.val()))
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  std::tan(expr.val()),
		  expr.dx(i)*
		    (value_type(1)+ std::tan(expr.val())* std::tan(expr.val())),
		  expr.fastAccessDx(i)*
		    (value_type(1)+ std::tan(expr.val())* std::tan(expr.val())))
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  std::acos(expr.val()),
		  -expr.dx(i)/ std::sqrt(value_type(1)-expr.val()*expr.val()),
		  -expr.fastAccessDx(i) /
		    std::sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  std::asin(expr.val()),
		  expr.dx(i)/ std::sqrt(value_type(1)-expr.val()*expr.val()),
		  expr.fastAccessDx(i) /
		    std::sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  std::atan(expr.val()),
		  expr.dx(i)/(value_type(1)+expr.val()*expr.val()),
		  expr.fastAccessDx(i)/(value_type(1)+expr.val()*expr.val()))
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  std::cosh(expr.val()),
		  expr.dx(i)* std::sinh(expr.val()),
		  expr.fastAccessDx(i)* std::sinh(expr.val()))
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  std::sinh(expr.val()),
		  expr.dx(i)* std::cosh(expr.val()),
		  expr.fastAccessDx(i)* std::cosh(expr.val()))
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  std::tanh(expr.val()),
		  expr.dx(i)/( std::cosh(expr.val())* std::cosh(expr.val())),
		  expr.fastAccessDx(i) / 
		    ( std::cosh(expr.val())* std::cosh(expr.val())))
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  std::acosh(expr.val()),
		  expr.dx(i)/ std::sqrt((expr.val()-value_type(1)) * 
				       (expr.val()+value_type(1))),
		  expr.fastAccessDx(i)/ std::sqrt((expr.val()-value_type(1)) * 
						 (expr.val()+value_type(1))))
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  std::asinh(expr.val()),
		  expr.dx(i)/ std::sqrt(value_type(1)+expr.val()*expr.val()),
		  expr.fastAccessDx(i)/ std::sqrt(value_type(1)+
						 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  std::atanh(expr.val()),
		  expr.dx(i)/(value_type(1)-expr.val()*expr.val()),
		  expr.fastAccessDx(i)/(value_type(1)-
						 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(abs,
		  AbsOp, 
		  std::abs(expr.val()),
		  expr.val() >= 0 ? value_type(+expr.dx(i)) : 
		    value_type(-expr.dx(i)),
		  expr.val() >= 0 ? value_type(+expr.fastAccessDx(i)) : 
		    value_type(-expr.fastAccessDx(i)))
FAD_UNARYOP_MACRO(fabs,
		  FAbsOp, 
		  std::fabs(expr.val()),
		  expr.val() >= 0 ? value_type(+expr.dx(i)) : 
		    value_type(-expr.dx(i)),
		  expr.val() >= 0 ? value_type(+expr.fastAccessDx(i)) : 
		    value_type(-expr.fastAccessDx(i)))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,VALUE,DX,FASTACCESSDX,VAL_CONST_DX_1,VAL_CONST_DX_2,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {							\
  namespace Fad {							\
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
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      int size() const {						\
	int sz1 = expr1.size(), sz2 = expr2.size();			\
	return sz1 > sz2 ? sz1 : sz2;					\
      }									\
									\
      bool hasFastAccess() const {					\
	return expr1.hasFastAccess() && expr2.hasFastAccess();		\
      }									\
									\
      bool isPassive() const {						\
	return expr1.isPassive() && expr2.isPassive();			\
      }									\
									\
      bool updateValue() const {					\
	return expr1.updateValue() && expr2.updateValue();		\
      }									\
									\
      const value_type val() const {					\
	return VALUE;							\
      }									\
									\
      const value_type dx(int i) const {				\
	return DX;							\
      }									\
									\
      const value_type fastAccessDx(int i) const {			\
	return FASTACCESSDX;						\
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
      Expr(const ExprT1& expr1_, const ConstT& c_) :			\
	expr1(expr1_), c(c_) {}						\
									\
      int size() const {						\
	return expr1.size();						\
      }									\
									\
      bool hasFastAccess() const {					\
	return expr1.hasFastAccess();					\
      }									\
									\
      bool isPassive() const {						\
	return expr1.isPassive();					\
      }									\
									\
      bool updateValue() const { return expr1.updateValue(); }		\
									\
      const value_type val() const {					\
	return VAL_CONST_DX_2;						\
      }									\
									\
      const value_type dx(int i) const {				\
	return CONST_DX_2;						\
      }									\
									\
      const value_type fastAccessDx(int i) const {			\
	return CONST_FASTACCESSDX_2;					\
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
      Expr(const ConstT& c_, const ExprT2& expr2_) :			\
	c(c_), expr2(expr2_) {}						\
									\
      int size() const {						\
	return expr2.size();						\
      }									\
									\
      bool hasFastAccess() const {					\
	return expr2.hasFastAccess();					\
      }									\
									\
      bool isPassive() const {						\
	return expr2.isPassive();					\
      }									\
									\
      bool updateValue() const { return expr2.updateValue(); }		\
									\
      const value_type val() const {					\
	return VAL_CONST_DX_1;						\
      }									\
									\
      const value_type dx(int i) const {				\
	return CONST_DX_1;						\
      }									\
									\
      const value_type fastAccessDx(int i) const {			\
	return CONST_FASTACCESSDX_1;					\
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


FAD_BINARYOP_MACRO(operator+,
		   AdditionOp, 
		   expr1.val() + expr2.val(),
		   expr1.dx(i) + expr2.dx(i),
		   expr1.fastAccessDx(i) + expr2.fastAccessDx(i),
		   c + expr2.val(),
		   expr1.val() + c,
		   expr2.dx(i),
		   expr1.dx(i),
		   expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator-,
		   SubtractionOp, 
		   expr1.val() - expr2.val(),
		   expr1.dx(i) - expr2.dx(i),
		   expr1.fastAccessDx(i) - expr2.fastAccessDx(i),
		   c - expr2.val(),
		   expr1.val() - c,
		   -expr2.dx(i),
		   expr1.dx(i),
		   -expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i))
// FAD_BINARYOP_MACRO(operator*,
// 		   MultiplicationOp, 
// 		   expr1.val() * expr2.val(),
// 		   expr1.val()*expr2.dx(i) + expr1.dx(i)*expr2.val(),
// 		   expr1.val()*expr2.fastAccessDx(i) + 
// 		     expr1.fastAccessDx(i)*expr2.val(),
// 		   c * expr2.val(),
// 		   expr1.val() * c,
// 		   c*expr2.dx(i),
// 		   expr1.dx(i)*c,
// 		   c*expr2.fastAccessDx(i),
// 		   expr1.fastAccessDx(i)*c)
FAD_BINARYOP_MACRO(operator/,
		   DivisionOp, 
		   expr1.val() / expr2.val(),
		   (expr1.dx(i)*expr2.val() - expr2.dx(i)*expr1.val()) /
		     (expr2.val()*expr2.val()),
		   (expr1.fastAccessDx(i)*expr2.val() - 
		      expr2.fastAccessDx(i)*expr1.val()) /
		      (expr2.val()*expr2.val()),
		   c / expr2.val(),
		   expr1.val() / c,
		   -expr2.dx(i)*c / (expr2.val()*expr2.val()),
		   expr1.dx(i)/c,
		   -expr2.fastAccessDx(i)*c / (expr2.val()*expr2.val()),
		   expr1.fastAccessDx(i)/c)
FAD_BINARYOP_MACRO(atan2,
		   Atan2Op,
		   std::atan2(expr1.val(), expr2.val()),
		   (expr2.val()*expr1.dx(i) - expr1.val()*expr2.dx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   (expr2.val()*expr1.fastAccessDx(i) - expr1.val()*expr2.fastAccessDx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   std::atan2(c, expr2.val()),
		   std::atan2(expr1.val(), c),
		   (-c*expr2.dx(i)) / (c*c + expr2.val()*expr2.val()),
		   (c*expr1.dx(i))/ (expr1.val()*expr1.val() + c*c),
		   (-c*expr2.fastAccessDx(i))/ (c*c + expr2.val()*expr2.val()),
		   (c*expr1.fastAccessDx(i))/ (expr1.val()*expr1.val() + c*c))
FAD_BINARYOP_MACRO(pow,
		   PowerOp,
		   std::pow(expr1.val(), expr2.val()),
		   expr1.val() == value_type(0) ? value_type(0) : value_type((expr2.dx(i)*std::log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*std::pow(expr1.val(),expr2.val())),
		   expr1.val() == value_type(0) ? value_type(0.0) : value_type((expr2.fastAccessDx(i)*std::log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*std::pow(expr1.val(),expr2.val())),
		   std::pow(c, expr2.val()),
		   std::pow(expr1.val(), c),
		   c == value_type(0) ? value_type(0) : value_type(expr2.dx(i)*std::log(c)*std::pow(c,expr2.val())),
		   expr1.val() == value_type(0) ? value_type(0.0) : value_type(c*expr1.dx(i)/expr1.val()*std::pow(expr1.val(),c)),
		   c == value_type(0) ? value_type(0) : value_type(expr2.fastAccessDx(i)*std::log(c)*std::pow(c,expr2.val())),
		   expr1.val() == value_type(0) ? value_type(0.0) : value_type(c*expr1.fastAccessDx(i)/expr1.val()*std::pow(expr1.val(),c)))
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   std::max(expr1.val(), expr2.val()),
                   expr1.val() >= expr2.val() ? expr1.dx(i) : expr2.dx(i),
                   expr1.val() >= expr2.val() ? expr1.fastAccessDx(i) : 
                                                expr2.fastAccessDx(i),
		   std::max(c, expr2.val()),
		   std::max(expr1.val(), c),
                   c >= expr2.val() ? value_type(0) : expr2.dx(i),
                   expr1.val() >= c ? expr1.dx(i) : value_type(0),
                   c >= expr2.val() ? value_type(0) : expr2.fastAccessDx(i),
                   expr1.val() >= c ? expr1.fastAccessDx(i) : value_type(0))
FAD_BINARYOP_MACRO(min,
                   MinOp,
                   std::min(expr1.val(), expr2.val()),
                   expr1.val() <= expr2.val() ? expr1.dx(i) : expr2.dx(i),
                   expr1.val() <= expr2.val() ? expr1.fastAccessDx(i) : 
                                                expr2.fastAccessDx(i),
		   std::min(c, expr2.val()),
		   std::min(expr1.val(), c),
                   c <= expr2.val() ? value_type(0) : expr2.dx(i),
                   expr1.val() <= c ? expr1.dx(i) : value_type(0),
                   c <= expr2.val() ? value_type(0) : expr2.fastAccessDx(i),
                   expr1.val() <= c ? expr1.fastAccessDx(i) : value_type(0))


#undef FAD_BINARYOP_MACRO

namespace Sacado {							
  namespace Fad {							
									
    template <typename ExprT1, typename ExprT2>				
    class MultiplicationOp {};							
									
    template <typename T1, typename T2>					
    class Expr< MultiplicationOp< Expr<T1>, Expr<T2> > > {				
									
    public:								
									
      typedef Expr<T1> ExprT1;						
      typedef Expr<T2> ExprT2;						
      typedef typename ExprT1::value_type value_type_1;			
      typedef typename ExprT2::value_type value_type_2;			
      typedef typename Sacado::Promote<value_type_1,			
				       value_type_2>::type value_type;  
									
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		
	expr1(expr1_), expr2(expr2_) {}					
									
      int size() const {						
	int sz1 = expr1.size(), sz2 = expr2.size();			
	return sz1 > sz2 ? sz1 : sz2;					
      }									
									
      bool hasFastAccess() const {					
	return expr1.hasFastAccess() && expr2.hasFastAccess();		
      }									
									
      bool isPassive() const {						
	return expr1.isPassive() && expr2.isPassive();			
      }

      bool updateValue() const { 
	return expr1.updateValue() && expr2.updateValue(); 
      }

      const value_type val() const {					
	return expr1.val()*expr2.val();
      }									
									
      const value_type dx(int i) const {
	if (expr1.size() > 0 && expr2.size() > 0)
	  return expr1.val()*expr2.dx(i) + expr1.dx(i)*expr2.val();
	else if (expr1.size() > 0)
	  return expr1.dx(i)*expr2.val();
	else
	  return expr1.val()*expr2.dx(i);
      }									
									
      const value_type fastAccessDx(int i) const {			
	return expr1.val()*expr2.fastAccessDx(i) + 
	  expr1.fastAccessDx(i)*expr2.val();
      }									
      									
    protected:								
									
      const ExprT1& expr1;						
      const ExprT2& expr2;						
									
    };									
									
    template <typename T1>						
    class Expr< MultiplicationOp< Expr<T1>, typename Expr<T1>::value_type> > {	
									
    public:								
									
      typedef Expr<T1> ExprT1;						
      typedef typename ExprT1::value_type value_type;			
      typedef typename ExprT1::value_type ConstT;			
									
      Expr(const ExprT1& expr1_, const ConstT& c_) :			
	expr1(expr1_), c(c_) {}						
									
      int size() const {						
	return expr1.size();						
      }									
									
      bool hasFastAccess() const {					
	return expr1.hasFastAccess();					
      }									
									
      bool isPassive() const {						
	return expr1.isPassive();					
      }

      bool updateValue() const { return expr1.updateValue(); }
									
      const value_type val() const {					
	return expr1.val()*c;						
      }									
									
      const value_type dx(int i) const {				
	return expr1.dx(i)*c;						
      }									
									
      const value_type fastAccessDx(int i) const {			
	return expr1.fastAccessDx(i)*c;					
      }									
									
    protected:								
									
      const ExprT1& expr1;						
      const ConstT& c;							
    };									
									
    template <typename T2>						
    class Expr< MultiplicationOp< typename Expr<T2>::value_type, Expr<T2> > > {	
									
    public:								
									
      typedef Expr<T2> ExprT2;						
      typedef typename ExprT2::value_type value_type;			
      typedef typename ExprT2::value_type ConstT;			
									
      Expr(const ConstT& c_, const ExprT2& expr2_) :			
	c(c_), expr2(expr2_) {}						
									
      int size() const {						
	return expr2.size();						
      }									
									
      bool hasFastAccess() const {					
	return expr2.hasFastAccess();					
      }									
									
      bool isPassive() const {						
	return expr2.isPassive();					
      }

      bool updateValue() const { return expr2.updateValue(); }
									
      const value_type val() const {					
	return c*expr2.val();						
      }									
									
      const value_type dx(int i) const {				
	return c*expr2.dx(i);						
      }									
									
      const value_type fastAccessDx(int i) const {			
	return c*expr2.fastAccessDx(i);					
      }									
      									
    protected:								
									
      const ConstT& c;							
      const ExprT2& expr2;						
    };									
									
    template <typename T1, typename T2>					
    inline Expr< MultiplicationOp< Expr<T1>, Expr<T2> > >				
    operator* (const Expr<T1>& expr1, const Expr<T2>& expr2)		
    {									
      typedef MultiplicationOp< Expr<T1>, Expr<T2> > expr_t;				
    									
      return Expr<expr_t>(expr1, expr2);				
    }									
									
    template <typename T>						
    inline Expr< MultiplicationOp< Expr<T>, Expr<T> > >				
    operator* (const Expr<T>& expr1, const Expr<T>& expr2)			
    {									
      typedef MultiplicationOp< Expr<T>, Expr<T> > expr_t;				
    									
      return Expr<expr_t>(expr1, expr2);				
    }									
									
    template <typename T>						
    inline Expr< MultiplicationOp< typename Expr<T>::value_type, Expr<T> > >		
    operator* (const typename Expr<T>::value_type& c,			
	    const Expr<T>& expr)					
    {									
      typedef typename Expr<T>::value_type ConstT;			
      typedef MultiplicationOp< ConstT, Expr<T> > expr_t;				
									
      return Expr<expr_t>(c, expr);					
    }									
									
    template <typename T>						
    inline Expr< MultiplicationOp< Expr<T>, typename Expr<T>::value_type > >		
    operator* (const Expr<T>& expr,					
	    const typename Expr<T>::value_type& c)			
    {									
      typedef typename Expr<T>::value_type ConstT;			
      typedef MultiplicationOp< Expr<T>, ConstT > expr_t;				
									
      return Expr<expr_t>(expr, c);					
    }									
  }									
}

//-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace Fad {							\
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

  namespace Fad {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace Fad

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace Fad {

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
  namespace Fad {							\
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

  namespace Fad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os << x.val() << " [";
      
      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

  } // namespace Fad

} // namespace Sacado


#endif // SACADO_FAD_OPS_HPP
