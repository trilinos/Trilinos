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

#ifndef SACADO_ELRFAD_EXPRESSION_HPP
#define SACADO_ELRFAD_EXPRESSION_HPP

#include "Sacado_Traits.hpp"

namespace Sacado {

  namespace ELRFad {

    //! Base template specification for %ExprPromote
    /*!
     * The %ExprPromote classes provide a mechanism for computing the 
     * promoted expression-type of a binary operation.
     */
    template <typename A, typename B> struct ExprPromote {};

    //! Specialization of %ExprPromote for a single type
    template <typename A> struct ExprPromote<A,A> {
      typedef A type;
    };

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all Fad expression
     * template classes.
     */
    template <typename ExprT> 
    class Expr {

    public:

      //! Typename of values
      typedef typename ExprT::value_type value_type;

      //! Typename of base-expressions
      typedef typename ExprT::base_expr_type base_expr_type;

      //! Number of arguments
      static const int num_args = ExprT::num_args;

      //! Constructor with given expression \c expr
      explicit Expr(const ExprT& expr) : expr_(expr) {}

      //! Return size of derivative array of expression
      int size() const {return expr_.size();}

      //! Return value of expression
      value_type val() const { return expr_.val();}

      //! Return partials w.r.t. arguments
      void computePartials(const value_type& bar, 
			   value_type partials[]) const { 
	expr_.computePartials(bar, partials); }

      //! Return tangent component \c i of arguments
      void getTangents(int i, value_type dots[]) const { 
	expr_.getTangents(i, dots); }

      //! Return tangent component \c i of argument \c Arg
      template <int Arg>
      value_type getTangent(int i) const {
	return expr_.template getTangent<Arg>(i);
      }

      //! Return whether argument is active
      template <int Arg>
      bool isActive() const { return expr_.template isActive<Arg>(); }
      
    protected:

      //! Disallow default constructor
      Expr() {}

      //! Expression
      ExprT expr_;

    }; // class Expr

    //! Constant expression template
    /*!
     * This template class represents a constant expression.
     */
    template <typename ConstT> 
    class ConstExpr {

    public:

      //! Typename of argument values
      typedef ConstT value_type;

      //! Typename of base-expressions
      typedef ConstT base_expr_type;

      //! Number of arguments
      static const int num_args = 0;

      //! Constructor
      ConstExpr(const ConstT& constant) : constant_(constant) {}

      //! Return size of the derivative array of the operation
      int size() const { return 0; }

      //! Return value of operation
      value_type val() const { return constant_; }

      //! Return partials w.r.t. arguments
      void computePartials(const value_type& bar, 
			   value_type partials[]) const {}

      //! Rturn tangent component \c i of arguments
      void getTangents(int i, value_type dots[]) const {}

      //! Return tangent component \c i of argument \c Arg
      template <int Arg>
      value_type getTangent(int i) const { return 0.0; }

      //! Return whether argument is active
      template <int Arg>
      bool isActive() const { return false; }

    protected:
      
      //! The constant
      const ConstT& constant_;

    }; // class ConstExpr

    //! Unary expression template
    /*!
     * This template class represents a unary operation of the form
     * op(a) where a is the argument of type \c ExprT and op is the 
     * operation represented by type \c Op. The operation is evaluated by the 
     * static methods Op::computeValue() and Op::computeDx().
     */
    template <typename ExprT, template<typename> class Op> 
    class UnaryExpr {

    public:

      //! Typename of argument value
      typedef typename ExprT::value_type value_type;

      //! Typename of base-expressions
      typedef typename ExprT::base_expr_type base_expr_type;

      //! Typename of op
      typedef Op<ExprT> OpT;

      //! Number of arguments
      static const int num_args = ExprT::num_args;

      //! Constructor
      UnaryExpr(const ExprT& expr) : expr_(expr) {}

      //! Return size of the derivative array of the operation
      int size() const { return expr_.size(); }

      //! Return value of operation
      value_type val() const { return OpT::computeValue(expr_); }

      //! Return partials w.r.t. arguments
      void computePartials(const value_type& bar, 
			   value_type partials[]) const { 
	expr_.computePartials(OpT::computeAdjoint(bar, expr_), partials); }

      //! Rturn tangent component \c i of arguments
      void getTangents(int i, value_type dots[]) const { 
	expr_.getTangents(i, dots); }

      //! Return tangent component \c i of argument \c Arg
      template <int Arg>
      value_type getTangent(int i) const {
	return expr_.template getTangent<Arg>(i);
      }

      //! Return whether argument is active
      template <int Arg>
      bool isActive() const { return expr_.template isActive<Arg>(); }

    protected:
      
      //! Left argument
      const ExprT& expr_;

    }; // class UnaryExpr

    // The Sun compiler has difficulty with class partial specialization and
    // template templates, which the original BinaryExpr classes below use.
    // However the only significant difference between the specializations
    // of BinaryExpr for constant arguments is removing the reference
    // for the corresponding data member.  The type functions below allow
    // us to do this without specializing BinaryExpr.  We could also 
    // remove the template templates, but that drastically increases
    // compile times due to the increased length of template argument lists.
    template <typename T> struct ExprConstRef {
      typedef const T& type;
    };
    template <typename T> struct ExprConstRef< ConstExpr<T> > {
      typedef const ConstExpr<T> type;
    };

    //! Binary expression template
    /*!
     * This template class represents a binary operation of the form
     * op(a1,a2) where a1 is the left argument of type \c ExprT1, r is 
     * the right argument of type \c ExprT2, and op is the operation 
     * represented by type \c Op. The operation is evaluated by the static 
     * methods Op::computeValue() and Op::computeDx().
     */
    template <typename ExprT1, typename ExprT2, 
	      template<typename,typename> class Op> 
    class BinaryExpr {
		
    public:

      //! Typename of the first argument value
      typedef typename ExprT1::value_type value_type_1;

      //! Typename of the second argument value
      typedef typename ExprT2::value_type value_type_2;

      //! Typename of the expression values
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      //! Typename of base-expression of first argument
      typedef typename ExprT1::base_expr_type base_expr_type_1;

      //! Typename of base-expression of second argument
      typedef typename ExprT2::base_expr_type base_expr_type_2;

      //! Typename of base-expression
      typedef typename ExprPromote<base_expr_type_1,
				   base_expr_type_2>::type base_expr_type;

      //! Number of arguments of left arg
      static const int num_args1 = ExprT1::num_args;

      //! Number of arguments of right arg
      static const int num_args2 = ExprT2::num_args;

      //! Number of arguments
      static const int num_args = num_args1 + num_args2;

      typedef Op<ExprT1,ExprT2> OpT;

      //! Constructor
      BinaryExpr(const ExprT1& expr1, const ExprT2& expr2) : 
	expr1_(expr1), expr2_(expr2) {}

      //! Return size of the derivative array of the operation
      int size() const {
	int sz1 = expr1_.size(), sz2 = expr2_.size(); 
	return sz1 > sz2 ? sz1 : sz2;
      }

      //! Return partials w.r.t. arguments
      void computePartials(const value_type& bar, 
			   value_type partials[]) const { 
	if (num_args1 > 0)
	  expr1_.computePartials(OpT::computeLeftAdjoint(bar, expr1_, expr2_), 
				 partials);
	if (num_args2 > 0)
	  expr2_.computePartials(OpT::computeRightAdjoint(bar, expr1_, expr2_),
				 partials+num_args1);
      }

      //! Rturn tangent component \c i of arguments
      void getTangents(int i, value_type dots[]) const { 
	expr1_.getTangents(i, dots); 
	expr2_.getTangents(i, dots+num_args1); }

      //! Return value of operation
      value_type val() const { 
	return OpT::computeValue(expr1_,expr2_); }

      //! Return tangent component \c i of argument \c Arg
      template <int Arg>
      value_type getTangent(int i) const {
	if (Arg < num_args1)
	  return expr1_.template getTangent<Arg>(i);
	else
	  return expr2_.template getTangent<Arg-num_args1>(i);
      }

      //! Return whether argument is active
      template <int Arg>
      bool isActive() const { 
	if (Arg < num_args1)
	  return expr1_.template isActive<Arg>(); 
	else
	  return expr2_.template isActive<Arg-num_args1>();
      }

    protected:
      
      //! Left argument
      typename ExprConstRef<ExprT1>::type expr1_;

      //! Right argument
      typename ExprConstRef<ExprT2>::type expr2_;

    }; // class BinaryExpr

  } // namespace ELRFad

} // namespace Sacado

#endif // SACADO_ELRFAD_EXPRESSION_HPP
