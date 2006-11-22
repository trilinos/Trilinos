// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_EXPRESSION_HPP
#define SACADO_CACHEFAD_EXPRESSION_HPP

#include "Sacado_Traits.hpp"

namespace Sacado {

  namespace CacheFad {

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

      //! Constructor with given expression \c expr
      explicit Expr(const ExprT& expr) : expr_(expr) {}

      //! Return size of derivative array of expression
      int size() const {return expr_.size();}
      
      //! Return if expression has fast access
      bool hasFastAccess() const { return expr_.hasFastAccess();}

      //! Return value of expression
      value_type val() const { return expr_.val();}

      //! Return derivative component \c i of expression
      value_type dx(int i) const { return expr_.dx(i);}

      //! Rturn derivative component \c i of expression
      value_type fastAccessDx(int i) const { return expr_.fastAccessDx(i);}
      
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

      //! Constructor
      ConstExpr(const ConstT& constant) : constant_(constant) {}

      //! Return size of the derivative array of the operation
      int size() const { return 0; }
      
      //! Return if operation has fast access
      bool hasFastAccess() const { return 1; }

      //! Return value of operation
      value_type val() const { return constant_; }

      //! Return derivative component \c i of operation
      value_type dx(int i) const { return value_type(0); }
      
      //! Return derivative component \c i of operation
      value_type fastAccessDx(int i) const { return value_type(0); }

    protected:
      
      //! The constant
      const ConstT& constant_;

    }; // class ConstExpr

    //! Unary expression template
    /*!
     * This template class represents a unary operation of the form
     * op(a) where a is the argument of type \c ExprT and op is the 
     * operation represented by type \c Op. The operation is evaluated by the 
     * non-static methods Op::computeValue() and Op::computeDx().
     *
     * It is assumed Op::computeValue() will cache its result for later
     * Op::computeDx() calls.
     */
    template <typename ExprT, template<typename> class Op> 
    class UnaryExpr {

    public:

      //! Typename of argument value
      typedef typename ExprT::value_type value_type;

      //! Constructor
      UnaryExpr(const ExprT& expr) : expr_(expr), op_(expr) {}

      //! Return size of the derivative array of the operation
      int size() const { return expr_.size(); }
      
      //! Return if operation has fast access
      bool hasFastAccess() const { return expr_.hasFastAccess(); }

      //! Return value of operation
      value_type val() const { return op_.computeValue(expr_); }

      //! Return derivative component \c i of operation
      value_type dx(int i) const { return op_.computeDx(i,expr_); }
      
      //! Return derivative component \c i of operation
      value_type fastAccessDx(int i) const { 
	return op_.computeFastAccessDx(i,expr_); 
      }

    protected:
      
      //! Left argument
      const ExprT& expr_;

      //! Operator
      Op<ExprT> op_;

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
     * represented by type \c Op. The operation is evaluated by the non-static 
     * methods Op::computeValue() and Op::computeDx().
     *
     * It is assumed Op::computeValue() will cache its result for later
     * Op::computeDx() calls.
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

      //! Constructor
      BinaryExpr(const ExprT1& expr1, const ExprT2& expr2) : 
	expr1_(expr1), expr2_(expr2), op_(expr1,expr2) {}

      //! Return size of the derivative array of the operation
      int size() const {
	int sz1 = expr1_.size(), sz2 = expr2_.size(); 
	return sz1 > sz2 ? sz1 : sz2;
      }
      
      //! Return if operation has fast access
      bool hasFastAccess() const { 
	return expr1_.hasFastAccess() && expr2_.hasFastAccess();}

      //! Return value of operation
      value_type val() const { 
	return op_.computeValue(expr1_,expr2_); }

      //! Return derivative component \c i of operation
      value_type dx(int i) const { 
	return op_.computeDx(i,expr1_,expr2_); }
      
      //! Return derivative component \c i of operation
      value_type fastAccessDx(int i) const { 
	return op_.computeFastAccessDx(i,expr1_,expr2_); 
      }

    protected:
      
      //! Left argument
      typename ExprConstRef<ExprT1>::type expr1_;

      //! Right argument
      typename ExprConstRef<ExprT2>::type expr2_;

      //! Operator
      Op<ExprT1,ExprT2> op_;

    }; // class BinaryExpr

  } // namespace CacheFad

} // namespace Sacado

#endif // SACADO_CACHEFAD_EXPRESSION_HPP
