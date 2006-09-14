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

#ifndef SACADO_CACHEFAD_DFAD_HPP
#define SACADO_CACHEFAD_DFAD_HPP

#include <valarray>

#include "Sacado_ConfigDefs.h"
#include "Sacado_CacheFad_Expression.hpp"

// forward decalarations
namespace Sacado {
  namespace CacheFad {
    template <class ExprT> class UnaryPlusOp;
    template <class ExprT> class UnaryMinusOp;
  }
}

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace CacheFad {

    //! Forward-mode AD class using dynamic memory allocation
    /*!
     * This class provides the implementation of the Fad object required
     * for expression templating.  Class DFad provides the complete
     * user inteface.
     */
    template <typename T> 
    class DFadImplementation {

    public:

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      DFadImplementation() : val_( T(0)), dx_() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DFadImplementation(const T & x) : val_(x), dx_() {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DFadImplementation(const int sz, const T & x) : val_(x), dx_(T(0),sz) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DFadImplementation(const int sz, const int i, const T & x) : 
	val_(x), dx_(T(0),sz) { 
	dx_[i]=1.; 
      }

      //! Copy constructor
      DFadImplementation(const DFadImplementation& x) : 
	val_(x.val_), dx_(x.dx_) {}

      //! Destructor
      ~DFadImplementation() {}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const T& val() const { return val_;}

      //! Returns value
      T& val() { return val_;}

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      //! Returns number of derivative components
      int size() const { return dx_.size();}

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return dx_.size()!=0;}

      //! Returns derivative array
      const std::valarray<T>& dx() const { return dx_;}

      //! Returns derivative component \c i with bounds checking
      const T dx(int i) const { T tmp= dx_.size()? dx_[i]:T(0); return tmp;}

      //! Returns derivative component \c i with bounds checking
      T dx(int i) { T tmp= dx_.size()? dx_[i]:T(0); return tmp;}
    
      //! Returns derivative component \c i without bounds checking
      T& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      T fastAccessDx(int i) const { return dx_[i];}
    
      //@}

    protected:

      //! Value
      T val_;

      //! Derivative array
      std::valarray<T> dx_;

    }; // class DFadImplementation

    //! DFad expression template specialization
    /*!
     * This template class represents a simple DFad expression.
     */
    template <typename T> 
    class Expr< DFadImplementation<T> > : 
      public DFadImplementation<T> {

    public:

      //! Default constructor
      Expr() : DFadImplementation<T>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Expr(const T & x) : DFadImplementation<T>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const T & x) : DFadImplementation<T>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const T & x) : 
	DFadImplementation<T>(sz,i,x) {}

      //! Copy constructor
      Expr(const Expr& x) : DFadImplementation<T>(x) {}

    }; // class Expr< DFadImplementation<T> >

    //! Forward-mode AD class using dynamic memory allocation
    /*!
     * This class provides the user interface of the Fad object.  Class
     * DFadImplementation provides the implementation.
     */
    template <typename T>
    class DFad : public Expr< DFadImplementation<T> > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DFad() : Expr< DFadImplementation<T> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DFad(const T & x) : Expr< DFadImplementation<T> >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DFad(const int sz, const T & x) : 
	Expr< DFadImplementation<T> >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DFad(const int sz, const int i, const T & x) : 
	Expr< DFadImplementation<T> >(sz,i,x) {}

      //! Copy constructor
      DFad(const DFad& x) : Expr< DFadImplementation<T> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DFad(const Expr<S>& x);

      //! Set %DFad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the 
       * DFadImplementation(const int sz, const int i, const T & x) 
       * constructor.
       */
      void diff(const int ith, const int n);

      //@}

      //! Destructor
      ~DFad() {}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      DFad<T>& operator=(const T& val);

      //! Assignment with DFad right-hand-side
      DFad<T>& operator=(const DFad<T>& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> DFad<T>& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      inline Expr< UnaryExpr< DFad<T>, UnaryPlusOp > >
      operator + () const {
	typedef UnaryExpr< DFad<T>, UnaryPlusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Unary-minus operator
      inline Expr< UnaryExpr< DFad<T>, UnaryMinusOp > >
      operator - () const {
	typedef UnaryExpr< DFad<T>, UnaryMinusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Addition-assignment operator with constant right-hand-side
      DFad<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      DFad<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      DFad<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      DFad<T>& operator /= (const T& x);

      //! Addition-assignment operator with Fad right-hand-side
      template <typename S> DFad<T>& operator += (const S& x);

      //! Subtraction-assignment operator with Fad right-hand-side
      template <typename S> DFad<T>& operator -= (const S& x);
  
      //! Multiplication-assignment operator with Fad right-hand-side
      template <typename S> DFad<T>& operator *= (const S& x);

      //! Division-assignment operator with Fad right-hand-side
      template <typename S> DFad<T>& operator /= (const S& x);

      //@}

    }; // class DFad<T>

  } // namespace CacheFad

} // namespace Sacado

#include "Sacado_CacheFad_DFadTraits.hpp"
#include "Sacado_CacheFad_DFadImp.hpp"
#include "Sacado_CacheFad_Ops.hpp"

#endif // SACADO_CACHEFAD_DFAD_HPP
