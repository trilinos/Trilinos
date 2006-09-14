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

#ifndef SACADO_FAD_SFAD_HPP
#define SACADO_FAD_SFAD_HPP

#include <valarray>

#include "Sacado_ConfigDefs.h"
#include "Sacado_Fad_Expression.hpp"

// forward decalarations
namespace Sacado {
  namespace Fad {
    template <class ExprT> class UnaryPlusOp;
    template <class ExprT> class UnaryMinusOp;
  }
}

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace Fad {

    //! Forward-mode AD class using static memory allocation
    /*!
     * This class provides the implementation of the Fad object required
     * for expression templating.  Class SFad provides the complete
     * user inteface.
     */
    template <typename T, int Num> 
    class SFadImplementation {

    public:

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      SFadImplementation() : val_( T(0)) { 
	memset(dx_,0,Num*sizeof(T)); }

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      SFadImplementation(const T & x) : val_(x) { 
	memset(dx_,0,Num*sizeof(T)); }

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SFadImplementation(const int sz, const T & x) : val_(x) {
	memset(dx_,0,Num*sizeof(T));}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SFadImplementation(const int sz, const int i, const T & x) : val_(x) {
	memset(dx_,0,Num*sizeof(T));
	dx_[i]=1.; 
      }

      //! Copy constructor
      SFadImplementation(const SFadImplementation& x) : val_(x.val_) {
	for (int i=0; i<Num; i++)
	  dx_[i] = x.dx_[i];
      }

      //! Destructor
      ~SFadImplementation() {}

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
      int size() const { return Num;}

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return true;}

      //! Returns derivative array
      const T* dx() const { return &dx_[0]; }

      //! Returns derivative component \c i with bounds checking
      const T dx(int i) const { return dx_[i]; }

      //! Returns derivative component \c i with bounds checking
      T dx(int i) { return dx_[i]; }
    
      //! Returns derivative component \c i without bounds checking
      T& fastAccessDx(int i) { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      T fastAccessDx(int i) const { return dx_[i]; }
    
      //@}

    protected:

      //! Value
      T val_;

      //! Derivative array
      T dx_[Num];

    }; // class SFadImplementation

    //! SFad expression template specialization
    /*!
     * This template class represents a simple SFad expression.
     */
    template <typename T, int Num> 
    class Expr< SFadImplementation<T,Num> > : 
      public SFadImplementation<T,Num> {

    public:

      //! Default constructor
      Expr() : SFadImplementation<T,Num>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Expr(const T & x) : SFadImplementation<T,Num>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const T & x) : SFadImplementation<T,Num>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const T & x) : 
	SFadImplementation<T,Num>(sz,i,x) {}

      //! Copy constructor
      Expr(const Expr& x) : SFadImplementation<T,Num>(x) {}

    }; // class Expr< SFadImplementation<T,Num> >

    //! Forward-mode AD class using static memory allocation
    /*!
     * This class provides the user interface of the Fad object.  Class
     * SFadImplementation provides the implementation.
     *
     * Current this class does no checking to make sure the size argument
     * of the constructor matches the template, does not check whether
     * arguments have the right dimensions, and seems to be pretty inefficient
     * when dealing with constants.
     */
    template <typename T, int Num>
    class SFad : public Expr< SFadImplementation<T,Num> > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      SFad() : Expr< SFadImplementation<T,Num> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      SFad(const T & x) : Expr< SFadImplementation<T,Num> >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SFad(const int sz, const T & x) : 
	Expr< SFadImplementation<T,Num> >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SFad(const int sz, const int i, const T & x) : 
	Expr< SFadImplementation<T,Num> >(sz,i,x) {}

      //! Copy constructor
      SFad(const SFad& x) : Expr< SFadImplementation<T,Num> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> SFad(const Expr<S>& x);

      //! Set %SFad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the 
       * SFadImplementation(const int sz, const int i, const T & x) 
       * constructor.
       */
      void diff(const int ith, const int n);

      //@}

      //! Destructor
      ~SFad() {}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      SFad<T,Num>& operator=(const T& val);

      //! Assignment with SFad right-hand-side
      SFad<T,Num>& operator=(const SFad<T,Num>& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> SFad<T,Num>& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      inline Expr< UnaryExpr< SFad<T,Num>, UnaryPlusOp > >
      operator + () const {
	typedef UnaryExpr< SFad<T,Num>, UnaryPlusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Unary-minus operator
      inline Expr< UnaryExpr< SFad<T,Num>, UnaryMinusOp > >
      operator - () const {
	typedef UnaryExpr< SFad<T,Num>, UnaryMinusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Addition-assignment operator with constant right-hand-side
      SFad<T,Num>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      SFad<T,Num>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      SFad<T,Num>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      SFad<T,Num>& operator /= (const T& x);

      //! Addition-assignment operator with Fad right-hand-side
      template <typename S> SFad<T,Num>& operator += (const S& x);

      //! Subtraction-assignment operator with Fad right-hand-side
      template <typename S> SFad<T,Num>& operator -= (const S& x);
  
      //! Multiplication-assignment operator with Fad right-hand-side
      template <typename S> SFad<T,Num>& operator *= (const S& x);

      //! Division-assignment operator with Fad right-hand-side
      template <typename S> SFad<T,Num>& operator /= (const S& x);

      //@}

    }; // class SFad<T,Num>

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_SFadTraits.hpp"
#include "Sacado_Fad_SFadImp.hpp"
#include "Sacado_Fad_Ops.hpp"

#endif // SACADO_FAD_SFAD_HPP
