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

#ifndef SACADO_TAYLOR_DTAYLOR_HPP
#define SACADO_TAYLOR_DTAYLOR_HPP

#include <valarray>

#include "Sacado_ConfigDefs.h"
#include "Sacado_Taylor_Expression.hpp"

// forward decalarations
namespace Sacado {
  namespace Taylor {
    template <class ExprT> class UnaryPlusOp;
    template <class ExprT> class UnaryMinusOp;
  }
}

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace Taylor {

    //! Forward-mode AD class using dynamic memory allocation
    /*!
     * This class provides the implementation of the Taylor object required
     * for expression templating.  Class DTaylor provides the complete
     * user inteface.
     */
    template <typename T> 
    class DTaylorImplementation {

    public:

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      DTaylorImplementation() : coeff_(T(0),1) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      DTaylorImplementation(const T& x) : coeff_(x,1) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      DTaylorImplementation(unsigned int d, const T & x) : coeff_(T(0),d+1) {
	coeff_[0] = x;
      }

      //! Copy constructor
      DTaylorImplementation(const DTaylorImplementation& x) : 
	coeff_(x.coeff_) {}

      //! Destructor
      ~DTaylorImplementation() {}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const T& val() const { return coeff_[0];}

      //! Returns value
      T& val() { return coeff_[0];}

      //@}

      /*!
       * @name Taylor coefficient accessor methods
       */
      //@{

      //! Returns degree of polynomial
      unsigned int degree() const { return coeff_.size()-1;}

      //! Returns true if polynomial has degree >= d
      bool hasFastAccess(unsigned int d) const { return coeff_.size()>=d+1;}

      //! Returns Taylor coefficient array
      const std::valarray<T>& coeff() const { return coeff_;}

      //! Returns degree \c i term with bounds checking
      const T coeff(unsigned int i) const { 
	T tmp= i<coeff_.size() ? coeff_[i]:T(0); return tmp;}

      //! Returns degree \c i term with bounds checking
      T coeff(unsigned int i) { 
	T tmp= i<coeff_.size() ? coeff_[i]:T(0); return tmp;}
    
      //! Returns degree \c i term without bounds checking
      T& fastAccessCoeff(unsigned int i) { return coeff_[i];}

      //! Returns degree \c i term without bounds checking
      T fastAccessCoeff(unsigned int i) const { return coeff_[i];}

      //! Allocate coefficient cache
      void allocateCache(unsigned int d) const {}
    
      //@}

    protected:

      //! Resize coefficient array to new size
      void resizeCoeffs(unsigned int dnew) {
	std::valarray<T> tmp = coeff_;
	std::slice s(0,coeff_.size(),1);
	coeff_.resize(dnew+1,T(0));
	coeff_[s] = tmp;
      }

    protected:

      //! Taylor polynomial coefficients
      std::valarray<T> coeff_;

    }; // class DTaylorImplementation

    //! DTaylor expression template specialization
    /*!
     * This template class represents a simple DTaylor expression.
     */
    template <typename T> 
    class Expr< DTaylorImplementation<T> > : 
      public DTaylorImplementation<T> {

    public:

      //! Default constructor
      Expr() : DTaylorImplementation<T>() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      Expr(const T & x) : DTaylorImplementation<T>(x) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      Expr(unsigned int d, const T & x) : DTaylorImplementation<T>(d,x) {}

      //! Copy constructor
      Expr(const Expr& x) : DTaylorImplementation<T>(x) {}

    }; // class Expr< DTaylorImplementation<T> >

    //! Forward-mode AD class using dynamic memory allocation
    /*!
     * This class provides the user interface of the Taylor object.  Class
     * DTaylorImplementation provides the implementation.
     */
    template <typename T>
    class DTaylor : public Expr< DTaylorImplementation<T> > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      DTaylor() : Expr< DTaylorImplementation<T> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      DTaylor(const T & x) : Expr< DTaylorImplementation<T> >(x) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      DTaylor(unsigned int d, const T & x) : 
	Expr< DTaylorImplementation<T> >(d,x) {}

      //! Copy constructor
      DTaylor(const DTaylor& x) : Expr< DTaylorImplementation<T> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DTaylor(const Expr<S>& x);

      //@}

      //! Destructor
      ~DTaylor() {}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      DTaylor<T>& operator=(const T& val);

      //! Assignment with DTaylor right-hand-side
      DTaylor<T>& operator=(const DTaylor<T>& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> DTaylor<T>& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      inline Expr< UnaryExpr< DTaylor<T>, UnaryPlusOp > >
      operator + () const {
	typedef UnaryExpr< DTaylor<T>, UnaryPlusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Unary-minus operator
      inline Expr< UnaryExpr< DTaylor<T>, UnaryMinusOp > >
      operator - () const {
	typedef UnaryExpr< DTaylor<T>, UnaryMinusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Addition-assignment operator with constant right-hand-side
      DTaylor<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      DTaylor<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      DTaylor<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      DTaylor<T>& operator /= (const T& x);

      //! Addition-assignment operator with Taylor right-hand-side
      template <typename S> DTaylor<T>& operator += (const S& x);

      //! Subtraction-assignment operator with Taylor right-hand-side
      template <typename S> DTaylor<T>& operator -= (const S& x);
  
      //! Multiplication-assignment operator with Taylor right-hand-side
      template <typename S> DTaylor<T>& operator *= (const S& x);

      //! Division-assignment operator with Taylor right-hand-side
      template <typename S> DTaylor<T>& operator /= (const S& x);

      //@}

    }; // class DTaylor<T>

  } // namespace Taylor

} // namespace Sacado

#include "Sacado_Taylor_DTaylorTraits.hpp"
#include "Sacado_Taylor_DTaylorImp.hpp"
#include "Sacado_Taylor_Ops.hpp"

#endif // SACADO_TAYLOR_DTAYLOR_HPP
