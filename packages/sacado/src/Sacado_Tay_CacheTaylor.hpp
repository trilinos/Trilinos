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
// @HEADER

#ifndef SACADO_TAY_CACHETAYLOR_HPP
#define SACADO_TAY_CACHETAYLOR_HPP

#include <valarray>

#include "Sacado_Tay_CacheTaylorExpr.hpp"

// forward decalarations
namespace Sacado {
  namespace Tay {
    template <class ExprT> class UnaryPlusOp;
    template <class ExprT> class UnaryMinusOp;
  }
}

namespace Sacado {

  namespace Tay {

    //! Taylor polynomial class using caching expression templates
    /*!
     * This class provides the implementation of the Taylor object required
     * for expression templating.  Class CacheTaylor provides the complete
     * user inteface.
     */
    template <typename T> 
    class CacheTaylorImplementation {

    public:

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      CacheTaylorImplementation() : coeff_(T(0.),1) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      CacheTaylorImplementation(const T& x) : coeff_(x,1) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      CacheTaylorImplementation(unsigned int d, const T & x) : 
	coeff_(T(0.),d+1) {
	coeff_[0] = x;
      }

      //! Copy constructor
      CacheTaylorImplementation(const CacheTaylorImplementation& x) : 
	coeff_(x.coeff_) {}

      //! Destructor
      ~CacheTaylorImplementation() {}

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
	T tmp= i<coeff_.size() ? coeff_[i]:T(0.); return tmp;}

      //! Returns degree \c i term with bounds checking
      T coeff(unsigned int i) { 
	T tmp= i<coeff_.size() ? coeff_[i]:T(0.); return tmp;}
    
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
	coeff_.resize(dnew+1,T(0.));
	coeff_[s] = tmp;
      }

    protected:

      //! Taylor polynomial coefficients
      std::valarray<T> coeff_;

    }; // class CacheTaylorImplementation

    //! CacheTaylor expression template specialization
    /*!
     * This template class represents a simple CacheTaylor expression.
     */
    template <typename T> 
    class Expr< CacheTaylorImplementation<T> > : 
      public CacheTaylorImplementation<T> {

    public:

      //! Default constructor
      Expr() : CacheTaylorImplementation<T>() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      Expr(const T & x) : CacheTaylorImplementation<T>(x) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      Expr(unsigned int d, const T & x) : CacheTaylorImplementation<T>(d,x) {}

      //! Copy constructor
      Expr(const Expr& x) : CacheTaylorImplementation<T>(x) {}

    }; // class Expr< CacheTaylorImplementation<T> >

    //! Forward-mode AD class using dynamic memory allocation
    /*!
     * This class provides the user interface of the Taylor object.  Class
     * CacheTaylorImplementation provides the implementation.
     */
    template <typename T>
    class CacheTaylor : public Expr< CacheTaylorImplementation<T> > {

    public:

      //! Turn CacheTaylor into a meta-function class usable with mpl::apply
      template <typename U> 
      struct apply {
	typedef CacheTaylor<U> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      CacheTaylor() : Expr< CacheTaylorImplementation<T> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      CacheTaylor(const T & x) : Expr< CacheTaylorImplementation<T> >(x) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      CacheTaylor(unsigned int d, const T & x) : 
	Expr< CacheTaylorImplementation<T> >(d,x) {}

      //! Copy constructor
      CacheTaylor(const CacheTaylor& x) : Expr< CacheTaylorImplementation<T> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> CacheTaylor(const Expr<S>& x);

      //@}

      //! Destructor
      ~CacheTaylor() {}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      CacheTaylor<T>& operator=(const T& v);

      //! Assignment with CacheTaylor right-hand-side
      CacheTaylor<T>& operator=(const CacheTaylor<T>& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> CacheTaylor<T>& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      inline Expr< UnaryExpr< CacheTaylor<T>, UnaryPlusOp > >
      operator + () const {
	typedef UnaryExpr< CacheTaylor<T>, UnaryPlusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Unary-minus operator
      inline Expr< UnaryExpr< CacheTaylor<T>, UnaryMinusOp > >
      operator - () const {
	typedef UnaryExpr< CacheTaylor<T>, UnaryMinusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Addition-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator /= (const T& x);

      //! Addition-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator += (const S& x);

      //! Subtraction-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator -= (const S& x);
  
      //! Multiplication-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator *= (const S& x);

      //! Division-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator /= (const S& x);

      //@}

    }; // class CacheTaylor<T>

  } // namespace Tay

} // namespace Sacado

#include "Sacado_Tay_CacheTaylorTraits.hpp"
#include "Sacado_Tay_CacheTaylorImp.hpp"
#include "Sacado_Tay_CacheTaylorOps.hpp"

#endif // SACADO_TAYLOR_CACHETAYLOR_HPP
