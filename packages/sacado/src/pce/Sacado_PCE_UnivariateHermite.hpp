// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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

#ifndef SACADO_PCE_UNIVARIATEHERMITE_HPP
#define SACADO_PCE_UNIVARIATEHERMITE_HPP

#ifdef HAVE_SACADO_STOKHOS

#include "Sacado_Handle.hpp"
#include "Stokhos_HermiteEBasis.hpp"
#include "Stokhos_HermiteExpansion.hpp"

#include <cmath>
#include <algorithm>	// for std::min and std::max
#include <ostream>	// for std::ostream

namespace Sacado {

  //! Namespace for polynomial chaos expansion classes
  namespace PCE {

    //! Hermite polynomial chaos expansion class
    /*!
     * Uses a handle and a "copy-on-write" strategy for efficient copying, but
     * no expression templating.
     */
    template <typename T> 
    class UnivariateHermite {
    public:

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      UnivariateHermite();

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      UnivariateHermite(const T& x);

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      UnivariateHermite(unsigned int d, const T & x);

      //! Constructor with degree d
      /*!
       * Initializes all components to zero
       */
      UnivariateHermite(unsigned int d);

      //! Copy constructor
      UnivariateHermite(const UnivariateHermite& x);

      //! Destructor
      ~UnivariateHermite();

      //! Resize polynomial to degree d
      /*!
       * Coefficients are preserved if \c keep_coeffs is \c true, otherwise 
       * all coefficients are reset to zero.
       */
      void resize(unsigned int d, bool keep_coeffs);

      //! Reserve space for a degree d polynomial
      /*!
       * Coefficients are preserved.
       */
      void reserve(unsigned int d);

      //! Prepare polynomial for writing 
      /*!
       * This method prepares the polynomial for writing through coeff() and 
       * fastAccessCoeff() member functions.  It ensures the handle for the
       * %Hermite coefficients is not shared among any other %Hermite polynomial
       * objects.  If the handle is not shared it does nothing, so there
       * is no cost in calling this method in this case.  If the handle is 
       * shared and this method is not called, any changes to the coefficients
       * by coeff() or fastAccessCoeff() may change other polynomial objects.
       */
      void copyForWrite() { th.makeOwnCopy(); }

      //! Initialize expansion
      /*! 
       * Intializes static expansion data.
       */
      static void initExpansion(unsigned int d);

      //! Write coefficients in standard basis
      Stokhos::StandardPoly<T> toStandardBasis() const;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      UnivariateHermite<T>& operator=(const T& val);

      //! Assignment with UnivariateHermite right-hand-side
      UnivariateHermite<T>& operator=(const UnivariateHermite<T>& x);

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const T& val() const { return th->coeff_[0];}

      //! Returns value
      T& val() { return th->coeff_[0];}

      //@}

      /*!
       * @name Hermite coefficient accessor methods
       */
      //@{

      //! Returns degree of polynomial
      unsigned int degree() const { return th->degree();}

      //! Returns true if polynomial has degree >= d
      bool hasFastAccess(unsigned int d) const { return th->degree()>=d;}

      //! Returns Hermite coefficient array
      const T* coeff() const { return th->coeff();}

      //! Returns Hermite coefficient array
      T* coeff() { return th->coeff();}

      //! Returns degree \c i term with bounds checking
      T coeff(unsigned int i) const { 
	T tmp= i<=th->degree() ? (*th)[i]:T(0.); return tmp;}
    
      //! Returns degree \c i term without bounds checking
      T& fastAccessCoeff(unsigned int i) { return (*th)[i];}

      //! Returns degree \c i term without bounds checking
      T fastAccessCoeff(unsigned int i) const { return (*th)[i];}
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      UnivariateHermite<T> operator + () const;

      //! Unary-minus operator
      UnivariateHermite<T> operator - () const;

      //! Addition-assignment operator with constant right-hand-side
      UnivariateHermite<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      UnivariateHermite<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      UnivariateHermite<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      UnivariateHermite<T>& operator /= (const T& x);

      //! Addition-assignment operator with Hermite right-hand-side
      UnivariateHermite<T>& operator += (const UnivariateHermite<T>& x);

      //! Subtraction-assignment operator with Hermite right-hand-side
      UnivariateHermite<T>& operator -= (const UnivariateHermite<T>& x);
  
      //! Multiplication-assignment operator with Hermite right-hand-side
      UnivariateHermite<T>& operator *= (const UnivariateHermite<T>& x);

      //! Division-assignment operator with Hermite right-hand-side
      UnivariateHermite<T>& operator /= (const UnivariateHermite<T>& x);

      //@}

      //! Get underlying Hermite polynomial
      const Stokhos::HermitePoly<T>& getHermitePoly() const { return *th; }

      //! Get underlying Hermite polynomial
      Stokhos::HermitePoly<T>& getHermitePoly() { return *th; }

    public:

      //! Expansion type
      typedef Stokhos::HermiteExpansion<T, Stokhos::HermiteEBasis<T> > he_type;

      static he_type expansion;

    protected:

      Sacado::Handle< Stokhos::HermitePoly<T> > th;

    }; // class Hermite

    // Operations
    template <typename T> UnivariateHermite<T> 
    operator+(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator+(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator+(const UnivariateHermite<T>& a, const T& b);

    template <typename T> UnivariateHermite<T> 
    operator-(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator-(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator-(const UnivariateHermite<T>& a, const T& b);

    template <typename T> UnivariateHermite<T> 
    operator*(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator*(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator*(const UnivariateHermite<T>& a, const T& b);

    template <typename T> UnivariateHermite<T> 
    operator/(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator/(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    operator/(const UnivariateHermite<T>& a, const T& b);

    template <typename T> UnivariateHermite<T> 
    exp(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    log(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    log10(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    sqrt(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    pow(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    pow(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    pow(const UnivariateHermite<T>& a, const T& b);

    template <typename T> UnivariateHermite<T> 
    cos(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    sin(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    tan(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    cosh(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T>
    sinh(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    tanh(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    acos(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    asin(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    atan(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    atan2(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    atan2(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    atan2(const UnivariateHermite<T>& a, const T& b);

    template <typename T> UnivariateHermite<T> 
    acosh(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    asinh(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    atanh(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    abs(const UnivariateHermite<T>& a);
    
    template <typename T> UnivariateHermite<T> 
    fabs(const UnivariateHermite<T>& a);

    template <typename T> UnivariateHermite<T> 
    max(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    max(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    max(const UnivariateHermite<T>& a, const T& b);

    template <typename T> UnivariateHermite<T> 
    min(const UnivariateHermite<T>& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    min(const T& a, const UnivariateHermite<T>& b);

    template <typename T> UnivariateHermite<T> 
    min(const UnivariateHermite<T>& a, const T& b);

    template <typename T> bool operator==(const UnivariateHermite<T>& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator==(const T& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator==(const UnivariateHermite<T>& a,
					  const T& b);
    template <typename T> bool operator!=(const UnivariateHermite<T>& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator!=(const T& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator!=(const UnivariateHermite<T>& a,
					  const T& b);
    template <typename T> bool operator<=(const UnivariateHermite<T>& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator<=(const T& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator<=(const UnivariateHermite<T>& a,
					  const T& b);
    template <typename T> bool operator>=(const UnivariateHermite<T>& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator>=(const T& a,
					  const UnivariateHermite<T>& b);
    template <typename T> bool operator>=(const UnivariateHermite<T>& a,
					  const T& b);
    template <typename T> bool operator<(const UnivariateHermite<T>& a,
					 const UnivariateHermite<T>& b);
    template <typename T> bool operator<(const T& a,
					 const UnivariateHermite<T>& b);
    template <typename T> bool operator<(const UnivariateHermite<T>& a,
					 const T& b);
    template <typename T> bool operator>(const UnivariateHermite<T>& a,
					 const UnivariateHermite<T>& b);
    template <typename T> bool operator>(const T& a,
					 const UnivariateHermite<T>& b);
    template <typename T> bool operator>(const UnivariateHermite<T>& a,
					 const T& b);

    template <typename T> std::ostream& 
    operator << (std::ostream& os, const UnivariateHermite<T>& a);

  } // namespace PCE

} // namespace Sacado

#include "Sacado_PCE_UnivariateHermiteTraits.hpp"
#include "Sacado_PCE_UnivariateHermiteImp.hpp"

#endif // HAVE_SACADO_STOKHOS

#endif // SACADO_PCE_UNIVARIATEHERMITE_HPP
