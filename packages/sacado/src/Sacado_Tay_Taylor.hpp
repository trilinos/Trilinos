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

#ifndef SACADO_TAY_TAYLOR_HPP
#define SACADO_TAY_TAYLOR_HPP

#include "Sacado_Handle.hpp"
#include <cmath>
#include <algorithm>	// for std::min and std::max
#include <ostream>	// for std::ostream

namespace Sacado {

  //! Namespace for Taylor polynomial AD classes
  namespace Tay {

    //! Taylor polynomial class
    /*!
     * Uses a handle and a "copy-on-write" strategy for efficient copying, but
     * no expression templating.
     */
    template <typename T> 
    class Taylor {
    public:

      //! Turn Taylor into a meta-function class usable with mpl::apply
      template <typename U> 
      struct apply {
	typedef Taylor<U> type;
      };

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      Taylor();

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      Taylor(const T& x);

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      Taylor(unsigned int d, const T & x);

      //! Constructor with degree d
      /*!
       * Initializes all components to zero
       */
      Taylor(unsigned int d);

      //! Copy constructor
      Taylor(const Taylor& x);

      //! Destructor
      ~Taylor();

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
       * %Taylor coefficients is not shared among any other %Taylor polynomial
       * objects.  If the handle is not shared it does nothing, so there
       * is no cost in calling this method in this case.  If the handle is 
       * shared and this method is not called, any changes to the coefficients
       * by coeff() or fastAccessCoeff() may change other polynomial objects.
       */
      void copyForWrite() { th.makeOwnCopy(); }

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      Taylor<T>& operator=(const T& val);

      //! Assignment with Taylor right-hand-side
      Taylor<T>& operator=(const Taylor<T>& x);

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
       * @name Taylor coefficient accessor methods
       */
      //@{

      //! Returns degree of polynomial
      unsigned int degree() const { return th->deg_;}

      //! Returns true if polynomial has degree >= d
      bool hasFastAccess(unsigned int d) const { return th->deg_>=d;}

      //! Returns Taylor coefficient array
      const T* coeff() const { return th->coeff_;}

      //! Returns Taylor coefficient array
      T* coeff() { return th->coeff_;}

      //! Returns degree \c i term with bounds checking
      T coeff(unsigned int i) const { 
	T tmp= i<=th->deg_ ? th->coeff_[i]:T(0.); return tmp;}
    
      //! Returns degree \c i term without bounds checking
      T& fastAccessCoeff(unsigned int i) { return th->coeff_[i];}

      //! Returns degree \c i term without bounds checking
      T fastAccessCoeff(unsigned int i) const { return th->coeff_[i];}
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      Taylor<T> operator + () const;

      //! Unary-minus operator
      Taylor<T> operator - () const;

      //! Addition-assignment operator with constant right-hand-side
      Taylor<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      Taylor<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      Taylor<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      Taylor<T>& operator /= (const T& x);

      //! Addition-assignment operator with Taylor right-hand-side
      Taylor<T>& operator += (const Taylor<T>& x);

      //! Subtraction-assignment operator with Taylor right-hand-side
      Taylor<T>& operator -= (const Taylor<T>& x);
  
      //! Multiplication-assignment operator with Taylor right-hand-side
      Taylor<T>& operator *= (const Taylor<T>& x);

      //! Division-assignment operator with Taylor right-hand-side
      Taylor<T>& operator /= (const Taylor<T>& x);

      //@}

    protected:

      //! Return length of array
      unsigned int length() const { return th->len_; }

      //! Resize coefficient array to new size
      void resizeCoeffs(unsigned int len);

    protected:

      struct TaylorData {

	//! Taylor polynomial coefficients
	T* coeff_;

	//! Degree of polynomial
	unsigned int deg_;

	//! Length of allocated polynomial array
	unsigned int len_;

	//! Default constructor
	TaylorData();

	//! Constructor with supplied value \c x
	TaylorData(const T& x);

	//! Constructor with degree d and value \c x
	TaylorData(unsigned int d, const T & x);

	//! Constructor with degree d
	TaylorData(unsigned int d);

	//! Constructor with degree d and length l
	TaylorData(unsigned int d, unsigned int l);

	//! Copy constructor
	TaylorData(const TaylorData& x);

	//! Destructor
	~TaylorData();

	//! Assignment operator
	TaylorData& operator=(const TaylorData& x);

      };

      Sacado::Handle<TaylorData> th;

    }; // class Taylor

    // Operations
    template <typename T> Taylor<T> operator+(const Taylor<T>& a,
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator+(const T& a, 
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator+(const Taylor<T>& a, 
					      const T& b);
    template <typename T> Taylor<T> operator-(const Taylor<T>& a, 
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator-(const T& a, 
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator-(const Taylor<T>& a, 
					      const T& b);
    template <typename T> Taylor<T> operator*(const Taylor<T>& a, 
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator*(const T& a, 
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator*(const Taylor<T>& a, 
					      const T& b);
    template <typename T> Taylor<T> operator/(const Taylor<T>& a, 
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator/(const T& a, 
					      const Taylor<T>& b);
    template <typename T> Taylor<T> operator/(const Taylor<T>& a, 
					      const T& b);
    template <typename T> Taylor<T> exp(const Taylor<T>& a);
    template <typename T> Taylor<T> log(const Taylor<T>& a);
    template <typename T> Taylor<T> log10(const Taylor<T>& a);
    template <typename T> Taylor<T> sqrt(const Taylor<T>& a);
    template <typename T> Taylor<T> pow(const Taylor<T>& a, 
					const Taylor<T>& b);
    template <typename T> Taylor<T> pow(const T& a,
					const Taylor<T>& b);
    template <typename T> Taylor<T> pow(const Taylor<T>& a,
					const T& b);
    template <typename T> void sincos(const Taylor<T>& a,
				      Taylor<T>& s, Taylor<T>& c);
    template <typename T> Taylor<T> cos(const Taylor<T>& a);
    template <typename T> Taylor<T> sin(const Taylor<T>& a);
    template <typename T> Taylor<T> tan(const Taylor<T>& a);
    template <typename T> void sinhcosh(const Taylor<T>& a,
					Taylor<T>& s, Taylor<T>& c);
    template <typename T> Taylor<T> cosh(const Taylor<T>& a);
    template <typename T> Taylor<T> sinh(const Taylor<T>& a);
    template <typename T> Taylor<T> tanh(const Taylor<T>& a);
    template <typename T> Taylor<T> quad(const T& c0,
					 const Taylor<T>& a,
					 const Taylor<T>& b);
    template <typename T> Taylor<T> acos(const Taylor<T>& a);
    template <typename T> Taylor<T> asin(const Taylor<T>& a);
    template <typename T> Taylor<T> atan(const Taylor<T>& a);
    template <typename T> Taylor<T> atan2(const Taylor<T>& a,
					  const Taylor<T>& b);
    template <typename T> Taylor<T> atan2(const T& a,
					  const Taylor<T>& b);
    template <typename T> Taylor<T> atan2(const Taylor<T>& a,
					  const T& b);
    template <typename T> Taylor<T> acosh(const Taylor<T>& a);
    template <typename T> Taylor<T> asinh(const Taylor<T>& a);
    template <typename T> Taylor<T> atanh(const Taylor<T>& a);
    template <typename T> Taylor<T> abs(const Taylor<T>& a);
    template <typename T> Taylor<T> fabs(const Taylor<T>& a);
    template <typename T> Taylor<T> max(const Taylor<T>& a,
					const Taylor<T>& b);
    template <typename T> Taylor<T> max(const T& a,
					const Taylor<T>& b);
    template <typename T> Taylor<T> max(const Taylor<T>& a,
					const T& b);
    template <typename T> Taylor<T> min(const Taylor<T>& a,
					const Taylor<T>& b);
    template <typename T> Taylor<T> min(const T& a,
					const Taylor<T>& b);
    template <typename T> Taylor<T> min(const Taylor<T>& a,
					const T& b);
    template <typename T> bool operator==(const Taylor<T>& a,
					  const Taylor<T>& b);
    template <typename T> bool operator==(const T& a,
					  const Taylor<T>& b);
    template <typename T> bool operator==(const Taylor<T>& a,
					  const T& b);
    template <typename T> bool operator!=(const Taylor<T>& a,
					  const Taylor<T>& b);
    template <typename T> bool operator!=(const T& a,
					  const Taylor<T>& b);
    template <typename T> bool operator!=(const Taylor<T>& a,
					  const T& b);
    template <typename T> bool operator<=(const Taylor<T>& a,
					  const Taylor<T>& b);
    template <typename T> bool operator<=(const T& a,
					  const Taylor<T>& b);
    template <typename T> bool operator<=(const Taylor<T>& a,
					  const T& b);
    template <typename T> bool operator>=(const Taylor<T>& a,
					  const Taylor<T>& b);
    template <typename T> bool operator>=(const T& a,
					  const Taylor<T>& b);
    template <typename T> bool operator>=(const Taylor<T>& a,
					  const T& b);
    template <typename T> bool operator<(const Taylor<T>& a,
					 const Taylor<T>& b);
    template <typename T> bool operator<(const T& a,
					 const Taylor<T>& b);
    template <typename T> bool operator<(const Taylor<T>& a,
					 const T& b);
    template <typename T> bool operator>(const Taylor<T>& a,
					 const Taylor<T>& b);
    template <typename T> bool operator>(const T& a,
					 const Taylor<T>& b);
    template <typename T> bool operator>(const Taylor<T>& a,
					 const T& b);
    template <typename T> std::ostream& operator << (std::ostream& os,
						     const Taylor<T>& a);

  } // namespace Tay

} // namespace Sacado

#include "Sacado_Tay_TaylorTraits.hpp"
#include "Sacado_Tay_TaylorImp.hpp"

#endif // SACADO_TAY_TAYLOR_HPP
