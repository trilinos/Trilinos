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

#ifndef SACADO_PCE_HERMITE_HPP
#define SACADO_PCE_HERMITE_HPP

#include "Sacado_Handle.hpp"
#include "Sacado_PCE_HermiteBasis.hpp"
#include "Sacado_PCE_HermiteEBasis.hpp"
#include "Sacado_PCE_UnitHermiteBasis.hpp"
#include "Sacado_PCE_Workspace.hpp"

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
    class Hermite {
    public:

      //! Turn Hermite into a meta-function class usable with mpl::apply
      template <typename U> 
      struct apply {
	typedef Hermite<U> type;
      };

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      Hermite();

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      Hermite(const T& x);

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      Hermite(unsigned int d, const T & x);

      //! Constructor with degree d
      /*!
       * Initializes all components to zero
       */
      Hermite(unsigned int d);

      //! Copy constructor
      Hermite(const Hermite& x);

      //! Destructor
      ~Hermite();

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

      //! Initialize workspace
      /*! 
       * Intializes static workspace data.
       */
      static void initWorkspace(unsigned int d);

      //! Write coefficients in standard basis
      StandardPoly<T> toStandardBasis() const;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      Hermite<T>& operator=(const T& val);

      //! Assignment with Hermite right-hand-side
      Hermite<T>& operator=(const Hermite<T>& x);

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
      unsigned int degree() const { return th->deg_;}

      //! Returns true if polynomial has degree >= d
      bool hasFastAccess(unsigned int d) const { return th->deg_>=d;}

      //! Returns Hermite coefficient array
      const T* coeff() const { return th->coeff_;}

      //! Returns Hermite coefficient array
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
      Hermite<T> operator + () const;

      //! Unary-minus operator
      Hermite<T> operator - () const;

      //! Addition-assignment operator with constant right-hand-side
      Hermite<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      Hermite<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      Hermite<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      Hermite<T>& operator /= (const T& x);

      //! Addition-assignment operator with Hermite right-hand-side
      Hermite<T>& operator += (const Hermite<T>& x);

      //! Subtraction-assignment operator with Hermite right-hand-side
      Hermite<T>& operator -= (const Hermite<T>& x);
  
      //! Multiplication-assignment operator with Hermite right-hand-side
      Hermite<T>& operator *= (const Hermite<T>& x);

      //! Division-assignment operator with Hermite right-hand-side
      Hermite<T>& operator /= (const Hermite<T>& x);

      //@}

    protected:

      //! Return length of array
      unsigned int length() const { return th->len_; }

      //! Resize coefficient array to new size
      void resizeCoeffs(unsigned int len);

    protected:

      struct HermiteData {

	//! Hermite polynomial coefficients
	T* coeff_;

	//! Degree of polynomial
	unsigned int deg_;

	//! Length of allocated polynomial array
	unsigned int len_;

	//! Default constructor
	HermiteData();

	//! Constructor with supplied value \c x
	HermiteData(const T& x);

	//! Constructor with degree d and value \c x
	HermiteData(unsigned int d, const T & x);

	//! Constructor with degree d
	HermiteData(unsigned int d);

	//! Constructor with degree d and length l
	HermiteData(unsigned int d, unsigned int l);

	//! Copy constructor
	HermiteData(const HermiteData& x);

	//! Destructor
	~HermiteData();

	//! Assignment operator
	HermiteData& operator=(const HermiteData& x);

      };

      Sacado::Handle<HermiteData> th;

    public:

      //! Workspace type
      typedef Workspace< HermiteEBasis<T> > ws_type;

      //! Workspace
      static ws_type workspace;

    }; // class Hermite

    // Operations
    template <typename T> Hermite<T> operator+(const Hermite<T>& a,
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator+(const T& a, 
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator+(const Hermite<T>& a, 
					      const T& b);
    template <typename T> Hermite<T> operator-(const Hermite<T>& a, 
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator-(const T& a, 
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator-(const Hermite<T>& a, 
					      const T& b);
    template <typename T> Hermite<T> operator*(const Hermite<T>& a, 
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator*(const T& a, 
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator*(const Hermite<T>& a, 
					      const T& b);
    template <typename T> Hermite<T> operator/(const Hermite<T>& a, 
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator/(const T& a, 
					      const Hermite<T>& b);
    template <typename T> Hermite<T> operator/(const Hermite<T>& a, 
					      const T& b);
    template <typename T> Hermite<T> exp(const Hermite<T>& a);
    template <typename T> Hermite<T> log(const Hermite<T>& a);
    template <typename T> Hermite<T> log10(const Hermite<T>& a);
    template <typename T> Hermite<T> sqrt(const Hermite<T>& a);
    template <typename T> Hermite<T> pow(const Hermite<T>& a, 
					const Hermite<T>& b);
    template <typename T> Hermite<T> pow(const T& a,
					const Hermite<T>& b);
    template <typename T> Hermite<T> pow(const Hermite<T>& a,
					const T& b);
    template <typename T> void sincos(const Hermite<T>& a,
				      Hermite<T>& s, Hermite<T>& c);
    template <typename T> Hermite<T> cos(const Hermite<T>& a);
    template <typename T> Hermite<T> sin(const Hermite<T>& a);
    template <typename T> Hermite<T> tan(const Hermite<T>& a);
    template <typename T> void sinhcosh(const Hermite<T>& a,
					Hermite<T>& s, Hermite<T>& c);
    template <typename T> Hermite<T> cosh(const Hermite<T>& a);
    template <typename T> Hermite<T> sinh(const Hermite<T>& a);
    template <typename T> Hermite<T> tanh(const Hermite<T>& a);
    template <typename T, typename OpT> Hermite<T> quad(const OpT& quad_func,
							const Hermite<T>& a,
							const Hermite<T>& b);
    template <typename T> Hermite<T> acos(const Hermite<T>& a);
    template <typename T> Hermite<T> asin(const Hermite<T>& a);
    template <typename T> Hermite<T> atan(const Hermite<T>& a);
    template <typename T> Hermite<T> atan2(const Hermite<T>& a,
					  const Hermite<T>& b);
    template <typename T> Hermite<T> atan2(const T& a,
					  const Hermite<T>& b);
    template <typename T> Hermite<T> atan2(const Hermite<T>& a,
					  const T& b);
    template <typename T> Hermite<T> acosh(const Hermite<T>& a);
    template <typename T> Hermite<T> asinh(const Hermite<T>& a);
    template <typename T> Hermite<T> atanh(const Hermite<T>& a);
    template <typename T> Hermite<T> abs(const Hermite<T>& a);
    template <typename T> Hermite<T> fabs(const Hermite<T>& a);
//     template <typename T> Hermite<T> deriv(const Hermite<T>& a);
    template <typename T> Hermite<T> max(const Hermite<T>& a,
					const Hermite<T>& b);
    template <typename T> Hermite<T> max(const T& a,
					const Hermite<T>& b);
    template <typename T> Hermite<T> max(const Hermite<T>& a,
					const T& b);
    template <typename T> Hermite<T> min(const Hermite<T>& a,
					const Hermite<T>& b);
    template <typename T> Hermite<T> min(const T& a,
					const Hermite<T>& b);
    template <typename T> Hermite<T> min(const Hermite<T>& a,
					const T& b);
    template <typename T> bool operator==(const Hermite<T>& a,
					  const Hermite<T>& b);
    template <typename T> bool operator==(const T& a,
					  const Hermite<T>& b);
    template <typename T> bool operator==(const Hermite<T>& a,
					  const T& b);
    template <typename T> bool operator!=(const Hermite<T>& a,
					  const Hermite<T>& b);
    template <typename T> bool operator!=(const T& a,
					  const Hermite<T>& b);
    template <typename T> bool operator!=(const Hermite<T>& a,
					  const T& b);
    template <typename T> bool operator<=(const Hermite<T>& a,
					  const Hermite<T>& b);
    template <typename T> bool operator<=(const T& a,
					  const Hermite<T>& b);
    template <typename T> bool operator<=(const Hermite<T>& a,
					  const T& b);
    template <typename T> bool operator>=(const Hermite<T>& a,
					  const Hermite<T>& b);
    template <typename T> bool operator>=(const T& a,
					  const Hermite<T>& b);
    template <typename T> bool operator>=(const Hermite<T>& a,
					  const T& b);
    template <typename T> bool operator<(const Hermite<T>& a,
					 const Hermite<T>& b);
    template <typename T> bool operator<(const T& a,
					 const Hermite<T>& b);
    template <typename T> bool operator<(const Hermite<T>& a,
					 const T& b);
    template <typename T> bool operator>(const Hermite<T>& a,
					 const Hermite<T>& b);
    template <typename T> bool operator>(const T& a,
					 const Hermite<T>& b);
    template <typename T> bool operator>(const Hermite<T>& a,
					 const T& b);
    template <typename T> std::ostream& operator << (std::ostream& os,
						     const Hermite<T>& a);

  } // namespace PCE

} // namespace Sacado

#include "Sacado_PCE_HermiteTraits.hpp"
#include "Sacado_PCE_HermiteImp.hpp"

#endif // SACADO_PCE_HERMITE_HPP
