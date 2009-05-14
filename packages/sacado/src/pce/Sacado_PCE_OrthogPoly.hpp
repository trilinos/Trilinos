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

#ifndef SACADO_PCE_ORTHOGPOLY_HPP
#define SACADO_PCE_ORTHOGPOLY_HPP

#include "Sacado_ConfigDefs.h"

#ifdef HAVE_SACADO_STOKHOS

#include "Teuchos_RCP.hpp"

#include "Sacado_Handle.hpp"

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"

#include <cmath>
#include <algorithm>	// for std::min and std::max
#include <ostream>	// for std::ostream

namespace Sacado {

  //! Namespace for polynomial chaos expansion classes
  namespace PCE {

    //! Generalized polynomial chaos expansion class
    /*!
     * Uses a handle and a "copy-on-write" strategy for efficient copying, but
     * no expression templating.
     */
    template <typename T> 
    class OrthogPoly {
    public:

      //! Turn OrthogPoly into a meta-function class usable with mpl::apply
      template <typename U> 
      struct apply {
	typedef OrthogPoly<U> type;
      };

      //! Typename of values
      typedef T value_type;

      //! Typename of ordinals
      typedef int ordinal_type;

      //! Basis type
      typedef Stokhos::OrthogPolyBasis<ordinal_type,T> basis_type;

      //! Expansion type
      typedef Stokhos::OrthogPolyExpansion<ordinal_type,T> expansion_type;

      //! Default constructor
      OrthogPoly();

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      OrthogPoly(const value_type& x);

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of size \c sz
       */
      OrthogPoly(ordinal_type sz, const value_type& x);

      //! Constructor with size \c sz
      /*!
       * Initializes all components to zero
       */
      OrthogPoly(ordinal_type sz);

      //! Copy constructor
      OrthogPoly(const OrthogPoly& x);

      //! Destructor
      ~OrthogPoly();

      //! Resize polynomial to size \c sz
      /*!
       * Coefficients are preserved.
       */
      void resize(ordinal_type sz);

      //! Reserve space for polynomial of size \c sz
      /*!
       * Coefficients are preserved.
       */
      void reserve(ordinal_type sz);

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
      static void initExpansion(const Teuchos::RCP<expansion_type>& e);

      //! Write coefficients in standard basis
      Stokhos::Polynomial<value_type> toStandardBasis() const;

      //! Evaluate polynomial approximation at a point
      value_type evaluate(const std::vector<value_type>& point) const;

      //! Evaluate polynomial approximation at a point with given basis values
      value_type evaluate(const std::vector<value_type>& point,
                          const std::vector<value_type>& bvals) const;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      OrthogPoly<T>& operator=(const value_type& val);

      //! Assignment with OrthogPoly right-hand-side
      OrthogPoly<T>& operator=(const OrthogPoly<T>& x);

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const value_type& val() const { return (*th)[0]; }

      //! Returns value
      value_type& val() { return (*th)[0]; }

      //@}

      /*!
       * @name Hermite coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      ordinal_type size() const { return th->size();}

      //! Returns true if polynomial has size >= sz
      bool hasFastAccess(ordinal_type sz) const { return th->size()>=sz;}

      //! Returns Hermite coefficient array
      const value_type* coeff() const { return th->coeff();}

      //! Returns Hermite coefficient array
      value_type* coeff() { return th->coeff();}

      //! Returns degree \c i term with bounds checking
      value_type coeff(ordinal_type i) const { 
	value_type tmp= i<th->size() ? (*th)[i]:value_type(0.); return tmp;}
    
      //! Returns degree \c i term without bounds checking
      value_type& fastAccessCoeff(ordinal_type i) { return (*th)[i];}

      //! Returns degree \c i term without bounds checking
      value_type fastAccessCoeff(ordinal_type i) const { return (*th)[i];}
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      OrthogPoly<T> operator + () const;

      //! Unary-minus operator
      OrthogPoly<T> operator - () const;

      //! Addition-assignment operator with constant right-hand-side
      OrthogPoly<T>& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      OrthogPoly<T>& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      OrthogPoly<T>& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      OrthogPoly<T>& operator /= (const value_type& x);

      //! Addition-assignment operator with Hermite right-hand-side
      OrthogPoly<T>& operator += (const OrthogPoly<T>& x);

      //! Subtraction-assignment operator with Hermite right-hand-side
      OrthogPoly<T>& operator -= (const OrthogPoly<T>& x);
  
      //! Multiplication-assignment operator with Hermite right-hand-side
      OrthogPoly<T>& operator *= (const OrthogPoly<T>& x);

      //! Division-assignment operator with Hermite right-hand-side
      OrthogPoly<T>& operator /= (const OrthogPoly<T>& x);

      //@}

      //! Get underlying Hermite polynomial
      const Stokhos::OrthogPolyApprox<int,value_type>& getOrthogPolyApprox() const 
      { return *th; }

      //! Get underlying Hermite polynomial
      Stokhos::OrthogPolyApprox<int,value_type>& getOrthogPolyApprox() { 
	return *th; }

    public:

      //! Expansion type
      static Teuchos::RCP<expansion_type> expansion;

    protected:

      Sacado::Handle< Stokhos::OrthogPolyApprox<int,value_type> > th;

    }; // class Hermite

    // Operations
    template <typename T> OrthogPoly<T> 
    operator+(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator+(const typename OrthogPoly<T>::value_type& a, 
	      const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator+(const OrthogPoly<T>& a, 
	      const typename OrthogPoly<T>::value_type& b);

    template <typename T> OrthogPoly<T> 
    operator-(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator-(const typename OrthogPoly<T>::value_type& a, 
	      const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator-(const OrthogPoly<T>& a, 
	      const typename OrthogPoly<T>::value_type& b);

    template <typename T> OrthogPoly<T> 
    operator*(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator*(const typename OrthogPoly<T>::value_type& a, 
	      const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator*(const OrthogPoly<T>& a, 
	      const typename OrthogPoly<T>::value_type& b);

    template <typename T> OrthogPoly<T> 
    operator/(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator/(const typename OrthogPoly<T>::value_type& a, 
	      const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    operator/(const OrthogPoly<T>& a, 
	      const typename OrthogPoly<T>::value_type& b);

    template <typename T> OrthogPoly<T> 
    exp(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    log(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    log10(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    sqrt(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    pow(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    pow(const T& a, 
	const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    pow(const OrthogPoly<T>& a, 
	const T& b);

    template <typename T> OrthogPoly<T> 
    cos(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    sin(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    tan(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    cosh(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T>
    sinh(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    tanh(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    acos(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    asin(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    atan(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    atan2(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    atan2(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    atan2(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b);

    template <typename T> OrthogPoly<T> 
    acosh(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    asinh(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    atanh(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    abs(const OrthogPoly<T>& a);
    
    template <typename T> OrthogPoly<T> 
    fabs(const OrthogPoly<T>& a);

    template <typename T> OrthogPoly<T> 
    max(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    max(const typename OrthogPoly<T>::value_type& a, 
	const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    max(const OrthogPoly<T>& a, 
	const typename OrthogPoly<T>::value_type& b);

    template <typename T> OrthogPoly<T> 
    min(const OrthogPoly<T>& a, const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    min(const typename OrthogPoly<T>::value_type& a, 
	const OrthogPoly<T>& b);

    template <typename T> OrthogPoly<T> 
    min(const OrthogPoly<T>& a, 
	const typename OrthogPoly<T>::value_type& b);

    template <typename T> bool 
    operator==(const OrthogPoly<T>& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator==(const typename OrthogPoly<T>::value_type& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator==(const OrthogPoly<T>& a,
	       const typename OrthogPoly<T>::value_type& b);

    template <typename T> bool 
    operator!=(const OrthogPoly<T>& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator!=(const typename OrthogPoly<T>::value_type& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator!=(const OrthogPoly<T>& a,
	       const typename OrthogPoly<T>::value_type& b);

    template <typename T> bool 
    operator<=(const OrthogPoly<T>& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator<=(const typename OrthogPoly<T>::value_type& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator<=(const OrthogPoly<T>& a,
	       const typename OrthogPoly<T>::value_type& b);

    template <typename T> bool 
    operator>=(const OrthogPoly<T>& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator>=(const typename OrthogPoly<T>::value_type& a,
	       const OrthogPoly<T>& b);

    template <typename T> bool 
    operator>=(const OrthogPoly<T>& a,
	       const typename OrthogPoly<T>::value_type& b);

    template <typename T> bool 
    operator<(const OrthogPoly<T>& a,
	      const OrthogPoly<T>& b);

    template <typename T> bool 
    operator<(const typename OrthogPoly<T>::value_type& a,
	      const OrthogPoly<T>& b);

    template <typename T> bool 
    operator<(const OrthogPoly<T>& a,
	      const typename OrthogPoly<T>::value_type& b);

    template <typename T> bool 
    operator>(const OrthogPoly<T>& a,
	      const OrthogPoly<T>& b);

    template <typename T> bool 
    operator>(const typename OrthogPoly<T>::value_type& a,
	      const OrthogPoly<T>& b);

    template <typename T> bool 
    operator>(const OrthogPoly<T>& a,
	      const typename OrthogPoly<T>::value_type& b);

    template <typename T> std::ostream& 
    operator << (std::ostream& os, const OrthogPoly<T>& a);

  } // namespace PCE

} // namespace Sacado

#include "Sacado_PCE_OrthogPolyTraits.hpp"
#include "Sacado_PCE_OrthogPolyImp.hpp"

#endif // HAVE_SACADO_STOKHOS

#endif // SACADO_PCE_UNIVARIATEHERMITE_HPP
