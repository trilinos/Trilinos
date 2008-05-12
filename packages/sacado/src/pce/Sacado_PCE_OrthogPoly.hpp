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
    template <typename BasisT> 
    class OrthogPoly {
    public:

      //! Typename of values
      typedef typename BasisT::value_type value_type;

      //! Basis type
      typedef BasisT basis_type;

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
      OrthogPoly(unsigned int sz, const value_type& x);

      //! Constructor with size \c sz
      /*!
       * Initializes all components to zero
       */
      OrthogPoly(unsigned int sz);

      //! Copy constructor
      OrthogPoly(const OrthogPoly& x);

      //! Destructor
      ~OrthogPoly();

      //! Resize polynomial to size \c sz
      /*!
       * Coefficients are preserved.
       */
      void resize(unsigned int sz);

      //! Reserve space for polynomial of size \c sz
      /*!
       * Coefficients are preserved.
       */
      void reserve(unsigned int sz);

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
      static void initExpansion(
	const Teuchos::RCP< Stokhos::OrthogPolyExpansion<BasisT> >& e);

      //! Write coefficients in standard basis
      Stokhos::Polynomial<value_type> toStandardBasis() const;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      OrthogPoly<BasisT>& operator=(const value_type& val);

      //! Assignment with OrthogPoly right-hand-side
      OrthogPoly<BasisT>& operator=(const OrthogPoly<BasisT>& x);

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const value_type& val() const { return th->coeff_[0];}

      //! Returns value
      value_type& val() { return th->coeff_[0];}

      //@}

      /*!
       * @name Hermite coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      unsigned int size() const { return th->size();}

      //! Returns true if polynomial has size >= sz
      bool hasFastAccess(unsigned int sz) const { return th->size()>=sz;}

      //! Returns Hermite coefficient array
      const value_type* coeff() const { return th->coeff();}

      //! Returns Hermite coefficient array
      value_type* coeff() { return th->coeff();}

      //! Returns degree \c i term with bounds checking
      value_type coeff(unsigned int i) const { 
	value_type tmp= i<th->size() ? (*th)[i]:value_type(0.); return tmp;}
    
      //! Returns degree \c i term without bounds checking
      value_type& fastAccessCoeff(unsigned int i) { return (*th)[i];}

      //! Returns degree \c i term without bounds checking
      value_type fastAccessCoeff(unsigned int i) const { return (*th)[i];}
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      OrthogPoly<BasisT> operator + () const;

      //! Unary-minus operator
      OrthogPoly<BasisT> operator - () const;

      //! Addition-assignment operator with constant right-hand-side
      OrthogPoly<BasisT>& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      OrthogPoly<BasisT>& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      OrthogPoly<BasisT>& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      OrthogPoly<BasisT>& operator /= (const value_type& x);

      //! Addition-assignment operator with Hermite right-hand-side
      OrthogPoly<BasisT>& operator += (const OrthogPoly<BasisT>& x);

      //! Subtraction-assignment operator with Hermite right-hand-side
      OrthogPoly<BasisT>& operator -= (const OrthogPoly<BasisT>& x);
  
      //! Multiplication-assignment operator with Hermite right-hand-side
      OrthogPoly<BasisT>& operator *= (const OrthogPoly<BasisT>& x);

      //! Division-assignment operator with Hermite right-hand-side
      OrthogPoly<BasisT>& operator /= (const OrthogPoly<BasisT>& x);

      //@}

      //! Get underlying Hermite polynomial
      const Stokhos::OrthogPolyApprox<value_type>& getOrthogPolyApprox() const 
      { return *th; }

      //! Get underlying Hermite polynomial
      Stokhos::OrthogPolyApprox<value_type>& getOrthogPolyApprox() { 
	return *th; }

    public:

      //! Expansion type
      static Teuchos::RCP< Stokhos::OrthogPolyExpansion<BasisT> > expansion;

    protected:

      Sacado::Handle< Stokhos::OrthogPolyApprox<value_type> > th;

    }; // class Hermite

    // Operations
    template <typename BasisT> OrthogPoly<BasisT> 
    operator+(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator+(const typename OrthogPoly<BasisT>::value_type& a, 
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator+(const OrthogPoly<BasisT>& a, 
	      const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator-(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator-(const typename OrthogPoly<BasisT>::value_type& a, 
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator-(const OrthogPoly<BasisT>& a, 
	      const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator*(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator*(const typename OrthogPoly<BasisT>::value_type& a, 
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator*(const OrthogPoly<BasisT>& a, 
	      const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator/(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator/(const typename OrthogPoly<BasisT>::value_type& a, 
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    operator/(const OrthogPoly<BasisT>& a, 
	      const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    exp(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    log(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    log10(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    sqrt(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    pow(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    pow(const typename OrthogPoly<BasisT>::value_type& a, 
	const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    pow(const OrthogPoly<BasisT>& a, 
	const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    cos(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    sin(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    tan(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    cosh(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT>
    sinh(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    tanh(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    acos(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    asin(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    atan(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    atan2(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    atan2(const typename OrthogPoly<BasisT>::value_type& a, 
	  const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    atan2(const OrthogPoly<BasisT>& a, 
	  const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    acosh(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    asinh(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    atanh(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    abs(const OrthogPoly<BasisT>& a);
    
    template <typename BasisT> OrthogPoly<BasisT> 
    fabs(const OrthogPoly<BasisT>& a);

    template <typename BasisT> OrthogPoly<BasisT> 
    max(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    max(const typename OrthogPoly<BasisT>::value_type& a, 
	const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    max(const OrthogPoly<BasisT>& a, 
	const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    min(const OrthogPoly<BasisT>& a, const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    min(const typename OrthogPoly<BasisT>::value_type& a, 
	const OrthogPoly<BasisT>& b);

    template <typename BasisT> OrthogPoly<BasisT> 
    min(const OrthogPoly<BasisT>& a, 
	const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> bool 
    operator==(const OrthogPoly<BasisT>& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator==(const typename OrthogPoly<BasisT>::value_type& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator==(const OrthogPoly<BasisT>& a,
	       const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> bool 
    operator!=(const OrthogPoly<BasisT>& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator!=(const typename OrthogPoly<BasisT>::value_type& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator!=(const OrthogPoly<BasisT>& a,
	       const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> bool 
    operator<=(const OrthogPoly<BasisT>& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator<=(const typename OrthogPoly<BasisT>::value_type& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator<=(const OrthogPoly<BasisT>& a,
	       const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> bool 
    operator>=(const OrthogPoly<BasisT>& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator>=(const typename OrthogPoly<BasisT>::value_type& a,
	       const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator>=(const OrthogPoly<BasisT>& a,
	       const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> bool 
    operator<(const OrthogPoly<BasisT>& a,
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator<(const typename OrthogPoly<BasisT>::value_type& a,
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator<(const OrthogPoly<BasisT>& a,
	      const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> bool 
    operator>(const OrthogPoly<BasisT>& a,
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator>(const typename OrthogPoly<BasisT>::value_type& a,
	      const OrthogPoly<BasisT>& b);

    template <typename BasisT> bool 
    operator>(const OrthogPoly<BasisT>& a,
	      const typename OrthogPoly<BasisT>::value_type& b);

    template <typename BasisT> std::ostream& 
    operator << (std::ostream& os, const OrthogPoly<BasisT>& a);

  } // namespace PCE

} // namespace Sacado

#include "Sacado_PCE_OrthogPolyTraits.hpp"
#include "Sacado_PCE_OrthogPolyImp.hpp"

#endif // HAVE_SACADO_STOKHOS

#endif // SACADO_PCE_UNIVARIATEHERMITE_HPP
