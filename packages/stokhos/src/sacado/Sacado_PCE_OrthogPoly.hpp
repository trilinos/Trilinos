// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PCE_ORTHOGPOLY_HPP
#define SACADO_PCE_ORTHOGPOLY_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include "Teuchos_RCP.hpp"

#include "Sacado_Traits.hpp"
#include "Sacado_Handle.hpp"
#include "Sacado_mpl_apply.hpp"

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
    template <typename T, typename Storage > 
    class OrthogPoly {
    public:

      //! Typename of values
      typedef T value_type;

      //! Typename of ordinals
      typedef int ordinal_type;

      //! Typename of storage class
      typedef Storage storage_type;

      //! Basis type
      typedef Stokhos::OrthogPolyBasis<ordinal_type,T> basis_type;

      //! Expansion type
      typedef Stokhos::OrthogPolyExpansion<ordinal_type,T,Storage> expansion_type;

      //! Stokhos approximation type
      typedef Stokhos::OrthogPolyApprox<ordinal_type,T,Storage> approx_type;

      typedef typename approx_type::pointer pointer;
      typedef typename approx_type::const_pointer const_pointer;
      typedef typename approx_type::reference reference;
      typedef typename approx_type::const_reference const_reference;

      //! Turn OrthogPoly into a meta-function class usable with mpl::apply
      template <typename S> 
      struct apply {
	typedef typename Sacado::mpl::apply<Storage,ordinal_type,S>::type storage_type;
	typedef OrthogPoly<S,storage_type> type;
      };

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      OrthogPoly();

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      OrthogPoly(const value_type& x);

      //! Constructor with expansion \c expansion (General case)
      /*!
       * Creates array of correct size and initializes coeffiencts to 0.
       */
      OrthogPoly(const Teuchos::RCP<expansion_type>& expansion);

      //! Constructor with expansion \c expansion and specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      OrthogPoly(const Teuchos::RCP<expansion_type>& expansion,
		 ordinal_type sz);

      //! Copy constructor
      OrthogPoly(const OrthogPoly& x);

      //! Destructor
      ~OrthogPoly();

      //! Initialize coefficients to value
      void init(const T& v) { th->init(v); }

      //! Initialize coefficients to an array of values
      void init(const T* v) { th->init(v); }

      //! Initialize coefficients from an OrthogPoly with different storage
      template <typename S>
      void init(const OrthogPoly<T,S>& v) { th->init(v.getOrthogPolyApprox()); }

      //! Load coefficients to an array of values
      void load(T* v) { th->load(v); }

      //! Load coefficients into an OrthogPoly with different storage
      template <typename S>
      void load(OrthogPoly<T,S>& v) { th->load(v.getOrthogPolyApprox()); }

      //! Reset expansion
      /*!
       * May change size of array.  Coefficients are preserved.  
       */
      void reset(const Teuchos::RCP<expansion_type>& expansion);

      //! Reset expansion and size
      /*!
       * Coefficients are preserved.  
       */
      void reset(const Teuchos::RCP<expansion_type>& expansion,
		 ordinal_type sz);

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

      //! Evaluate polynomial approximation at a point
      value_type evaluate(const Teuchos::Array<value_type>& point) const;

      //! Evaluate polynomial approximation at a point with given basis values
      value_type evaluate(const Teuchos::Array<value_type>& point,
                          const Teuchos::Array<value_type>& bvals) const;

      //! Compute mean of expansion
      value_type mean() const {return th->mean(); }

      //! Compute standard deviation of expansion
      value_type standard_deviation() const { return th->standard_deviation(); }

      //! Compute the two-norm of expansion
      value_type two_norm() const { return th->two_norm(); }

      //! Compute the squared two-norm of expansion
      value_type two_norm_squared() const { return th->two_norm_squared(); }

      //! Compute the L2 inner product of 2 PCEs
      value_type inner_product(const OrthogPoly& b) const { 
	return th->inner_product(b.getOrthogPolyApprox()); }

      //! Print approximation in basis
      std::ostream& print(std::ostream& os) const { return th->print(os); }

      //! Returns whether two PCE objects have the same values
      bool isEqualTo(const OrthogPoly& x) const;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      OrthogPoly<T,Storage>& operator=(const value_type& val);

      //! Assignment with OrthogPoly right-hand-side
      OrthogPoly<T,Storage>& operator=(const OrthogPoly<T,Storage>& x);

      //@}

      /*!
       * Accessor methods
       */
      //@{

      //! Get basis
      Teuchos::RCP<const basis_type> basis() const { return th->basis(); }

      //! Get expansion
      Teuchos::RCP<expansion_type> expansion() const { return expansion_; }

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const_reference val() const { return (*th)[0]; }

      //! Returns value
      reference val() { return (*th)[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      ordinal_type size() const { return th->size();}

      //! Returns true if polynomial has size >= sz
      bool hasFastAccess(ordinal_type sz) const { return th->size()>=sz;}

      //! Returns Hermite coefficient array
      const_pointer coeff() const { return th->coeff();}

      //! Returns Hermite coefficient array
      pointer coeff() { return th->coeff();}

      //! Returns degree \c i term with bounds checking
      value_type coeff(ordinal_type i) const { 
	value_type tmp= i<th->size() ? (*th)[i]:value_type(0.); return tmp;}
    
      //! Returns degree \c i term without bounds checking
      reference fastAccessCoeff(ordinal_type i) { return (*th)[i];}

      //! Returns degree \c i term without bounds checking
      value_type fastAccessCoeff(ordinal_type i) const { return (*th)[i];}

      //! Get coefficient term for given dimension and order
      reference term(ordinal_type dimension, ordinal_type order) {
	return th->term(dimension, order); }

      //! Get coefficient term for given dimension and order
      const_reference term(ordinal_type dimension, ordinal_type order) const {
	return th->term(dimension, order); }

      //! Get orders for a given term
      Teuchos::Array<ordinal_type> order(ordinal_type term) const {
	return th->order(term); }
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      OrthogPoly<T,Storage> operator + () const;

      //! Unary-minus operator
      OrthogPoly<T,Storage> operator - () const;

      //! Addition-assignment operator with constant right-hand-side
      OrthogPoly<T,Storage>& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      OrthogPoly<T,Storage>& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      OrthogPoly<T,Storage>& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      OrthogPoly<T,Storage>& operator /= (const value_type& x);

      //! Addition-assignment operator with Hermite right-hand-side
      OrthogPoly<T,Storage>& operator += (const OrthogPoly<T,Storage>& x);

      //! Subtraction-assignment operator with Hermite right-hand-side
      OrthogPoly<T,Storage>& operator -= (const OrthogPoly<T,Storage>& x);
  
      //! Multiplication-assignment operator with Hermite right-hand-side
      OrthogPoly<T,Storage>& operator *= (const OrthogPoly<T,Storage>& x);

      //! Division-assignment operator with Hermite right-hand-side
      OrthogPoly<T,Storage>& operator /= (const OrthogPoly<T,Storage>& x);

      //@}

      //! Get underlying Stokhos::OrthogPolyApprox
      const approx_type& getOrthogPolyApprox() const { return *th; }

      //! Get underlying Stokhos::OrthogPolyApprox
      approx_type& getOrthogPolyApprox() { return *th; }

    protected:

      //! Expansion class
      Teuchos::RCP<expansion_type> expansion_;

      //! Static constant expansion class for constants
      static Teuchos::RCP<expansion_type> const_expansion_;

      Sacado::Handle< Stokhos::OrthogPolyApprox<int,value_type,Storage> > th;

    }; // class Hermite

    // Operations
    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator+(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator+(const typename OrthogPoly<T,Storage>::value_type& a, 
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator+(const OrthogPoly<T,Storage>& a, 
	      const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator-(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator-(const typename OrthogPoly<T,Storage>::value_type& a, 
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator-(const OrthogPoly<T,Storage>& a, 
	      const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator*(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator*(const typename OrthogPoly<T,Storage>::value_type& a, 
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator*(const OrthogPoly<T,Storage>& a, 
	      const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator/(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator/(const typename OrthogPoly<T,Storage>::value_type& a, 
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    operator/(const OrthogPoly<T,Storage>& a, 
	      const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    exp(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    log(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> void
    log(OrthogPoly<T,Storage>& c, const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    log10(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    sqrt(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    pow(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    pow(const T& a, 
	const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    pow(const OrthogPoly<T,Storage>& a, 
	const T& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    cos(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    sin(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    tan(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    cosh(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage>
    sinh(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    tanh(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    acos(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    asin(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    atan(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    atan2(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    atan2(const typename OrthogPoly<T,Storage>::value_type& a, 
	  const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    atan2(const OrthogPoly<T,Storage>& a, 
	  const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    acosh(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    asinh(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    atanh(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    abs(const OrthogPoly<T,Storage>& a);
    
    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    fabs(const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    max(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    max(const typename OrthogPoly<T,Storage>::value_type& a, 
	const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    max(const OrthogPoly<T,Storage>& a, 
	const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    min(const OrthogPoly<T,Storage>& a, const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    min(const typename OrthogPoly<T,Storage>::value_type& a, 
	const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> OrthogPoly<T,Storage> 
    min(const OrthogPoly<T,Storage>& a, 
	const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> bool 
    operator==(const OrthogPoly<T,Storage>& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator==(const typename OrthogPoly<T,Storage>::value_type& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator==(const OrthogPoly<T,Storage>& a,
	       const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> bool 
    operator!=(const OrthogPoly<T,Storage>& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator!=(const typename OrthogPoly<T,Storage>::value_type& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator!=(const OrthogPoly<T,Storage>& a,
	       const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> bool 
    operator<=(const OrthogPoly<T,Storage>& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator<=(const typename OrthogPoly<T,Storage>::value_type& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator<=(const OrthogPoly<T,Storage>& a,
	       const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> bool 
    operator>=(const OrthogPoly<T,Storage>& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator>=(const typename OrthogPoly<T,Storage>::value_type& a,
	       const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator>=(const OrthogPoly<T,Storage>& a,
	       const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> bool 
    operator<(const OrthogPoly<T,Storage>& a,
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator<(const typename OrthogPoly<T,Storage>::value_type& a,
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator<(const OrthogPoly<T,Storage>& a,
	      const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> bool 
    operator>(const OrthogPoly<T,Storage>& a,
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator>(const typename OrthogPoly<T,Storage>::value_type& a,
	      const OrthogPoly<T,Storage>& b);

    template <typename T, typename Storage> bool 
    operator>(const OrthogPoly<T,Storage>& a,
	      const typename OrthogPoly<T,Storage>::value_type& b);

    template <typename T, typename Storage> std::ostream& 
    operator << (std::ostream& os, const OrthogPoly<T,Storage>& a);

    template <typename T, typename Storage> std::istream& 
    operator >> (std::istream& os, OrthogPoly<T,Storage>& a);

  } // namespace PCE

} // namespace Sacado

#include "Sacado_PCE_OrthogPolyTraits.hpp"
#include "Sacado_PCE_OrthogPolyImp.hpp"

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_PCE_ORTHOGPOLY_HPP
