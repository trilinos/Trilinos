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

#ifndef SACADO_MP_VECTOR_HPP
#define SACADO_MP_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_dummy_arg.hpp"

#include <ostream>	// for std::ostream

namespace Sacado {

  //! Namespace for multipoint classes
  namespace MP {

    //! Wrapper for a generic expression template
    /*!
     * This class is used to limit the overload set for building up 
     * expressions.  Each expression object should derive from this
     * using CRTP:
     *
     * \code
     * class T : public Expr<T> { ... };
     * \endcode
     *
     * In this case the default implementation here should be correct for
     * any expression class.  If not, an expression class is free to change
     * the implementation through partial specialization.
     */
    template <typename T> class Expr {
    public:

      //! Typename of derived object, returned by derived()
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>
       */
      typedef T derived_type;

      //! Return derived object
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>.  This will only compile if this infact the case.
       */
      const derived_type& derived() const {
	return static_cast<const derived_type&>(*this);
      }

    };

    //! Vectorized evaluation class
    template <typename T, typename Storage> 
    class Vector : public Expr< Vector<T,Storage> > {
    public:

      //! Typename of values
      typedef T value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<T>::type scalar_type;

      //! Typename of ordinals
      typedef int ordinal_type;

      //! Typename of storage class
      typedef Storage storage_type;

      typedef typename storage_type::pointer pointer;
      typedef typename storage_type::const_pointer const_pointer;
      typedef typename storage_type::reference reference;
      typedef typename storage_type::const_reference const_reference;

      //! Turn Vector into a meta-function class usable with mpl::apply
      template <typename S> 
      struct apply {
	typedef typename Sacado::mpl::apply<Storage,ordinal_type,S>::type storage_type;
	typedef Vector<S,storage_type> type;
      };

      //! Number of arguments
      static const int num_args = 1;

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      Vector();

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      Vector(const value_type& x);

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      Vector(ordinal_type sz, const value_type& x);

      //! Copy constructor
      Vector(const Vector& x);

      //! Copy constructor from any Expression object
      template <typename S> Vector(const Expr<S>& x);

      //! Destructor
      ~Vector() {}

      //! Initialize coefficients to value
      void init(const T& v) { s.init(v); }

      //! Initialize coefficients to an array of values
      void init(const T* v) { s.init(v); }

      //! Initialize coefficients from an Vector with different storage
      template <typename S>
      void init(const Vector<T,S>& v) { 
	s.init(v.s.coeff(), v.s.size()); 
      }

      //! Load coefficients to an array of values
      void load(T* v) { s.load(v); }

      //! Load coefficients into an Vector with different storage
      template <typename S>
      void load(Vector<T,S>& v) { s.load(v.s.coeff()); }

      //! Reset size
      /*!
       * Coefficients are preserved.  
       */
      void reset(ordinal_type sz);

      //! Prepare vector for writing 
      /*!
       * This method prepares the vector for writing through coeff() and 
       * fastAccessCoeff() member functions.  It ensures the handle for the
       * coefficients is not shared among any other vector.  
       * If the handle is not shared it does nothing, so there
       * is no cost in calling this method in this case.  If the handle is 
       * shared and this method is not called, any changes to the coefficients
       * by coeff() or fastAccessCoeff() may change other vector objects.
       */
      void copyForWrite() {  }

      //! Returns whether two ETV objects have the same values
      template <typename S>
      bool isEqualTo(const Expr<S>& xx) const {
	const typename Expr<S>::derived_type& x = xx.derived();
	typedef IsEqual<value_type> IE;
	if (x.size() != this->size()) return false;
	bool eq = true;
	for (int i=0; i<this->size(); i++)
	  eq = eq && IE::eval(x.coeff(i), this->coeff(i));
	return eq;
      }

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      Vector& operator=(const value_type& val);

      //! Assignment with Vector right-hand-side
      Vector& operator=(const Vector& x);

      //! Assignment with any expression right-hand-side
      template <typename S> 
      Vector& operator=(const Expr<S>& x);

      //@}

      /*!
       * Accessor methods
       */

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const_reference val() const { return s[0]; }

      //! Returns value
      reference val() { return s[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      ordinal_type size() const { return s.size();}

      //! Returns true if polynomial has size >= sz
      bool hasFastAccess(ordinal_type sz) const { return s.size()>=sz;}

      //! Returns Hermite coefficient array
      const_pointer coeff() const { return s.coeff();}

      //! Returns Hermite coefficient array
      pointer coeff() { return s.coeff();}

      //! Returns degree \c i term with bounds checking
      value_type coeff(ordinal_type i) const { 
	return i<s.size() ? s[i] : s[0]; }
    
      //! Returns degree \c i term without bounds checking
      reference fastAccessCoeff(ordinal_type i) { return s[i];}

      //! Returns degree \c i term without bounds checking
      value_type fastAccessCoeff(ordinal_type i) const { return s[i];}
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      Vector& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      Vector& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      Vector& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      Vector& operator /= (const value_type& x);

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator += (const Expr<S>& x) {
	*this = *this + x;
	return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator -= (const Expr<S>& x) {
	*this = *this - x;
	return *this;
      }
  
      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator *= (const Expr<S>& x) {
	*this = *this * x;
	return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator /= (const Expr<S>& x) {
	*this = *this / x;
	return *this;
      }

      //@}

      std::string name() const { return "x"; }

    protected:

      Storage s;

    }; // class Vector

  } // namespace MP

} // namespace Sacado

#include "Sacado_MP_VectorTraits.hpp"
#include "Sacado_MP_VectorImp.hpp"
#include "Sacado_MP_VectorOps.hpp"

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_MP_VECTOR_HPP
