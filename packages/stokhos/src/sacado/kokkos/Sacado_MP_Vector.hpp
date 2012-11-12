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

#include <ostream>	// for std::ostream

#include "KokkosArray_Macros.hpp"

#include "Sacado_MP_VectorTraits.hpp"
#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_mpl_range_c.hpp"
#include "Stokhos_mpl_for_each.hpp"

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
    template <typename T, typename node> 
    class Expr {
    public:

      //! Node type
      typedef node node_type;

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
      KOKKOSARRAY_INLINE_FUNCTION
      const derived_type& derived() const {
	return static_cast<const derived_type&>(*this);
      }

    };

    //! Vectorized evaluation class
    template <typename Storage, typename Node> 
    class Vector : public Expr< Vector<Storage,Node>, Node > {
    public:

      //! Typename of storage class
      typedef Storage storage_type;

       //! Node type
      typedef Node node_type;

      typedef typename storage_type::value_type value_type;
      typedef typename storage_type::ordinal_type ordinal_type;
      typedef typename storage_type::pointer pointer;
      typedef typename storage_type::const_pointer const_pointer;
      typedef typename storage_type::reference reference;
      typedef typename storage_type::const_reference const_reference;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Turn Vector into a meta-function class usable with mpl::apply
      template <typename S> 
      struct apply {
	typedef typename Sacado::mpl::apply<Storage,ordinal_type,S>::type new_storage_type;
	typedef Vector<new_storage_type,node_type> type;
      };

      //! Number of arguments
      static const int num_args = 1;

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      KOKKOSARRAY_INLINE_FUNCTION
      Vector() : s(1) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      KOKKOSARRAY_INLINE_FUNCTION
      Vector(const value_type& x) : s(1) { s.init(x); }

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      KOKKOSARRAY_INLINE_FUNCTION
      Vector(ordinal_type sz, const value_type& x) : s(sz,x) {}

      //! Copy constructor
      KOKKOSARRAY_INLINE_FUNCTION
      Vector(const Vector& x) : s(x.s) {}

      //! Copy constructor from any Expression object
      template <typename S> 
      KOKKOSARRAY_INLINE_FUNCTION
      Vector(const Expr<S,node_type>& xx) : 
	s(xx.derived().size()) {
	typedef typename Expr<S,node_type>::derived_type expr_type;
	const expr_type& x = xx.derived();
	
	if (storage_type::is_static) {
	  typedef Sacado::mpl::range_c< int, 0, 
					storage_type::static_size > range_type;
	  StaticOp<expr_type> op(s,x);
	  Stokhos::mpl::for_each< range_type, node_type > f(op);
	}
	else if (x.hasFastAccess(s.size())) {
	  for (ordinal_type i=0; i<s.size(); i++)
	    s[i] = x.fastAccessCoeff(i);
	}
	else {
	  for (ordinal_type i=0; i<s.size(); i++)
	    s[i] = x.coeff(i);
	}
      }

      //! Destructor
      KOKKOSARRAY_INLINE_FUNCTION
      ~Vector() {}

      //! Initialize coefficients to value
      KOKKOSARRAY_INLINE_FUNCTION
      void init(const value_type& v) { s.init(v); }

      //! Initialize coefficients to an array of values
      KOKKOSARRAY_INLINE_FUNCTION
      void init(const value_type* v) { s.init(v); }

      //! Initialize coefficients from an Vector with different storage
      template <typename S, typename N>
      KOKKOSARRAY_INLINE_FUNCTION
      void init(const Vector<S,N>& v) { 
	s.init(v.s.coeff(), v.s.size()); 
      }

      //! Load coefficients to an array of values
      KOKKOSARRAY_INLINE_FUNCTION
      void load(value_type* v) { s.load(v); }

      //! Load coefficients into an Vector with different storage
      template <typename S, typename N>
      KOKKOSARRAY_INLINE_FUNCTION
      void load(Vector<S,N>& v) { s.load(v.s.coeff()); }

      //! Reset size
      /*!
       * Coefficients are preserved.  
       */
      KOKKOSARRAY_INLINE_FUNCTION
      void reset(ordinal_type sz_new) {
	ordinal_type sz = this->size();
	s.resize(sz_new);
	if (sz == 1 && sz_new > sz)
	  for (ordinal_type i=1; i<sz_new; i++)
	    s[i] = s[0];
      }

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
      KOKKOSARRAY_INLINE_FUNCTION
      void copyForWrite() {  }

      //! Returns whether two ETV objects have the same values
      template <typename S>
      KOKKOSARRAY_INLINE_FUNCTION
      bool isEqualTo(const Expr<S,node_type>& xx) const {
	const typename Expr<S,node_type>::derived_type& x = xx.derived();
	typedef IsEqual<value_type> IE;
	if (x.size() != this->size()) return false;
	bool eq = true;
	for (ordinal_type i=0; i<this->size(); i++)
	  eq = eq && IE::eval(x.coeff(i), this->coeff(i));
	return eq;
      }

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator=(const value_type& x) {
	s.init(x);
	return *this;
      }

      //! Assignment with Vector right-hand-side
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator=(const Vector& x) {
	s = x.s;
	return *this;
      }

      //! Assignment with any expression right-hand-side
      template <typename S> 
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator=(const Expr<S,node_type>& xx) {
	typedef typename Expr<S,node_type>::derived_type expr_type;
	const expr_type& x = xx.derived();
	
	this->reset(x.size());
	if (storage_type::is_static) {
	  typedef Sacado::mpl::range_c< int, 0, 
					storage_type::static_size > range_type;
	  StaticOp<expr_type > op(s,x);
	  Stokhos::mpl::for_each< range_type, node_type > f(op);
	}
	else if (x.hasFastAccess(s.size())) {
	  for (ordinal_type i=0; i<s.size(); i++)
	    s[i] = x.fastAccessCoeff(i);
	}
	else {
	  for (ordinal_type i=0; i<s.size(); i++)
	    s[i] = x.coeff(i);
	}
	return *this;
      }

      //@}

      /*!
       * Accessor methods
       */

      //! Returns storage object
      KOKKOSARRAY_INLINE_FUNCTION
      const storage_type& storage() const { return s; }

      //! Returns storage object
      KOKKOSARRAY_INLINE_FUNCTION
      storage_type& storage() { return s; }

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      KOKKOSARRAY_INLINE_FUNCTION
      const_reference val() const { return s[0]; }

      //! Returns value
      KOKKOSARRAY_INLINE_FUNCTION
      reference val() { return s[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      KOKKOSARRAY_INLINE_FUNCTION
      ordinal_type size() const { return s.size();}

      //! Returns true if polynomial has size >= sz
      KOKKOSARRAY_INLINE_FUNCTION
      bool hasFastAccess(ordinal_type sz) const { return s.size()>=sz;}

      //! Returns Hermite coefficient array
      KOKKOSARRAY_INLINE_FUNCTION
      const_pointer coeff() const { return s.coeff();}

      //! Returns Hermite coefficient array
      KOKKOSARRAY_INLINE_FUNCTION
      pointer coeff() { return s.coeff();}

      //! Returns degree \c i term with bounds checking
      KOKKOSARRAY_INLINE_FUNCTION
      value_type coeff(ordinal_type i) const { 
	return i<s.size() ? s[i] : s[0]; }
    
      //! Returns degree \c i term without bounds checking
      KOKKOSARRAY_INLINE_FUNCTION
      reference fastAccessCoeff(ordinal_type i) { return s[i];}

      //! Returns degree \c i term without bounds checking
      KOKKOSARRAY_INLINE_FUNCTION
      value_type fastAccessCoeff(ordinal_type i) const { return s[i];}

      template <int i>
      KOKKOSARRAY_INLINE_FUNCTION
      value_type getCoeff() const { return s.template getCoeff<i>(); }

      template <int i>
      KOKKOSARRAY_INLINE_FUNCTION
      reference getCoeff() { return s.template getCoeff<i>(); }
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator += (const value_type& x) {
	for (ordinal_type i=0; i<s.size(); i++)
	  s[i] += x;
	return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator -= (const value_type& x) {
	for (ordinal_type i=0; i<s.size(); i++)
	  s[i] -= x;
	return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator *= (const value_type& x) {
	for (ordinal_type i=0; i<s.size(); i++)
	  s[i] *= x;
	return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator /= (const value_type& x) {
	for (ordinal_type i=0; i<s.size(); i++)
	  s[i] /= x;
	return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S> 
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator += (const Expr<S,node_type>& x) {
	*this = *this + x;
	return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S> 
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator -= (const Expr<S,node_type>& x) {
	*this = *this - x;
	return *this;
      }
  
      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S> 
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator *= (const Expr<S,node_type>& x) {
	*this = *this * x;
	return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S> 
      KOKKOSARRAY_INLINE_FUNCTION
      Vector& operator /= (const Expr<S,node_type>& x) {
	*this = *this / x;
	return *this;
      }

      //@}

      KOKKOSARRAY_INLINE_FUNCTION
      std::string name() const { return "x"; }

    protected:

      Storage s;

      template <typename expr_type>
      struct StaticOp {
	storage_type& s;
	const expr_type& x;

	KOKKOSARRAY_INLINE_FUNCTION
	StaticOp(storage_type& s_, const expr_type& x_) : s(s_), x(x_) {}

	template <typename ArgT>
	KOKKOSARRAY_INLINE_FUNCTION
	void operator() (ArgT arg) const {
	  const int Arg = ArgT::value;
	  s.template getCoeff<Arg>() = x.template getCoeff<Arg>();
	}

      };

    }; // class Vector

    //! Type for storing nodes in expression graph
    /*!
     * Since expression nodes are returned by value in the overloaded
     * operators, we can't store them by reference in general.
     */
    template <typename T> struct const_expr_ref {
      typedef const T type;
    };

    //! Type for storing nodes in expression graph
    /*!
     * Specialization for leaf-nodes, which can be stored by reference
     * since they are an argument to the expression.
     */
    template <typename S, typename N> struct const_expr_ref< Vector<S,N> > {
      typedef const Vector<S,N>& type;
    };

    template <typename Storage, typename Node>
    KOKKOSARRAY_INLINE_FUNCTION
    std::ostream& 
    operator << (std::ostream& os, 
		 const Vector<Storage,Node>& a)
    {
      typedef typename Vector<Storage,Node>::ordinal_type ordinal_type;
      
      os << "[ ";
      
      for (ordinal_type i=0; i<a.size(); i++) {
	os << a.coeff(i) << " ";
      }

      os << "]\n";
      return os;
    }

  } // namespace MP

} // namespace Sacado

#include "Sacado_MP_Vector_ops.hpp"

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_MP_VECTOR_HPP
