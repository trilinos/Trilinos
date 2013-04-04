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

#ifndef SACADO_ETV_VECTOR2_HPP
#define SACADO_ETV_VECTOR2_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_dummy_arg.hpp"

#include <ostream>      // for std::ostream

namespace Sacado {

  //! Namespace for expression templated vector classes
  namespace ETV {

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all expression
     * template classes.
     */
    //template <typename ExprT> class Expr {};

    //! Vectorized evaluation class
    template <typename T, typename Storage >
    class Vector2Impl {
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


      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      Vector2Impl();

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      Vector2Impl(const value_type& x);

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      Vector2Impl(ordinal_type sz, const value_type& x);

      //! Copy constructor
      Vector2Impl(const Vector2Impl& x);

      //! Copy constructor from any Expression object
      template <typename S> Vector2Impl(const Expr<S>& x);

      //! Destructor
      ~Vector2Impl() {}

      //! Initialize coefficients to value
      void init(const T& v) { s.init(v); }

      //! Initialize coefficients to an array of values
      void init(const T* v) { s.init(v); }

      //! Initialize coefficients from an Vector2Impl with different storage
      template <typename S>
      void init(const Vector2Impl<T,S>& v) {
        s.init(v.s.coeff(), v.s.size());
      }

      //! Load coefficients to an array of values
      void load(T* v) { s.load(v); }

      //! Load coefficients into an Vector2Impl with different storage
      template <typename S>
      void load(Vector2Impl<T,S>& v) { s.load(v.s.coeff()); }

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
      bool isEqualTo(const Expr<S>& x) const {
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
      Vector2Impl& operator=(const value_type& val);

      //! Assignment with Vector2Impl right-hand-side
      Vector2Impl& operator=(const Vector2Impl& x);

      //! Assignment with any expression right-hand-side
      template <typename S>
      Vector2Impl& operator=(const Expr<S>& x);

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
      Vector2Impl& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      Vector2Impl& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      Vector2Impl& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      Vector2Impl& operator /= (const value_type& x);

      //@}

    protected:

      Storage s;

    }; // class Vector2Impl

    //! Vector2Impl expression template specialization
    /*!
     * This template class represents a simple Vector2Impl expression and
     * mixes-in the Vector2Impl interface and the expression template
     * interface.
     */
    template <typename T, typename Storage>
    class Expr< Vector2Impl<T,Storage> > :
        public Vector2Impl<T,Storage> {

    public:

      //! Typename of values
      typedef typename Vector2Impl<T,Storage>::value_type value_type;
      typedef typename Vector2Impl<T,Storage>::scalar_type scalar_type;
      typedef typename Vector2Impl<T,Storage>::storage_type storage_type;
      typedef typename Vector2Impl<T,Storage>::const_reference const_reference;

      //! Number of arguments
      static const int num_args = 1;

      //! Default constructor
      Expr() :
        Vector2Impl<T,Storage>() {}

      //! Constructor with supplied value \c x
      Expr(const T & x) :
        Vector2Impl<T,Storage>(x) {}

      //! Constructor with specified size \c sz
      Expr(typename Vector2Impl<T,Storage>::ordinal_type sz, const T& x) :
        Vector2Impl<T,Storage>(sz,x) {}

      //! Copy constructor
      Expr(const Expr& x) :
        Vector2Impl<T,Storage>(static_cast<const Vector2Impl<T,Storage>&>(x)) {}

      //! Copy constructor
      Expr(const Vector2Impl<T,Storage>& x) :
        Vector2Impl<T,Storage>(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Expr(const Expr<S>& x) :
        Vector2Impl<T,Storage>(x) {}

      //! Destructor
      ~Expr() {}

      std::string name() const { return "x"; }

    }; // class Expr<Vector2Impl>

  } // namespace ETV

} // namespace Sacado

#include "Sacado_ETV_ExpressionTraits.hpp"
#include "Sacado_ETV_Vector2Traits.hpp"
#include "Sacado_ETV_Vector2Imp.hpp"
#include "Sacado_ETV_VectorOps.hpp"

namespace Sacado {

  namespace ETV {

    //! Generalized polynomial chaos expansion class
    template <typename T, typename Storage>
    class Vector2 : public Expr< Vector2Impl<T,Storage> > {
    public:

      //! Typename of values
      typedef typename Vector2Impl<T,Storage>::value_type value_type;

      //! Typename of scalars
      typedef typename Vector2Impl<T,Storage>::scalar_type scalar_type;

      //! Typename of ordinals
      typedef typename Vector2Impl<T,Storage>::ordinal_type ordinal_type;

      //! Typename of storage class
      typedef typename Vector2Impl<T,Storage>::storage_type storage_type;

      typedef typename Vector2Impl<T,Storage>::pointer pointer;
      typedef typename Vector2Impl<T,Storage>::const_pointer const_pointer;
      typedef typename Vector2Impl<T,Storage>::reference reference;
      typedef typename Vector2Impl<T,Storage>::const_reference const_reference;

      //! Turn Vector2 into a meta-function class usable with mpl::apply
      template <typename S>
      struct apply {
        typedef typename Sacado::mpl::apply<Storage,ordinal_type,S>::type storage_type;
        typedef Vector2<S,storage_type> type;
      };

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      Vector2() :
        Expr< Vector2Impl<T,Storage> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      Vector2(const value_type& x) :
        Expr< Vector2Impl<T,Storage> >(x) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      Vector2(const typename dummy<value_type,scalar_type>::type& x) :
        Expr< Vector2Impl<T,Storage> >(value_type(x)) {}

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      explicit Vector2(ordinal_type sz, const value_type& x) :
        Expr< Vector2Impl<T,Storage> >(sz,x) {}

      //! Copy constructor
      Vector2(const Vector2& x) :
        Expr< Vector2Impl<T,Storage> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Vector2(const Expr<S>& x) :
        Expr< Vector2Impl<T,Storage> >(x) {}

      //! Destructor
      ~Vector2() {}

      //! Assignment operator with constant right-hand-side
      Vector2& operator=(const value_type& val) {
        Vector2Impl<T,Storage>::operator=(val);
        return *this;
      }

      //! Assignment operator with constant right-hand-side
      Vector2& operator=(const typename dummy<value_type,scalar_type>::type& val) {
        Vector2Impl<T,Storage>::operator=(value_type(val));
        return *this;
      }

      //! Assignment with Vector2 right-hand-side
      Vector2& operator=(const Vector2& x) {
        Vector2Impl<T,Storage>::operator=(static_cast<const Vector2Impl<T,Storage>&>(x));
        return *this;
      }

      //! Assignment with Expr< Vector2Impl > right-hand-side
      Vector2& operator=(const Expr< Vector2Impl<T,Storage> >& x) {
        Vector2Impl<T,Storage>::operator=(static_cast<const Vector2Impl<T,Storage>&>(x));
        return *this;
      }

      //! Assignment with any expression right-hand-side
      template <typename S>
      Vector2& operator=(const Expr<S>& x) {
        Vector2Impl<T,Storage>::operator=(x);
        return *this;
      }

      //@{

      //! Addition-assignment operator with constant right-hand-side
      Vector2& operator += (const value_type& x) {
        Vector2Impl<T,Storage>::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      Vector2& operator -= (const value_type& x) {
        Vector2Impl<T,Storage>::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      Vector2& operator *= (const value_type& x) {
        Vector2Impl<T,Storage>::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      Vector2& operator /= (const value_type& x) {
        Vector2Impl<T,Storage>::operator/=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      Vector2& operator /= (const typename dummy<value_type,scalar_type>::type& x) {
        Vector2Impl<T,Storage>::operator/=(value_type(x));
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      Vector2& operator += (const typename dummy<value_type,scalar_type>::type& x) {
        Vector2Impl<T,Storage>::operator+=(value_type(x));
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      Vector2& operator -= (const typename dummy<value_type,scalar_type>::type& x) {
        Vector2Impl<T,Storage>::operator-=(value_type(x));
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      Vector2& operator *= (const typename dummy<value_type,scalar_type>::type& x) {
        Vector2Impl<T,Storage>::operator*=(value_type(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      Vector2& operator += (const Expr<S>& x) {
        *this = *this + x;
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      Vector2& operator -= (const Expr<S>& x) {
        *this = *this - x;
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      Vector2& operator *= (const Expr<S>& x) {
        *this = *this * x;
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      Vector2& operator /= (const Expr<S>& x) {
        *this = *this / x;
        return *this;
      }

      //@}

    }; // class Vector2

  } // namespace ETV

} // namespace Sacado

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_ETV_VECTOR_HPP
