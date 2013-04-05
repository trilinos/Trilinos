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

#ifndef SACADO_ETV_VECTOR_HPP
#define SACADO_ETV_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include "Sacado_Traits.hpp"
#include "Sacado_Handle.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_dummy_arg.hpp"

#include <ostream>      // for std::ostream

#ifdef HAVE_STOKHOS_THRUST
#include "thrust/tuple.h"
#endif

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

namespace Sacado {

  //! Namespace for expression templated vector classes
  namespace ETV {

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all expression
     * template classes.
     */
    template <typename ExprT> class Expr {};

    //! Vectorized evaluation class
    template <typename T, typename Storage >
    class VectorImpl {
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
      VectorImpl();

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      VectorImpl(const value_type& x);

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      VectorImpl(ordinal_type sz, const value_type& x);

      //! Copy constructor
      VectorImpl(const VectorImpl& x);

      //! Copy constructor from any Expression object
      template <typename S> VectorImpl(const Expr<S>& x);

      //! Destructor
      ~VectorImpl() {}

      //! Initialize coefficients to value
      void init(const T& v) { th_->init(v); }

      //! Initialize coefficients to an array of values
      void init(const T* v) { th_->init(v); }

      //! Initialize coefficients from an VectorImpl with different storage
      template <typename S>
      void init(const VectorImpl<T,S>& v) {
        th_->init(v.th_->coeff(), v.th_->size());
      }

      //! Load coefficients to an array of values
      void load(T* v) { th_->load(v); }

      //! Load coefficients into an VectorImpl with different storage
      template <typename S>
      void load(VectorImpl<T,S>& v) { th_->load(v.th_->coeff()); }

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
      void copyForWrite() { th_.makeOwnCopy(); }

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
      VectorImpl& operator=(const value_type& val);

      //! Assignment with VectorImpl right-hand-side
      VectorImpl& operator=(const VectorImpl& x);

      //! Assignment with any expression right-hand-side
      template <typename S>
      VectorImpl& operator=(const Expr<S>& x);

      //@}

      /*!
       * Accessor methods
       */

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const_reference val() const { return (*th_)[0]; }

      //! Returns value
      reference val() { return (*th_)[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      ordinal_type size() const { return th_->size();}

      //! Returns true if polynomial has size >= sz
      bool hasFastAccess(ordinal_type sz) const { return th_->size()>=sz;}

      //! Returns Hermite coefficient array
      const_pointer coeff() const { return th_->coeff();}

      //! Returns Hermite coefficient array
      pointer coeff() { return th_->coeff();}

      //! Returns degree \c i term with bounds checking
      value_type coeff(ordinal_type i) const {
        return i<th_->size() ? (*th_)[i] : (*th_)[0]; }

      //! Returns degree \c i term without bounds checking
      reference fastAccessCoeff(ordinal_type i) { return (*th_)[i];}

      //! Returns degree \c i term without bounds checking
      value_type fastAccessCoeff(ordinal_type i) const { return (*th_)[i];}

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      VectorImpl& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      VectorImpl& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      VectorImpl& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      VectorImpl& operator /= (const value_type& x);

      //@}

    protected:

      //! Handle to underlying vector
      Sacado::Handle<storage_type> th_;

    }; // class VectorImpl

    //! VectorImpl expression template specialization
    /*!
     * This template class represents a simple VectorImpl expression and
     * mixes-in the VectorImpl interface and the expression template
     * interface.
     */
    template <typename T, typename Storage>
    class Expr< VectorImpl<T,Storage> > :
public VectorImpl<T,Storage> {

    public:

      //! Typename of values
      typedef typename VectorImpl<T,Storage>::value_type value_type;
      typedef typename VectorImpl<T,Storage>::scalar_type scalar_type;
      typedef typename VectorImpl<T,Storage>::storage_type storage_type;
      typedef typename VectorImpl<T,Storage>::const_reference const_reference;

      //! Number of arguments
      static const int num_args = 1;

      //! Default constructor
      Expr() :
        VectorImpl<T,Storage>() {}

      //! Constructor with supplied value \c x
      Expr(const T & x) :
        VectorImpl<T,Storage>(x) {}

      //! Constructor with specified size \c sz
      Expr(typename VectorImpl<T,Storage>::ordinal_type sz, const T& x) :
        VectorImpl<T,Storage>(sz,x) {}

      //! Copy constructor
      Expr(const Expr& x) :
        VectorImpl<T,Storage>(static_cast<const VectorImpl<T,Storage>&>(x)) {}

      //! Copy constructor
      Expr(const VectorImpl<T,Storage>& x) :
        VectorImpl<T,Storage>(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Expr(const Expr<S>& x) :
        VectorImpl<T,Storage>(x) {}

      //! Destructor
      ~Expr() {}

      std::string name() const { return "x"; }

    }; // class Expr<VectorImpl>

  } // namespace ETV

} // namespace Sacado

#include "Sacado_ETV_ExpressionTraits.hpp"
#include "Sacado_ETV_VectorTraits.hpp"
#include "Sacado_ETV_VectorImp.hpp"
#include "Sacado_ETV_VectorOps.hpp"

namespace Sacado {

  namespace ETV {

    //! Generalized polynomial chaos expansion class
    template <typename T, typename Storage>
    class Vector : public Expr< VectorImpl<T,Storage> > {
    public:

      //! Typename of values
      typedef typename VectorImpl<T,Storage>::value_type value_type;

      //! Typename of scalars
      typedef typename VectorImpl<T,Storage>::scalar_type scalar_type;

      //! Typename of ordinals
      typedef typename VectorImpl<T,Storage>::ordinal_type ordinal_type;

      //! Typename of storage class
      typedef typename VectorImpl<T,Storage>::storage_type storage_type;

      typedef typename VectorImpl<T,Storage>::pointer pointer;
      typedef typename VectorImpl<T,Storage>::const_pointer const_pointer;
      typedef typename VectorImpl<T,Storage>::reference reference;
      typedef typename VectorImpl<T,Storage>::const_reference const_reference;

      //! Turn Vector into a meta-function class usable with mpl::apply
      template <typename S>
      struct apply {
        typedef typename Sacado::mpl::apply<Storage,ordinal_type,S>::type storage_type;
        typedef Vector<S,storage_type> type;
      };

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      Vector() :
        Expr< VectorImpl<T,Storage> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      Vector(const value_type& x) :
        Expr< VectorImpl<T,Storage> >(x) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      Vector(const typename dummy<value_type,scalar_type>::type& x) :
        Expr< VectorImpl<T,Storage> >(value_type(x)) {}

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      Vector(ordinal_type sz, const value_type& x) :
        Expr< VectorImpl<T,Storage> >(sz,x) {}

      //! Copy constructor
      Vector(const Vector& x) :
        Expr< VectorImpl<T,Storage> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Vector(const Expr<S>& x) :
        Expr< VectorImpl<T,Storage> >(x) {}

      //! Destructor
      ~Vector() {}

      //! Assignment operator with constant right-hand-side
      Vector& operator=(const value_type& val) {
        VectorImpl<T,Storage>::operator=(val);
        return *this;
      }

      //! Assignment operator with constant right-hand-side
      Vector& operator=(const typename dummy<value_type,scalar_type>::type& val) {
        VectorImpl<T,Storage>::operator=(value_type(val));
        return *this;
      }

      //! Assignment with Vector right-hand-side
      Vector& operator=(const Vector& x) {
        VectorImpl<T,Storage>::operator=(static_cast<const VectorImpl<T,Storage>&>(x));
        return *this;
      }

      //! Assignment with Expr< VectorImpl > right-hand-side
      Vector& operator=(const Expr< VectorImpl<T,Storage> >& x) {
        VectorImpl<T,Storage>::operator=(static_cast<const VectorImpl<T,Storage>&>(x));
        return *this;
      }

      //! Assignment with any expression right-hand-side
      template <typename S>
      Vector& operator=(const Expr<S>& x) {
        VectorImpl<T,Storage>::operator=(x);
        return *this;
      }

      //@{

      //! Addition-assignment operator with constant right-hand-side
      Vector& operator += (const value_type& x) {
        VectorImpl<T,Storage>::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      Vector& operator -= (const value_type& x) {
        VectorImpl<T,Storage>::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      Vector& operator *= (const value_type& x) {
        VectorImpl<T,Storage>::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      Vector& operator /= (const value_type& x) {
        VectorImpl<T,Storage>::operator/=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      Vector& operator /= (const typename dummy<value_type,scalar_type>::type& x) {
        VectorImpl<T,Storage>::operator/=(value_type(x));
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      Vector& operator += (const typename dummy<value_type,scalar_type>::type& x) {
        VectorImpl<T,Storage>::operator+=(value_type(x));
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      Vector& operator -= (const typename dummy<value_type,scalar_type>::type& x) {
        VectorImpl<T,Storage>::operator-=(value_type(x));
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      Vector& operator *= (const typename dummy<value_type,scalar_type>::type& x) {
        VectorImpl<T,Storage>::operator*=(value_type(x));
        return *this;
      }

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

    }; // class Vector

  } // namespace ETV

} // namespace Sacado

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_ETV_VECTOR_HPP
