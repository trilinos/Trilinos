// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_SIMPLEFAD_HPP
#define SACADO_FAD_SIMPLEFAD_HPP

#include "Sacado_Fad_SimpleFadTraits.hpp"
#include "Sacado_Fad_GeneralFad.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"

namespace Sacado {

  namespace Fad {

    /*!
     * \brief Forward-mode AD class using dynamic memory allocation but no
     * expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::Fad::GeneralFad.
     */
    template <typename ValueT>
    class SimpleFad : public GeneralFad<ValueT,DynamicStorage<ValueT> > {

    public:

      //! Base classes
      typedef DynamicStorage<ValueT> StorageType;
      typedef GeneralFad<ValueT,StorageType> GeneralFadType;

      //! Typename of values
      typedef typename GeneralFadType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename GeneralFadType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SimpleFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef SimpleFad<T> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      SimpleFad() :
        GeneralFadType() {}

      //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      SimpleFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        GeneralFadType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SimpleFad(const int sz, const ValueT& x, const DerivInit zero_out = InitDerivArray) :
        GeneralFadType(sz,x,zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SimpleFad(const int sz, const int i, const ValueT & x) :
        GeneralFadType(sz,i,x) {}

      //! Copy constructor
      SimpleFad(const SimpleFad& x) :
        GeneralFadType(x) {}

      //! Tangent copy constructor
      SimpleFad(const SimpleFad& x, const ValueT& v, const ValueT& partial) :
        GeneralFadType(x.size(), v) {
        for (int i=0; i<this->size(); i++)
          this->fastAccessDx(i) = x.fastAccessDx(i)*partial;
      }

      //@}

      //! Destructor
      ~SimpleFad() {}

      //! Returns whether two Fad objects have the same values
      bool isEqualTo(const SimpleFad& x) const {
        typedef IsEqual<value_type> IE;
        if (x.size() != this->size()) return false;
        bool eq = IE::eval(x.val(), this->val());
        for (int i=0; i<this->size(); i++)
          eq = eq && IE::eval(x.dx(i), this->dx(i));
        return eq;
      }

      //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(SimpleFad&) operator=(const S& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with SimpleFad right-hand-side
      SimpleFad& operator=(const SimpleFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SimpleFad&) operator += (const S& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SimpleFad&) operator -= (const S& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SimpleFad&) operator *= (const S& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SimpleFad&) operator /= (const S& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with SimpleFad right-hand-side
      SACADO_INLINE_FUNCTION
      SimpleFad& operator += (const SimpleFad& x) {
        GeneralFadType::operator+=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with SimpleFad right-hand-side
      SACADO_INLINE_FUNCTION
      SimpleFad& operator -= (const SimpleFad& x) {
        GeneralFadType::operator-=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with SimpleFad right-hand-side
      SACADO_INLINE_FUNCTION
      SimpleFad& operator *= (const SimpleFad& x) {
        GeneralFadType::operator*=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Division-assignment operator with SimpleFad right-hand-side
      SACADO_INLINE_FUNCTION
      SimpleFad& operator /= (const SimpleFad& x) {
        GeneralFadType::operator/=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

    }; // class SimpleFad<ValueT>

  } // namespace Fad

} // namespace Sacado

// Include elementary operation overloads
#include "Sacado_Fad_SimpleFadOps.hpp"

#endif // SACADO_FAD_SIMPLEFAD_HPP
