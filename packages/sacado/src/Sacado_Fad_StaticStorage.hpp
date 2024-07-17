// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_STATICSTORAGE_HPP
#define SACADO_FAD_STATICSTORAGE_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_StaticStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T, int N>
    using StaticStorage = Exp::StaticStorage<T,N>;

  }
}

#else

#include "Sacado_ConfigDefs.h"
#include "Sacado_StaticArrayTraits.hpp"

namespace Sacado {

  namespace Fad {

    //! Derivative array storage class using static memory allocation
    /*!
     * This class uses a statically allocated array whose dimension is fixed
     * by the template parameter \c Num.
     */
    template <typename T, int Num>
    class StaticStorage {

    public:

       typedef T value_type;

      //! Default constructor
      template <typename S>
      SACADO_INLINE_FUNCTION
      StaticStorage(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        val_(x), sz_(0) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      StaticStorage(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        val_(x), sz_(sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz > Num)
          throw "StaticStorage::StaticStorage() Error:  Supplied derivative dimension exceeds maximum length.";
#endif
        if (zero_out == InitDerivArray)
          ss_array<T>::zero(dx_, sz_);
      }

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      StaticStorage(const StaticStorage& x) :
        val_(x.val_), sz_(x.sz_) {
        //ss_array<T>::copy(x.dx_, dx_, sz_);
        for (int i=0; i<sz_; i++)
          dx_[i] = x.dx_[i];
      }

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~StaticStorage() {}

      //! Assignment
      SACADO_INLINE_FUNCTION
      StaticStorage& operator=(const StaticStorage& x) {
        if (this != &x) {
          val_ = x.val_;
          sz_ = x.sz_;
          //ss_array<T>::copy(x.dx_, dx_, sz_);
          for (int i=0; i<sz_; i++)
            dx_[i] = x.dx_[i];
        }
        return *this;
      }

      //! Returns number of derivative components
      SACADO_INLINE_FUNCTION
      int size() const { return sz_;}

      //! Returns array length
      SACADO_INLINE_FUNCTION
      int length() const { return Num; }

      //! Resize the derivative array to sz
      SACADO_INLINE_FUNCTION
      void resize(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz > Num)
          throw "StaticStorage::resize() Error:  Supplied derivative dimension exceeds maximum length.";
#endif
        sz_ = sz;
      }

      //! Resize the derivative array to sz
      /*!
       * This method doest not preserve any existing derivative components but
       * sets any that are added to zero.
       */
      SACADO_INLINE_FUNCTION
      void resizeAndZero(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz > Num)
          throw "StaticStorage::resize() Error:  Supplied derivative dimension exceeds maximum length.";
#endif
        if (sz > sz_)
          ss_array<T>::zero(dx_+sz_, sz-sz_);
        sz_ = sz;
      }

      //! Expand derivative array to size sz
      /*!
       * This method preserves any existing derivative components and
       * sets any that are added to zero.
       */
      SACADO_INLINE_FUNCTION
      void expand(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz > Num)
          throw "StaticStorage::resize() Error:  Supplied derivative dimension exceeds maximum length.";
#endif
        if (sz > sz_)
          ss_array<T>::zero(dx_+sz_, sz-sz_);
        sz_ = sz;
      }


      //! Zero out derivative array
      SACADO_INLINE_FUNCTION
      void zero() { ss_array<T>::zero(dx_, sz_); }

      //! Returns value
      SACADO_INLINE_FUNCTION
      const T& val() const { return val_; }

      //! Returns value
      SACADO_INLINE_FUNCTION
      T& val() { return val_; }

      //! Returns derivative array
      SACADO_INLINE_FUNCTION
      const T* dx() const { return dx_;}

      //! Returns derivative component \c i with bounds checking
      SACADO_INLINE_FUNCTION
      T dx(int i) const { return sz_ ? dx_[i] : T(0.); }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      T& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      const T& fastAccessDx(int i) const { return dx_[i];}

    protected:

      //! Value
      T val_;

      //! Derivative array
      T dx_[Num];

      //! Size of derivative array
      int sz_;

    }; // class StaticStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_STATICSTORAGE_HPP
