// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_DYNAMICSTORAGE_HPP
#define SACADO_FAD_DYNAMICSTORAGE_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_DynamicStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T, typename U = T>
    using DynamicStorage = Exp::DynamicStorage<T,U>;

  }
}

#else

#include "Sacado_Traits.hpp"
#include "Sacado_DynamicArrayTraits.hpp"

namespace Sacado {

  namespace Fad {

    //! Derivative array storage class using dynamic memory allocation
    template <typename T, typename U = T>
    class DynamicStorage {

    public:

      typedef T value_type;

      //! Default constructor
      template <typename S>
      SACADO_INLINE_FUNCTION
      DynamicStorage(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        val_(x), sz_(0), len_(0), dx_(NULL) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      DynamicStorage(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        val_(x), sz_(sz), len_(sz) {
        if (zero_out == InitDerivArray)
          dx_ = ds_array<U>::get_and_fill(sz_);
        else
          dx_ = ds_array<U>::get(sz_);
      }

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      DynamicStorage(const DynamicStorage& x) :
        val_(x.val_), sz_(x.sz_), len_(x.sz_) {
        dx_ = ds_array<U>::get_and_fill(x.dx_, sz_);
      }

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~DynamicStorage() {
        if (len_ != 0)
          ds_array<U>::destroy_and_release(dx_, len_);
      }

      //! Assignment
      SACADO_INLINE_FUNCTION
      DynamicStorage& operator=(const DynamicStorage& x) {
        if (this != &x) {
          val_ = x.val_;
          if (sz_ != x.sz_) {
            sz_ = x.sz_;
            if (x.sz_ > len_) {
              if (len_ != 0)
                ds_array<U>::destroy_and_release(dx_, len_);
              len_ = x.sz_;
              dx_ = ds_array<U>::get_and_fill(x.dx_, sz_);
            }
            else
              ds_array<U>::copy(x.dx_, dx_, sz_);
          }
          else
            ds_array<U>::copy(x.dx_, dx_, sz_);
        }
        return *this;
      }

      //! Returns number of derivative components
      SACADO_INLINE_FUNCTION
      int size() const { return sz_;}

      //! Returns array length
      SACADO_INLINE_FUNCTION
      int length() const { return len_; }

      //! Resize the derivative array to sz
      /*!
       * Note:  This does not necessarily preserve derivative components.
       */
      SACADO_INLINE_FUNCTION
      void resize(int sz) {
        if (sz > len_) {
          if (len_ != 0)
            ds_array<U>::destroy_and_release(dx_, len_);
          dx_ = ds_array<U>::get_and_fill(sz);
          len_ = sz;
        }
        sz_ = sz;
      }

      //! Resize the derivative array to sz
      /*!
       * This method doest not preserve any existing derivative components but
       * sets any that are added to zero.
       */
      SACADO_INLINE_FUNCTION
      void resizeAndZero(int sz) {
        if (sz > len_) {
          if (len_ != 0)
            ds_array<U>::destroy_and_release(dx_, len_);
          dx_ = ds_array<U>::get_and_fill(sz);
          len_ = sz;
        }
        else if (sz > sz_)
          ds_array<U>::zero(dx_+sz_, sz-sz_);
        sz_ = sz;
      }

      //! Expand derivative array to size sz
      /*!
       * This method preserves any existing derivative components and
       * sets any that are added to zero.
       */
      SACADO_INLINE_FUNCTION
      void expand(int sz) {
        if (sz > len_) {
          U* dx_new = ds_array<U>::get_and_fill(sz);
          ds_array<U>::copy(dx_, dx_new, sz_);
          if (len_ > 0)
            ds_array<U>::destroy_and_release(dx_, len_);
          dx_ = dx_new;
          len_ = sz;
        }
        else if (sz > sz_)
          ds_array<U>::zero(dx_+sz_, sz-sz_);
        sz_ = sz;
      }

      //! Zero out derivative array
      SACADO_INLINE_FUNCTION
      void zero() {
        ds_array<U>::zero(dx_, sz_);
      }

      //! Returns value
      SACADO_INLINE_FUNCTION
      const T& val() const { return val_; }

      //! Returns value
      SACADO_INLINE_FUNCTION
      T& val() { return val_; }

      //! Returns derivative array
      SACADO_INLINE_FUNCTION
      const U* dx() const { return dx_;}

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(__CUDA_ARCH__)

      //! Returns derivative component \c i with bounds checking
      SACADO_INLINE_FUNCTION
      U dx(int i) const { return sz_ ? dx_[i*blockDim.x] : U(0.); }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      U& fastAccessDx(int i) { return dx_[i*blockDim.x];}

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      const U& fastAccessDx(int i) const { return dx_[i*blockDim.x];}

#else

      //! Returns derivative component \c i with bounds checking
      SACADO_INLINE_FUNCTION
      U dx(int i) const { return sz_ ? dx_[i] : U(0.); }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      U& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      const U& fastAccessDx(int i) const { return dx_[i];}

#endif

    protected:

      //! Value
      T val_;

      //! Derivative array size
      int sz_;

      //! Derivative array length
      int len_;

      //! Derivative array
      U* dx_;

    }; // class DynamicStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_DYNAMICSTORAGE_HPP
