// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_DYNAMICSTORAGE_HPP
#define SACADO_FAD_EXP_DYNAMICSTORAGE_HPP

#include <type_traits>
#include <utility>

#include "Sacado_ConfigDefs.h"
#include "Sacado_DynamicArrayTraits.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    //! Derivative array storage class using dynamic memory allocation
    template <typename T, typename U = T>
    class DynamicStorage {

    public:

      typedef typename std::remove_cv<T>::type value_type;
      static constexpr bool is_statically_sized = false;
      static constexpr int static_size = 0;
      static constexpr bool is_view = false;

      //! Turn DynamicStorage into a meta-function class usable with mpl::apply
      template <typename TT, typename UU = TT>
      struct apply {
        typedef DynamicStorage<TT,UU> type;
      };

      //! Replace static derivative length
      template <int N>
      struct apply_N {
        typedef DynamicStorage<T,U> type;
      };

      //! Default constructor
      SACADO_INLINE_FUNCTION
      DynamicStorage() :
        val_(), sz_(0), len_(0), dx_(nullptr) {}

      //! Constructor with value
      SACADO_INLINE_FUNCTION
      DynamicStorage(const T & x) :
        val_(x), sz_(0), len_(0), dx_(nullptr) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      DynamicStorage(const int sz, const T & x,
                     const DerivInit zero_out = InitDerivArray) :
        val_(x), sz_(sz), len_(sz) {
        if (zero_out == InitDerivArray)
          dx_ = ds_array<U>::get_and_fill(sz_);
        else
          dx_ = ds_array<U>::get(sz_);
      }

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SACADO_INLINE_FUNCTION
      DynamicStorage(const int sz, const int i, const value_type & x) :
        DynamicStorage(sz, x, InitDerivArray) {
        dx_[i]=1.;
      }

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      DynamicStorage(const DynamicStorage& x) :
        val_(x.val_), sz_(x.sz_), len_(x.sz_) {
        dx_ = ds_array<U>::get_and_fill(x.dx_, sz_);
      }

      //! Move constructor
      SACADO_INLINE_FUNCTION
      DynamicStorage(DynamicStorage&& x) :
        val_(std::move(x.val_)), sz_(x.sz_), len_(x.len_), dx_(x.dx_) {
        x.sz_ = 0;
        x.len_ = 0;
        x.dx_ = nullptr;
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

      //! Move assignment
      SACADO_INLINE_FUNCTION
      DynamicStorage& operator=(DynamicStorage&& x) {
        if (this != &x) {
          if (len_ != 0)
            ds_array<U>::destroy_and_release(dx_, len_);
          val_ = std::move(x.val_);
          sz_ = x.sz_; x.sz_ = 0;
          len_ = x.len_; x.len_ = 0;
          dx_ = x.dx_; x.dx_ = nullptr;
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

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )

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

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXP_DYNAMICSTORAGE_HPP
