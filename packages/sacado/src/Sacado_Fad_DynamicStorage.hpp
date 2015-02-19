// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_DYNAMICSTORAGE_HPP
#define SACADO_FAD_DYNAMICSTORAGE_HPP

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
      KOKKOS_INLINE_FUNCTION
      DynamicStorage(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        val_(x), sz_(0), len_(0), dx_(NULL) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      DynamicStorage(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        val_(x), sz_(sz), len_(sz) {
        if (zero_out == InitDerivArray)
          dx_ = ds_array<U>::get_and_fill(sz_);
        else
          dx_ = ds_array<U>::get(sz_);
      }

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      DynamicStorage(const DynamicStorage& x) :
        val_(x.val_), sz_(x.sz_), len_(x.sz_) {
        dx_ = ds_array<U>::get_and_fill(x.dx_, sz_);
      }

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~DynamicStorage() {
        if (len_ != 0)
          ds_array<U>::destroy_and_release(dx_, len_);
      }

      //! Assignment
      KOKKOS_INLINE_FUNCTION
      DynamicStorage& operator=(const DynamicStorage& x) {
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

        return *this;
      }

      //! Returns number of derivative components
      KOKKOS_INLINE_FUNCTION
      int size() const { return sz_;}

      //! Returns array length
      KOKKOS_INLINE_FUNCTION
      int length() const { return len_; }

      //! Resize the derivative array to sz
      /*!
       * Note:  This does not necessarily preserve derivative components.
       */
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
      void zero() {
        ds_array<U>::zero(dx_, sz_);
      }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const T& val() const { return val_; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      T& val() { return val_; }

      //! Returns derivative array
      KOKKOS_INLINE_FUNCTION
      const U* dx() const { return dx_;}

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      U dx(int i) const { return sz_ ? dx_[i] : U(0.); }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      U& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const U& fastAccessDx(int i) const { return dx_[i];}

    private:

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

#endif // SACADO_FAD_DYNAMICSTORAGE_HPP
