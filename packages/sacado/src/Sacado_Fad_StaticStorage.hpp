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

#ifndef SACADO_FAD_STATICSTORAGE_HPP
#define SACADO_FAD_STATICSTORAGE_HPP

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
      KOKKOS_INLINE_FUNCTION
      StaticStorage(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        val_(x), sz_(0) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      StaticStorage(const int sz, const T & x, const bool zero_out = true) :
        val_(x), sz_(sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz > Num)
          throw "StaticStorage::StaticStorage() Error:  Supplied derivative dimension exceeds maximum length.";
#endif
        if (zero_out)
          ss_array<T>::zero(dx_, sz_);
      }

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      StaticStorage(const StaticStorage& x) :
        val_(x.val_), sz_(x.sz_) {
        //ss_array<T>::copy(x.dx_, dx_, sz_);
        for (int i=0; i<sz_; i++)
          dx_[i] = x.dx_[i];
      }

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~StaticStorage() {}

      //! Assignment
      KOKKOS_INLINE_FUNCTION
      StaticStorage& operator=(const StaticStorage& x) {
        val_ = x.val_;
        sz_ = x.sz_;
        //ss_array<T>::copy(x.dx_, dx_, sz_);
        for (int i=0; i<sz_; i++)
          dx_[i] = x.dx_[i];
        return *this;
      }

      //! Returns number of derivative components
      KOKKOS_INLINE_FUNCTION
      int size() const { return sz_;}

      //! Returns array length
      KOKKOS_INLINE_FUNCTION
      int length() const { return Num; }

      //! Resize the derivative array to sz
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
      void zero() { ss_array<T>::zero(dx_, sz_); }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const T& val() const { return val_; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      T& val() { return val_; }

      //! Returns derivative array
      KOKKOS_INLINE_FUNCTION
      const T* dx() const { return dx_;}

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      T dx(int i) const { return sz_ ? dx_[i] : T(0.); }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      T& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& fastAccessDx(int i) const { return dx_[i];}

    private:

      //! Value
      T val_;

      //! Derivative array
      T dx_[Num];

      //! Size of derivative array
      int sz_;

    }; // class StaticStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_STATICSTORAGE_HPP
