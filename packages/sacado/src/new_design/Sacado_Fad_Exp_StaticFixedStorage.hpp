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

#ifndef SACADO_FAD_EXP_STATICFIXEDSTORAGE_HPP
#define SACADO_FAD_EXP_STATICFIXEDSTORAGE_HPP

#include <type_traits>

#include "Sacado_ConfigDefs.h"
#include "Sacado_StaticArrayTraits.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    //! Derivative array storage class using static, fixed memory allocation
    /*!
     * This class uses a statically allocated array whose dimension is fixed
     * by the template parameter \c Num.  The dimension cannot be resized.
     */
    template <typename T, int Num>
    class StaticFixedStorage {

    public:

      typedef typename std::remove_cv<T>::type value_type;
      static constexpr bool is_statically_sized = true;
      static constexpr int static_size = Num;

      //! Turn StaticFixedStorage into a meta-function class usable with mpl::apply
      template <typename TT>
      struct apply {
        typedef StaticFixedStorage<TT,Num> type;
      };

      //! Replace static derivative length
      template <int N>
      struct apply_N {
        typedef StaticFixedStorage<T,N> type;
      };

      //! Default constructor
      StaticFixedStorage() = default;

      //! Constructor with value
      KOKKOS_INLINE_FUNCTION
      StaticFixedStorage(const T & x) :
        val_(x) {
        ss_array<T>::zero(dx_, Num);
      }

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      StaticFixedStorage(const int sz, const T & x, const DerivInit zero_out) :
        val_(x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "StaticFixedStorage::StaticFixedStorage() Error:  Supplied derivative dimension does not equal static length.";
#endif
        if (zero_out == InitDerivArray)
          ss_array<T>::zero(dx_, Num);
      }

      //! Copy constructor
      /*!
       * Can't make this " = default" because of scalar types that don't
       * define a const copy consturctor (like Rad).  We also can't leave it
       * and let it be implicitly generated because of KOKKOS_INLINE_FUNCTION.
       */
      KOKKOS_INLINE_FUNCTION
      StaticFixedStorage(const StaticFixedStorage& x) :
        val_(x.val_) {
         for (int i=0; i<Num; i++)
           dx_[i] = x.dx_[i];
      }

      //! Destructor
      ~StaticFixedStorage() = default;

      //! Assignment
      /*!
       * Can't make this " = default" because of scalar types that don't
       * define a const operator= (like Rad).  We also can't leave it
       * and let it be implicitly generated because of KOKKOS_INLINE_FUNCTION.
       */
      KOKKOS_INLINE_FUNCTION
      StaticFixedStorage& operator=(const StaticFixedStorage& x) {
        if (this != &x) {
          val_ = x.val_;
          for (int i=0; i<Num; i++)
            dx_[i] = x.dx_[i];
        }
        return *this;
      }

      //! Returns number of derivative components
      KOKKOS_INLINE_FUNCTION
      constexpr int size() const { return Num; }

      //! Returns array length
      KOKKOS_INLINE_FUNCTION
      constexpr int length() const { return Num; }

      //! Resize the derivative array to sz
      KOKKOS_INLINE_FUNCTION
      void resize(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != 0 && sz != Num)
          throw "StaticFixedStorage::resize() Error:  Cannot resize fixed storage length.";
#endif
        // Because we don't track a "used" length and can't set the length to 0,
        // we need to instead zero out derivative components if the length
        // requested is 0
        if (sz == 0)
          ss_array<T>::zero(dx_, Num);
      }

      //! Resize the derivative array to sz
      KOKKOS_INLINE_FUNCTION
      void resizeAndZero(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != 0 && sz != Num)
          throw "StaticFixedStorage::resize() Error:  Cannot resize fixed storage length.";
#endif
        ss_array<T>::zero(dx_, Num);
      }

      //! Expand derivative array to size sz
      KOKKOS_INLINE_FUNCTION
      void expand(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "StaticFixedStorage::expand() Error:  Cannot resize fixed storage length.";
#endif
      }


      //! Zero out derivative array
      KOKKOS_INLINE_FUNCTION
      void zero() { ss_array<T>::zero(dx_, Num); }

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
      T dx(int i) const { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      T& fastAccessDx(int i) { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& fastAccessDx(int i) const { return dx_[i]; }

    protected:

      //! Value
      T val_;

      //! Derivative array
      T dx_[Num];

    }; // class StaticFixedStorage

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXP_STATICFIXEDSTORAGE_HPP
