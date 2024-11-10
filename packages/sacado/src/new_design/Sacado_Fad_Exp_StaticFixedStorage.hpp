// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_STATICFIXEDSTORAGE_HPP
#define SACADO_FAD_EXP_STATICFIXEDSTORAGE_HPP

#include <type_traits>
#include <utility>

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
      static constexpr bool is_view = false;

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
#ifdef SACADO_SFAD_INIT_DEFAULT_CONSTRUCTOR
      SACADO_INLINE_FUNCTION
      StaticFixedStorage() :
        val_(T(0.0)) {
        ss_array<T>::zero(dx_, Num);
      }
#else
      SACADO_DEFAULTED_FUNCTION
      StaticFixedStorage() = default;
#endif

      //! Constructor with value
      SACADO_INLINE_FUNCTION
      StaticFixedStorage(const T & x) :
        val_(x) {
        ss_array<T>::zero(dx_, Num);
      }

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      StaticFixedStorage(const int sz, const T & x,
                         const DerivInit zero_out = InitDerivArray) :
        val_(x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ ) && !defined(__HIP_DEVICE_COMPILE__)
        if (sz != Num)
          throw "StaticFixedStorage::StaticFixedStorage() Error:  Supplied derivative dimension does not equal static length.";
#endif
        if (zero_out == InitDerivArray)
          ss_array<T>::zero(dx_, Num);
      }

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SACADO_INLINE_FUNCTION
      StaticFixedStorage(const int sz, const int i, const value_type & x) :
        StaticFixedStorage(sz, x, InitDerivArray) {
        dx_[i]=1.;
      }

      //! Copy constructor
      /*!
       * Can't make this " = default" because of scalar types that don't
       * define a const copy consturctor (like Rad).  We also can't leave it
       * and let it be implicitly generated because of SACADO_INLINE_FUNCTION.
       */
      SACADO_INLINE_FUNCTION
      StaticFixedStorage(const StaticFixedStorage& x) :
        val_(x.val_) {
         for (int i=0; i<Num; i++)
           dx_[i] = x.dx_[i];
      }

      //! Move constructor
      SACADO_INLINE_FUNCTION
      StaticFixedStorage(StaticFixedStorage&& x) :
        val_(std::move(x.val_)) {
         for (int i=0; i<Num; i++)
           dx_[i] = std::move(x.dx_[i]);
      }

      //! Destructor
      SACADO_DEFAULTED_FUNCTION
      ~StaticFixedStorage() = default;

      //! Assignment
      /*!
       * Can't make this " = default" because of scalar types that don't
       * define a const operator= (like Rad).  We also can't leave it
       * and let it be implicitly generated because of SACADO_INLINE_FUNCTION.
       */
      SACADO_INLINE_FUNCTION
      StaticFixedStorage& operator=(const StaticFixedStorage& x) {
        if (this != &x) {
          val_ = x.val_;
          for (int i=0; i<Num; i++)
            dx_[i] = x.dx_[i];
        }
        return *this;
      }

      //! Move assignment
      SACADO_INLINE_FUNCTION
      StaticFixedStorage& operator=(StaticFixedStorage&& x) {
        if (this != &x) {
          val_ = std::move(x.val_);
          for (int i=0; i<Num; i++)
            dx_[i] = std::move(x.dx_[i]);
        }
        return *this;
      }

      //! Returns number of derivative components
      SACADO_INLINE_FUNCTION
      constexpr int size() const { return Num; }

      //! Returns array length
      SACADO_INLINE_FUNCTION
      constexpr int length() const { return Num; }

      //! Resize the derivative array to sz
      SACADO_INLINE_FUNCTION
      void resize(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ ) && !defined(__HIP_DEVICE_COMPILE__)
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
      SACADO_INLINE_FUNCTION
      void resizeAndZero(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ ) && !defined(__HIP_DEVICE_COMPILE__)
        if (sz != 0 && sz != Num)
          throw "StaticFixedStorage::resize() Error:  Cannot resize fixed storage length.";
#endif
        ss_array<T>::zero(dx_, Num);
      }

      //! Expand derivative array to size sz
      SACADO_INLINE_FUNCTION
      void expand(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ ) && !defined(__HIP_DEVICE_COMPILE__)
        if (sz != Num)
          throw "StaticFixedStorage::expand() Error:  Cannot resize fixed storage length.";
#endif
      }


      //! Zero out derivative array
      SACADO_INLINE_FUNCTION
      void zero() { ss_array<T>::zero(dx_, Num); }

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
      T dx(int i) const { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      T& fastAccessDx(int i) { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
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
