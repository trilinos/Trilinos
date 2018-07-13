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

#ifndef SACADO_FAD_EXP_VIEWSTORAGE_HPP
#define SACADO_FAD_EXP_VIEWSTORAGE_HPP

#include <type_traits>
#include <memory>

#include "Sacado_DynamicArrayTraits.hpp"
#include "Sacado_mpl_integral_nonzero_constant.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    // Class representing a pointer to ViewFad so that &ViewFad is supported
    template <typename T, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr;

    /*!
     * \brief Derivative array storage class that is a view into a contiguous
     * memory allocation.  It does not provide proper value semantics and
     * thus should not be used in a general-purpose scalar type.
     */
    template <typename T, unsigned static_length, unsigned static_stride, typename U>
    class ViewStorage {

    private:

      // Enumerated flag so logic is evaluated at compile-time
      enum { stride_one = 1 == static_stride };

    public:

      typedef typename std::remove_cv<T>::type value_type;
      static constexpr bool is_statically_sized = (static_length > 0);
      static constexpr int static_size = static_length;
      typedef U base_fad_type;

      //! Turn ViewStorage into a meta-function class usable with mpl::apply
      template <typename TT>
      struct apply {
        typedef ViewStorage<TT,static_length,static_stride,U> type;
      };

      //! Replace static derivative length
      template <int N>
      struct apply_N {
        typedef ViewStorage<T,static_length,static_stride,U> type;
      };

      //! Default constructor (needed to satisfy interface)
      KOKKOS_INLINE_FUNCTION
      ViewStorage() :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor with value (needed to satisfy interface)
      KOKKOS_INLINE_FUNCTION
      ViewStorage(const T & x) :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor with size \c sz (needed to satisfy interface)
      KOKKOS_INLINE_FUNCTION
      ViewStorage(const int sz, const T & x, const DerivInit zero_out) :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor
      KOKKOS_INLINE_FUNCTION
      ViewStorage(T* v, const int arg_size = 0, const int arg_stride = 0) :
        sz_(arg_size), stride_(arg_stride), val_(v+sz_.value*stride_.value), dx_(v) {}

      //! Constructor
      KOKKOS_INLINE_FUNCTION
      ViewStorage(T* arg_dx, T* arg_val, const int arg_size = 0,
                  const int arg_stride = 0) :
        sz_(arg_size), stride_(arg_stride), val_(arg_val), dx_(arg_dx) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      ViewStorage(const ViewStorage& x) :
        sz_(x.sz_), stride_(x.stride_), val_(x.val_), dx_(x.dx_) {}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~ViewStorage() {}

      //! Assignment
      KOKKOS_INLINE_FUNCTION
      ViewStorage& operator=(const ViewStorage& x) {
        if (this != std::addressof(x)) {
          *val_ = *x.val_;
          if (stride_one)
            for (int i=0; i<sz_.value; ++i)
              dx_[i] = x.dx_[i];
          else
            for (int i=0; i<sz_.value; ++i)
              dx_[i*stride_.value] = x.dx_[i*x.stride_.value];
        }
        return *this;
      }

      //! Returns number of derivative components
      KOKKOS_INLINE_FUNCTION
      constexpr int size() const { return sz_.value;}

      //! Returns array length
      KOKKOS_INLINE_FUNCTION
      constexpr int length() const { return sz_.value; }

      //! Resize the derivative array to sz
      /*!
       * Since we can't actually resize, we check for resizes to zero,
       * which signify assigning a constant.  Thus we zero out the derivative
       * components.
       */
      KOKKOS_INLINE_FUNCTION
      void resize(int sz) {
        if (sz == 0) ds_array<T>::strided_zero(dx_, stride_.value, sz_.value);
      }

      //! Resize the derivative array to sz
      /*!
       * We don't do anything here as this is used in the context of resizing
       * the derivative array to zero and then back to some size > 0.  Instead
       * we zero out components when it is resized to zero above.
       */
      KOKKOS_INLINE_FUNCTION
      void resizeAndZero(int sz) {}

      //! Expand derivative array to size sz
      KOKKOS_INLINE_FUNCTION
      void expand(int sz) {}

      //! Zero out derivative array
      KOKKOS_INLINE_FUNCTION
      void zero() {
        ds_array<T>::strided_zero(dx_, stride_.value, sz_.value);
      }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const T& val() const { return *val_; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      T& val() { return *val_; }

      //! Returns derivative array
      KOKKOS_INLINE_FUNCTION
      const T* dx() const { return dx_;}

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      T dx(int i) const {
        return sz_.value ? dx_[ stride_one ? i : i * stride_.value ] : T(0.);
      }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      T& fastAccessDx(int i) {
        return dx_[ stride_one ? i : i * stride_.value ];
      }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& fastAccessDx(int i) const {
        return dx_[ stride_one ? i : i * stride_.value ];
      }

      //! Overload of addressof operator
      KOKKOS_INLINE_FUNCTION
      ViewFadPtr<T,static_length,static_stride,U> operator&() const {
        return ViewFadPtr<T,static_length,static_stride,U>(
          this->dx_, this->val_, this->sz_.value, this->stride_.value);
      }

    protected:

      //! Derivative array size
      const mpl::integral_nonzero_constant< int, static_length > sz_;

      //! Derivative array stride
      const mpl::integral_nonzero_constant< int, static_stride > stride_;

      //! Value
      T *val_;

      //! Derivative array
      T *dx_;

    }; // class ViewStorage

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXP_VIEWSTORAGE_HPP
