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

#ifndef SACADO_FAD_VIEWSTORAGE_HPP
#define SACADO_FAD_VIEWSTORAGE_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_DynamicArrayTraits.hpp"
#include "Sacado_mpl_integral_nonzero_constant.hpp"

namespace Sacado {

  namespace Fad {

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

      typedef T value_type;

      //! Default constructor (needed to satisfy interface)
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      ViewStorage(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor with size \c sz (needed to satisfy interface)
      KOKKOS_INLINE_FUNCTION
      ViewStorage(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor
      KOKKOS_INLINE_FUNCTION
      ViewStorage(T* v, const int arg_size = 0, const int arg_stride = 0) :
        sz_(arg_size), stride_(arg_stride), val_(v+sz_.value*stride_.value), dx_(v) {}

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
        if (this != &x) {
          *val_ = *x.val_;
          if (stride_one)
            for (int i=0; i<sz_.value; ++i)
              dx_[i] = x.dx_[i];
          else
            for (int i=0; i<sz_.value; ++i)
              dx_[i*stride_.value] = x.dx_[i*stride_.value];
        }
        return *this;
      }

      //! Returns number of derivative components
      KOKKOS_INLINE_FUNCTION
      int size() const { return sz_.value;}

      //! Returns array length
      KOKKOS_INLINE_FUNCTION
      int length() const { return sz_.value; }

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

    private:

      //! Derivative array size
      const mpl::integral_nonzero_constant< int, static_length > sz_;

      //! Derivative array stride
      const mpl::integral_nonzero_constant< int, static_stride > stride_;

      //! Value
      T *val_;

      //! Derivative array
      T *dx_;

    }; // class ViewStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_VIEWSTORAGE_HPP
