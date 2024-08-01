// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_VIEWSTORAGE_HPP
#define SACADO_FAD_VIEWSTORAGE_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_ViewStorage.hpp"

namespace Sacado {
  namespace Fad {

    template <typename T, unsigned static_length, unsigned static_stride,
              typename U>
    using ViewStorage = Exp::ViewStorage<T,static_length,static_stride,U>;

  }
}

#else

#include "Sacado_Traits.hpp"
#include "Sacado_DynamicArrayTraits.hpp"
#include "Sacado_mpl_integral_nonzero_constant.hpp"

namespace Sacado {

  namespace Fad {

#ifndef SACADO_FAD_DERIV_LOOP
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && defined(__CUDA_ARCH__)
#define SACADO_FAD_DERIV_LOOP(I,SZ) for (int I=threadIdx.x; I<SZ; I+=blockDim.x)
#else
#define SACADO_FAD_DERIV_LOOP(I,SZ) for (int I=0; I<SZ; ++I)
#endif
#endif

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
      SACADO_INLINE_FUNCTION
      ViewStorage(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor with size \c sz (needed to satisfy interface)
      SACADO_INLINE_FUNCTION
      ViewStorage(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor
      SACADO_INLINE_FUNCTION
      ViewStorage(T* v, const int arg_size = 0, const int arg_stride = 0) :
        sz_(arg_size), stride_(arg_stride), val_(v+sz_.value*stride_.value), dx_(v) {}

      //! Constructor
      SACADO_INLINE_FUNCTION
      ViewStorage(T* arg_dx, T* arg_val, const int arg_size = 0,
                  const int arg_stride = 0) :
        sz_(arg_size), stride_(arg_stride), val_(arg_val), dx_(arg_dx) {}

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      ViewStorage(const ViewStorage& x) :
        sz_(x.sz_), stride_(x.stride_), val_(x.val_), dx_(x.dx_) {}

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~ViewStorage() {}

      //! Assignment
      SACADO_INLINE_FUNCTION
      ViewStorage& operator=(const ViewStorage& x) {
        if (this != &x) {
          *val_ = *x.val_;
          if (stride_one)
            //for (int i=0; i<sz_.value; ++i)
            SACADO_FAD_DERIV_LOOP(i,sz_.value)
              dx_[i] = x.dx_[i];
          else
            //for (int i=0; i<sz_.value; ++i)
            SACADO_FAD_DERIV_LOOP(i,sz_.value)
              dx_[i*stride_.value] = x.dx_[i*x.stride_.value];
        }
        return *this;
      }

      //! Returns number of derivative components
      SACADO_INLINE_FUNCTION
      int size() const { return sz_.value;}

      //! Returns array length
      SACADO_INLINE_FUNCTION
      int length() const { return sz_.value; }

      //! Resize the derivative array to sz
      /*!
       * Since we can't actually resize, we check for resizes to zero,
       * which signify assigning a constant.  Thus we zero out the derivative
       * components.
       */
      SACADO_INLINE_FUNCTION
      void resize(int sz) {
        if (sz == 0) ds_array<T>::strided_zero(dx_, stride_.value, sz_.value);
      }

      //! Resize the derivative array to sz
      /*!
       * We don't do anything here as this is used in the context of resizing
       * the derivative array to zero and then back to some size > 0.  Instead
       * we zero out components when it is resized to zero above.
       */
      SACADO_INLINE_FUNCTION
      void resizeAndZero(int /* sz */) {}

      //! Expand derivative array to size sz
      SACADO_INLINE_FUNCTION
      void expand(int sz) {}

      //! Zero out derivative array
      SACADO_INLINE_FUNCTION
      void zero() {
        ds_array<T>::strided_zero(dx_, stride_.value, sz_.value);
      }

      //! Returns value
      SACADO_INLINE_FUNCTION
      const T& val() const { return *val_; }

      //! Returns value
      SACADO_INLINE_FUNCTION
      T& val() { return *val_; }

      //! Returns derivative array
      SACADO_INLINE_FUNCTION
      const T* dx() const { return dx_;}

      //! Returns derivative component \c i with bounds checking
      SACADO_INLINE_FUNCTION
      T dx(int i) const {
        return unsigned(sz_.value) > 0 ? dx_[ unsigned(stride_one) == 1 ? i : i * stride_.value ] : T(0.);
      }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      T& fastAccessDx(int i) {
        return dx_[ unsigned(stride_one) == 1 ? i : i * stride_.value ];
      }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      const T& fastAccessDx(int i) const {
        return dx_[ unsigned(stride_one) == 1 ? i : i * stride_.value ];
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

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_VIEWSTORAGE_HPP
