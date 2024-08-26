// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_VIEWSTORAGE_HPP
#define SACADO_FAD_EXP_VIEWSTORAGE_HPP

#include <type_traits>
#include <utility>
#include <memory>

#include "Sacado_DynamicArrayTraits.hpp"
#include "Sacado_mpl_integral_nonzero_constant.hpp"
#include "Sacado_mpl_apply.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

#ifndef SACADO_FAD_DERIV_LOOP
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD) && !defined(SACADO_DISABLE_CUDA_IN_KOKKOS) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
#define SACADO_FAD_DERIV_LOOP(I,SZ) for (int I=threadIdx.x; I<SZ; I+=blockDim.x)
#else
#define SACADO_FAD_DERIV_LOOP(I,SZ) for (int I=0; I<SZ; ++I)
#endif
#endif

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
      static constexpr bool stride_one = (1 == static_stride);

    public:

      typedef typename std::remove_cv<T>::type value_type;
      static constexpr bool is_statically_sized = (static_length > 0);
      static constexpr int static_size = static_length;
      static constexpr bool is_view = true;
      typedef U base_fad_type;

      //! Turn ViewStorage into a meta-function class usable with mpl::apply
      template <typename TT>
      struct apply {
        typedef typename std::remove_cv<TT>::type non_const_TT;
        typedef typename BaseExprType<non_const_TT>::type base_expr_TT;
        typedef typename mpl::apply<U,base_expr_TT>::type UU;
        typedef ViewStorage<TT,static_length,static_stride,UU> type;
      };

      //! Replace static derivative length
      template <int N>
      struct apply_N {
        typedef ViewStorage<T,static_length,static_stride,U> type;
      };

      //! Default constructor (needed to satisfy interface)
      SACADO_INLINE_FUNCTION
      ViewStorage() :
        sz_(0), stride_(0), val_(0), dx_(0) {}

      //! Constructor
      SACADO_INLINE_FUNCTION
      ViewStorage(T* v, const int arg_size = 0, const int arg_stride = 0) :
        sz_(arg_size), stride_(arg_stride), val_(v+sz_.value*static_cast<int>(stride_.value)), dx_(v) {}

      //! Constructor
      SACADO_INLINE_FUNCTION
      ViewStorage(T* arg_dx, T* arg_val, const int arg_size = 0,
                  const int arg_stride = 0) :
        sz_(arg_size), stride_(arg_stride), val_(arg_val), dx_(arg_dx) {}

      //! Copy constructor
      /*!
        Allow, e.g., ViewStorage<double,...>(ViewStorage<double,...>),
        ViewStorage<const double,...>(ViewStorage<double,...>), and
        ViewStorage<const double,...>(ViewStorage<const double,...>) but not
        ViewStorage<double,...>(ViewStorage<const double,...>).
      */
      template <typename TT>
      SACADO_INLINE_FUNCTION
      ViewStorage(
        const ViewStorage<TT,static_length,static_stride,U>& x,
        typename std::enable_if< std::is_same<TT,T>::value ||
                                 std::is_same<const TT,T>::value >::type* = 0)
        : sz_(x.sz_), stride_(x.stride_), val_(x.val_), dx_(x.dx_) {}

      // Move does not make sense for this storage since it is always tied to
      // some preallocated data.  Don't define move constructor so compiler will
      // always fall-back to copy

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~ViewStorage() {}

      //! Assignment
      SACADO_INLINE_FUNCTION
      ViewStorage& operator=(const ViewStorage& x) {
        // Can't use std::addressof() on the GPU, so this is equivalent
        // according to cppreference.com
        if (this != reinterpret_cast<ViewStorage*>(
              &const_cast<char&>(
                reinterpret_cast<const volatile char&>(x)))) {
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

      // Move does not make sense for this storage since it is always tied to
      // some preallocated data.  Don't define move assignment so compiler will
      // always fall-back to copy

      //! Returns number of derivative components
      SACADO_INLINE_FUNCTION
      constexpr int size() const { return sz_.value;}

      //! Returns array length
      SACADO_INLINE_FUNCTION
      constexpr int length() const { return sz_.value; }

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
      void resizeAndZero(int sz) {}

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
        return unsigned(sz_.value) ? dx_[ stride_one ? i : i * stride_.value ] : T(0.);
      }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      T& fastAccessDx(int i) {
        return dx_[ stride_one ? i : i * stride_.value ];
      }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      const T& fastAccessDx(int i) const {
        return dx_[ stride_one ? i : i * stride_.value ];
      }

      //! Overload of addressof operator
      SACADO_INLINE_FUNCTION
      ViewFadPtr<T,static_length,static_stride,U> operator&() const {
        return ViewFadPtr<T,static_length,static_stride,U>(
          this->dx_, this->val_, this->sz_.value, this->stride_.value);
      }

    protected:

      template <typename, unsigned, unsigned, typename>
      friend class ViewStorage;

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
