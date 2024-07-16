// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_mpl_apply.hpp"
#include "Sacado_mpl_if.hpp"
#include "Sacado_mpl_type_wrap.hpp"

#if defined(HAVE_SACADO_KOKKOS)
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Error.hpp"
#endif

namespace Sacado {

  namespace FAD_NS {

    // Class representing a pointer to ViewFad so that &ViewFad is supported
    template <typename T, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr;

    /*!
     * \brief Forward-mode AD class using dynamic memory allocation and
     * expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::FAD_NS::GeneralFad.
     */
    template <typename ValueT, unsigned length, unsigned stride,
              typename BaseFadT>
    class ViewFad :
      public Expr< GeneralFad<ValueT,Fad::ViewStorage<ValueT,length,stride,BaseFadT> > > {

    public:

      //! Base classes
      typedef Fad::ViewStorage<ValueT,length,stride,BaseFadT> StorageType;
      typedef GeneralFad<ValueT,StorageType> GeneralFadType;
      typedef Expr<GeneralFadType> ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Fad type view is based on
      typedef BaseFadT base_fad_type;

      //! Turn ViewFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        // BaseExprType<T>::type is T if T is not a view
        typedef typename BaseExprType<T>::type T_for_base;
        typedef typename mpl::apply<base_fad_type,T_for_base>::type new_base_fad_type;
        typedef ViewFad<T,length,stride,new_base_fad_type> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      // ViewFad cannot be created with the usual constructors, so we remove
      // them here.

      //! Constructor with supplied storage \c s
      SACADO_INLINE_FUNCTION
      ViewFad(const StorageType& s) :
        ExprType(s) {}

      //! View-specific constructor
      SACADO_INLINE_FUNCTION
      ViewFad(ValueT* v, const int arg_size = 0, const int arg_stride = 0) :
        ExprType( StorageType(v,arg_size,arg_stride) ) {}

      //! View-specific constructor
      SACADO_INLINE_FUNCTION
      ViewFad(ValueT* dx_ptr, ValueT* val_ptr, const int arg_size = 0,
              const int arg_stride = 0) :
        ExprType( StorageType(dx_ptr,val_ptr,arg_size,arg_stride) ) {}

      //@}

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~ViewFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator=(const S& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with ViewFad right-hand-side
      SACADO_INLINE_FUNCTION
      ViewFad& operator=(const ViewFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator=(const Expr<S>& x)
      {
        GeneralFadType::operator=(x);
        return *this;
      }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator += (const S& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator -= (const S& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator *= (const S& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator /= (const S& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with ViewFad right-hand-side
      SACADO_INLINE_FUNCTION
      ViewFad& operator += (const ViewFad& x) {
        GeneralFadType::operator+=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with ViewFad right-hand-side
      SACADO_INLINE_FUNCTION
      ViewFad& operator -= (const ViewFad& x) {
        GeneralFadType::operator-=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with ViewFad right-hand-side
      SACADO_INLINE_FUNCTION
      ViewFad& operator *= (const ViewFad& x) {
        GeneralFadType::operator*=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Division-assignment operator with ViewFad right-hand-side
      SACADO_INLINE_FUNCTION
      ViewFad& operator /= (const ViewFad& x) {
        GeneralFadType::operator/=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator += (const Expr<S>& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator -= (const Expr<S>& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator *= (const Expr<S>& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator /= (const Expr<S>& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      // Overload of addressof operator
      SACADO_INLINE_FUNCTION
      ViewFadPtr<ValueT,length,stride,BaseFadT> operator&() const {
        return ViewFadPtr<ValueT,length,stride,BaseFadT>(
          this->dx_, this->val_, this->sz_.value, this->stride_.value);
      }

      //@}

    }; // class ViewFad<ValueT>

    // Class representing a pointer to ViewFad so that &ViewFad is supported
    template <typename T, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr : public ViewFad<T,sl,ss,U> {
    public:

      // Storage type base class
      typedef ViewFad<T,sl,ss,U> view_fad_type;

      // Bring in constructors
      using view_fad_type::view_fad_type;

      // Add overload of dereference operator
      SACADO_INLINE_FUNCTION
      view_fad_type* operator->() { return this; }

      // Add overload of dereference operator
      SACADO_INLINE_FUNCTION
      view_fad_type& operator*() { *this; }
    };

#if defined(HAVE_SACADO_KOKKOS)
    // Overload of Kokkos::atomic_add for ViewFad types.
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x) {
      using Kokkos::atomic_add;

      x.cache();

      const int xsz = x.size();
      const int sz = dst->size();

      // We currently cannot handle resizing since that would need to be
      // done atomically.
      if (xsz > sz)
        Kokkos::abort(
          "Sacado error: Fad resize within atomic_add() not supported!");

      if (xsz != sz && sz > 0 && xsz > 0)
        Kokkos::abort(
          "Sacado error: Fad assignment of incompatiable sizes!");


      if (sz > 0 && xsz > 0) {
        SACADO_FAD_DERIV_LOOP(i,sz)
          atomic_add(&(dst->fastAccessDx(i)), x.fastAccessDx(i));
      }
      SACADO_FAD_THREAD_SINGLE
        atomic_add(&(dst->val()), x.val());
    }
#endif

    template <typename T, unsigned l, unsigned s, typename U>
    struct BaseExpr< GeneralFad<T,Fad::ViewStorage<T,l,s,U> > > {
      //typedef ViewFad<T,l,s,U> type;
      typedef U type;
    };

    template <typename T, unsigned l, unsigned s, typename U>
    struct ExprLevel< ViewFad<T,l,s,U> > {
      static const unsigned value =
        ExprLevel< typename ViewFad<T,l,s,U>::value_type >::value + 1;
    };

    template <typename T, unsigned l, unsigned s, typename U>
    struct IsFadExpr< ViewFad<T,l,s,U> > {
      static const bool value = true;
    };

  } // namespace Fad

  template <typename T, unsigned l, unsigned s, typename U>
  struct IsView< Sacado::FAD_NS::ViewFad<T,l,s,U> > {
    static const bool value = true;
  };

  template <typename T, unsigned l, unsigned s, typename U>
  struct IsFad< FAD_NS::ViewFad<T,l,s,U> > {
    static const bool value = true;
  };

  template <typename T, unsigned l, unsigned s, typename U>
  struct IsExpr< FAD_NS::ViewFad<T,l,s,U> > {
    static const bool value = true;
  };

  template <typename T, unsigned l, unsigned s, typename U>
  struct BaseExprType< FAD_NS::ViewFad<T,l,s,U> > {
    typedef typename FAD_NS::ViewFad<T,l,s,U>::base_expr_type type;
  };

} // namespace Sacado
