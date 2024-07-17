// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#if defined(HAVE_SACADO_KOKKOS)
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Error.hpp"
#endif

namespace Sacado {

  namespace FAD_NS {

    // Forward declaration
    template <typename T, int Num>
    class StaticStorage;

    /*!
     * \brief Forward-mode AD class using static memory allocation
     * with long arrays and expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The largest size
     * of the derivative array is fixed by the template parameter \c Num
     * while the actual size used is set by the \c sz argument to the
     * constructor or the \c n argument to diff().  The user
     * interface is provided by Sacado::FAD_NS::GeneralFad.
     */
    template <typename ValueT, int Num>
    class SLFad :
      public Expr< GeneralFad<ValueT,Fad::StaticStorage<ValueT,Num> > > {

    public:

      //! Base classes
      typedef Fad::StaticStorage<ValueT,Num> StorageType;
      typedef GeneralFad<ValueT,StorageType> GeneralFadType;
      typedef Expr<GeneralFadType> ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SLFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef SLFad<T,Num> type;
      };

      //! Replace static derivative length
      /*! For SLFad, N is treated as the static dimension, so we don't change
       *  the array length.
       */
      template <int N>
      struct apply_N {
        typedef SLFad<ValueT,Num> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      SACADO_INLINE_FUNCTION
      SLFad() :
        ExprType() {}

      //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      SACADO_INLINE_FUNCTION
      SLFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      SLFad(const int sz, const ValueT & x, const DerivInit zero_out = InitDerivArray) :
        ExprType(sz,x,zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SACADO_INLINE_FUNCTION
      SLFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      SLFad(const SLFad& x) :
        ExprType(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      SACADO_INLINE_FUNCTION
      SLFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~SLFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator=(const S& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with SLFad right-hand-side
      SACADO_INLINE_FUNCTION
      SLFad& operator=(const SLFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator=(const Expr<S>& x)
      {
        GeneralFadType::operator=(x);
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator += (const S& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator -= (const S& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator *= (const S& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator /= (const S& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with SLFad right-hand-side
      SACADO_INLINE_FUNCTION
      SLFad& operator += (const SLFad& x) {
        GeneralFadType::operator+=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with SLFad right-hand-side
      SACADO_INLINE_FUNCTION
      SLFad& operator -= (const SLFad& x) {
        GeneralFadType::operator-=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with SLFad right-hand-side
      SACADO_INLINE_FUNCTION
      SLFad& operator *= (const SLFad& x) {
        GeneralFadType::operator*=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Division-assignment operator with SLFad right-hand-side
      SACADO_INLINE_FUNCTION
      SLFad& operator /= (const SLFad& x) {
        GeneralFadType::operator/=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator += (const Expr<S>& x) {
        GeneralFadType::operator+=(x);
        //*this = *this + x;
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator -= (const Expr<S>& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator *= (const Expr<S>& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator /= (const Expr<S>& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

    }; // class SLFad<ValueT,Num>

    template <typename T, int N>
    struct BaseExpr< GeneralFad<T,Fad::StaticStorage<T,N> > > {
      typedef SLFad<typename GeneralFad<T,Fad::StaticStorage<T,N> >::value_type,N> type;
    };

    template <typename T, int N>
    struct ExprLevel< SLFad<T,N> > {
      static const unsigned value =
        ExprLevel< typename SLFad<T,N>::value_type >::value + 1;
    };

    template <typename T, int N>
    struct IsFadExpr< SLFad<T,N> > {
      static const bool value = true;
    };

  } // namespace Fad

  template <typename T, int N>
  struct IsFad< FAD_NS::SLFad<T,N> > {
    static const bool value = true;
  };

  template <typename T, int N>
  struct IsExpr< FAD_NS::SLFad<T,N> > {
    static const bool value = true;
  };

  template <typename T, int N>
  struct BaseExprType< FAD_NS::SLFad<T,N> > {
    typedef typename FAD_NS::SLFad<T,N>::base_expr_type type;
  };

  template <typename T,unsigned,unsigned> struct ViewFadType;
  namespace FAD_NS {
    template <typename,unsigned,unsigned,typename> class ViewFad;
  }

  //! The View Fad type associated with this type
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< Sacado::FAD_NS::SLFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< ValueType,length,stride,Sacado::FAD_NS::SLFad<ValueType,N> > type;
  };

  //! The View Fad type associated with this type
  /*!
   * Do not include the const in the base expr type.
   */
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< const Sacado::FAD_NS::SLFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< const ValueType,length,stride,Sacado::FAD_NS::SLFad<ValueType,N > > type;
  };

} // namespace Sacado

#if defined(HAVE_SACADO_KOKKOS)

//-------------------------- Atomic Operators -----------------------

namespace Sacado {

  namespace FAD_NS {

    // Overload of Kokkos::atomic_add for Fad types.
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SLFad<T,N>* dst, const SLFad<T,N>& x) {
      using Kokkos::atomic_add;

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

  } // namespace Fad

} // namespace Sacado

#endif // HAVE_SACADO_KOKKOS
