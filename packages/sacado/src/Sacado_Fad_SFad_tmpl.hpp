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

    //! Forward-mode AD class using static memory allocation
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The size
     * of the derivative array is fixed by the template parameter \c Num.
     */
    template <typename ValueT, int Num>
    class SFad :
      public Expr< SFadExprTag<ValueT,Num > > {

    public:

      //! Base classes
      typedef Expr< SFadExprTag<ValueT,Num > > ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef SFad<T,Num> type;
      };

      //! Replace static derivative length
      template <int N>
      struct apply_N {
        typedef SFad<ValueT,N> type;
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
      SFad() :
        ExprType() {}

      //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      SACADO_INLINE_FUNCTION
      SFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      SFad(const int sz, const ValueT & x, const DerivInit zero_out = InitDerivArray) :
        ExprType(sz,x,zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SACADO_INLINE_FUNCTION
      SFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      SFad(const SFad& x) :
        ExprType(static_cast<const ExprType&>(x)) {}

      //! Copy constructor from any Expression object
      template <typename S>
      SACADO_INLINE_FUNCTION
      SFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~SFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator=(const S& v) {
        ExprType::operator=(v);
        return *this;
      }

      //! Assignment operator with SFad right-hand-side
      SACADO_INLINE_FUNCTION
      SFad& operator=(const SFad& x) {
        ExprType::operator=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SFad&) operator=(const Expr<S>& x)
      {
        ExprType::operator=(x);
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator += (const S& x) {
        ExprType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator -= (const S& x) {
        ExprType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator *= (const S& x) {
        ExprType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator /= (const S& x) {
        ExprType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with SFad right-hand-side
      SACADO_INLINE_FUNCTION
      SFad& operator += (const SFad& x) {
        ExprType::operator+=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with SFad right-hand-side
      SACADO_INLINE_FUNCTION
      SFad& operator -= (const SFad& x) {
        ExprType::operator-=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with SFad right-hand-side
      SACADO_INLINE_FUNCTION
      SFad& operator *= (const SFad& x) {
        ExprType::operator*=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Division-assignment operator with SFad right-hand-side
      SACADO_INLINE_FUNCTION
      SFad& operator /= (const SFad& x) {
        ExprType::operator/=(static_cast<const ExprType&>(x));
        return *this;
      }

       //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SFad&) operator += (const Expr<S>& x) {
        ExprType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SFad&) operator -= (const Expr<S>& x) {
        ExprType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SFad&) operator *= (const Expr<S>& x) {
        ExprType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SFad&) operator /= (const Expr<S>& x) {
        ExprType::operator/=(x);
        return *this;
      }

    }; // class SFad<ValueT,Num>

    template <typename T, int Num>
    std::ostream& operator << (std::ostream& os,
                               const Expr< SFadExprTag<T,Num> >& x) {
      os << x.val() << " [";

      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

    template <typename T, int N>
    struct ExprLevel< SFad<T,N> > {
      static const unsigned value =
        ExprLevel< typename SFad<T,N>::value_type >::value + 1;
    };

    template <typename T, int N>
    struct IsFadExpr< SFad<T,N> > {
      static const bool value = true;
    };

  } // namespace Fad

  template <typename T, int N>
  struct IsFad< FAD_NS::SFad<T,N> > {
    static const bool value = true;
  };

  template <typename T, int N>
  struct IsExpr< FAD_NS::SFad<T,N> > {
    static const bool value = true;
  };

  template <typename T, int N>
  struct BaseExprType< FAD_NS::SFad<T,N> > {
    typedef typename FAD_NS::SFad<T,N>::base_expr_type type;
  };

  template <typename T,unsigned,unsigned> struct ViewFadType;
  namespace FAD_NS {
    template <typename,unsigned,unsigned,typename> class ViewFad;
  }

  //! The View Fad type associated with this type
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< Sacado::FAD_NS::SFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< ValueType,length,stride,Sacado::FAD_NS::SFad< ValueType,N> > type;
};

  //! The View Fad type associated with this type
  /*!
   * Do not include the const in the base expr type.
   */
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< const Sacado::FAD_NS::SFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< const ValueType,length,stride,Sacado::FAD_NS::SFad<ValueType,N> > type;
  };

} // namespace Sacado

#if defined(HAVE_SACADO_KOKKOS)

//-------------------------- Atomic Operators -----------------------

namespace Sacado {

  namespace FAD_NS {

    // Overload of Kokkos::atomic_add for Fad types.
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SFad<T,N>* dst, const SFad<T,N>& x) {
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
