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

#include "Sacado_mpl_apply.hpp"
#include "Sacado_mpl_if.hpp"
#include "Sacado_mpl_type_wrap.hpp"

namespace Sacado {

  namespace FAD_NS {

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
      KOKKOS_INLINE_FUNCTION
      ViewFad(const StorageType& s) :
        ExprType(s) {}

      //! View-specific constructor
      KOKKOS_INLINE_FUNCTION
      ViewFad(ValueT* v, const int arg_size = 0, const int arg_stride = 0) :
        ExprType( StorageType(v,arg_size,arg_stride) ) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~ViewFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator=(const S& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with ViewFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator=(const ViewFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator += (const S& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator -= (const S& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator *= (const S& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(ViewFad&) operator /= (const S& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with ViewFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator += (const ViewFad& x) {
        GeneralFadType::operator+=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with ViewFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator -= (const ViewFad& x) {
        GeneralFadType::operator-=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with ViewFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator *= (const ViewFad& x) {
        GeneralFadType::operator*=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Division-assignment operator with ViewFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator /= (const ViewFad& x) {
        GeneralFadType::operator/=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator += (const Expr<S>& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator -= (const Expr<S>& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator *= (const Expr<S>& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(ViewFad&) operator /= (const Expr<S>& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //@}

    }; // class ViewFad<ValueT>

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
  struct IsExpr< FAD_NS::ViewFad<T,l,s,U> > {
    static const bool value = true;
  };

  template <typename T, unsigned l, unsigned s, typename U>
  struct BaseExprType< FAD_NS::ViewFad<T,l,s,U> > {
    typedef typename FAD_NS::ViewFad<T,l,s,U>::base_expr_type type;
  };

} // namespace Sacado
