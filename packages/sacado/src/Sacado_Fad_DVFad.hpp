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

#ifndef SACADO_FAD_DVFAD_HPP
#define SACADO_FAD_DVFAD_HPP

#include "Sacado_Fad_GeneralFadExpr.hpp"
#include "Sacado_Fad_DVFadTraits.hpp"
#include "Sacado_Fad_VectorDynamicStorage.hpp"

namespace Sacado {

  namespace Fad {

    /*!
     * \brief Forward-mode AD class using dynamic memory allocation and
     * expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::Fad::GeneralVFad.
     */
    template <typename ValueT>
    class DVFad :
      public Expr< GeneralFad<ValueT,VectorDynamicStorage<ValueT> > > {

    public:

      //! Base classes
      typedef VectorDynamicStorage<ValueT> StorageType;
      typedef GeneralFad<ValueT,StorageType> GeneralFadType;
      typedef Expr<GeneralFadType> ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn DVFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef DVFad<T> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DVFad() :
        ExprType() {}

     //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      DVFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DVFad(const int sz, const ValueT& x) :
        ExprType(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DVFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Constructor with supplied memory
      /*!
       * Initializes value to point to \c x and derivative array to point
       * to\c dx.  Derivative array is zero'd out if \c zero_out is true.
       */
      DVFad(const int sz, ValueT* x, ValueT* dx, bool zero_out = false) :
        ExprType(sz,x,dx,zero_out) {}

      //! Constructor with supplied memory and index \c i
      /*!
       * Initializes value to point to \c x and derivative array to point
       * to\c dx.  Initializes derivative array row \c i of the identity matrix,
       * i.e., sets derivative component \c i to 1 and all other's to zero.
       */
      DVFad(const int sz, const int i, ValueT* x, ValueT* dx) :
        ExprType(sz,i,x,dx) {}

      //! Copy constructor
      DVFad(const DVFad& x) :
        ExprType(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      DVFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      ~DVFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(DVFad&) operator=(const S& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with DVFad right-hand-side
      DVFad& operator=(const DVFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(DVFad&) operator=(const Expr<S>& x)
      {
        GeneralFadType::operator=(x);
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DVFad&) operator += (const S& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DVFad&) operator -= (const S& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DVFad&) operator *= (const S& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DVFad&) operator /= (const S& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with DVFad right-hand-side
      DVFad& operator += (const DVFad& x) {
        GeneralFadType::operator+=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with DVFad right-hand-side
      DVFad& operator -= (const DVFad& x) {
        GeneralFadType::operator-=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with DVFad right-hand-side
      DVFad& operator *= (const DVFad& x) {
        GeneralFadType::operator*=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Division-assignment operator with DVFad right-hand-side
      DVFad& operator /= (const DVFad& x) {
        GeneralFadType::operator/=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DVFad&) operator += (const Expr<S>& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DVFad&) operator -= (const Expr<S>& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DVFad&) operator *= (const Expr<S>& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DVFad&) operator /= (const Expr<S>& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

    }; // class DVFad<ValueT>

    template <typename T>
    struct BaseExpr< GeneralFad<T,Fad::VectorDynamicStorage<T> > > {
      typedef DVFad< typename GeneralFad<T,Fad::VectorDynamicStorage<T> >::value_type > type;
    };

    template <typename T>
    struct ExprLevel< DVFad<T> > {
      static const unsigned value =
        ExprLevel< typename DFad<T>::value_type >::value + 1;
    };

    template <typename T>
    struct IsFadExpr< DVFad<T> > {
      static const bool value = true;
    };

  } // namespace Fad

  template <typename T>
  struct IsExpr< Fad::DVFad<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct BaseExprType< Fad::DVFad<T> > {
    typedef typename Fad::DVFad<T>::base_expr_type type;
  };

  //! The View Fad type associated with this type
  template< class ValueType, unsigned length, unsigned stride >
  struct ViewFadType< Sacado::Fad::DVFad< ValueType >, length, stride > {
    typedef Sacado::Fad::ViewFad< ValueType,length,stride,Sacado::Fad::DVFad<ValueType> > type;
  };

  //! The View Fad type associated with this type
  template< class ValueType, unsigned length, unsigned stride >
  struct ViewFadType< const Sacado::Fad::DVFad< ValueType >, length, stride > {
    typedef Sacado::Fad::ViewFad< const ValueType,length,stride,Sacado::Fad::DVFad<const ValueType> > type;
  };

} // namespace Sacado

#endif // SACADO_FAD_DFAD_HPP
